"""Check joins, selection, published estimates and the exact five submitted table files."""
import csv
import gzip
import hashlib
import math
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from tempfile import TemporaryDirectory

from scripts.test_set_panel import counts

ROOT = Path(__file__).resolve().parent


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def require(test, label):
    if not test:
        raise ValueError(label)


def close(value, expected, label, tolerance=2e-7):
    require(math.isclose(float(value), float(expected), abs_tol=tolerance, rel_tol=1e-9),
            f"{label}: observed {value}, expected {expected}")


def sequence_mean(records, field):
    by_seq = defaultdict(list)
    for record in records:
        by_seq[record["sequence_id"]].append(float(record[field]))
    return mean(mean(values) for values in by_seq.values())


def check_tables():
    table_dir = ROOT / "tables"
    with TemporaryDirectory() as stage:
        subprocess.run([sys.executable, "-B", str(ROOT / "scripts/export_tables.py"),
                        "--output-dir", stage], check=True)
        manifest = (table_dir / "SHA256SUMS").read_text().splitlines()
        sizes = [12, 310, 192, 252, 156]
        require(len(manifest) == len(sizes), "Expected five submitted tables")
        for line, size in zip(manifest, sizes):
            digest, name = line.split()
            source = (table_dir / name).read_bytes()
            require(hashlib.sha256(source).hexdigest() == digest, name + " differs from submission")
            require((Path(stage) / name).read_bytes() == source, name + " table export differs")
            require(len(rows(table_dir / name)) == size, name + " row count")
    return {n: rows(next(table_dir.glob(f"S{n}_*.tsv"))) for n in range(1, 6)}


def case_observations():
    cases = {r["case_id"]: r for r in rows(ROOT / "data/cases.tsv") if r["common_six"] == "1"}
    settings = {r["method"] for r in rows(ROOT / "data/predictor_settings.tsv")}
    groups = defaultdict(list)
    for r in rows(ROOT / "scores/candidates.tsv.gz"):
        case = cases[r["case_id"]]
        require(r["method"] in settings and
                all(r[k] == case[k] for k in ("in_B", "group_id", "sequence_id")), r["candidate_id"])
        groups[r["method"], r["case_id"]].append(r)
    require(len(groups) == 1962 and len(settings) == 6 and len(cases) == 327, "Input group coverage")
    inputs = rows(ROOT / "scores/inputs.tsv")
    require(len(inputs) == len(groups) == len({(r["method"], r["case_id"]) for r in inputs}),
            "Input table must match each method and case once")
    selections = {}
    selected_sets, input_sets = defaultdict(list), defaultdict(list)
    for r in inputs:
        key = r["method"], r["case_id"]
        draws = groups[key]
        require(sorted(int(x["candidate_index"]) for x in draws) == list(range(5)), str(key))
        picks = [x for x in draws if x["is_native_top1"] == "1"]
        require(len(picks) == 1 and picks[0]["candidate_index"] == r["top1_candidate_index"], str(key))
        pick = picks[0]
        require(all(r[k] == pick[k] for k in ("in_B", "group_id", "sequence_id")), str(key))
        close(pick["native_score"], max(float(x["native_score"]) for x in draws), str(key), 1e-12)
        close(pick["native_score"], r["score"], str(key))
        metric = "glycan_protein_frame_macro_RMSD_A"
        errors = [float(x[metric]) for x in draws]
        for col, value in (("top1_rmsd", float(pick[metric])), ("mean5_rmsd", mean(errors)),
                           ("best5_rmsd", min(errors))):
            close(value, r[col], f"{key} {col}")
        selections[key] = pick
        for test_set in (("A", "B") if r["in_B"] == "1" else ("A",)):
            selected_sets[test_set, key[0]].append(pick)
            input_sets[test_set, key[0]].append(r)

    by_length = defaultdict(list)
    site_keys = set()
    site_counts = Counter()
    for r in rows(ROOT / "scores/selected_sites.tsv.gz"):
        key = r["method"], r["case_id"]
        require(r["candidate_id"] == selections[key]["candidate_id"], "Site/candidate mismatch")
        site = (*key, r["glycan_id"])
        require(site not in site_keys, "Duplicate site " + str(site))
        site_keys.add(site)
        site_counts[key] += 1
        case = cases[key[1]]
        for name in (("A", "B") if case["in_B"] == "1" else ("A",)):
            by_length[name, key[0], r["node_count"]].append(r)
    require(len(site_keys) == 5304 and all(site_counts[key] == int(cases[key[1]]["n_sites"])
            for key in selections), "Site coverage")
    return selected_sets, input_sets, by_length


def check_published(tables, selected, inputs, by_length):
    columns = {f"glycan_{paper}_aligned_RMSD_A": f"glycan_{score}_frame_macro_RMSD_A"
               for paper, score in (("self", "self"), ("root_ring", "root"),
                                    ("local_protein", "site"), ("whole_protein", "protein"))}
    columns.update({k: k for k in ("protein_CA_RMSD_A", "protein_CA_lDDT",
                                   "contact_precision_4.5A", "contact_recall_4.5A", "contact_F1_4.5A")})
    for r in tables[1]:
        sample = selected[r["test_set"], r["method"]]
        require(len(sample) == int(r["n_inputs"]) and
                len({x["sequence_id"] for x in sample}) == int(r["n_sequences"]), str(r))
        for target, source in columns.items():
            close(sequence_mean(sample, source), r[target], f"S1 {r['method']} {target}")
    for r in tables[2]:
        if r["metric"] == "chemical_extended_local_pass_fraction":
            sample = selected[r["test_set"], r["method"]]
            close(100 * sequence_mean(sample, r["metric"]), r["value"], "S2 " + r["method"])
    for r in tables[3]:
        if r["outcome"] == "chem_pass_geo_pass":
            sample = selected[r["test_set"], r["method"]]
            col = "joint_extended_local_" + format(float(r["positioning_cutoff_A"]), ".1f") + "A"
            close(100 * sequence_mean(sample, col), r["percent"], "S3 " + r["method"])
    for r in tables[4]:
        if r["analysis"] == "by_length":
            sample = by_length[r["test_set"], r["method"], r["length"]]
            require(len(sample) == int(r["n_sites"]), "S4 site count")
            within_case = defaultdict(list)
            for site in sample:
                within_case[site["sequence_id"], site["case_id"]].append(float(site["self_frame_rmsd_A"]))
            case_means = [{"sequence_id": seq, "value": mean(values)}
                          for (seq, case), values in within_case.items()]
            close(sequence_mean(case_means, "value"), r["value_A"], "S4 " + r["length"])
    choices = {"native_top1": "top1_rmsd", "random_expected": "mean5_rmsd",
               "full_candidate_oracle": "best5_rmsd"}
    for r in tables[5]:
        if r["analysis"] == "candidate_choice" and r["measure"] == "mean_positioning_RMSD":
            col = choices[r["stratum"]]
            close(sequence_mean(inputs[r["test_set"], r["method"]], col), r["value"], "S5 " + col)


def main():
    subprocess.run([sys.executable, "-B", str(ROOT / "data/validate.py")], check=True)
    subprocess.run([sys.executable, "-B", str(ROOT / "data/reference/validate.py")], check=True)
    fig1 = counts()
    for name, expected in (("A", (327, 884, 484, 299, 101)), ("B", (57, 171, 78, 63, 30))):
        actual = tuple(fig1[name][k] for k in ("Inputs", "Sites", "Monosaccharide", "Linear", "Branched"))
        require(actual == expected, "Figure 1B " + name)
    tables = check_tables()
    check_published(tables, *case_observations())
    print("Verified 9,810 candidates, 1,962 native selections and 5,304 selected-site records.")
    print("Figure 1B, A/B score estimates and all five submitted TSV files passed.")


if __name__ == "__main__":
    main()
