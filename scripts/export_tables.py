"""Export the five manuscript tables from the verified benchmark and figure data."""
import argparse
import csv
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
AGG = ROOT / "scores/aggregates"
METHODS = {"af3": "AlphaFold 3", "boltz2": "Boltz-2", "chai1": "Chai-1",
           "protenix": "Protenix", "rfaa": "RFAA", "sweetfold": "SweetFold"}
SETS = ("A", "B")


def read(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def method(name):
    if name in METHODS.values() or name == "Experiment":
        return name
    return METHODS[name.lower().replace("-", "")]


def percent(value):
    if value in ("NA", "", None):
        return "NA"
    return f"{Decimal(value) * 100:.9f}".rstrip("0").rstrip(".")


def write(folder, name, columns, rows):
    if not rows:
        raise ValueError("Empty supplementary table: " + name)
    with (folder / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns.split(), delimiter="\t",
                                lineterminator="\n", restval="NA")
        writer.writeheader()
        writer.writerows(rows)


def structure(folder):
    metrics = {
        "protein_CA_lDDT": "protein_CA_lDDT",
        "glycan_self_frame_macro_RMSD_A": "glycan_self_aligned_RMSD_A",
        "glycan_root_frame_macro_RMSD_A": "glycan_root_ring_aligned_RMSD_A",
        "glycan_site_frame_macro_RMSD_A": "glycan_local_protein_aligned_RMSD_A",
        "glycan_protein_frame_macro_RMSD_A": "glycan_whole_protein_aligned_RMSD_A",
        "contact_precision_4.5A": "contact_precision_4.5A",
        "contact_recall_4.5A": "contact_recall_4.5A",
        "contact_F1_4.5A": "contact_F1_4.5A",
    }
    estimates = defaultdict(dict)
    for r in read(AGG / "performance.tsv"):
        if r["analysis"] == "method_mean" and r["metric"] in metrics:
            row = estimates[r["cohort"], method(r["method"])]
            row.update(test_set=r["cohort"], method=method(r["method"]), n_inputs=r["cases"],
                       n_sequences=r["queries"], n_sites=r["sites"])
            row[metrics[r["metric"]]] = r["value"]
    for r in read(AGG / "candidate_statistics.tsv"):
        if r["endpoint"] == "protein_CA_RMSD_A":
            estimates[r["cohort"], method(r["method"])]["protein_CA_RMSD_A"] = r["value"]
    for r in read(AGG / "backbone_correlations.tsv"):
        estimates[r["test_set"], method(r["method"])].update(
            backbone_lddt_vs_local_glycan_RMSD_rho=r["weighted_spearman_rho"],
            rho_ci95_low=r["ci95_low"], rho_ci95_high=r["ci95_high"],
            n_lddt_ge_0p90=r["n_inputs_lddt_ge_0p90"],
            n_lddt_ge_0p90_and_local_RMSD_gt_5A=r["n_inputs_lddt_ge_0p90_local_rmsd_gt_5"])
    for r in read(AGG / "alternative_references.tsv"):
        if r["row_type"] == "method_summary" and r["frame"] == "protein":
            estimates[r["cohort"], method(r["method"])].update(
                n_multireference_inputs=r["cases"],
                mean_positioning_RMSD_reduction_A=r["sequence_weighted_mean_A"],
                max_single_input_RMSD_reduction_A=r["maximum_input_reduction_A"])
    rows = [estimates[c, m] for c in SETS for m in METHODS.values()]
    for r in rows:
        assert set(metrics.values()) <= r.keys() and "protein_CA_RMSD_A" in r
        assert "backbone_lddt_vs_local_glycan_RMSD_rho" in r
        assert ("n_multireference_inputs" in r) == (r["test_set"] == "A")
    write(folder, "S1_structure_and_contacts.tsv",
          "test_set method n_inputs n_sequences n_sites protein_CA_RMSD_A " +
          " ".join(metrics.values()) + " backbone_lddt_vs_local_glycan_RMSD_rho rho_ci95_low "
          "rho_ci95_high n_lddt_ge_0p90 n_lddt_ge_0p90_and_local_RMSD_gt_5A n_multireference_inputs "
          "mean_positioning_RMSD_reduction_A max_single_input_RMSD_reduction_A", rows)


def chemistry(folder):
    rows = []
    for r in read(AGG / "bond_angle_and_site_statistics.tsv"):
        rows.append(dict(test_set="A", method=r["method"], subset=r["group"], metric=r["metric"],
                         statistic=r["statistic"], value=r["value"], unit=r["unit"], n=r["n"],
                         observation_unit=r["counted"], n_inputs=r["n_inputs"], n_sequences=r["n_sequences"]))
    for r in read(AGG / "attachment_lengths.tsv"):
        if r["panel"] == "A":
            rows.append(dict(test_set=r["test_set"], method=r["method"], subset="all sites",
                             metric=r["metric"], statistic="count", unit="sites", value=r["value"],
                             n=r["raw_denominator_sites"], observation_unit="sites",
                             n_inputs=r["n_inputs"], n_sequences=r["n_protein_sequences"]))
    failure_types = {"fail_root", "fail_internal", "fail_sugar_stereo",
                     "fail_clash", "fail_sugar_stereo_and_clash"}
    units = {"root_supported": "protein–sugar bonds", "internal_supported": "sugar–sugar bonds",
             "tree_supported": "oligosaccharide sites", "anomer_correct": "anomeric centers"}
    for r in read(AGG / "chemical_checks.tsv"):
        if r["selection"] != "native_top1" or r["subset"] != "all":
            continue
        metric = r["metric"]
        if metric not in failure_types and not (r["cohort"] == "B" and metric in units):
            continue
        rows.append(dict(test_set=r["cohort"], method=method(r["method"]),
                         subset="sites with a chemistry or positioning error" if metric in failure_types else "All",
                         metric=metric, statistic="sequence-weighted mean", unit="%", value=percent(r["value"]),
                         ci95_low=percent(r["lo"]), ci95_high=percent(r["hi"]),
                         n=r["n"], observation_unit=units.get(metric, "sites"),
                         n_inputs=r["n_cases"], n_sequences=r["n_queries"], raw_numerator=r["raw_sum"]))
    for r in read(AGG / "performance.tsv"):
        if r["cohort"] == "B" and r["analysis"] == "method_mean" and r["metric"] == "chemical_extended_local_pass_fraction":
            rows.append(dict(test_set="B", method=method(r["method"]), subset="All", metric=r["metric"],
                             statistic="sequence-weighted mean", unit="%", value=percent(r["value"]),
                             n=r["sites"], observation_unit="glycosylation sites", n_inputs=r["cases"],
                             n_sequences=r["queries"]))
    for r in read(AGG / "ring_categories.tsv"):
        if r["cohort"] == "B":
            rows.append(dict(test_set="B", method=method(r["method"]), subset="All", metric="ring_agreement",
                             statistic="sequence-weighted mean", unit="%", value=percent(r["value"]),
                             n=r["comparison_n"], observation_unit="classifiable rings", n_inputs=r["case_n"],
                             n_sequences=r["query_n"], raw_numerator=r["match_n"]))
    for r in read(AGG / "torsion_medians.tsv"):
        rows.append(dict(test_set="A", method=method(r["method"]), subset=r["link_type"],
                         metric=r["angle"] + "_absolute_periodic_error", statistic="median",
                         value=r["median_absolute_error_deg"], unit="degrees", n=r["valid_edges"],
                         observation_unit="linkages"))
    for r in read(AGG / "linkage_omega.tsv"):
        if r["selection"] == "native_top1" and r["metric"] == "omega_abs":
            rows.append(dict(test_set=r["cohort"], method=method(r["method"]), subset=r["signature"],
                             metric="omega_absolute_periodic_error", statistic="sequence-weighted mean",
                             value=r["value"], unit="degrees", ci95_low=r["lo"], ci95_high=r["hi"],
                             n=r["n"], observation_unit="linkages", n_inputs=r["n_cases"], n_sequences=r["n_queries"]))
    write(folder, "S2_site_chemistry.tsv",
          "test_set method subset metric statistic value unit ci95_low ci95_high n observation_unit "
          "n_inputs n_sequences raw_numerator", rows)


def joint_outcomes(folder):
    groups = {"chem_pass_geo_pass", "chem_pass_geo_fail", "chem_fail_geo_pass", "chem_fail_geo_fail"}
    rows = []
    for r in read(AGG / "joint_outcomes.tsv"):
        if r["chemistry_profile"] == "extended_local" and r["outcome"] in groups:
            rows.append(dict(test_set=r["cohort"], method=method(r["method"]),
                             positioning_cutoff_A=r["protein_frame_rmsd_cutoff_A"], outcome=r["outcome"],
                             percent=percent(r["value"]), ci95_low_percent=percent(r["ci95_low"]),
                             ci95_high_percent=percent(r["ci95_high"]), observed_sites=r["raw_numerator"],
                             evaluated_sites=r["raw_denominator"], n_inputs=r["cases"], n_sequences=r["queries"]))
    assert len(rows) == len(SETS) * len(METHODS) * 4 * len(groups)
    write(folder, "S3_joint_chemistry_and_placement.tsv",
          "test_set method positioning_cutoff_A outcome percent ci95_low_percent ci95_high_percent "
          "observed_sites evaluated_sites n_inputs n_sequences", rows)


def length_and_branching(folder):
    rows = []
    for r in read(AGG / "length.tsv"):
        if r["endpoint"] == "self_frame_rmsd_A":
            rows.append(dict(test_set=r["cohort"], method=method(r["method"]), analysis="by_length",
                             length=r["length"], topology="all", metric=r["endpoint"], value_A=r["value"],
                             ci95_low_A=r["ci_low"], ci95_high_A=r["ci_high"], n_sites=r["n_sites"],
                             n_inputs=r["n_cases"], n_sequences=r["n_queries"]))
    for r in read(AGG / "topology.tsv"):
        branch = (r["unit"] == "site" and r["metric"] == "self_frame_rmsd_A"
                  and r["analysis"] in {"unadjusted_all_oligosaccharides", "length_standardized"})
        if branch:
            rows.append(dict(test_set=r["cohort"], method=method(r["method"]), analysis=r["analysis"],
                             length=r["length"], topology=r["statistic"],
                             metric=r["metric"], value_A=r["value"], ci95_low_A=r["ci95_low"],
                             ci95_high_A=r["ci95_high"], n_sites=r["sites"], n_inputs=r["cases"],
                             n_sequences=r["queries"], linear_sites=r["linear_sites"],
                             branched_sites=r["branched_sites"], length_weights=r["length_weights"]))
    write(folder, "S4_glycan_length_and_branching.tsv",
          "test_set method analysis length topology metric value_A ci95_low_A ci95_high_A n_sites "
          "n_inputs n_sequences linear_sites branched_sites length_weights", rows)


def selection(folder):
    rows = []
    for r in read(AGG / "rank_quartiles.tsv"):
        if r["endpoint"] == "rmsd":
            rows.append(dict(test_set=r["cohort"], method=method(r["method"]), analysis="score_quartile",
                             measure="mean_positioning_RMSD", stratum="Q" + r["bin"], value=r["value"],
                             unit="angstrom", ci95_low=r["ci_low"], ci95_high=r["ci_high"],
                             min_score=r["min_score"], max_score=r["max_score"],
                             n_inputs=r["n_cases"], n_sequences=r["n_queries"]))
    choices = {"native_rmsd": "native_top1", "mean5_rmsd": "random_expected", "best5_rmsd": "full_candidate_oracle"}
    published = {(r["cohort"], method(r["method"]), r["selection"]): float(r["value"])
                 for r in read(AGG / "selection_checks.tsv")
                 if r["analysis"] == "candidate_selection" and r["target"] == "glycan_protein_frame_macro_RMSD_A"}
    for r in read(AGG / "candidate_statistics.tsv"):
        ep = r["endpoint"]
        if ep in choices:
            reference = published[r["cohort"], method(r["method"]), choices[ep]]
            assert abs(float(r["value"]) - reference) < 1e-8
            label, measure, stratum, unit = "candidate_choice", "mean_positioning_RMSD", choices[ep], "angstrom"
        elif ep in {"dispersion_rho", "dispersion_partial_rho"}:
            label, measure, stratum, unit = "candidate_agreement", ep, "all", "correlation"
        else:
            continue
        rows.append(dict(test_set=r["cohort"], method=method(r["method"]), analysis=label,
                         measure=measure, stratum=stratum, value=r["value"], unit=unit,
                         ci95_low=r["ci_low"], ci95_high=r["ci_high"], n_inputs=r["n_cases"],
                         n_sequences=r["n_queries"]))
    cases = defaultdict(list)
    for r in read(ROOT / "scores/inputs.tsv"):
        cases["A", r["method"]].append(r)
        assert r["in_B"] in {"0", "1"}
        if r["in_B"] == "1":
            cases["B", r["method"]].append(r)
    for r in read(AGG / "score_spread.tsv"):
        block = cases[r["cohort"], method(r["method"])]
        assert len(block) == int(r["n_inputs"])
        base = dict(test_set=r["cohort"], method=method(r["method"]), stratum="all",
                    n_inputs=r["n_inputs"], n_sequences=r["n_sequences"])
        for metric, unit in (("median_within_input_range", "native score"), ("between_input_top1_iqr", "native score"),
                             ("median_range_over_between_input_top1_iqr", "ratio")):
            rows.append(dict(base, analysis="score_spread", measure=metric, value=r[metric], unit=unit))
        better = sum(float(x["best5_rmsd"]) < float(x["top1_rmsd"]) for x in block)
        rows.append(dict(base, analysis="candidate_choice", measure="inputs_with_better_candidate",
                         value=100 * better / len(block), unit="%", raw_numerator=better))
    write(folder, "S5_scores_selection_and_agreement.tsv",
          "test_set method analysis measure stratum value unit ci95_low ci95_high min_score "
          "max_score n_inputs n_sequences raw_numerator", rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "build/tables")
    destination = parser.parse_args().output_dir
    destination.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=destination) as tmp:
        stage = Path(tmp)
        for build in (structure, chemistry, joint_outcomes, length_and_branching, selection):
            build(stage)
        for table in sorted(stage.glob("*.tsv")):
            table.replace(destination / table.name)


if __name__ == "__main__":
    main()
