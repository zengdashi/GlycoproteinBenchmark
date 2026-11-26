#!/usr/bin/env python3
"""Score one mapped benchmark prediction, or verify the released manuscript example."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import fmean

from chemistry import chemical_scores
from coordinates import read_mapping
from kernels import geometry_metrics, physics_metrics
from kernels.ensemble_metrics import evaluate_tm
from topology import build_topology

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EXAMPLE = HERE / "example"
EXAMPLE_CASE = "GLYCO_38627a5f10e0e98f"
EXAMPLE_ID = "chai1__" + EXAMPLE_CASE + "__c001"


def load_case(case_id, catalog):
    data = json.loads(Path(catalog).read_text())
    case = next((case for case in data["cases"] if case["case_id"] == case_id), None)
    if case is None:
        raise ValueError(f"Case {case_id} missing from catalog")
    return case, build_topology(case, data["ccd_components"])


def average(values, require_all=True):
    values = list(values)
    finite = [v for v in values if v is not None and math.isfinite(v)]
    return fmean(finite) if finite and (not require_all or len(finite) == len(values)) else None


def verify_reference_mask(case_id, observed):
    path = REPO / "data/reference/records" / (case_id + ".tsv.gz")
    if not path.is_file():
        raise ValueError(f"No primary reference for {case_id}. Supply --alternative-reference for a different reference.")
    with gzip.open(path, "rt") as f:
        expected = {r["canonical_atom"] for r in csv.DictReader(f, delimiter="\t") if r["x"] != "."}
    actual = set(observed["atoms"])
    if actual != expected:
        raise ValueError(f"Reference mask mismatch: {len(expected - actual)} missing and {len(actual - expected)} extra keys")


def score(case, ref, observed, pred, tmalign=None):
    ref.update(observed)
    # Frozen physics caches by reference; include both the coordinate and canonical masks.
    mask_digest = hashlib.sha256(json.dumps(observed["atoms"], sort_keys=True).encode()).hexdigest()
    ref["reference_source"] = {"reference_only_sha256": mask_digest}
    geometry = geometry_metrics.evaluate_geometry(ref, pred, case)
    neighbors = physics_metrics._prepare(ref)
    contacts = physics_metrics.contacts(ref, pred, neighbors)
    clashes = physics_metrics.clashes(ref, pred, neighbors)
    chemistry = chemical_scores(ref, pred, geometry, clashes)
    ca = geometry["protein"]
    summary = {"protein_CA_RMSD_A": ca["ca_global"]["rmsd_A"], "protein_CA_lDDT": ca["ca_lddt"]["score"]}
    for frame in ["protein", "site", "root", "self"]:
        summary[f"glycan_{frame}_frame_macro_RMSD_A"] = average(
            x["frames"][frame]["rmsd_A"] for x in geometry["sites"].values())
    for metric in ["precision", "recall", "F1"]:
        summary[f"contact_{metric}_4.5A"] = average(
            (s["by_threshold_A"]["4.5"][metric] for s in contacts["sites"].values()), require_all=False)
    if chemistry["status"] == "scored":
        chem_sites = list(chemistry["sites"].values())
        summary["chemical_extended_local_pass_fraction"] = average(
            None if s["extended_status"] == "unknown" else float(s["extended_status"] == "pass") for s in chem_sites)
        for cutoff in (1., 2., 3., 5.):
            summary[f"joint_extended_local_{cutoff}A"] = average(s["joint"][str(cutoff)] for s in chem_sites)
    result = {"case_id": case["case_id"], "summary": summary, "geometry": geometry,
              "contacts": contacts, "clashes": clashes, "chemistry": chemistry}
    if tmalign is not None:
        result["fixed_correspondence_TM_score"] = evaluate_tm(ref, pred, executable=tmalign)
    return result


def close(name, value, expected, tolerance):
    if value is None or not math.isclose(value, float(expected), rel_tol=0., abs_tol=tolerance):
        raise ValueError(f"Score mismatch {name}: {value!r} versus {expected!r}")


def verify_example(result):
    if result["chemistry"]["status"] != "scored":
        raise RuntimeError("The full example check requires RDKit. Install the repository requirements.")
    tests = 0
    with gzip.open(REPO / "scores/candidates.tsv.gz", "rt") as f:
        published = next(r for r in csv.DictReader(f, delimiter="\t") if r["candidate_id"] == EXAMPLE_ID)
    for key, val in result["summary"].items():
        close(key, val, published[key], 1e-7)
        tests += 1
    # Original Figure 4d framing values were exported with six decimal places.
    frames = {"protein": 0.546103145549088, "site": 0.401107, "root": 0.184264, "self": 0.114920}
    for frame, expected in frames.items():
        val = result["geometry"]["sites"]["G002"]["frames"][frame]["rmsd_A"]
        close("G002 " + frame, val, expected, 5.1e-7)
        tests += 1
    with gzip.open(REPO / "scores/selected_sites.tsv.gz", "rt") as f:
        sites = [row for row in csv.DictReader(f, delimiter="\t") if row["candidate_id"] == EXAMPLE_ID]
    for row in sites:
        gid = row["glycan_id"]
        site = result["geometry"]["sites"][gid]
        root = next(r for r in result["geometry"]["linkages"]
                    if r["glycan_id"] == gid and r["kind"] == "protein_glycan")
        bond, angle = root["bond_geometry"], root["root_attachment_angle"]
        values = {"self_frame_rmsd_A": site["frames"]["self"]["rmsd_A"],
                  "glycan_global": site["frames"]["protein"]["rmsd_A"],
                  "glycan_local": site["frames"]["site"]["rmsd_A"],
                  "sidechain_local": result["geometry"]["protein"]["sites"][gid]["sidechain_site_frame"]["rmsd_A"],
                  "root_length_ref": bond["reference"], "root_length_pred": bond["prediction"],
                  "root_length_abs": bond["absolute_error"], "root_angle_abs": angle["absolute_error"],
                  "root_angle_ref": angle["reference"], "root_angle_pred": angle["prediction"]}
        for t in root["torsions"]:
            if t["name"] in {"phi","psi"}:
                values["root_" + t["name"] + "_abs"] = t["absolute_error"]
        for field, val in values.items():
            close(gid + " " + field, val, row[field], 5.1e-5)
            tests += 1
    return tests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--example", action="store_true", help="Recompute the Chai-1 example and check the score tables.")
    parser.add_argument("--case", help="Case ID in data/glycan_connections.json.")
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--prediction", type=Path)
    parser.add_argument("--mapping", type=Path, help="TSV with canonical_atom, reference_atom_id and prediction_atom_id.")
    parser.add_argument("--alternative-reference", action="store_true",
                        help="Score a different reference mask; the published scores use the primary mask.")
    parser.add_argument("--catalog", type=Path, default=REPO / "data/glycan_connections.json")
    parser.add_argument("--tmalign", type=Path, help="Optionally compute the fixed-correspondence TM-score.")
    parser.add_argument("--out", type=Path, help="Write all per-site, per-atom and diagnostic values as JSON (optionally .gz).")
    args = parser.parse_args()
    if args.example:
        if any((args.case, args.reference, args.prediction, args.mapping)):
            parser.error("--example and custom coordinate arguments are mutually exclusive")
        args.case = EXAMPLE_CASE
        args.reference = EXAMPLE / "9N3S.cif"
        args.prediction = EXAMPLE / "Chai1_c001.cif"
        args.mapping = EXAMPLE / "atom_map.tsv"
    if not all((args.case, args.reference, args.prediction, args.mapping)):
        parser.error("Specify --example or --case, --reference, --prediction and --mapping together.")
    case, topology = load_case(args.case, args.catalog)
    ref, pred = read_mapping(args.mapping, args.reference, args.prediction, topology)
    if not args.alternative_reference:
        verify_reference_mask(args.case, ref)
    result = score(case, topology, ref, pred, args.tmalign)
    result["reference_mask"] = "user_supplied_alternative" if args.alternative_reference else "benchmark_primary"
    if args.example:
        result["example_checks_passed"] = verify_example(result)
    print(json.dumps({k: v for k, v in result.items() if k in {"case_id","summary","example_checks_passed"}}
                     | {"chemistry_status": result["chemistry"]["status"]}, indent=2, allow_nan=False))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        opener = gzip.open if args.out.suffix == ".gz" else open
        with opener(args.out, "wt") as f:
            json.dump(result, f, separators=(",",":"), allow_nan=False)
            f.write("\n")


if __name__ == "__main__":
    main()
