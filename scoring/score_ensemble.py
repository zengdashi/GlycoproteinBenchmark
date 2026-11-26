#!/usr/bin/env python3
"""Score the ten candidate pairs from one method and input on the input-defined mask."""
import argparse
import csv
import json
import math
from pathlib import Path
from statistics import median

from coordinates import load_atoms
from kernels.ensemble_metrics import evaluate_candidate_pairs
from score import REPO, load_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True)
    parser.add_argument("--manifest", required=True, type=Path,
                        help="TSV with candidate_index, prediction, mapping. Paths relative to manifest.")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    _, topology = load_case(args.case, REPO / "data/glycan_connections.json")
    with args.manifest.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    indices = [int(r["candidate_index"]) for r in rows]
    if len(indices) != 5 or set(indices) != set(range(5)):
        parser.error("One row for each candidate_index 0,1,2,3,4 is required.")
    views = []
    for r in rows:
        source = load_atoms(args.manifest.parent / r["prediction"])
        coordinates = {}
        used = set()
        with (args.manifest.parent / r["mapping"]).open() as f:
            for atom in csv.DictReader(f, delimiter="\t"):
                key, serial = atom["canonical_atom"], atom["prediction_atom_id"].strip()
                if serial in {"", ".", "?"}:
                    continue
                if key in coordinates or serial in used or serial not in source:
                    raise ValueError(f"Missing or ambiguous prediction atom {key} in {r['prediction']}")
                xyz = source[serial][0]
                if not all(math.isfinite(x) for x in xyz):
                    raise ValueError(f"Nonfinite {key} in {r['prediction']}")
                coordinates[key] = xyz
                used.add(serial)
        views.append({"candidate_index": int(r["candidate_index"]), "atoms": coordinates})
    out = evaluate_candidate_pairs(topology, sorted(views, key=lambda v:v["candidate_index"]))
    out["case_id"] = args.case
    values = [p["global_protein_frame_macro_RMSD_A"] for p in out["pairs"]]
    out["median_pairwise_glycan_RMSD_protein_frame_A"] = (
        median(values) if all(v is not None for v in values) else None)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, separators=(",",":"), allow_nan=False) + "\n")
    print(json.dumps({"status": out["status"], "pairs": out["pair_count"],
                      "median_pairwise_glycan_RMSD_protein_frame_A": out["median_pairwise_glycan_RMSD_protein_frame_A"]}))


if __name__ == "__main__":
    main()
