#!/usr/bin/env python3
"""Join the anonymous reference map to a native model file on explicit atom identity."""
import argparse
import csv
import gzip
import json
from pathlib import Path

from kernels.native_names import canonical_mapping, read_structure

HERE = Path(__file__).resolve().parent
METHODS = ("af3", "boltz2", "chai1", "protenix", "rfaa", "sweetfold")
AMINO_ACIDS = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
                       "ARNDCQEGHILKMFPSTWYV"))


def native_map(method, case, coordinate_file, override=None):
    if override:
        rule = json.loads(Path(override).read_text())
    else:
        if method == "rfaa":
            raise ValueError("RFAA requires --native-map with the prediction's runtime_atom_map.json")
        with gzip.open(HERE / "kernels/native_ports.json.gz", "rt") as f:
            rule = json.load(f)[case["case_id"]][method]
    if method == "rfaa" and rule.get("schema") != "rfaa-anonymous-native-atom-map-v1":
        raise ValueError("RFAA requires the recorded runtime map, including native atomized token order")
    atoms = read_structure(coordinate_file)[0]
    assignment, errors, proxy, _ = canonical_mapping(method, case, rule, atoms)
    if errors:
        raise ValueError(f"Native mapping failed ({len(errors)} errors); first: {errors[0]}")
    result = {}
    for i, key in assignment.items():
        a = atoms[i]
        name = "/".join(key)
        if name in result:
            raise ValueError(f"Duplicate native atom mapping {name}")
        if key[0] == "protein" and key[-1] == "CA":
            pos = int(key[2])
            if not 1 <= pos <= len(case["protein"]["sequence"]):
                raise ValueError(f"Canonical position outside sequence at {pos}")
            if pos not in proxy and AMINO_ACIDS.get(a["resname"]) != case["protein"]["sequence"][pos - 1]:
                raise ValueError(f"Native and input amino-acid identity disagree at {pos}")
        result[name] = a["serial"]
    return result


def join(reference_map, prediction, output):
    opener = gzip.open if str(reference_map).endswith(".gz") else open
    with opener(reference_map, "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    reference = {r["canonical_atom"]: r["reference_atom_id"] for r in rows}
    if len(reference) != len(rows):
        raise ValueError("Duplicate canonical keys in reference map")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["canonical_atom", "reference_atom_id", "prediction_atom_id"])
        for key in sorted(set(reference) | set(prediction)):
            writer.writerow((key, reference.get(key, "."), prediction.get(key, ".")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True)
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--prediction", required=True, type=Path, help="Uncompressed native PDB/mmCIF")
    parser.add_argument("--reference-map", required=True, type=Path, help="From data/reference/export.py")
    parser.add_argument("--native-map", type=Path, help="The actual input/runtime map; required for RFAA")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads((HERE.parent / "data/glycan_connections.json").read_text())
    case = next(c for c in data["cases"] if c["case_id"] == args.case)
    prediction = native_map(args.method, case, args.prediction, args.native_map)
    join(args.reference_map, prediction, args.out)
    print(f"Mapped {len(prediction)} heavy atoms to {args.out.name}")


if __name__ == "__main__":
    main()
