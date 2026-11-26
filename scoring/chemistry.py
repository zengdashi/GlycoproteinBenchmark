"""Compute the manuscript's site chemistry and joint pass from the mapped model."""
import gzip
import json
import math
from pathlib import Path

from kernels import chemical_metrics as frozen

HERE = Path(__file__).resolve().parent


def assign_products(reference, prediction, case_products, templates):
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    centers, coverage = [], {}
    for gid, site in reference["sites"].items():
        if gid not in case_products:
            coverage[gid] = "product_template_unavailable"
            continue
        token, node_map = case_products[gid]
        template = templates[token]
        keys = []
        for origin, _ in template["atoms"]:
            unit, name = origin.split(":", 1)
            if unit == "AA":
                key = f"protein/P1/{site['position']}/{name}"
            else:
                key = f"glycan/{gid}/{node_map[unit]}/{name}"
            keys.append(key)
        if any(key not in prediction["atoms"] for key in keys):
            coverage[gid] = "missing_product_atom"
            continue
        mol = Chem.MolFromMolBlock(template["mol_block"], sanitize=True, removeHs=False)
        if mol is None or mol.GetNumAtoms() != len(keys):
            raise ValueError(f"Invalid stored product graph {token}")
        Chem.RemoveStereochemistry(mol)
        mol.RemoveAllConformers()
        conf = Chem.Conformer(len(keys))
        conf.Set3D(True)
        for i, key in enumerate(keys):
            conf.SetAtomPosition(i, Point3D(*prediction["atoms"][key]))
        mol.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        for i, (_, expected) in enumerate(template["atoms"]):
            if expected in {"R", "S"}:
                atom = mol.GetAtomWithIdx(i)
                observed = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
                centers.append({"atom": keys[i], "expected_product_CIP": expected,
                                "CIP_from_output_coordinates": observed, "pass": expected == observed})
        coverage[gid] = "scored"
    return centers, coverage


def chemical_scores(reference, prediction, geometry, clashes):
    inv = frozen.inventory(reference, prediction)
    try:
        from rdkit import rdBase
    except ImportError:
        return {"status": "RDKit_unavailable", "inventory": inv}
    with gzip.open(HERE / "kernels/reacted_products.json.gz", "rt") as f:
        products = json.load(f)
    centers, coverage = assign_products(reference, prediction,
                                        products["by_case"][reference["case_id"]], products["templates"])
    linkage_bonds, internal_bonds = [], []
    coords = prediction["atoms"]
    for (a, b, _), kind in zip(reference["bonds"], reference["bond_kinds"]):
        if kind == "protein_sequence":
            continue
        distance = math.dist(coords[a], coords[b]) if a in coords and b in coords else None
        row = {"atom_1": a, "atom_2": b, "distance_A": distance,
               "broad_anomaly": distance is None or distance < .9 or distance > 2.2, "kind": kind}
        if kind == "within_sugar":
            internal_bonds.append(row)
        else:
            row["within_fixed_routine_range"] = distance is not None and 1.1 <= distance <= 1.9
            linkage_bonds.append(row)
    report = {"stereochemistry": {"status": "scored", "centers": centers},
              "glycosidic_geometry": {"bonds": linkage_bonds},
              "within_monosaccharide_geometry": {"bond_distances": internal_bonds}}
    local = frozen.local_v1(reference, report, inv)["sites"]
    parent = frozen.convert(frozen.protein_stereochemistry(reference, prediction), reference["protein_sequence"])
    reference["protein_standard_parent_stereocenters"] = [{"atom": r["atom"]} for r in parent["centers"]]
    result = {}
    for gid, value in local.items():
        root = frozen._root_parent_stereo(reference, parent, reference["sites"][gid]["position"])
        clash = clashes.get("sites", {}).get(gid, {}).get("severe_gt_0_6_A")
        gates = {"legacy_local": value["status"],
                 "inventory_global": "fail" if inv["global_hard_fail"] else "pass",
                 "severe_clash": "unknown" if clash is None else "fail" if clash else "pass",
                 "root_parent_stereo": root["status"]}
        if coverage[gid] == "product_template_unavailable" and gates["legacy_local"] == "pass":
            gates["legacy_local"] = "unknown"
        extended = frozen._three_value_and(list(gates.values()))
        frame = geometry["sites"][gid]["frames"]["protein"]
        hard = inv["global_hard_fail"] or frame["status"] == "model_missing_or_nonfinite"
        joint = {str(t): frozen._joint_bit(extended, frame.get("rmsd_A"), t, hard) for t in (1., 2., 3., 5.)}
        result[gid] = {"legacy": value, "extended_components": gates, "extended_status": extended,
                       "joint": joint, "product_coverage": coverage[gid]}
    return {"status": "scored", "rdkit_version": rdBase.rdkitVersion, "inventory": inv, "sites": result,
            "parent_protein_stereo": parent}
