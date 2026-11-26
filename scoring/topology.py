"""Build the benchmark's expected heavy-atom scoring graph from a release case."""
import json
from pathlib import Path

from kernels.geometry_metrics import _site_data, _node_data
from kernels.torsions import registered_torsions

HERE = Path(__file__).resolve().parent


def build_topology(case, components):
    templates = json.loads((HERE / "kernels/topology_templates.json").read_text())
    sequence = case["protein"]["sequence"]
    elements, bonds, kinds = {}, [], []

    def bond(a, b, order, kind):
        bonds.append([a, b, order])
        kinds.append(kind)

    for pos, residue in enumerate(sequence, 1):
        template = templates["amino_acids"][residue]
        prefix = f"protein/P1/{pos}/"
        for name, element in template["atoms"].items():
            elements[prefix + name] = {"element": element}
        for a, b, order in template["bonds"]:
            bond(prefix + a, prefix + b, order, "protein_sequence")
        if pos > 1:
            bond(f"protein/P1/{pos-1}/C", prefix + "N", "SING", "protein_sequence")

    nodes = _node_data({}, case)
    for glycan in case["glycans"]:
        gid = glycan["glycan_id"]
        for node in glycan["nodes"]:
            nid, ccd = node["node_id"], node["ccd_id"]
            prefix = f"glycan/{gid}/{nid}/"
            expected = set(node["expected_reacted_heavy_atom_names"])
            comp = components[ccd]
            for atom in comp["atoms"]:
                if atom["atom_name"] in expected:
                    elements[prefix + atom["atom_name"]] = {"element": atom["element"]}
            for edge in comp["bonds"]:
                a, b = edge["atom_1"], edge["atom_2"]
                if a in expected and b in expected:
                    bond(prefix + a, prefix + b, edge["bond_order"], "within_sugar")
            rings = templates["sugar_rings"][ccd]
            nodes[nid]["ring_atom_keys"] = [prefix + a for a in rings[0]] if rings else []
            nodes[nid]["ring_registry_status"] = (
                "single_5_or_6_member_oxygen_ring" if len(rings) == 1 else
                "multiple_rings" if rings else "unsupported_ring")
        for edge in glycan["bonds"]:
            a = f"glycan/{gid}/{edge['node_1']}/{edge['atom_1']}"
            b = f"glycan/{gid}/{edge['node_2']}/{edge['atom_2']}"
            bond(a, b, edge["bond_order"], "glycosidic")
        at = glycan["attachment"]
        a = f"protein/P1/{at['protein_residue_index']}/{at['protein_atom']}"
        b = f"glycan/{gid}/{at['glycan_node_id']}/{at['glycan_atom']}"
        bond(a, b, at["bond_order"], "attachment")
    if case.get("additional_covalent_bonds"):
        raise ValueError("This release registers extra covalent edges only when explicitly supported.")
    return {"case_id": case["case_id"], "protein_sequence": sequence,
            "expected_atom_meta": elements, "bonds": bonds, "bond_kinds": kinds,
            "nodes": nodes, "sites": _site_data({}, case),
            "torsion_registry": registered_torsions(bonds, kinds, elements, nodes)}
