#!/usr/bin/env python3
"""Export one anonymous experimental reference and its scoring-map template."""
import argparse
import csv
import gzip
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
AA3 = dict(zip("ARNDCQEGHILKMFPSTWYV", (
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"
).split()))


def read_rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def case_graphs():
    annotation = json.loads((HERE.parent / "glycan_connections.json").read_text(encoding="utf-8"))
    return {item["case_id"]: item for item in annotation["cases"]}


def linkage_keys(case):
    for glycan in case["glycans"]:
        gid = glycan["glycan_id"]
        root = glycan["attachment"]
        protein = f"protein/P1/{root['protein_residue_index']}/{root['protein_atom']}"
        sugar = f"glycan/{gid}/{root['glycan_node_id']}/{root['glycan_atom']}"
        yield protein, sugar, root["bond_order"]
        for edge in glycan["bonds"]:
            left = f"glycan/{gid}/{edge['node_1']}/{edge['atom_1']}"
            right = f"glycan/{gid}/{edge['node_2']}/{edge['atom_2']}"
            yield left, right, edge["bond_order"]


def cif_loop(stream, category, fields, payload):
    stream.write("#\nloop_\n")
    stream.writelines(f"_{category}.{f}\n" for f in fields)
    for record in payload:
        # Released IDs and CCD/atom names are single, unquoted CIF tokens.
        values = [str(v) for v in record]
        if any(not v or any(c.isspace() for c in v) or v[0] in "'\"_;#$[]"
               for v in values):
            raise ValueError(f"Unsupported CIF token in {category}: {values}")
        stream.write(" ".join(values) + "\n")


def write_reference(case, records, coordinate_stream, mapping_stream):
    case_id = case["case_id"]
    sequence = case["protein"]["sequence"]
    nodes = {n["node_id"]: n["ccd_id"] for g in case["glycans"] for n in g["nodes"]}
    entities = {ccd: str(i) for i, ccd in enumerate(sorted(set(nodes.values())), 2)}
    chain_entities = {"P1": "1", **{node: entities[ccd] for node, ccd in nodes.items()}}
    coordinate_stream.write(f"data_{case_id}\n_entry.id {case_id}\n")
    cif_loop(coordinate_stream, "entity", ["id", "type"], [
        ("1", "polymer"), *((num, "non-polymer") for num in entities.values())
    ])
    coordinate_stream.write("#\n_entity_poly.entity_id 1\n_entity_poly.type 'polypeptide(L)'\n")
    cif_loop(coordinate_stream, "entity_poly_seq", ["entity_id", "num", "mon_id", "hetero"],
             (("1", p, AA3[letter], "n") for p, letter in enumerate(sequence, 1)))
    cif_loop(coordinate_stream, "struct_asym", ["id", "entity_id"], chain_entities.items())

    atom_fields = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
                   "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
                   "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
                   "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id",
                   "auth_atom_id", "pdbx_PDB_model_num"]
    writer = csv.writer(mapping_stream, delimiter="\t", lineterminator="\n")
    writer.writerow(("canonical_atom", "reference_atom_id", "prediction_atom_id"))
    present = {}

    def atom_rows():
        for row in records:
            key = row["canonical_atom"]
            if row["x"] == ".":
                writer.writerow((key, ".", "."))
                continue
            kind, group, unit, atom = key.split("/")
            protein = kind == "protein"
            chain = group if protein else unit
            pos = unit if protein else "."
            auth_pos = unit if protein else "1"
            ccd = row["ccd"]
            atom_id = len(present) + 1
            present[key] = (chain, ccd, pos, atom, auth_pos)
            writer.writerow((key, atom_id, "."))
            yield ("ATOM" if protein else "HETATM", atom_id, row["element"], atom, ".",
                   ccd, chain, chain_entities[chain], pos, "?", row["x"], row["y"],
                   row["z"], row["occupancy"], row["b_iso"], auth_pos, ccd, chain, atom, 1)

    cif_loop(coordinate_stream, "atom_site", atom_fields, atom_rows())
    endpoints = ["label_asym_id", "label_comp_id", "label_seq_id", "label_atom_id", "auth_seq_id"]
    conn_fields = ["id", "conn_type_id", *("ptnr1_"+s for s in endpoints),
                   *("ptnr2_"+s for s in endpoints), "pdbx_value_order"]
    conn_rows = []
    for i, (a, b, order) in enumerate(linkage_keys(case), 1):
        if a in present and b in present:
            conn_rows.append((f"link{i}", "covale", *present[a], *present[b], order))
    if conn_rows:
        cif_loop(coordinate_stream, "struct_conn", conn_fields, conn_rows)
    coordinate_stream.write("#\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_id", help="Standardized GLYCO identifier in index.tsv")
    parser.add_argument("--out", required=True, type=Path, help="Directory for reference.cif and atom_map.tsv")
    args = parser.parse_args()
    known = {r["case_id"] for r in read_rows(HERE / "index.tsv")}
    if args.case_id not in known:
        parser.error(f"No scored primary reference for {args.case_id}; see index.tsv")
    data = read_rows(HERE / "records" / (args.case_id + ".tsv.gz"))
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "reference.cif").open("w", encoding="utf-8") as structure:
        with (args.out / "atom_map.tsv").open("w", encoding="utf-8", newline="") as mapping:
            write_reference(case_graphs()[args.case_id], data, structure, mapping)


if __name__ == "__main__":
    main()
