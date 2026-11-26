#!/usr/bin/env python3
"""Check released reference coverage, observed masks and frozen linkage angles."""
import hashlib
import math
from collections import defaultdict

from export import AA3, HERE, case_graphs, linkage_keys, read_rows


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dihedral(points):
    a, b, c, d = points
    axis = subtract(c, b)
    axis = tuple(v / math.sqrt(sum(x * x for x in axis)) for v in axis)

    def projected(vector):
        dot = sum(x * y for x, y in zip(vector, axis))
        return tuple(x - dot * y for x, y in zip(vector, axis))

    v, w = projected(subtract(a, b)), projected(subtract(d, c))
    cross = (axis[1]*v[2] - axis[2]*v[1], axis[2]*v[0] - axis[0]*v[2],
             axis[0]*v[1] - axis[1]*v[0])
    return math.degrees(math.atan2(sum(x*y for x, y in zip(cross, w)),
                                   sum(x*y for x, y in zip(v, w))))


def check_reference(index, case, checks):
    cid = case["case_id"]
    path = HERE / "records" / (cid + ".tsv.gz")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == index["record_sha256"], cid
    records = read_rows(path)
    sequence = case["protein"]["sequence"]
    sugars = {n["node_id"]: (g["glycan_id"], n) for g in case["glycans"] for n in g["nodes"]}
    expected, observed = set(), {}
    counts = defaultdict(int)
    for row in records:
        key, ccd = row["canonical_atom"], row["ccd"]
        kind, group, unit, atom = key.split("/")
        assert key not in expected and row["element"] not in {"H", "D"}, (cid, key)
        expected.add(key)
        if kind == "protein":
            assert group == "P1" and 0 < int(unit) <= len(sequence)
            assert AA3[sequence[int(unit) - 1]] == ccd
        else:
            glycan_id, node = sugars[unit]
            assert kind == "glycan" and glycan_id == group and node["ccd_id"] == ccd
            assert atom in node["expected_reacted_heavy_atom_names"]
            assert (row["x"] != ".") == (atom in node["observed_atom_names"])
        if row["x"] == ".":
            assert all(row[k] == "." for k in ("y", "z", "occupancy", "b_iso"))
        else:
            xyz = tuple(float(row[k]) for k in ("x", "y", "z"))
            assert all(math.isfinite(x) for x in xyz) and float(row["occupancy"]) > 0
            observed[key] = xyz
            counts[kind] += 1
            if kind == "protein" and atom == "CA":
                counts["CA"] += 1
    for node_id, (gid, node) in sugars.items():
        assert all(f"glycan/{gid}/{node_id}/{n}" in expected for n in node["expected_reacted_heavy_atom_names"])
    for a, b, order in linkage_keys(case):
        assert a in expected and b in expected and order == "SING"
    mask = hashlib.sha256(("\n".join(sorted(observed)) + "\n").encode()).hexdigest()
    assert mask == index["observed_atom_mask_sha256"], cid
    assert int(index["expected_atoms"]) == len(expected)
    for name, field in (("protein", "observed_protein_atoms"), ("CA", "observed_protein_CA"),
                        ("glycan", "observed_glycan_atoms")):
        assert counts[name] == int(index[field])
    assert len(checks) == 2 * int(index["frozen_angle_pairs"])
    for t in checks:
        angle = dihedral([observed[t[f"atom_{n}"]] for n in range(1, 5)])
        error = abs((angle - float(t["archived_reference_degrees"]) + 180.) % 360. - 180.)
        assert error < 1e-8, (cid, t["edge_id"], t["torsion"], error)
    return len(observed)


def main():
    catalog = {r["case_id"]: r for r in read_rows(HERE.parent / "cases.tsv")}
    index = {r["case_id"]: r for r in read_rows(HERE / "index.tsv")}
    assert set(index) == {cid for cid, r in catalog.items() if int(r["n_methods_complete"]) > 0}
    assert set(index) == {p.name.split(".")[0] for p in (HERE / "records").glob("*.tsv.gz")}
    checks = defaultdict(list)
    for r in read_rows(HERE / "torsion_checks.tsv.gz"):
        checks[r["case_id"]].append(r)
    assert set(checks) == set(index)
    graphs = case_graphs()
    observed = 0
    for cid, row in index.items():
        assert all(row[k] == catalog[cid][k] for k in ("in_B", "common_six"))
        observed += check_reference(row, graphs[cid], checks[cid])
    print(f"Validated {len(index)} primary references, {observed:,} observed atoms, "
          f"and {sum(map(len, checks.values())):,} archived experimental torsions.")


if __name__ == "__main__":
    main()
