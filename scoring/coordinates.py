"""Resolve explicit file-local atom IDs. Coordinates never determine an atom match."""
import csv
import gzip
import math
from pathlib import Path


def coordinate_rows(path):
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    suffix = Path(path.stem).suffix if path.suffix.lower() == ".gz" else path.suffix
    with opener(path, "rt") as handle:
        if suffix.lower() in {".cif", ".mmcif"}:
            from Bio.PDB.MMCIF2Dict import MMCIF2Dict
            d = MMCIF2Dict(handle)
            ids = d["_atom_site.id"]
            models = d.get("_atom_site.pdbx_PDB_model_num", ["1"] * len(ids))
            for i, serial in enumerate(ids):
                if models[i] == models[0]:
                    yield serial, [float(d["_atom_site.Cartn_" + v][i]) for v in "xyz"], d["_atom_site.type_symbol"][i]
        elif suffix.lower() in {".pdb", ".ent"}:
            for line in handle:
                if line.startswith("ENDMDL"):
                    break
                if line.startswith(("ATOM  ", "HETATM")):
                    yield line[6:11].strip(), [float(line[a:b]) for a, b in ((30,38),(38,46),(46,54))], line[76:78].strip()
        else:
            raise ValueError("Coordinate input requires PDB or mmCIF (optionally .gz).")


def load_atoms(path):
    atoms = {}
    for serial, xyz, element in coordinate_rows(path):
        if serial in atoms:
            raise ValueError(f"Duplicate atom serial {serial} in {Path(path).name}")
        atoms[serial] = (xyz, element.upper())
    if not atoms:
        raise ValueError(f"No atoms in {Path(path).name}")
    return atoms


def read_mapping(path, reference, prediction, topology):
    sources = [load_atoms(reference), load_atoms(prediction)]
    views = [{"atoms": {}, "atom_meta": {}}, {"atoms": {}, "atom_meta": {}}]
    seen_atoms, used = set(), [set(), set()]
    meta = topology["expected_atom_meta"]
    sequence_length = len(topology["protein_sequence"])
    optional = f"protein/P1/{sequence_length}/OXT"
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        required = {"canonical_atom", "reference_atom_id", "prediction_atom_id"}
        if not required <= set(rows.fieldnames or []):
            raise ValueError(f"Mapping needs columns {sorted(required)}")
        for row in rows:
            key = row["canonical_atom"]
            parts = key.split("/")
            if key in seen_atoms or len(parts) != 4:
                raise ValueError(f"Invalid or duplicate canonical atom {key}")
            if key not in meta and key != optional:
                if parts[0] != "glycan" or parts[2] not in topology["nodes"] or topology["nodes"][parts[2]]["glycan_id"] != parts[1]:
                    raise ValueError(f"Unregistered canonical atom {key}")
            seen_atoms.add(key)
            for i, field in enumerate(("reference_atom_id", "prediction_atom_id")):
                serial = row[field].strip()
                if serial in {"", ".", "?"}:
                    continue
                if serial in used[i] or serial not in sources[i]:
                    raise ValueError(f"Unknown or reused {field} {serial} for {key}")
                used[i].add(serial)
                xyz, el = sources[i][serial]
                if not all(math.isfinite(x) for x in xyz):
                    raise ValueError(f"Nonfinite mapped coordinates for {key}")
                if i == 0 and key not in meta:
                    raise ValueError(f"Reference key {key} outside the fixed scoring inventory")
                if i == 0 and el and el != meta[key]["element"].upper():
                    raise ValueError(f"Reference elemental identity disagrees at {key}")
                views[i]["atoms"][key] = xyz
                views[i]["atom_meta"][key] = {"element": el}
    if not all(v["atoms"] for v in views):
        raise ValueError("Each coordinate file must contain at least one mapped atom.")
    return views
