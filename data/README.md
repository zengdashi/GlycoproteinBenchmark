# Benchmark Dataset Schema (`glycoproteins.json`)

[`glycoproteins.json`](glycoproteins.json) defines the 327 nonredundant single-chain glycoprotein inputs (884 glycosylation sites; 1,831 resolved monosaccharides; 947 inter-sugar glycosidic bonds) completed by all six evaluated predictors.

---

## Record Hierarchy and Field Definitions

The JSON file is a top-level list of **327 chain records**.

### 1. Protein Chain Record

| Field | Type | Description |
| --- | --- | --- |
| `pdb_id` | `string` | Four-character Worldwide Protein Data Bank (wwPDB) entry identifier (uppercase). |
| `protein_chain` | `string` | mmCIF strand identifier (`label_asym_id`) of the glycosylated protein chain. |
| `protein_auth_chain` | `string` | Author-assigned chain identifier (`auth_asym_id`) in the deposited mmCIF structure. |
| `test_sets` | `list[string]` | Cohort membership: `["A"]` (270 inputs) or `["A", "B"]` (57 later-released nested inputs). |
| `sequence` | `string` | Full canonical one-letter amino-acid sequence (`< 2,000` residues; non-carbohydrate PTMs `HYP` and `OCS` standardized to parent `P` and `C`). |
| `glycosylation_sites` | `list[object]` | All experimentally resolved covalent glycosylation sites on this protein chain (1 to 17 sites per chain). |

### 2. Glycosylation Site Object (`glycosylation_sites[]`)

| Field | Type | Description |
| --- | --- | --- |
| `label_seq_id` | `int` | 1-based residue index in `sequence` (`_atom_site.label_seq_id`). |
| `auth_seq_id` | `string` | Author residue numbering (including any insertion code) in the deposited mmCIF. |
| `residue` | `string` | Three-letter parent residue code of the modified amino acid (`ASN`, `SER`, or `THR`). |
| `link_type` | `string` | Covalent attachment class: `"N"` (*N*-linked, 844 sites in A) or `"O"` (*O*-linked, 40 sites in A). |
| `protein_atom` | `string` | Bonding heavy-atom name on the modified amino-acid side chain (e.g., `ND2`, `OG`, `OG1`). |
| `root_sugar` | `string` | Identifier (`id`) of the root monosaccharide covalently attached to `protein_atom`. |
| `root_atom` | `string` | Bonding heavy-atom name on the root monosaccharide (e.g., `C1`). |
| `sugars` | `list[object]` | Experimentally resolved monosaccharides belonging to the glycan tree at this site. |
| `bonds` | `list[object]` | Directed list of inter-sugar glycosidic covalent bonds (`[]` for single-monosaccharide sites). |

### 3. Monosaccharide Object (`sugars[]`)

| Field | Type | Description |
| --- | --- | --- |
| `id` | `string` | Unique monosaccharide identifier within the PDB entry formatted as `label_asym_id:auth_seq_id`. |
| `ccd` | `string` | PDB Chemical Component Dictionary (CCD) monosaccharide code (16 distinct codes across the dataset: `A1H1B`, `A1H1F`, `A2G`, `BGC`, `BMA`, `FUC`, `FUL`, `GAL`, `GCU`, `GLA`, `GLC`, `MAN`, `NAG`, `RAM`, `XYP`, `XYS`). |
| `auth_chain` | `string` | Author chain identifier (`auth_asym_id`) carrying the carbohydrate residue in mmCIF. |
| `branch_num` | `int` *(optional)* | Monosaccharide index (`_pdbx_branch_scheme.num`) when modeled as an mmCIF branched carbohydrate entity; omitted for non-polymer single-sugar HETATM records. |

### 4. Inter-Sugar Bond Object (`bonds[]`)

| Field | Type | Description |
| --- | --- | --- |
| `sugar_1` | `string` | Monosaccharide `id` of the first bonded sugar residue. |
| `atom_1` | `string` | Heavy-atom name on `sugar_1` participating in the glycosidic bond (e.g., `O4`, `O6`, `O3`, `O2`). |
| `sugar_2` | `string` | Monosaccharide `id` of the second bonded sugar residue. |
| `atom_2` | `string` | Heavy-atom name on `sugar_2` participating in the glycosidic bond (typically `C1`). |

---

## Minimal Python Usage Example

```python
import json
from pathlib import Path

records = json.loads(Path("data/glycoproteins.json").read_text(encoding="utf-8"))

set_b = [r for r in records if "B" in r["test_sets"]]
print(f"Set A chains: {len(records)}, Set B chains: {len(set_b)}")

for site in records[0]["glycosylation_sites"]:
    glycan_len = len(site["sugars"])
    print(
        f"{records[0]['pdb_id']} chain {records[0]['protein_chain']} "
        f"{site['residue']}{site['label_seq_id']} ({site['link_type']}-linked): "
        f"{glycan_len} sugar(s), root={site['root_sugar']}"
    )
```
