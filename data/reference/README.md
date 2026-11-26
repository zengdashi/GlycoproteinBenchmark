# Experimental reference coordinates

Each entry contains the primary experimental protein chain and its observed covalently attached sugars, labeled with the same residue positions, glycan IDs and sugar-node IDs as the model input. Use these coordinates for evaluation after prediction.

| Scope | Set A | Nested Set B |
| --- | ---: | ---: |
| Input graphs in the release | 348 | 59 |
| Protein length <2,000 residues | 337 | 58 |
| At least one method returned five candidates and a scored primary reference is available here | 336 | 57 |
| All six methods returned five candidates | 327 | 57 |

The nine additional scored Set A inputs outside the common-six comparison have results from five methods. The input graphs for all 348 cases remain available in the parent directory. Each covered case has one primary reference. Alternate experimental copies, crystal neighbors and assembly mates fall outside this coordinate package.

## Export a structure for scoring

From the repository root, run

```sh
python data/reference/export.py GLYCO_03ac15f95b2abb85 --out build/example
```

- `reference.cif` contains the selected observed atoms in the experimental coordinate frame, with occupancy and experimental B factors where reported. The protein is chain P1. Each sugar node has its own chain with its anonymous node ID and residue number 1. Standard mmCIF `_struct_conn` records identify protein–sugar and sugar–sugar links with both endpoint coordinates present. The complete input connectivity, including links with an unobserved endpoint, is defined in `data/glycan_connections.json`.
- `atom_map.tsv` gives the canonical atom key and the corresponding `_atom_site.id` in the exported reference. The column `prediction_atom_id` is initialized to `.`. Fill it using atom IDs from the specific prediction and that method’s input-to-output correspondence. Reference serials and prediction serials have independent numbering. Keep rows with missing reference or prediction coordinates. Examples of canonical names are `protein/P1/118/ND2` and `glycan/G001/G001_R001/C1`.

The exporter and the reference validator use the Python standard library. The export directory may be outside the repository. Each invocation writes the two named files at that location.

## Records and atom selection

`records/<case_id>.tsv.gz` holds one row per expected heavy atom of the standardized input. Columns are the canonical key, standardized CCD residue name, element, Cartesian coordinates in ångströms, occupancy and isotropic B factor in Å². A period denotes a missing value. Missing atoms have no coordinates. All 336 records together contain 1,477,664 expected atom positions and 1,128,393 observed positions, including 1,101,294 protein atoms and 27,099 glycan atoms.

The primary reference was reconstructed from the same deposited mmCIF, source chain and sugar-node correspondence used in frozen benchmark scoring. For each residue, the selection uses the first experimental model, positive occupancies and a single compatible alternate conformation. The scoring rule prioritizes heavy-atom coverage, then summed occupancy, with a fixed tie-break. Protein coordinates map to matching heavy-atom names in the standardized parent residue. The terminal OXT and removed modification-specific atoms are excluded. Sugars use the frozen experimental mask for the reacted component. Coordinates for missing atoms are never imputed.

`index.tsv` lists set flags, per-case atom counts and checksums. The observed mask digest is SHA-256 of lexicographically sorted canonical atom keys, one key per line with a final newline. The record digest is calculated on the distributed gzip bytes. The column `alternate_view_mask_available` marks the 105 cases whose primary observed-mask digest was also recorded in the separate analysis of alternate experimental copies.

## Validate the release

```sh
python data/reference/validate.py
```

This checks all record digests, observed-atom masks, glycan atom names, link endpoints and the 4,210 experimentally measured φ and ψ values for 2,105 linkages in `torsion_checks.tsv.gz`. Each angle has four ordered canonical atom IDs and the corresponding archived reference value in degrees. Every scored case contributes at least one angle pair. The acceptance tolerance is 10⁻⁸ degrees after angular wraparound.

During packaging, deposited coordinate-file and canonical-source fingerprints were verified against the private primary-reference assignment. All 336 reconstructed protein masks, alternate-conformer selections and counts of observed glycan atoms agreed with the archived scoring tables. All 27,099 sugar coordinates were also matched by canonical identity against the independently stored source-chain annotation. The 105 available full observed-mask digests agreed with the alternate-copy analysis. As an independent round-trip test, BioPython parsed every exported mmCIF. All 1,128,393 coordinates matched the released records through the emitted atom maps, and all 2,105 observed intercomponent bonds resolved to the two intended canonical atoms. Original PDB and author-chain metadata are held separately.
