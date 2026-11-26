# Score a glycoprotein structure

Run commands from the repository root after installing `requirements.txt`. These are CPU calculations. The single-model command reports the same 14 coordinate-derived input-level fields as `scores/candidates.tsv.gz`, with optional detailed site, bond, ring and clash outputs. The five-model command reports the input-level dispersion used in Figure 5D.

## Reproduce a manuscript case

```bash
python scoring/score.py --example --out build/example_scores.json.gz
```

The example contains the deposited 9N3S mmCIF, the original Chai-1 model (candidate 001), and the complete atom mapping. Both glycosylation sites on the benchmark chain are scored. The candidate-level protein-aligned glycan RMSD is **2.013110 Å**. At the second site depicted in Figure 4D, this RMSD is **0.546103 Å**. The site-averaged noncovalent contact F1 is **0.876190** and one of the two sites passes extended chemistry. The example checks 42 values against the released candidate/site scores and the reported Figure 4D frames. It reads both original coordinate files, without using the saved scores during scoring. The PDB structure is identified in the manuscript. The other primary reference files in the release use anonymous identifiers.

## Score another benchmark prediction

Export the anonymous experimental structure and its fixed atom inventory, map the native prediction, then score. For the bundled example, this route reproduces the raw experimental mmCIF scores.

```bash
case_id=GLYCO_38627a5f10e0e98f
pred=scoring/example/Chai1_c001.cif
python data/reference/export.py "$case_id" --out build/reference
python scoring/map_prediction.py --method chai1 --case "$case_id" \
    --prediction "$pred" --reference-map build/reference/atom_map.tsv \
    --out build/matched_atoms.tsv
python scoring/score.py --case "$case_id" --reference build/reference/reference.cif \
    --prediction "$pred" --mapping build/matched_atoms.tsv --out build/scores.json.gz
```

`map_prediction.py` recognizes the frozen native output conventions for `af3`, `boltz2`, `chai1`, `protenix` and `sweetfold`. The compact name/port registry covers all 337 length-eligible inputs. It includes the per-atom ligand naming used in the Boltz-2 runs. An output produced with different chain names or a revised adapter needs the corresponding `--native-map` or an explicit mapping TSV. The mapper takes an uncompressed native mmCIF or PDB. The scorer also accepts gzip-compressed coordinates and maps.

For **RFAA**, add `--native-map path/to/runtime_atom_map.json` to the mapping command and use `--method rfaa`. Retain the runtime map written with that prediction. RFAA writes the atomized anchor and sugar atoms in native token order, with repeated atom names. The mapper checks their order, count and elements before restoring canonical identities. A planned token order is insufficient. The complete set of formal RFAA runtime maps belongs with the original native runs and is not distributed here.

For a different writer, provide a tab-delimited file with these columns.

| Column | Meaning |
|---|---|
| `canonical_atom` | Protein `protein/P1/<1-based sequence position>/<atom>`; sugar `glycan/<glycan_id>/<node_id>/<CCD atom>`, e.g. `glycan/G002/G002_R002/C1`. |
| `reference_atom_id` | Exact `_atom_site.id` in the selected mmCIF model or PDB atom serial. Use `.` for an unobserved reference atom. |
| `prediction_atom_id` | The corresponding serial in the native model, including predicted atoms that lack an experimental coordinate. Use `.` for missing predictions. |

Canonical keys and file-local serials must map one-to-one. For alternate locations the serial selects the conformer. Carry over the complete reference map and include every mapped heavy atom of the modeled entity, including any retained leaving atom. Rigid-body alignment is performed after mapping. The default checks the reference observation mask against the corresponding released primary reference. For a different experimental reference, `--alternative-reference` permits its own observation mask. Geometry routines retain missing predicted atoms in the fixed reference denominators.

## What the scores measure

| Endpoint | Coordinate calculation and outputs |
|---|---|
| Protein fold | Proper-rotation Cα RMSD and Cα-lDDT on the fixed observed protein mask. |
| Glycan conformation and placement | Separate whole-protein, site-neighborhood, root-ring and glycan-self fits, with equal weights for sites in the case-level RMSD. The JSON retains fitted atom sets, observed atom counts, transforms and individual site values. |
| Covalent chemistry | Expected protein–sugar and inter-sugar distances use the 1.1–1.9 Å interval. The within-sugar graph uses the broad 0.9–2.2 Å interval. Root angles, registered four-atom φ/ψ/ω errors, anchor side-chain χ angles and continuous ring puckering are returned separately. |
| Sugar and anchor stereochemistry | RDKit reassigns product stereocenters using the predicted coordinates. Eighty distinct reacted site graphs cover 1,045 sites in the 337 length-eligible cases. Parent-residue stereochemistry includes the frozen CCD isoleucine atom-name correction. |
| Contacts and clashes | Per-site nonbonded contact precision, recall and F1 at 4.0, 4.5 and 5.0 Å. Atom pairs separated by up to three covalent edges are excluded. Severe clashes have van der Waals overlap >0.6 Å; their denominator includes predicted heavy atoms at experimentally unresolved positions. |
| Extended chemistry and joint outcome | The extended gate combines the local chemical inventory, bond geometry and stereochemistry, global coordinate inventory, severe-clash status and anchor stereochemistry. For each threshold (1, 2, 3, 5 Å), joint success additionally requires sufficiently accurate protein-aligned placement. Unknown components preserve the corresponding missing-value status. |

The detailed JSON contains the atom identities and observed values that underlie each computation. The optional `--tmalign path/to/TMalign` activates the fixed-correspondence Cα TM-score adapter if an executable is installed. RDKit is required for the five chemical/joint fields and the full example comparison. Without RDKit, the custom scoring command returns the nine geometry/contact fields and an explicit chemistry status.

## Agreement of five predictions

Create a tab-separated manifest with `candidate_index`, `prediction` and `mapping` columns. Supply one row per index 0–4, all from the same method, case and seed. The two path columns are relative to the manifest. Generate each mapping as above.

```bash
python scoring/score_ensemble.py --case "$case_id" --manifest five_candidates.tsv \
    --out build/dispersion.json
```

The ensemble calculation uses the full input-defined protein Cα positions and reacted glycan atoms. It computes ten prediction–prediction comparisons and reports their **median** site-averaged, protein-aligned glycan RMSD. The experimental coordinates are unused. This is the `dispersion_protein` measure behind S5.

## Other manuscript endpoints

The readers can regenerate the figures and supplementary tables directly from the released scores. A single coordinate pair supplies neither native model confidence nor between-input summaries. The candidate and input tables carry the original score/rank values, and the reference-sensitivity table contains the aggregate over alternative experimental references. The coordinate CLI scores one specified reference at a time. The S2 ring-class agreement uses the independently run Privateer classifier. The geometric Cremer–Pople quantities returned here describe ring shapes; Privateer categories are supplied in the score tables. See [score sources](../scores/README.md) and [supplementary-table definitions](../tables/README.md).

## Implementation checks

The geometry and nonbonded kernels in `kernels/` are the ones used for the benchmark. The public wrappers rebuild their expected-atom and covalent-edge dictionaries from the release annotation and the standard parent-residue templates. Compact reagent graphs contain the fixed atom order and expected product CIP labels. The classes, keys, cutoff conventions and averaging rules used by the case-level scores remain those of the analysis.

During packaging, five archived formal predictions across three reference inputs reproduced all 70 candidate-level values in the score table (AlphaFold 3 twice, Boltz-2, Chai-1 and SweetFold). The largest absolute difference was below 6 × 10⁻¹⁵. The 42-value bundled example can be rerun without network access. See [verification details](VERIFICATION.md).
