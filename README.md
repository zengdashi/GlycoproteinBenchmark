# GlycanFrameBench

GlycanFrameBench evaluates how accurately six all-atom structure predictors place covalently attached glycans on proteins. It includes the two benchmark sets, atom-level glycan connectivity, experimental reference coordinates, candidate scores, and code to evaluate a prediction or regenerate the paper's statistical panels and Supplementary Tables S1–S5.

## Benchmark sets

| Selection | Set A | Set B |
| --- | ---: | ---: |
| Curated chemical inputs | 348 | 59 |
| Distinct protein sequences in those inputs | 269 | 49 |
| Inputs with protein length <2,000 residues | 337 | 58 |
| Inputs completed by all six methods | 327 | 57 |
| Distinct sequences in the six-method comparison | 256 | 47 |
| Glycosylation sites in the six-method comparison | 884 | 171 |

An input is one protein chain with its specified glycosylation sites and glycan structures. Inputs can share a protein sequence and differ in glycan chemistry. Set A uses PDB entries first released from 1 June 2023 through 25 October 2026. Set B was filtered separately from entries released from 1 April 2026 through 25 October 2026; it falls entirely within A. The manuscript's accession-level audit found matches to SweetFold's public processed archive for 172 of the 327 common inputs in A and none of the 57 in B. A public archive match measures possible exposure.

AlphaFold 3, Boltz-2, Chai-1, Protenix, RFAA and SweetFold each produced five candidates on every common input, giving 9,810 candidate structures scored in set A. Analysis of set B reuses the same predictions. [Prediction settings](data/predictor_settings.tsv) list the frozen versions, weights, recycles and ranking fields. All runs used seed 1, full-protein MSAs and disabled protein templates. The RFAA adapter accounted for atomized attachment residues when reading its full-sequence MSA.

## Get started

Run in the repository root. Dataset and table validation use the Python standard library. Plotting and coordinate scoring use the packages listed in [requirements.txt](requirements.txt). A clean installation was tested with Python 3.12.2, NumPy 2.5.5, pandas 3.1.1, Matplotlib 3.11.2, Biopython 1.88, RDKit 2026.9.2 and SciPy 1.18.1. The original figure-rendering environment is documented in the panel index.

```bash
python -m pip install -r requirements.txt
python verify.py
python scoring/score.py --example
python plot.py 5D
```

The example reads a real experimental structure, a Chai-1 candidate and their explicit atom correspondence, scores them and checks the results against the released observations.

| Task | Command or reference |
| --- | --- |
| Render the 17 statistical panels | `python plot.py all` |
| Render one panel to a chosen directory | `python plot.py 2B --out-dir my_figures` |
| Rebuild S1–S5 in `build/tables/` | `python scripts/export_tables.py` |
| Check plot statistics and selected bootstrap intervals against the tables | `python plotting/check_data.py` |
| Read both datasets or reconstruct covalent bonds | [Case IDs, FASTA, graph specification and worked example](data/README.md) |
| Export an experimental chain, sugars and the matching atom IDs | [Experimental reference coordinates](data/reference/README.md) |
| Score another prediction | [Coordinate mapping, metrics and worked example](scoring/README.md) |
| Inspect native selection and per-site results | [Score files and units](scores/README.md) |
| Find the data behind a figure or table | [Panel index](plotting/README.md) and [supplementary table index](tables/README.md) |

The default figure destination is `build/figures/`. Panel names are `1B`, `2A–2C`, `3A–3D`, `4A–4C`, `5A–5D` and `S1A–S1B`. Each panel has a vector PDF and a high-resolution PNG. Figure 1B re-typesets the counts in the test-set diagram from the case list and glycan graphs. The remaining 16 panels use the submission's statistical plots. Protein and glycan overlay illustrations are prepared separately.

## What can be reproduced

- Recover the protein sequence, attachment residue and atom, monosaccharide identities, and both endpoints of every sugar–sugar bond for each input. The two sets share one collection of sequences and graphs, with explicit membership lists.
- Read the scores of all five predictions per method and common input, join the native selection to its site observations, and recompute the principal numerical comparisons. The repository includes a worked coordinate example and primary experimental references for all 336 inputs with any completed predictions, including all 327 inputs in the main comparison. The alternative-reference sensitivity estimates are in S1; alternate experimental copies are separate from the primary coordinate collection. Bulk predictions, MSAs and model weights are separate resources.
- Rebuild the exact five submitted TSV files from the archived estimates. `verify.py` checks their submission checksums and independently recomputes structural, chemistry, joint pass, length-stratified and candidate-selection estimates from the observation-level data. See the scoring and plotting guides for executable-level coverage and for quantities supplied as frozen results.
