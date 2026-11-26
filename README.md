# GlycoproteinBenchmark

Benchmark of glycan conformation, placement and covalent geometry in predictions from AlphaFold 3, Boltz-2, Chai-1, Protenix, RFAA and SweetFold.

## Test sets

| | PDB entries | Protein chains | Glycosylation sites |
| --- | ---: | ---: | ---: |
| A | 214 | 327 | 884 |
| B (subset of A) | 35 | 57 | 171 |

In the [benchmark dataset](data/glycoproteins.json), each input is identified by its PDB ID and protein chain (`label_asym_id`). Each record includes the protein sequence, glycosylation sites, sugar CCD codes and atom pairs for every protein–sugar and sugar–sugar bond. Sites have both sequence and author numbering. Sugar residue IDs use `label_asym_id:auth_seq_id`.

## Statistical figures

```bash
python -m pip install -r requirements.txt
python plot.py all
python plot.py 3D
```

The plotting commands write vector PDF and PNG panels to `build/figures/`. Available panels are **1B, 2A–C, 3A–D, 4A–C, 5A–D and S1A–B**. Figure 1B reads the test sets directly. [Figure source tables](plotting/data/) are named by panel.
