# GlycoproteinBenchmark

Benchmark datasets, figure source data, and plotting code for **all-atom structure prediction of covalently glycosylated proteins**, accompanying the manuscript:

> **Benchmarking all-atom structure prediction of covalently glycosylated proteins**  
> Fuming Zeng\*, Yu Wang\*, Dongbo Bu & Shiwei Sun  
> \*Equal contribution · Correspondence: `dwsun@ict.ac.cn`

---

## Overview

`GlycoproteinBenchmark` evaluates six all-atom biomolecular structure predictors—**AlphaFold 3**, **Boltz-2**, **Chai-1**, **Protenix**, **RoseTTAFold All-Atom (RFAA)**, and **SweetFold**—on sequence-nonredundant single-chain glycoproteins from the Protein Data Bank (PDB) when glycosylation sites, monosaccharide identities, and atom-level covalent connectivity are provided as input.

The benchmark separates prediction performance across four hierarchical dimensions:
1. **Multi-frame spatial accuracy** — comparing internal glycan conformation (glycan self-alignment and root-ring alignment) against surface placement relative to the host protein (local and whole-protein Cα alignments) and noncovalent residue–monosaccharide contacts.
2. **Covalent linkage and local stereochemistry** — auditing protein–root-sugar and inter-sugar bond lengths, anomeric configurations, Privateer IUPAC ring conformations, attachment valence angles, and glycosidic torsions ($\phi$, $\psi$, and $(1{\to}6)$ $\omega$).
3. **Joint chemical–spatial validity and failure modes** — quantifying site-level pass rates that simultaneously require extended chemical plausibility and $\le 2\text{ \AA}$ whole-protein placement, alongside defect breakdowns (steric clashes, chirality inversions, and bond-length violations).
4. **Glycan complexity, candidate ranking, and spatial consistency** — controlling for chain length across linear and branched topologies, contrasting native score discrimination across inputs versus within five candidates per input, and evaluating inter-candidate glycan spatial agreement.

---

## Benchmark Test Sets

All 327 single-chain glycoprotein inputs completed by all six predictors (5 candidates per input; seed 1) are provided in [`data/glycoproteins.json`](data/glycoproteins.json).

| Metric | Test set A | Test set B (nested subset of A) |
| --- | ---: | ---: |
| **Release window** | 2023-06-01 – 2026-10-25 | 2026-04-01 – 2026-10-25 |
| **PDB entries** | 214 | 35 |
| **Protein chains (inputs)** | 327 | 57 |
| **Distinct protein sequences** | 256 | 47 |
| **Glycosylation sites** | **884** | **171** |
| ↳ *N*-linked / *O*-linked sites | 844 / 40 | 154 / 17 |
| ↳ Monosaccharide sites | 484 | 78 |
| ↳ Linear oligosaccharide sites | 299 | 63 |
| ↳ Branched oligosaccharide sites | 101 | 30 |
| **Total resolved monosaccharides** | 1,831 | 476 |
| **Inter-sugar glycosidic bonds** | 947 | 305 |

- **Sequence filtering:** Candidate entries were filtered with MMseqs2 ($<30\%$ sequence identity at $\ge 80\%$ bidirectional coverage) against pre-cutoff PDB protein sequences and candidate RFAA training sequences, then de-duplicated across clusters and restricted to single chains $<2{,}000$ residues.
- **Temporal holdout (Test set B):** All 57 inputs in subset B were released after SweetFold's reported data collection cutoff (2026-03-08) and have zero overlap with SweetFold's public processed-data manifests.

---

## Repository Structure

```text
GlycoproteinBenchmark/
├── README.md                  # Repository overview, quick start, and panel index
├── requirements.txt           # Python dependencies for reproducing figures
├── plot.py                    # Command-line entry point for rendering statistical panels
├── data/
│   ├── README.md              # Detailed JSON schema, numbering conventions, and loader example
│   └── glycoproteins.json     # Curated Test sets A & B (327 PDB chains, 884 sites, atom-level bonds)
└── plotting/
    ├── figure1.py             # Fig. 1B test-set composition table panel
    ├── figure2.py             # Fig. 2A–C backbone accuracy, 4-frame glycan RMSD, and lDDT correlation
    ├── figure3.py             # Fig. 3A–D & Fig. S1A–B bond lengths, stereochemistry, angles, and torsions
    ├── figure4.py             # Fig. 4A–C root-bond outliers, joint outcomes, and failure defects
    ├── figure5.py             # Fig. 5A–D chain length, score quartiles, 5-candidate choice, and spread
    ├── panels.py              # Panel registry and data-loading dispatcher
    ├── style.py               # Shared publication typography and method color palette
    └── data/
        ├── README.md          # Column dictionary and statistical weighting guide for all TSVs
        └── *.tsv              # 14 source data tables underlying Figs. 2–5 and Fig. S1
```

---

## Quick Start: Reproducing Statistical Figures

### 1. Install dependencies

Requires Python $\ge 3.9$ with `numpy`, `pandas`, and `matplotlib`:

```bash
python -m pip install -r requirements.txt
```

### 2. Render all statistical panels or individual panels

```bash
# Render all 17 statistical panels (1B, 2A–C, 3A–D, 4A–C, 5A–D, S1A–B)
python plot.py all

# Render a single panel (case-insensitive)
python plot.py 3D

# Customize output directory and PNG resolution
python plot.py 2B --out-dir build/figures --dpi 600
```

Each panel is exported as both an editable vector **PDF** (`build/figures/Fig<PANEL>.pdf`, TrueType `fonttype=42`) and a high-resolution **PNG** (`build/figures/Fig<PANEL>.png`, default 400 dpi for `1B` and 500 dpi for all other panels).

### 3. Panel-to-Source Mapping

| Panel | Description | Source File(s) |
| --- | --- | --- |
| **`1B`** | Benchmark test-set composition (inputs, sites, and topologies) | [`data/glycoproteins.json`](data/glycoproteins.json) |
| **`2A`** | Whole-chain protein Cα RMSD and Cα-lDDT on Sets A and B | [`plotting/data/fig2_structure.tsv`](plotting/data/fig2_structure.tsv) + `glycoproteins.json` |
| **`2B`** | Glycan heavy-atom RMSD across four rigid-body alignment frames | [`plotting/data/fig2_structure.tsv`](plotting/data/fig2_structure.tsv) + `glycoproteins.json` |
| **`2C`** | Protein Cα-lDDT versus local-protein-aligned glycan RMSD | [`plotting/data/fig2_structure.tsv`](plotting/data/fig2_structure.tsv) + `glycoproteins.json` |
| **`3A`** | Protein–sugar and inter-sugar covalent bond-length pass rates | [`plotting/data/fig3A_bond_lengths.tsv`](plotting/data/fig3A_bond_lengths.tsv) |
| **`3B`** | Anomeric chirality, Privateer ring class, and extended site chemistry | [`plotting/data/fig3B_stereochemistry.tsv`](plotting/data/fig3B_stereochemistry.tsv) |
| **`3C`** | *N*-linked and *O*-linked attachment valence angle distributions | [`plotting/data/fig3C_attachment_angles.tsv`](plotting/data/fig3C_attachment_angles.tsv) |
| **`3D`** | Attachment-angle error and Asn side-chain RMSD vs. disaccharide placement | [`plotting/data/fig3D_disaccharides.tsv`](plotting/data/fig3D_disaccharides.tsv) + `glycoproteins.json` |
| **`4A`** | Sites with out-of-range protein–root-sugar bond lengths (Sets A & B) | [`plotting/data/fig4A_root_bonds.tsv`](plotting/data/fig4A_root_bonds.tsv) |
| **`4B`** | Mutually exclusive site outcomes (chemistry × 2 Å placement) | [`plotting/data/fig4B_joint_outcomes.tsv`](plotting/data/fig4B_joint_outcomes.tsv) |
| **`4C`** | Prevalence of specific chemical defects among failed sites | [`plotting/data/fig4C_defects.tsv`](plotting/data/fig4C_defects.tsv) |
| **`5A`** | Self-aligned glycan RMSD stratified by observed monosaccharide count | [`plotting/data/fig5A_glycan_length.tsv`](plotting/data/fig5A_glycan_length.tsv) |
| **`5B`** | Whole-protein glycan placement RMSD across native score quartiles | [`plotting/data/fig5B_score_quartiles.tsv`](plotting/data/fig5B_score_quartiles.tsv) |
| **`5C`** | Five-candidate selection: native top-1 vs. random vs. reference-best | [`plotting/data/fig5C_candidate_choice.tsv`](plotting/data/fig5C_candidate_choice.tsv) |
| **`5D`** | Inter-candidate glycan spread vs. top-ranked error & partial Spearman $\rho$ | [`plotting/data/fig5D_AlphaFold3.tsv`](plotting/data/fig5D_AlphaFold3.tsv), [`fig5D_correlations.tsv`](plotting/data/fig5D_correlations.tsv) |
| **`S1A`**, **`S1B`** | Periodic $\phi$/$\psi$ torsion error distributions (Asn–GlcNAc & chitobiose) | [`plotting/data/figS1_torsions.tsv`](plotting/data/figS1_torsions.tsv) |

*(Note: Manuscript subpanels **1A**, **1C**, **2D** [PDB `9HYU`], and **4D** [PDB `9N3S`] are 3D structural illustrations rather than statistical plots.)*

---

## Benchmark Dataset Format (`data/glycoproteins.json`)

Each element in [`data/glycoproteins.json`](data/glycoproteins.json) represents one curated single-chain glycoprotein input identified by its PDB accession (`pdb_id`) and mmCIF protein chain (`protein_chain`, i.e., `label_asym_id`), containing the standardized parent protein sequence and all experimentally resolved covalently attached glycans:

```json
{
  "pdb_id": "23GF",
  "protein_chain": "A",
  "protein_auth_chain": "A",
  "test_sets": ["A"],
  "sequence": "VTRHRVTVLYNAPEDIGNHMRQNDTHLTVRGGSGVVLQQR...",
  "glycosylation_sites": [
    {
      "label_seq_id": 23,
      "auth_seq_id": "24",
      "residue": "ASN",
      "link_type": "N",
      "protein_atom": "ND2",
      "root_sugar": "C:1",
      "root_atom": "C1",
      "sugars": [
        {"id": "C:1", "ccd": "NAG", "auth_chain": "C", "branch_num": 1},
        {"id": "C:2", "ccd": "NAG", "auth_chain": "C", "branch_num": 2}
      ],
      "bonds": [
        {"sugar_1": "C:1", "atom_1": "O4", "sugar_2": "C:2", "atom_2": "C1"}
      ]
    }
  ]
}
```

See [`data/README.md`](data/README.md) for the complete field-by-field specification and [`plotting/data/README.md`](plotting/data/README.md) for column definitions of all figure source tables.
