# Figure Source Tables (`plotting/data/`)

This directory contains the 14 tab-separated value (`.tsv`) files read by [`plot.py`](../../plot.py) to generate all statistical panels in **Figures 2–5** and **Supplementary Figure S1** (note that **Figure 1B** is computed directly from [`../../data/glycoproteins.json`](../../data/glycoproteins.json)).

---

## Weighting & Statistical Conventions

- **Evaluated predictors (`method`):** `AlphaFold 3`, `Boltz-2`, `Chai-1`, `Protenix`, `RFAA`, and `SweetFold` (plus `Experiment` in `fig3C_attachment_angles.tsv`).
- **Sequence-weighted aggregation:** Unless otherwise noted, summary percentages, means, and Spearman rank correlations weight distinct protein sequences equally (averaging across sites within each input, then across inputs sharing an identical protein sequence, and finally across the 256 distinct protein sequences in Test set A or 47 in Test set B).
- **Confidence intervals (`ci95_low`, `ci95_high`):** Two-sided 95% intervals estimated via clustered bootstrap resampling (2,000 resamples) over 195 sequence-related PDB groups.

---

## File-by-File Schema

| File | Panel(s) | Rows | Observation Unit & Columns |
| --- | --- | ---: | --- |
| [`fig2_structure.tsv`](fig2_structure.tsv) | `2A`, `2B`, `2C` | 1,962 | **One row per predictor × Test set A input (`6 × 327`).** Columns: `pdb_id`, `protein_chain`, `method`, `protein_CA_RMSD_A` (whole-chain Cα RMSD, Å), `protein_CA_lDDT` (whole-chain Cα-lDDT), and site-mean glycan heavy-atom RMSDs (Å) across the four alignment frames: `glycan_self_frame_macro_RMSD_A` (glycan heavy-atom self-alignment), `glycan_root_frame_macro_RMSD_A` (root sugar ring), `glycan_site_frame_macro_RMSD_A` (local protein Cα), and `glycan_protein_frame_macro_RMSD_A` (whole-protein Cα). Joined with `data/glycoproteins.json` during plotting to obtain `sequence` and `in_set_b`. |
| [`fig3A_bond_lengths.tsv`](fig3A_bond_lengths.tsv) | `3A` | 18 | **Predictor × bond metric summary on Test set A.** Columns: `method`, `metric` (`root_supported` across 884 protein–sugar bonds, `internal_supported` across 947 inter-sugar bonds, `tree_supported` across 400 oligosaccharides), `n`, `percent` (sequence-weighted pass rate within 1.1–1.9 Å). |
| [`fig3B_stereochemistry.tsv`](fig3B_stereochemistry.tsv) | `3B` | 18 | **Predictor × stereochemical metric summary on Test set A.** Columns: `method`, `metric` (`anomer_correct` over 1,831 centers, `ring_agreement` over 1,827 Privateer-classifiable rings, `chemical_extended_local_pass_fraction` over 884 sites), `percent` (sequence-weighted %). |
| [`fig3C_attachment_angles.tsv`](fig3C_attachment_angles.tsv) | `3C` | 70 | **Site-level attachment valence angle percentiles on Test set A.** Columns: `method` (`Experiment` + 6 predictors), `linkage` (`N` for 844 Asn CG–ND2–C1 sites; `O` for 40 Ser/Thr CB–OG/OG1–C1 sites), `percentile` (`10`, `25`, `50`, `75`, `90`), `angle_deg` (degrees). |
| [`fig3D_disaccharides.tsv`](fig3D_disaccharides.tsv) | `3D` | 408 | **Matched site-level comparison across 204 Asn–GlcNAc disaccharide sites (`2 × 204` for `AlphaFold 3` and `SweetFold`).** Columns: `pdb_id`, `protein_chain`, `label_seq_id`, `method`, `angle_error_deg` (absolute attachment-angle error, degrees), `sidechain_RMSD_A` (local-protein-aligned Asn side-chain RMSD, Å), `local_glycan_RMSD_A` (local-protein-aligned disaccharide RMSD, Å). |
| [`fig4A_root_bonds.tsv`](fig4A_root_bonds.tsv) | `4A` | 12 | **Out-of-range protein–root-sugar bond counts.** Columns: `test_set` (`A` or `B`), `method`, `total_sites` (`884` or `171`), `sites_outside_range` (count of sites outside 1.1–1.9 Å). |
| [`fig4B_joint_outcomes.tsv`](fig4B_joint_outcomes.tsv) | `4B` | 24 | **Mutually exclusive site outcomes on Test set A (2 Å cutoff).** Columns: `method`, `outcome` (`chem_pass_geo_pass`, `chem_pass_geo_fail`, `chem_fail_geo_pass`, `chem_fail_geo_fail`), `percent` (sequence-weighted %, summing to 100% per method). |
| [`fig4C_defects.tsv`](fig4C_defects.tsv) | `4C` | 30 | **Defect prevalences among Test set A sites failing joint criteria.** Columns: `method`, `metric` (`fail_root`, `fail_internal`, `fail_sugar_stereo`, `fail_clash`, `fail_sugar_stereo_and_clash`), `percent` (sequence-weighted %, non-exclusive across columns). |
| [`fig5A_glycan_length.tsv`](fig5A_glycan_length.tsv) | `5A` | 84 | **Internal glycan error by observed monosaccharide count on Test set A (`6 × 14` observed lengths).** Columns: `method`, `length` (monosaccharides per site, 1–16), `n_sites` (site count), `glycan_RMSD_A` (sequence-weighted mean self-aligned RMSD, Å). |
| [`fig5B_score_quartiles.tsv`](fig5B_score_quartiles.tsv) | `5B` | 24 | **Placement error across native top-ranked score quartiles on Test set A.** Columns: `method`, `quartile` (`1` lowest to `4` highest), `placement_RMSD_A` (sequence-weighted mean whole-protein-aligned glycan RMSD, Å), `ci95_low`, `ci95_high`. |
| [`fig5C_candidate_choice.tsv`](fig5C_candidate_choice.tsv) | `5C` | 18 | **Five-candidate selection strategies on Test set A.** Columns: `method`, `choice` (`native_top1`, `random_expected`, `full_candidate_oracle`), `placement_RMSD_A` (sequence-weighted mean whole-protein-aligned glycan RMSD, Å). |
| [`fig5D_AlphaFold3.tsv`](fig5D_AlphaFold3.tsv) | `5D` (left) | 327 | **Per-input spread vs. error for AlphaFold 3 across all 327 Test set A inputs.** Columns: `pdb_id`, `protein_chain`, `candidate_spread_A` (median across 10 candidate pairs of site-mean whole-protein-aligned glycan RMSD, Å), `placement_RMSD_A` (top-ranked whole-protein-aligned glycan RMSD, Å). |
| [`fig5D_correlations.tsv`](fig5D_correlations.tsv) | `5D` (right) | 6 | **Partial rank correlations across Test set A inputs.** Columns: `method`, `partial_spearman_rho` (sequence-weighted partial Spearman $\rho$ between candidate spread and top-ranked placement error controlling for native score and mean glycan length), `ci95_low`, `ci95_high`. |
| [`figS1_torsions.tsv`](figS1_torsions.tsv) | `S1A`, `S1B` | 24 | **Periodic linkage torsion error histograms on Test set A.** Columns: `linkage` (`NAG|C1|ASN|ND2|root` with $n=803$; `NAG|C1|NAG|O4|interglycan` with $n=387$), `method`, `torsion` (`phi` or `psi`), `median_deg` (per-linkage median absolute error, degrees), and $10^\circ$ bin counts `n_0_10` through `n_170_180`. |
