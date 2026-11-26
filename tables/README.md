# Supplementary tables S1–S5

These UTF-8 TSVs contain the numerical summaries submitted with the manuscript. The point-level scatter data, box-plot observations and histogram counts are indexed in the [plotting guide](../plotting/README.md). The five tables can be regenerated with `python scripts/export_tables.py` from the repository root.

The common evaluation set comprises 327 inputs, 256 distinct sequences and 884 glycosylation sites for A. B contains 57 of those inputs, 47 sequences and 171 sites. A single input is a protein chain with its specified attached glycans. Structural and chemical measurements use each method’s native top-ranked candidate unless candidate selection or agreement is being evaluated.

For sequence-weighted means, the eligible site, bond, center or ring observations for a metric are averaged within each input, inputs sharing the same sequence are averaged, and distinct sequences receive equal weight. The observation counts give the number of eligible physical objects before weighting. The reported weighted percentages can therefore differ from 100 times the pooled numerator divided by the pooled denominator. Bootstrap intervals resample sequence-related groups, preserving the association between inputs within a group. The full A set contains 195 groups and B contains 32. `NA` marks an unavailable or inapplicable value.

## S1. Structure and noncovalent contacts

`S1_structure_and_contacts.tsv` contains 12 rows, one per method and test set, and supports Fig. 2, Table 1, and the placement RMSD columns of Table 2.

- `protein_CA_RMSD_A` is whole-chain Cα RMSD in ångströms. `protein_CA_lDDT` reports Cα-lDDT on a 0–1 scale. The four `glycan_*_aligned_RMSD_A` columns evaluate glycan heavy-atom RMSD after alignment of glycan heavy atoms, the root ring, local protein Cα atoms or whole-chain protein Cα atoms.
- `contact_precision_4.5A`, `contact_recall_4.5A` and `contact_F1_4.5A` summarize noncovalent protein residue–monosaccharide contacts. A contact requires a heavy-atom pair within 4.5 Å after excluding atom pairs separated by three or fewer covalent bonds. Site precision is undefined without predicted contacts. Site recall and F1 are undefined without reference contacts. Undefined observations are omitted in their respective means.
- `backbone_lddt_vs_local_glycan_RMSD_rho` and `rho_ci95_*` describe the sequence-weighted Spearman correlation and its 95% interval from 2,000 clustered bootstrap resamples. The two `n_lddt_*` columns count individual inputs with Cα-lDDT ≥0.90 and, within that group, local-protein-aligned glycan RMSD >5 Å.
- For 104 A inputs with multiple equivalent experimental references, the final three columns give the input count, sequence-weighted mean reduction and largest single-input reduction in whole-protein-aligned placement RMSD when the closest eligible reference replaces the primary reference. These columns are populated for A only.

## S2. Covalent linkage geometry and site chemistry

`S2_site_chemistry.tsv` contains 310 rows and supports Figs. 3 and 4A,C, Supplementary Fig. S1 and the chemistry pass columns of Table 2. Read `test_set`, `subset`, `metric`, `statistic` and `observation_unit` together. `n` gives the observed denominator and `raw_numerator`, when available, gives the corresponding unweighted count.

| Metric or subset | Observation and statistic |
|---|---|
| `root_supported`, `internal_supported` | Protein–root-sugar and inter-sugar bonds with the specified atom-pair distance within 1.1–1.9 Å, inclusive. Fig. 3A draws on 884 root bonds and 947 inter-sugar bonds in A. |
| `tree_supported` | Complete oligosaccharides of at least two sugars for which **every internal sugar–sugar bond** meets that distance criterion. There are 400 eligible glycans in A. Root attachment is assessed separately. |
| `anomer_correct`, `ring_agreement`, `chemical_extended_local_pass_fraction` | Correct specified anomeric configuration, agreement with the experimental Privateer ring class, and site-level extended chemistry pass, respectively. In A the denominators are 1,831 anomeric centers, 1,827 classifiable rings and 884 glycosylation sites. |
| `attachment_valence_angle` | Unweighted site percentiles of the observed **angle** in experiment or prediction for 844 N-linked and 40 O-linked sites in A. These rows report angles in degrees. The absolute prediction error in the matched-disaccharide analysis is `root_angle_abs`. |
| Asn–GlcNAc disaccharides | Local-protein-aligned side-chain, root-sugar and glycan RMSD, absolute attachment-angle error and weighted Spearman associations for AlphaFold 3 and SweetFold at 204 paired sites from 134 inputs and 115 sequences. |
| `root_bond_outside_1.1_to_1.9_A` | Raw count of sites outside the inclusive root-bond length window among all 884 sites in A or 171 in B. |
| `fail_*` | Sequence-weighted prevalence within sites failing chemistry, whole-protein-aligned glycan placement (RMSD >2 Å), or both. Each method has its own failed-site denominator. The defect classes overlap. |
| `phi_absolute_periodic_error`, `psi_absolute_periodic_error` | Unweighted linkage-level medians of absolute shortest periodic errors. Supplementary Fig. S1 shows the two most frequent linkages in A. |
| `omega_absolute_periodic_error` | Sequence-weighted **mean** absolute periodic errors for the exocyclic acceptor torsion O6–C6–C5–O5 in applicable (1→6) linkage classes, reported in degrees. |

The extended chemistry check includes chain identity and coordinate completeness, glycan and attachment-residue heavy-atom completeness, covalent bond lengths (0.9–2.2 Å within a sugar and 1.1–1.9 Å between linked residues), specified stereochemistry and severe steric clashes. A severe clash has van der Waals overlap >0.6 Å after excluding pairs within three covalent bonds. Bootstrap intervals, where supplied for defect or bond metrics and ω errors, use 2,000 resamples.

## S3. Joint chemical and spatial outcomes

`S3_joint_chemistry_and_placement.tsv` contains 192 rows and supports Fig. 4B and the joint-pass columns of Table 2. For each test set, method and `positioning_cutoff_A`, a site passes placement when its whole-protein-aligned glycan RMSD is **at or below** the cutoff. The four mutually exclusive `outcome` values are both pass, chemistry only, placement only and neither. Fig. 4B uses 2 Å. Sensitivity analyses also use 1, 3 and 5 Å.

`observed_sites` gives the raw count in an outcome and `evaluated_sites` is 884 in A or 171 in B. The four counts sum to that denominator. `percent` gives the sequence-weighted percentage and `ci95_low_percent`, `ci95_high_percent` its 95% interval from 5,000 clustered bootstrap resamples.

## S4. Glycan length and branching

`S4_glycan_length_and_branching.tsv` contains 252 rows and supports Fig. 5A and the linear–branched comparison. Length is the number of experimentally resolved monosaccharides at a site. A glycan is branched if at least one monosaccharide has two or more downstream sugars. `value_A` measures mean self-aligned glycan RMSD or a branched-minus-linear difference in ångströms.

| `analysis` | Population, weighting and intervals |
|---|---|
| `by_length` | Sites within each exact observed length and test set. Means use equal sequence weights. `ci95_low_A` and `ci95_high_A` come from 2,000 clustered bootstrap resamples when at least five sequence-related groups contribute. Fig. 5A plots A means with the raw number of sites beneath each length. |
| `unadjusted_all_oligosaccharides` | Linear and branched glycans of length ≥2 in each set, pooled across their observed lengths. Topology estimates and their difference use 5,000 clustered bootstrap resamples. |
| `length_standardized` | Each retained length has at least five sequence-related groups in each topology. Lengths 3 and 4 qualify in A. Their shared pooled **sequence-weighted** frequencies (0.80247 and 0.19753) weight the within-length estimates for each topology. No length qualified in B. These intervals also use 5,000 resamples. |

## S5. Candidate selection, scores and agreement

`S5_scores_selection_and_agreement.tsv` contains 156 rows and supports Fig. 5B–D. Placement error is whole-protein-aligned, site-mean glycan RMSD. All six methods generated five candidates for each evaluated input.

| `analysis` | Interpretation |
|---|---|
| `score_quartile` | Inputs are partitioned separately within each method and test set using sequence-weighted quartiles of the native top-ranked score. Q1 is lowest. `min_score` and `max_score` are the observed extremes within the quartile. `value` is its sequence-weighted mean placement error. |
| `candidate_choice` | `native_top1`, `random_expected` and `full_candidate_oracle` compare the native top candidate, the arithmetic mean of five placement errors and the smallest of the five errors within each input. The `inputs_with_better_candidate` row is a raw proportion: 100 × the number of inputs with a strictly smaller error among the five than the native top rank, divided by `n_inputs`. |
| `candidate_agreement` | Within each input, spread is the median over all ten candidate pairs of site-averaged glycan RMSD after pairwise whole-protein Cα alignment. The sequence-weighted correlation with top-ranked placement error is reported as `dispersion_rho`. `dispersion_partial_rho` adjusts the ranked variables for native score and mean input glycan length. |
| `score_spread` | Sequence-weighted median within-input five-candidate score range, sequence-weighted interquartile range of top-ranked scores across inputs, and their ratio. |

For quartile means, candidate-choice mean RMSD and correlations, `ci95_low` and `ci95_high` give the 2.5th and 97.5th percentiles of 2,000 clustered bootstrap resamples. The table reports the descriptive score-spread quantities and the raw fraction of inputs with a better candidate without intervals.
