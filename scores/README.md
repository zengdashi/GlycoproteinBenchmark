# Scores and analysis units

All observations use the paper's common six-method inputs. Set A contains 327 inputs. Set B selects 57 of those inputs with `in_B=1`. The primary key for an input is `case_id`. Distinct chemical inputs can share a `sequence_id`. The related-input `group_id` identifies the unit for clustered bootstrap resampling.

| File | One row per | Rows | Main use |
| --- | --- | ---: | --- |
| [candidates.tsv.gz](candidates.tsv.gz) | Predictor × chemical input × candidate | 9,810 | Evaluate each of the five draws and identify the model's selected candidate. |
| [inputs.tsv](inputs.tsv) | Predictor × chemical input | 1,962 | Reproduce native selection, compare five-candidate selection strategies, and relate inter-candidate dispersion to error. |
| [selected_sites.tsv.gz](selected_sites.tsv.gz) | Predictor × chemical input × glycosylation site of its selected candidate | 5,304 | Inspect local geometry, positioning error and failures at individual attachment sites. |
| [aggregates/](aggregates/) | Manuscript statistic or stratum | varies | Source estimates, denominators, reference sensitivity, pooled results and archived confidence bounds used by the S1–S5 exporter. |

## Candidate and input scores

`candidate_id` identifies one method output. `candidate_index` is zero-based; `is_native_top1=1` marks the candidate chosen by the original method ranking. Native score scales are method-specific. `protein_CA_RMSD_A` and `protein_CA_lDDT` summarize the protein backbone. The four glycan RMSDs evaluate the same resolved glycan atoms after fitting the glycan itself, the root sugar ring, the local protein around that site, or the whole protein. Candidate-level glycan RMSDs are means across its evaluable sites. Distances are expressed in ångströms.

The 4.5 Å contact scores compare predicted and experimental noncovalent protein-residue–monosaccharide contacts. Attachment is assessed separately. `chemical_extended_local_pass_fraction` is the proportion of sites passing the published local chemistry criteria. Each `joint_extended_local_*A` is the proportion of sites meeting both those criteria and the named whole-protein-frame glycan RMSD threshold.

In `inputs.tsv`, `score` is the selected candidate's native ranking score, and `score_quartile` is its within-method set-A score quartile. `top1_rmsd`, `mean5_rmsd` and `best5_rmsd` summarize the selected candidate, mean of the five, and minimum of the five, respectively. All three refer to whole-protein-frame glycan positioning. `dispersion_protein` is the median over the ten candidate pairs of their site-averaged glycan heavy-atom RMSD after pairwise whole-protein Cα superposition. Raw scores are written at the precision retained in each source export.

## Site observations

Join a site to the glycan graph by `case_id` and `glycan_id`. `root_ccd` names the root sugar. `node_count` counts monosaccharides. `link_type` specifies the protein attachment type. The local columns describe the root sugar, attachment side chain or whole glycan in the local protein frame. `glycan_global` uses the whole-protein frame. For the root bond and angle, `_pred` and `_ref` are predicted and experimental values, and `_abs` is their absolute difference. Lengths use ångströms and angles use degrees. Torsion errors use the shortest periodic difference.

`root_supported` and `tree_supported` summarize bond-distance support. Defect categories can overlap. A site can fail the joint chemistry-and-positioning criterion and have none of the listed chemical defects if its positioning alone fails. Keep the observation unit in [S2](../tables/S2_site_chemistry.tsv) when comparing rates for bonds, sugars and entire sites. Missing or inapplicable values are marked `NA` or empty according to the source column.

For a typical structural mean, average sites within an input, then inputs within each exact protein sequence, then across sequences. Consult each supplementary table for exceptions such as raw counts, medians or pooled histograms. [The table index](../tables/README.md) identifies the statistic and denominator. [Plot-specific observations](../plotting/source_data/) supply the matched disaccharide scatter and torsion histograms.
