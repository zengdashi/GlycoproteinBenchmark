# Coordinate-score verification

The packaged Figure 4D case is 9N3S chain B, scored against Chai-1 candidate 001. File hashes use SHA-256.

| File | SHA-256 |
|---|---|
| `example/9N3S.cif` | `d8d808b14c85e5aac69bec3d13ef865029ac2cfef638b0c49ad9869dae4c9ed4` |
| `example/Chai1_c001.cif` | `6a255120248cb104f8d760011f21e1a75de21f55f6b0e4e522018d64d46ca311` |
| `kernels/geometry_metrics.py` | `262c13caed28b918cbd45e8c07b027c8dc8d2d3df76538e485e482fafbe083a6` |
| `kernels/physics_metrics.py` | `2318a3f0645932cd6076228cd6406926690d0a301bbc72bb07f27c93d295f3c3` |

## Bundled test

`python scoring/score.py --example` tests the following values against files included in the release.

- **14 input-level values.** Compare the two protein metrics, four glycan frames, three nonbonded contact summaries, extended chemistry pass fraction, and four joint outcomes against `scores/candidates.tsv.gz` with absolute tolerance 10⁻⁷.
- **24 site-level values.** At each of the two distinct sites, compare three glycan frames, local anchor-side-chain RMSD, reference/predicted/absolute-error values for the root distance and bond angle, and absolute root φ and ψ errors against `scores/selected_sites.tsv.gz`. Rounded source fields use absolute tolerance 5.1 × 10⁻⁵.
- **Four displayed-site frames.** Compare G002 in the whole-protein, local-site, root-ring and glycan-self frames. The corresponding values are 0.546103, 0.401107, 0.184264 and 0.114920 Å. The original figure table used six decimal places.

Native mapping reproduced all 3,531 model atom IDs in the supplied TSV. The experimental mask has 2,923 observed heavy atoms and 362 protein Cα atoms. The G002 glycan contains 28 scored heavy atoms. Exporting the anonymous primary reference, remapping the model and rescoring returned identical values for all 14 input-level fields.

## Additional release checks

Five original benchmark predictions were independently mapped from their native files, then scored against the anonymous primary references. Each row covers all 14 coordinate-derived candidate fields in `scores/candidates.tsv.gz`.

| Method | Case | Candidate index | Largest absolute difference |
|---|---|---:|---:|
| AlphaFold 3 | GLYCO_ed0a8c4324b1bc39 | 0 | 5.78 × 10⁻¹⁵ |
| SweetFold | GLYCO_ed0a8c4324b1bc39 | 0 | 3.00 × 10⁻¹⁵ |
| AlphaFold 3 | GLYCO_efd2eb0fa029901d | 0 | 4.45 × 10⁻¹⁵ |
| Boltz-2 | GLYCO_38627a5f10e0e98f | 0 | 2.00 × 10⁻¹⁵ |
| Chai-1 | GLYCO_38627a5f10e0e98f | 1 | 8.89 × 10⁻¹⁶ |

Two additional native smoke-test predictions exercised the remaining readers. The input-specific RFAA runtime map resolved 3,412 atoms and the Protenix reader resolved 1,116 atoms. Both produced all 14 fields. These two checks establish reader execution; the smoke-test structures were generated separately from the formal predictions.

A consistency test supplied five copies of the same native file and map to `score_ensemble.py`. It returned ten pairs and a median of 1.59 × 10⁻¹⁴ Å (zero to numerical precision). S5 uses the median for five distinct model samples.
