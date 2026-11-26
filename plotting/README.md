# Reproducing the statistical panels

The repository entry point, `python plot.py all`, also generates the data-derived test-set panel (Fig. 1B), for 17 panels in total. To render only the 16 statistical result panels in Figs. 2–5 and S1, run the commands below from the repository root. Each command writes a vector PDF with editable labels and a 500 dpi PNG. The Python entry point locates its input tables relative to the source file. It can run offline on a CPU.

```bash
python plotting/plot.py 3D --out-dir results/figures
python plotting/plot.py all --out-dir results/figures
python plotting/check_data.py
```

Use `--dpi 150` for a smaller preview PNG. Output filenames follow the panel labels, for example `Fig3D.pdf` and `Fig3D.png`. The reference rendering used Python 3.9.12, NumPy 1.21.5, pandas 1.4.2 and Matplotlib 3.5.1 with its bundled DejaVu fonts. Other graphics-library versions can introduce small typographic differences.

## Panel-to-data index

In this index, S1–S5 refer to the TSVs in `tables/`. The default is each method’s native top-ranked structure. Panels 5C and 5D also use the five candidates generated for each input. Set A contains 327 inputs, 256 distinct sequences and 884 sites. Nested set B contains 57 inputs, 47 sequences and 171 sites. For sequence-weighted means, eligible observations are averaged within each input and inputs within each identical sequence. Each distinct sequence then contributes equally. Raw counts, site-angle percentiles and unweighted linkage-error histograms are marked below.

| Panel | Data read by the renderer | What is shown |
|---|---|---|
| 2A | S1 | Mean protein Cα RMSD and Cα lDDT for the two test sets. |
| 2B | `scores/candidates.tsv.gz`, rows with `is_native_top1=1` | Set A distributions of sequence means in four glycan alignment frames. Boxes show quartiles and medians, whiskers the 5th and 95th percentiles, and diamonds the means. |
| 2C | The same candidate rows; S1 | One point per input, relating backbone lDDT to local-protein-aligned glycan RMSD. The sequence-weighted log-linear fit is calculated from the points. The annotation gives the weighted Spearman correlation. |
| 3A | S2 | Sequence-weighted pass rates for 884 protein–sugar and 947 sugar–sugar bonds in the inclusive 1.1–1.9 Å window. The right axis tests all internal sugar–sugar bonds simultaneously in each of 400 oligosaccharides. |
| 3B | S2 | Sequence-weighted agreement in anomeric configuration (1,831 centers), Privateer ring class (1,827 rings) and extended site chemistry (884 sites). |
| 3C | S2 | Unweighted 10th, 25th, 50th, 75th and 90th percentiles of attachment angles at 844 N-linked and 40 O-linked sites. Shading marks empirical N- and O-specific reference bands. |
| 3D | `source_data/fig3_disaccharides.tsv`; S2 | AlphaFold 3 and SweetFold on 204 matched disaccharide sites per method. The site records provide the points and S2 supplies the sequence-weighted diamonds. |
| 4A | S2 | Raw site counts outside the protein–sugar bond-length window in sets A and B. |
| 4B | S3 | Sequence-weighted percentages of the 884 A sites in four mutually exclusive outcomes. Placement passes at whole-protein-aligned glycan RMSD ≤2 Å. |
| 4C | S2 | Sequence-weighted prevalence of overlapping defect classes among sites failing chemistry, placement, or both. Denominators are method-specific failed sites, 732–874 in A. |
| 5A | S4 | Mean glycan-aligned RMSD by observed chain length, with site counts below. Lengths represented by fewer than ten sites use individual markers. |
| 5B | S5 | Sequence-weighted mean whole-protein-aligned glycan RMSD by within-method, sequence-weighted native score quartile (Q1 lowest). Bars give 95% clustered bootstrap intervals. |
| 5C | S5 | Native top rank, expected random choice (mean of five), and reference-guided best of five within each input. |
| 5D | `scores/inputs.tsv`; S5 | AlphaFold 3 input-level spread versus placement. Spread is the median of ten pairwise, site-mean glycan RMSDs after pairwise whole-protein Cα alignment. Six sequence-weighted partial Spearman correlations adjust for native score and mean glycan length and are shown with 95% clustered bootstrap intervals. |
| S1A | `source_data/figS1_histograms.tsv`; S2 | Absolute periodic φ and ψ error for the Asn–GlcNAc protein attachment, 803 linkages per method. |
| S1B | The same histogram table; S2 | Absolute periodic φ and ψ error for GlcNAc β(1→4) GlcNAc, 387 linkages per method. |

## Two plot-specific source tables

- `fig3_disaccharides.tsv` contains one record per method and matched glycosylation site. The fields `candidate`, `case`, `site`, `seq` and `group` link observations to benchmark records. Attachment-angle values are in degrees. Side-chain, root-sugar and whole-glycan errors are in ångströms after local protein alignment. The supplied weights average sites within inputs and inputs within sequences.
- `figS1_histograms.tsv` contains observed integer counts in 18 bins of width 10° for each method, linkage and torsion. The renderer divides by the valid-edge count within that series. Median markers use the original unbinned medians in S2. The histograms record absolute periodic errors in the range 0–180°.

`check_data.py` connects these plotting inputs to the submitted tables. It recomputes the structural means and weighted rank correlations in S1 for both sets, checks the matched-site means and torsion histogram counts against S2, and independently reproduces all 24 quartile means in Fig. 5B and all six partial correlations in Fig. 5D. For both panels it recomputes the 95% clustered bootstrap intervals from the input-level observations using the original 2,000 resamples and seeds. The renderer reads these plotted estimates and intervals directly from S5. For release verification, each of the 16 data-generated PDFs was rasterized at 144 dpi and compared with the corresponding original panel embedded in the submission figures. All 16 pixel comparisons matched exactly in the reference environment.
