"""Cross-check the observation-level plotting files against deposited table values."""

from pathlib import Path
import zlib

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SOURCE = Path(__file__).resolve().parent / "source_data"


def read_table(path):
    return pd.read_csv(path, sep="\t")


def read_supplement(n):
    (path,) = (ROOT / "tables").glob(f"S{n}_*.tsv")
    return read_table(path)


def sequence_weights(frame):
    return 1.0 / frame.groupby("sequence_id").case_id.transform("size").to_numpy()


def bootstrap_group_counts(frame):
    groups, index = np.unique(frame.group_id, return_inverse=True)
    seed = 20261123 + zlib.crc32("|".join(groups).encode("utf-8"))
    generator = np.random.default_rng(seed)
    draws = generator.multinomial(len(groups), np.full(len(groups), 1 / len(groups)), size=2000)
    return draws, index


def rank_estimator(frame, x, y, controls=()):
    """Construct a weighted rank statistic; tied observations have the same midrank."""
    indices = [np.unique(frame[key], return_inverse=True)[1] for key in (x, y, *controls)]

    def estimate(weight):
        ranks = []
        for index in indices:
            mass = np.bincount(index, weights=weight)
            ranks.append((mass.cumsum() - 0.5 * mass)[index] / weight.sum())
        variables = np.column_stack(ranks[:2])
        design = np.column_stack([np.ones(len(frame)), *ranks[2:]])
        root_weight = np.sqrt(weight[:, None])
        coefficients = np.linalg.lstsq(design * root_weight, variables * root_weight, rcond=None)[0]
        residuals = variables - design @ coefficients
        covariance = (residuals * weight[:, None]).T @ residuals
        return covariance[0, 1] / np.sqrt(covariance[0, 0] * covariance[1, 1])

    return estimate


def bootstrap_mean(group):
    weights = sequence_weights(group)
    draws, index = bootstrap_group_counts(group)
    total = np.bincount(index, weights=weights * group.top1_rmsd)
    mass = np.bincount(index, weights=weights)
    means = (draws @ total) / (draws @ mass)
    return np.r_[total.sum() / mass.sum(), np.quantile(means, [0.025, 0.975])]


def check_backbone_and_frames():
    values = read_table(ROOT / "scores/candidates.tsv.gz").query("is_native_top1 == 1")
    source_to_summary = {
        "protein_CA_RMSD_A": "protein_CA_RMSD_A",
        "protein_CA_lDDT": "protein_CA_lDDT",
        "glycan_self_frame_macro_RMSD_A": "glycan_self_aligned_RMSD_A",
        "glycan_root_frame_macro_RMSD_A": "glycan_root_ring_aligned_RMSD_A",
        "glycan_site_frame_macro_RMSD_A": "glycan_local_protein_aligned_RMSD_A",
        "glycan_protein_frame_macro_RMSD_A": "glycan_whole_protein_aligned_RMSD_A",
    }
    for row in read_supplement(1).itertuples():
        group = values[values.method.eq(row.method)]
        if row.test_set == "B":
            group = group[group.in_B.eq(1)]
        per_sequence = group.groupby("sequence_id")[list(source_to_summary)].mean()
        assert len(per_sequence) == row.n_sequences and group.case_id.nunique() == row.n_inputs
        for source, target in source_to_summary.items():
            np.testing.assert_allclose(
                per_sequence[source].mean(), getattr(row, target), atol=1e-9, rtol=0
            )
        estimate = rank_estimator(group, "protein_CA_lDDT", "glycan_site_frame_macro_RMSD_A")
        rho = estimate(sequence_weights(group))
        np.testing.assert_allclose(
            rho, row.backbone_lddt_vs_local_glycan_RMSD_rho, atol=1e-8, rtol=0
        )
    print("Figure 2: both test sets, structural means and weighted correlations match S1.")


def check_chemistry():
    summary = read_supplement(2)
    sites = read_table(SOURCE / "fig3_disaccharides.tsv")
    selection = summary[
        (summary.subset == "Asn–GlcNAc disaccharides")
        & (summary.statistic == "sequence-weighted mean")
    ]
    for row in selection.itertuples():
        group = sites[sites.method.eq(row.method)]
        assert len(group) == row.n and group.seq.nunique() == row.n_sequences
        np.testing.assert_allclose(
            np.average(group[row.metric], weights=group.weight), row.value, atol=1e-7, rtol=0
        )
    bins = read_table(SOURCE / "figS1_histograms.tsv")
    methods = {
        "af3": "AlphaFold 3",
        "boltz2": "Boltz-2",
        "chai1": "Chai-1",
        "protenix": "Protenix",
        "rfaa": "RFAA",
        "sweetfold": "SweetFold",
    }
    for (link, model, angle), group in bins.groupby(["link_type", "method", "angle"]):
        (row,) = summary[
            (summary.subset == link)
            & (summary.method == methods[model])
            & (summary.metric == angle + "_absolute_periodic_error")
            & (summary.statistic == "median")
        ].itertuples()
        np.testing.assert_array_equal(group.left_degrees, np.arange(0, 180, 10))
        np.testing.assert_array_equal(group.right_degrees, np.arange(10, 190, 10))
        assert group["count"].sum() == row.n and group.valid_edges.eq(row.n).all()
        # With an even sample, the median may fall between the two central bins.
        indices = np.searchsorted(group["count"].cumsum(), [(row.n + 1) // 2, row.n // 2 + 1])
        assert (
            group.left_degrees.iloc[indices[0]] <= row.value <= group.right_degrees.iloc[indices[1]]
        )
    print("Figures 3D and S1: site-weighted means, torsion-bin counts and median bins match S2.")


def check_candidate_agreement():
    inputs = read_table(ROOT / "scores/inputs.tsv")
    table = read_supplement(5).query("test_set == 'A'")
    for row in table.query("analysis == 'score_quartile'").itertuples():
        quartile = int(row.stratum[1:])
        group = inputs[inputs.method.eq(row.method) & inputs.score_quartile.eq(quartile)]
        np.testing.assert_allclose(
            bootstrap_mean(group), [row.value, row.ci95_low, row.ci95_high], atol=1e-8, rtol=0
        )
    for row in table.query("measure == 'dispersion_partial_rho'").itertuples():
        group = inputs[inputs.method.eq(row.method)]
        assert len(group) == row.n_inputs and group.sequence_id.nunique() == row.n_sequences
        estimate = rank_estimator(
            group, "dispersion_protein", "top1_rmsd", ["score", "mean_length"]
        )
        weight = sequence_weights(group)
        draws, index = bootstrap_group_counts(group)
        bootstrap_rho = [estimate(weight * draw[index]) for draw in draws]
        observed = np.r_[estimate(weight), np.quantile(bootstrap_rho, [0.025, 0.975])]
        np.testing.assert_allclose(
            observed, [row.value, row.ci95_low, row.ci95_high], atol=1e-8, rtol=0
        )
    print(
        "Figures 5B/D: all 24 score quartiles, six partial correlations and their 95% bootstrap intervals match S5."
    )


if __name__ == "__main__":
    check_backbone_and_frames()
    check_chemistry()
    check_candidate_agreement()
