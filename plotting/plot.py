"""Render publication panels from the submitted tables and case/site-level scores."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figure2
import figure3
import figure4
import figure5
from style import METHODS, configure_style

ROOT = Path(__file__).resolve().parent.parent
DATA = Path(__file__).resolve().parent / "source_data"
TABLES = {
    1: "S1_structure_and_contacts.tsv",
    2: "S2_site_chemistry.tsv",
    3: "S3_joint_chemistry_and_placement.tsv",
    4: "S4_glycan_length_and_branching.tsv",
    5: "S5_scores_selection_and_agreement.tsv",
}
PANELS = (
    "2A",
    "2B",
    "2C",
    "3A",
    "3B",
    "3C",
    "3D",
    "4A",
    "4B",
    "4C",
    "5A",
    "5B",
    "5C",
    "5D",
    "S1A",
    "S1B",
)


def read(path):
    return pd.read_csv(path, sep="\t")


def supplement(number):
    return read(ROOT / "tables" / TABLES[number])


def trend_lines(inputs, summary):
    """Fit log(local glycan RMSD) on lDDT with each protein sequence carrying equal weight."""
    fits = []
    for method, block in inputs.groupby("method"):
        x = block.protein_CA_lDDT.to_numpy()
        y = block.glycan_site_frame_macro_RMSD_A.to_numpy()
        weights = (1 / block.groupby("sequence_id").case_id.transform("size")).to_numpy()
        design = np.column_stack([np.ones(len(x)), x]) * np.sqrt(weights)[:, None]
        intercept, slope = np.linalg.lstsq(design, np.log(y) * np.sqrt(weights), rcond=None)[0]
        order = np.argsort(x, kind="stable")
        cumulative = np.cumsum(weights[order])
        bounds = np.searchsorted(cumulative, np.array([0.025, 0.975]) * weights.sum(), side="left")
        low, high = x[order][bounds]
        fits.append(
            dict(
                method=METHODS[method],
                log_rmsd_intercept=intercept,
                log_rmsd_slope_per_lddt=slope,
                trend_xmin=low,
                trend_xmax=high,
            )
        )
    stats = summary[summary.test_set == "A"].rename(
        columns={
            "backbone_lddt_vs_local_glycan_RMSD_rho": "weighted_spearman_rho",
        }
    )
    return stats.merge(pd.DataFrame(fits), on="method", validate="one_to_one")


def render_panel(key):
    plt.rcdefaults()
    if key.startswith("2"):
        figure2.configure_style()
        if key == "2A":
            return figure2.backbone(supplement(1))
        inputs = read(ROOT / "scores" / "candidates.tsv.gz").query("is_native_top1 == 1").copy()
        inputs["method"] = inputs.method.map({name: key for key, name in METHODS.items()})
        if key == "2B":
            return figure2.frames(inputs)
        return figure2.assoc_panel(inputs, trend_lines(inputs, supplement(1)))

    configure_style()
    if key.startswith(("3", "S1")):
        chemistry = supplement(2).query("test_set == 'A'")
        if key == "3A":
            return figure3.panel_a(chemistry)
        if key == "3B":
            return figure3.panel_b(chemistry)
        if key == "3C":
            return figure3.panel_c(chemistry)
        if key == "3D":
            return figure3.panel_d(chemistry, read(DATA / "fig3_disaccharides.tsv"))
        return figure3.torsions(key[-1], read(DATA / "figS1_histograms.tsv"), chemistry)

    if key == "4A":
        roots = supplement(2).query("metric == 'root_bond_outside_1.1_to_1.9_A'")
        return figure4.a(roots)
    if key == "4B":
        outcomes = supplement(3).query("test_set == 'A' and positioning_cutoff_A == 2")
        return figure4.b(outcomes.rename(columns={"percent": "value", "outcome": "metric"}))
    if key == "4C":
        defects = supplement(2).query(
            "test_set == 'A' and subset == 'sites with a chemistry or positioning error'"
        )
        return figure4.c(defects)
    if key == "5A":
        length = supplement(4).query("test_set == 'A' and analysis == 'by_length'").copy()
        length["stratum"] = pd.to_numeric(length["length"])
        return figure5.a(length.rename(columns={"value_A": "value"}))

    scores = (
        supplement(5)
        .query("test_set == 'A'")
        .rename(columns={"ci95_low": "ci_low", "ci95_high": "ci_high"})
    )
    if key == "5B":
        quartiles = scores[scores.analysis == "score_quartile"].copy()
        quartiles["stratum"] = quartiles.stratum.str.removeprefix("Q").astype(int)
        return figure5.b(quartiles)
    if key == "5C":
        choice = scores.query(
            "analysis == 'candidate_choice' and measure == 'mean_positioning_RMSD'"
        ).copy()
        choice["metric"] = choice.stratum.map(
            {
                "native_top1": "native_rmsd",
                "random_expected": "mean5_rmsd",
                "full_candidate_oracle": "best5_rmsd",
            }
        )
        return figure5.c(choice)
    if key == "5D":
        agreement = scores[scores.measure == "dispersion_partial_rho"]
        return figure5.d(agreement, read(ROOT / "scores/inputs.tsv"))
    raise ValueError(f"Unknown panel {key}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "panel", choices=(*PANELS, "all"), help="Panel label or all statistical panels"
    )
    parser.add_argument("--out-dir", type=Path, default=Path("results/figures"))
    parser.add_argument("--dpi", type=int, default=500, help="PNG resolution (default 500)")
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for key in PANELS if args.panel == "all" else [args.panel]:
        fig = render_panel(key)
        try:
            stem = args.out_dir / f"Fig{key}"
            fig.savefig(stem.with_suffix(".pdf"))
            fig.savefig(stem.with_suffix(".png"), dpi=args.dpi)
            print(stem)
        finally:
            plt.close(fig)


if __name__ == "__main__":
    main()
