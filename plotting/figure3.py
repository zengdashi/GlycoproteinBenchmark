"""Figure 3 and Supplementary Figure S1. Chemistry and linkage geometry."""

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator
from style import METHODS, COLORS, INK, GRID, panel, style

SUPP_SIZE = (5.1, 5.1)
LINKS = {
    "NAG|C1|ASN|ND2|root": ("A", "Asn–GlcNAc attachment"),
    "NAG|C1|NAG|O4|interglycan": ("B", "GlcNAc β(1→4) GlcNAc"),
}
ANGLE_COLORS = {"phi": "#176d70", "psi": "#ac5939"}
D_METHODS = ["af3", "sweetfold"]


def panel_a(summary):
    f = panel("A", "Covalent bond lengths")
    rows = summary[summary.metric.isin(["root_supported", "internal_supported", "tree_supported"])]
    values = rows.set_index(["method", "metric"]).value
    counts = rows[["metric", "n"]].drop_duplicates().set_index("metric").n
    bonds = f.add_axes([0.206, 0.26, 0.326, 0.49])
    glycans = f.add_axes([0.634, 0.26, 0.326, 0.49])
    bonds.set_title("Per bond", fontsize=10.4, weight="bold", pad=7)
    glycans.set_title(
        f"Per glycan (n = {counts['tree_supported']:.0f})", fontsize=10.4, weight="bold", pad=7
    )
    bonds.set_xlabel("Bonds in range (%)", fontsize=9.1)
    glycans.set_xlabel("Glycans with all sugar–sugar\nbonds in range (%)", fontsize=9.1)
    for ax in (bonds, glycans):
        ax.set(
            yticks=range(len(METHODS)),
            ylim=(5.55, -0.55),
            xlim=(95.55, 100.43),
            xticks=[96, 98, 100],
        )
        style(ax, "x")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    bonds.set_yticklabels(METHODS.values())
    glycans.tick_params(axis="y", labelleft=False)
    for i, (m, name) in enumerate(METHODS.items()):
        color = COLORS[m]
        bonds.scatter(
            values.loc[name, "root_supported"],
            i - 0.13,
            marker="o",
            s=23,
            facecolor=color,
            edgecolor=color,
            linewidth=0.85,
            zorder=3,
        )
        bonds.scatter(
            values.loc[name, "internal_supported"],
            i + 0.13,
            marker="s",
            s=23,
            facecolor="white",
            edgecolor=color,
            linewidth=0.95,
            zorder=3,
        )
        glycans.scatter(
            values.loc[name, "tree_supported"],
            i,
            marker="D",
            s=24,
            facecolor=color,
            edgecolor=color,
            linewidth=0.85,
            zorder=3,
        )
    legend = [
        Line2D(
            [],
            [],
            marker="o",
            ls="none",
            color=INK,
            ms=4.5,
            label=f"Protein–sugar (n = {counts['root_supported']:.0f})",
        ),
        Line2D(
            [],
            [],
            marker="s",
            ls="none",
            markerfacecolor="white",
            markeredgecolor=INK,
            ms=4.5,
            label=f"Sugar–sugar (n = {counts['internal_supported']:.0f})",
        ),
    ]
    f.legend(
        handles=legend,
        loc="center",
        bbox_to_anchor=(0.301, 0.054),
        ncol=2,
        fontsize=8.1,
        handlelength=0.8,
        handletextpad=0.32,
        columnspacing=1,
    )
    return f


def panel_b(summary):
    f = panel("B", "Stereochemistry and ring conformation")
    ax = f.add_axes([0.205, 0.205, 0.76, 0.59])
    data = (
        summary[
            summary.metric.isin(
                ["anomer_correct", "ring_agreement", "chemical_extended_local_pass_fraction"]
            )
        ]
        .set_index(["method", "metric"])
        .value
    )
    cols = [
        ("anomer_correct", "Anomeric\nconfiguration", "#2f716e"),
        ("ring_agreement", "Ring class\nagreement", "#598bab"),
        ("chemical_extended_local_pass_fraction", "Site chemistry\npass rate", "#a37147"),
    ]
    for j, (metric, label, color) in enumerate(cols):
        for i, name in enumerate(METHODS.values()):
            val = data.loc[name, metric]
            origin = j - 0.37
            ax.hlines(i, origin, origin + 0.57, color=GRID, lw=2.7, zorder=1)
            ax.hlines(i, origin, origin + 0.57 * val / 100, color=color, lw=3.0, zorder=2)
            ax.text(origin + 0.62, i, f"{val:.1f}", ha="left", va="center", fontsize=7.6)
    ax.set(
        xlim=(-0.45, 2.63),
        ylim=(5.55, -0.55),
        yticks=range(6),
        yticklabels=METHODS.values(),
        xticks=range(3),
        xticklabels=[label for _, label, _ in cols],
    )
    ax.tick_params(axis="both", length=0, pad=4)
    ax.tick_params(axis="x", labelsize=8.7)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return f


def panel_c(summary):
    f = panel("C", "Protein–sugar attachment angles")
    data = (
        summary[summary.metric == "attachment_valence_angle"]
        .set_index(["subset", "method", "statistic"])
        .value
    )
    for j, (kind, win) in enumerate([("N", (110, 130)), ("O", (104, 124))]):
        ax = f.add_axes([0.215 + j * 0.405, 0.22, 0.33, 0.55])
        ax.axvspan(*win, color="#e9efdf", lw=0, zorder=0)
        for i, (m, name) in enumerate([("exp", "Experiment"), *METHODS.items()]):
            values = [data.loc[kind, name, f"site percentile {v}"] for v in [10, 25, 50, 75, 90]]
            color = INK if m == "exp" else COLORS[m]
            ax.hlines(i, values[0], values[4], color=color, lw=1.05, zorder=2)
            ax.hlines(i, values[1], values[3], color=color, lw=3.15, zorder=3)
            ax.scatter(values[2], i, color=color, s=17, zorder=4)
        ax.set_title(f"{kind}-linked", fontsize=10.4, weight="bold", pad=5)
        ax.set(
            xlim=(78, 181),
            xticks=[90, 120, 150, 180],
            ylim=(6.55, -0.60),
            yticks=range(7),
            yticklabels=["Experiment", *METHODS.values()] if not j else [],
        )
        style(ax, "x")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    f.text(0.57, 0.083, "Attachment bond angle (°)", ha="center", fontsize=10)
    return f


def panel_d(stats, raw):
    f = panel("D", "Attachment geometry and glycan placement")
    means = (
        stats[
            (stats.subset == "Asn–GlcNAc disaccharides")
            & (stats.statistic == "sequence-weighted mean")
        ]
        .set_index(["method", "metric"])
        .value
    )
    for j, (metric, label) in enumerate(
        [("root_angle_abs", "Bond-angle error (°)"), ("sidechain_local", "Side-chain RMSD (Å)")]
    ):
        ax = f.add_axes([0.156 + 0.421 * j, 0.28, 0.35, 0.49])
        for m, marker in [("sweetfold", "^"), ("af3", "o")]:
            frame = raw[raw.method == METHODS[m]]
            if m == "af3":
                ax.scatter(
                    frame[metric],
                    frame.glycan_local,
                    marker=marker,
                    s=8.2,
                    facecolor=COLORS[m],
                    edgecolors="none",
                    alpha=0.45,
                    zorder=3,
                )
            else:
                ax.scatter(
                    frame[metric],
                    frame.glycan_local,
                    marker=marker,
                    s=10.5,
                    facecolors="none",
                    edgecolors=COLORS[m],
                    linewidths=0.62,
                    alpha=0.62,
                    zorder=2,
                )
        for m in D_METHODS:
            ax.scatter(
                means.loc[METHODS[m], metric],
                means.loc[METHODS[m], "glycan_local"],
                marker="D",
                s=42,
                facecolor=COLORS[m],
                edgecolor="white",
                linewidth=0.85,
                zorder=5,
            )
        ax.set_yscale("log")
        ax.set(
            ylim=(0.32, 55),
            yticks=[0.5, 1, 3, 10, 30],
            yticklabels=["0.5", "1", "3", "10", "30"] if not j else [],
            xlabel=label,
        )
        if j:
            ax.set_xscale("log")
            ax.set(
                xlim=(0.15, 48),
                xticks=[0.3, 1, 3, 10, 30],
                xticklabels=["0.3", "1", "3", "10", "30"],
            )
            ax.xaxis.set_minor_locator(NullLocator())
        else:
            ax.set(xlim=(-1.5, 63), xticks=[0, 20, 40, 60], ylabel="Glycan RMSD (Å)")
        ax.yaxis.set_minor_locator(NullLocator())
        style(ax, "y")
        ax.tick_params(labelsize=8)
        ax.xaxis.label.set_size(8.8)
    handles = [
        Line2D([], [], ls="none", marker="o", color=COLORS["af3"], ms=5, label=METHODS["af3"]),
        Line2D(
            [],
            [],
            ls="none",
            marker="^",
            markerfacecolor="none",
            markeredgecolor=COLORS["sweetfold"],
            ms=5,
            label=METHODS["sweetfold"],
        ),
        Line2D([], [], ls="none", marker="D", color=INK, ms=4.7, label="Mean"),
    ]
    f.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(0.544, 0.081),
        ncol=3,
        handlelength=0.9,
        handletextpad=0.4,
        columnspacing=1.25,
        fontsize=8.7,
    )
    return f


def torsions(tag, table, summary):
    link, title = next((key, title) for key, (letter, title) in LINKS.items() if letter == tag)
    medians = summary[summary.statistic == "median"].set_index(["subset", "method", "metric"])
    f = panel(tag, title, size=SUPP_SIZE)
    ax = f.add_axes([0.228, 0.191, 0.69, 0.615])
    for i, (m, public) in enumerate(METHODS.items()):
        y = 5 - i
        ax.hlines(y, 0, 180, color="#c4ced1", lw=0.65)
        for angle, sign in [("phi", 1), ("psi", -1)]:
            df = table.query("link_type == @link and method == @m and angle == @angle").sort_values(
                "left_degrees"
            )
            row = medians.loc[link, public, angle + "_absolute_periodic_error"]
            n = int(row.n)
            assert len(df) == 18 and df["count"].sum() == n and (df.valid_edges == n).all()
            center = (df.left_degrees.to_numpy() + df.right_degrees.to_numpy()) / 2
            freq = df["count"].to_numpy() / n
            x = np.r_[0, center, 180]
            ax.fill_between(
                x,
                y,
                y + sign * 0.65 * np.r_[freq[0], freq, freq[-1]],
                color=ANGLE_COLORS[angle],
                alpha=0.83,
                lw=0,
                zorder=2,
            )
            ax.scatter(
                row.value,
                y + sign * 0.09,
                s=17,
                facecolors="white",
                edgecolors=ANGLE_COLORS[angle],
                linewidths=0.9,
                zorder=3,
            )
    ax.set(
        xlim=(0, 180),
        ylim=(-0.62, 5.62),
        xticks=[0, 45, 90, 135, 180],
        yticks=range(6),
        yticklabels=list(METHODS.values())[::-1],
        xlabel="Absolute circular error (°)",
    )
    style(ax, "x")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    f.legend(
        handles=[
            Line2D([], [], color=ANGLE_COLORS["phi"], lw=2.5, label="φ error"),
            Line2D([], [], color=ANGLE_COLORS["psi"], lw=2.5, label="ψ error"),
            Line2D(
                [],
                [],
                ls="none",
                marker="o",
                markerfacecolor="white",
                markeredgecolor=INK,
                ms=4,
                label="Median",
            ),
        ],
        loc="center",
        bbox_to_anchor=(0.56, 0.086),
        ncol=3,
        handlelength=1.2,
        columnspacing=1.35,
        fontsize=8.9,
    )
    return f
