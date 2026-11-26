"""Figure 2. Backbone accuracy and glycan positioning."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, to_hex
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator
from .style import METHODS, COLORS

N = list(METHODS.values())
C = list(COLORS.values())
INK = "#24363b"
GRID = "#e0e6e5"
OCHRE = "#a65e35"
SIZE = (5.1, 3.8)
GLYCAN_ALIGNMENT = "Glycan heavy-atom alignment"
LOCAL_ALIGNMENT = "Local protein Cα alignment"
FRAMES = [
    ("self", GLYCAN_ALIGNMENT),
    ("root", "Root sugar ring alignment"),
    ("site", LOCAL_ALIGNMENT),
    ("protein", "Whole-protein Cα alignment"),
]
METHOD_LABELS = N


def rmsd_column(frame):
    return "glycan_" + frame + "_frame_macro_RMSD_A"


def panel(letter, title):
    f = plt.figure(figsize=SIZE, facecolor="white")
    f.text(0.025, 0.97, letter, fontsize=16, weight="bold", va="top")
    f.text(0.104, 0.957, title, fontsize=12.3, weight="bold", va="top")
    return f


def light(color, fraction=0.65):
    c = np.array(to_rgb(color))
    return to_hex(c + fraction * (1 - c))


def tidy(ax, grid_axis=None):
    ax.set_axisbelow(True)
    ax.tick_params(length=2.6, width=0.65, pad=2.4)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID, lw=0.55)


def backbone(data):
    f = panel("A", "Protein backbone accuracy")
    # Each distinct protein sequence receives the same total weight.
    metrics = ["protein_CA_RMSD_A", "protein_CA_lDDT"]
    means = {}
    for cohort, subset in (("A", data), ("B", data[data.in_set_b])):
        by_sequence = subset.groupby(["method", "sequence"])[metrics].mean()
        means[cohort] = by_sequence.groupby("method").mean().reindex(N)
    for index, (metric, label, x, ticks) in enumerate(
        [
            ("protein_CA_RMSD_A", "Cα RMSD (Å)", 0.20, [4, 8, 12]),
            ("protein_CA_lDDT", "Cα lDDT", 0.63, [0.75, 0.80, 0.85, 0.90]),
        ]
    ):
        ax = f.add_axes([x, 0.215, 0.31, 0.575])
        for i in range(6):
            ax.plot(
                [means["A"].iloc[i][metric], means["B"].iloc[i][metric]],
                [i - 0.12, i + 0.12],
                color="#aab6b6",
                lw=0.85,
            )
        for key, marker, dy, col in [("A", "o", -0.12, INK), ("B", "D", 0.12, OCHRE)]:
            ax.scatter(
                means[key][metric],
                np.arange(6) + dy,
                color=col,
                edgecolors="white",
                linewidths=0.55,
                marker=marker,
                s=36,
                zorder=3,
            )
        ax.set(
            yticks=np.arange(6),
            yticklabels=METHOD_LABELS if not index else [],
            ylim=(5.5, -0.5),
            xticks=ticks,
            xlabel=label,
        )
        if index:
            ax.set_xlim(0.742, 0.908)
        else:
            ax.set_xlim(3.2, 13.2)
        tidy(ax, "x")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    f.legend(
        handles=[
            Line2D([], [], ls="none", marker=mk, color=co, label=lab, ms=5.4)
            for mk, co, lab in [("o", INK, "Set A"), ("D", OCHRE, "Set B")]
        ],
        loc="center",
        bbox_to_anchor=(0.55, 0.072),
        ncol=2,
        handlelength=0.9,
        columnspacing=1.75,
    )
    return f


def frames(data):
    f = panel("B", "Glycan RMSD by alignment frame")
    by_seq = data.groupby(["method", "sequence"])[[rmsd_column(k) for k, _ in FRAMES]].mean()
    for i, (frame, title) in enumerate(FRAMES):
        col, row = i % 2, i // 2
        ax = f.add_axes([0.175 + col * 0.422, 0.536 - row * 0.363, 0.363, 0.242])
        ax.set_title(title, loc="left", fontsize=8.5, fontweight="bold", pad=3)
        for j, (method, color) in enumerate(zip(N, C)):
            vals = by_seq.loc[method, rmsd_column(frame)].to_numpy()
            q5, q25, med, q75, q95 = np.quantile(vals, [0.05, 0.25, 0.50, 0.75, 0.95])
            stat = dict(
                whislo=q5,
                q1=q25,
                med=med,
                q3=q75,
                whishi=q95,
                fliers=vals[(vals < q5) | (vals > q95)],
                mean=np.mean(vals),
            )
            ax.bxp(
                [stat],
                positions=[j],
                vert=False,
                widths=0.58,
                showmeans=True,
                patch_artist=True,
                manage_ticks=False,
                boxprops=dict(facecolor=light(color, 0.55), edgecolor=color, linewidth=0.8),
                medianprops=dict(color=INK, linewidth=1.05),
                whiskerprops=dict(color=color, linewidth=0.8),
                capprops=dict(color=color, linewidth=0.8),
                flierprops=dict(
                    marker=".",
                    markersize=1.6,
                    markerfacecolor=light(color, 0.5),
                    markeredgewidth=0,
                    linestyle="none",
                ),
                meanprops=dict(
                    marker="D",
                    markersize=3.0,
                    markeredgecolor="white",
                    markeredgewidth=0.38,
                    markerfacecolor=INK,
                ),
            )
        ax.set_xscale("log")
        ax.set(
            xlim=(0.043, 140),
            ylim=(5.6, -0.6),
            yticks=range(6),
            yticklabels=METHOD_LABELS if col == 0 else [],
            xticks=[0.1, 1, 10, 100],
            xticklabels=["0.1", "1", "10", "100"] if row else [],
        )
        ax.xaxis.set_minor_locator(NullLocator())
        tidy(ax, "x")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", labelsize=7.8, length=0)
        ax.tick_params(axis="x", labelsize=8)
    f.text(0.56, 0.038, "Glycan RMSD (Å)", ha="center", fontsize=10)
    return f


def weighted_midrank(values, weights):
    _, inverse = np.unique(values, return_inverse=True)
    mass = np.bincount(inverse, weights=weights)
    return (mass.cumsum() - 0.5 * mass)[inverse] / weights.sum()


def association(data):
    x = data.protein_CA_lDDT.to_numpy()
    y = data[rmsd_column("site")].to_numpy()
    weights = 1.0 / data.groupby("sequence").protein_CA_lDDT.transform("size").to_numpy()
    root_weights = np.sqrt(weights)
    design = np.column_stack([np.ones(len(x)), x]) * root_weights[:, None]
    intercept, slope = np.linalg.lstsq(design, np.log(y) * root_weights, rcond=None)[0]
    order = np.argsort(x, kind="stable")
    bounds = np.searchsorted(weights[order].cumsum(), np.array([0.025, 0.975]) * weights.sum())
    x_range = x[order][bounds]

    ranks = np.vstack([weighted_midrank(x, weights), weighted_midrank(y, weights)])
    covariance = np.cov(ranks, aweights=weights, ddof=0)
    rho = covariance[0, 1] / np.sqrt(covariance[0, 0] * covariance[1, 1])
    return x_range, np.exp(intercept + x_range * slope), rho


def assoc_panel(data):
    f = panel("C", "Backbone accuracy and glycan placement")
    for i, (name, color) in enumerate(zip(N, C)):
        col, row = i % 3, i // 3
        ax = f.add_axes([0.139 + col * 0.283, 0.537 - row * 0.341, 0.241, 0.233])
        q = data[data.method == name]
        xx, yy, rho = association(q)
        ax.scatter(
            q.protein_CA_lDDT,
            q[rmsd_column("site")],
            s=5.7,
            color=color,
            alpha=0.39,
            edgecolors="none",
            zorder=1,
        )
        ax.axhline(5, color="#68777b", ls=(0, (1.2, 2.3)), lw=0.65, zorder=0)
        ax.axvline(0.9, color="#68777b", ls=(0, (1.2, 2.3)), lw=0.65, zorder=0)
        ax.plot(xx, yy, color="white", lw=2.4, zorder=2, solid_capstyle="round")
        ax.plot(xx, yy, color=INK, lw=1.4, zorder=3, solid_capstyle="round")
        ax.set_yscale("log")
        ax.set(
            xlim=(0.14, 1.01),
            ylim=(0.34, 51),
            xticks=[0.3, 0.6, 0.9],
            xticklabels=["0.3", "0.6", "0.9"] if row else [],
            yticks=[1, 5, 25],
            yticklabels=["1", "5", "25"] if col == 0 else [],
        )
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_title(METHOD_LABELS[i], fontsize=8.8, weight="bold", color=INK, pad=13.2)
        coefficient = f"ρ = {rho:.2f}".replace("-", "−")
        ax.text(
            0.5,
            1.025,
            coefficient,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            family="DejaVu Serif",
            weight="normal",
            fontsize=7.4,
            color="#536569",
        )
        tidy(ax)
        ax.tick_params(labelsize=7.8)
    f.text(0.533, 0.065, "Protein Cα lDDT", ha="center", fontsize=10)
    f.text(
        0.025,
        0.47,
        "Glycan RMSD (Å) after local protein Cα alignment",
        rotation=90,
        ha="center",
        va="center",
        fontsize=8.3,
    )
    return f


def configure_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "axes.edgecolor": "#96a2a2",
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.4,
            "ytick.labelsize": 9.5,
            "figure.facecolor": "white",
            "legend.fontsize": 9.7,
        }
    )
