"""Figure 5. Chain length, confidence and candidate agreement."""

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator
from style import METHODS, COLORS, GRID, INK, panel

LABELS = list(METHODS.values())
C = dict(zip(LABELS, COLORS.values()))
TITLES = {
    "A": "Glycan shape by chain length",
    "B": "Placement across score quartiles",
    "C": "Selection from five candidates",
    "D": "Candidate spread and placement",
}
MARKERS = ["o", "s", "^", "D", "v", "p"]


def styl(ax, grid=None):
    ax.tick_params(length=2.5, width=0.66, pad=2.5)
    ax.set_axisbelow(True)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.54)


def a(df):
    f = panel("A", TITLES["A"])
    ax = f.add_axes([0.148, 0.328, 0.816, 0.393])
    counts = f.add_axes([0.148, 0.161, 0.816, 0.107])
    lens = df
    ax.set(xlim=(0.5, 16.6), ylim=(0, 8.7), yticks=[0, 2, 4, 6, 8], ylabel="Glycan-fit RMSD (Å)")
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    styl(ax, "y")
    entries = []
    for i, m in enumerate(LABELS):
        d = lens[lens.method.eq(m)].sort_values("stratum")
        common = d[d.n_sites >= 10]
        rare = d[d.n_sites < 10]
        x = common.stratum.to_numpy(dtype=float)
        y = common.value.to_numpy(dtype=float)
        spec = dict(marker=MARKERS[i], ms=4.25, mfc=C[m], mec="white", mew=0.34, color=C[m])
        ax.plot(x, y, lw=1.18, **spec)
        xx = rare.stratum.to_numpy(dtype=float) + (i - 2.5) * 0.083
        ax.plot(xx, rare.value, ls="none", **spec)
        entries.append(Line2D([], [], lw=1.08, **spec, label=m))
    f.legend(
        handles=entries,
        loc="center",
        bbox_to_anchor=(0.54, 0.822),
        ncol=3,
        fontsize=8.4,
        columnspacing=1.19,
        handlelength=1.46,
        handletextpad=0.42,
        labelspacing=0.34,
    )
    d = lens[lens.method.eq("AlphaFold 3")].sort_values("stratum")
    xs = d.stratum.to_numpy(dtype=float)
    ns = d.n_sites.to_numpy(dtype=float)
    counts.bar(xs, ns - 0.65, bottom=0.65, width=0.58, color="#9babae", zorder=2)
    counts.set(
        xlim=ax.get_xlim(),
        ylim=(0.65, 650),
        yscale="log",
        yticks=[1, 100],
        yticklabels=["1", "100"],
        xticks=xs,
        xticklabels=[str(int(v)) for v in xs],
        ylabel="Sites",
    )
    styl(counts)
    counts.tick_params(axis="both", labelsize=7.4, length=2, pad=2)
    counts.yaxis.set_minor_locator(NullLocator())
    counts.yaxis.label.set_size(8.6)
    counts.yaxis.labelpad = 3.8
    f.text(0.55, 0.049, "Observed monosaccharides per site", va="center", ha="center", fontsize=9.4)
    return f


def b(df):
    f = panel("B", TITLES["B"])
    sub = df
    for i, m in enumerate(LABELS):
        data = sub[sub.method.eq(m)].sort_values("stratum")
        g = f.add_axes([0.145 + (i % 3) * 0.282, 0.536 - (i // 3) * 0.372, 0.231, 0.224])
        x = data.stratum.to_numpy(dtype=float)
        y = data.value.to_numpy(dtype=float)
        g.set(
            xlim=(0.65, 4.35),
            ylim=(0, 34.2),
            xticks=[1, 2, 3, 4],
            xticklabels=["Q1", "Q2", "Q3", "Q4"],
            yticks=[0, 15, 30],
        )
        styl(g, "y")
        g.tick_params(axis="x", labelbottom=i >= 3, labelsize=8.0, length=2.4, pad=2)
        g.tick_params(
            axis="y", labelleft=i % 3 == 0, labelsize=7.8, length=2.3 if i % 3 == 0 else 0
        )
        g.set_title(m, fontsize=9.4, color=C[m], weight="bold", pad=2.8)
        g.vlines(x, data.ci_low, data.ci_high, color=C[m], lw=1.1, alpha=0.67)
        g.plot(x, y, color=C[m], marker="o", ms=4.0, mew=0.33, mec="white", mfc=C[m], lw=1.22)

    f.text(
        0.038,
        0.455,
        "Protein-fit glycan RMSD (Å)",
        rotation=90,
        ha="center",
        va="center",
        fontsize=9.3,
    )
    f.text(0.546, 0.047, "Native score quartile", ha="center", va="center", fontsize=9.4)
    return f


def c(df):
    f = panel("C", TITLES["C"])
    sub = df
    ax = f.add_axes([0.222, 0.232, 0.737, 0.590])
    ax.set(
        yticks=range(6),
        yticklabels=LABELS,
        ylim=(5.5, -0.5),
        xlim=(5.9, 19.2),
        xticks=[6, 8, 10, 12, 14, 16, 18],
        xlabel="Protein-fit glycan RMSD (Å)",
    )
    styl(ax, "x")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=4)
    desc = [
        ("native_rmsd", "Native top rank", "o"),
        ("mean5_rmsd", "Random choice", "D"),
        ("best5_rmsd", "Reference best", "^"),
    ]
    for i, m in enumerate(LABELS):
        d = sub[sub.method.eq(m)].set_index("metric")
        sx = {key: d.loc[key, "value"] for key, _, _ in desc}
        ax.plot([min(sx.values()), max(sx.values())], [i, i], color=C[m], lw=1.7, alpha=0.68)
        for dy, (key, name, marker) in zip([-0.121, 0.121, 0], desc):
            x = sx[key]
            ax.plot(
                x,
                i + dy,
                marker=marker,
                mfc="white" if key == "mean5_rmsd" else C[m],
                mec=C[m],
                ls="none",
                mew=1.0 if key == "mean5_rmsd" else 0.6,
                ms=6.1 if key == "best5_rmsd" else 5.2,
            )
    handles = [
        Line2D(
            [],
            [],
            marker=k,
            ls="none",
            ms=5.8,
            mfc="white" if key == "mean5_rmsd" else INK,
            mec=INK,
            label=label,
        )
        for key, label, k in desc
    ]
    f.legend(
        handles=handles,
        ncol=3,
        loc="center",
        bbox_to_anchor=(0.53, 0.054),
        fontsize=8.4,
        columnspacing=0.91,
        handletextpad=0.3,
        handlelength=0.85,
    )
    return f


def d(df, inputs):
    f = panel("D", TITLES["D"])
    ax = f.add_axes([0.143, 0.237, 0.344, 0.504])
    corr = f.add_axes([0.705, 0.237, 0.249, 0.504])
    scatter = inputs[inputs.method.eq("AlphaFold 3")].sort_values("case_id")
    assert len(scatter) == 327
    for xy in [("x", (0.285, 73), [0.5, 1, 3, 10, 30]), ("y", (0.36, 115), [1, 3, 10, 30, 100])]:
        axis, limits, ticks = xy
        if axis == "x":
            ax.set_xscale("log")
            ax.set_xlim(limits)
            ax.set_xticks(ticks, [str(t) for t in ticks])
        else:
            ax.set_yscale("log")
            ax.set_ylim(limits)
            ax.set_yticks(ticks, [str(t) for t in ticks])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    styl(ax)
    ax.tick_params(labelsize=8, length=2.4, pad=2)
    ax.spines["bottom"].set_color("#7a9190")
    ax.spines["left"].set_color("#7a9190")
    ax.scatter(
        scatter.dispersion_protein,
        scatter.top1_rmsd,
        s=10.8,
        c=C["AlphaFold 3"],
        alpha=0.43,
        edgecolors="none",
        zorder=3,
    )
    ax.set(xlabel="Candidate spread (Å)", ylabel="Top-ranked glycan RMSD (Å)")
    ax.xaxis.label.set_size(9)
    ax.yaxis.label.set_size(9)
    ax.set_title("AlphaFold 3", fontsize=9.6, weight="bold", color=C["AlphaFold 3"], pad=5)
    corr.set(
        xlim=(-0.025, 0.605),
        xticks=[0, 0.3, 0.6],
        xticklabels=["0", "0.3", "0.6"],
        ylim=(5.55, -0.55),
        yticks=range(6),
        yticklabels=LABELS,
    )
    styl(corr, "x")
    corr.spines["left"].set_visible(False)
    corr.tick_params(axis="y", length=0, pad=2.4, labelsize=8.3)
    corr.tick_params(axis="x", labelsize=8.0, pad=2)
    corr.set_xlabel("Partial Spearman ρ", fontsize=9, labelpad=4.8)
    v = df.set_index("method")
    for i, m in enumerate(LABELS):
        line = v.loc[m]
        corr.hlines(i, line.ci_low, line.ci_high, color=C[m], lw=1.40)
        corr.plot(line.value, i, marker="o", ms=5.3, mfc=C[m], mec="white", mew=0.42)
    return f
