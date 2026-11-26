"""Figure 4. Chemical errors and positioning failures."""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex
from matplotlib.patches import Rectangle, Patch
from style import METHODS, INK, GRID, panel

OUTCOMES = {
    "chem_pass_geo_pass": ("Both pass", "#246d67"),
    "chem_pass_geo_fail": ("Chemistry only", "#accec2"),
    "chem_fail_geo_pass": ("Placement only", "#e8ae6a"),
    "chem_fail_geo_fail": ("Neither passes", "#aa5a56"),
}
METRICS = {
    "fail_root": "Root\nbond",
    "fail_internal": "Internal\nbond",
    "fail_sugar_stereo": "Sugar\nchirality",
    "fail_clash": "Clash",
    "fail_sugar_stereo_and_clash": "Chirality\n+ clash",
}
TITLES = {
    "A": "Out-of-range protein–sugar bonds",
    "B": "Chemistry and glycan placement",
    "C": "Defects among failed sites (%)",
}
HEAT = LinearSegmentedColormap.from_list(
    "defects", ["#f0f3f2", "#d0e3df", "#98c0b9", "#578a82", "#235954"]
)


def methods(ax):
    ax.set_yticks(range(6), METHODS.values())
    ax.set_ylim(5.58, -0.58)
    ax.tick_params(axis="y", length=0, pad=4.5, labelsize=9)
    ax.spines["left"].set_visible(False)


def get(data, method, **query):
    rows = data[data.method == method]
    for k, v in query.items():
        rows = rows[rows[k] == v]
    assert len(rows) == 1
    return rows.iloc[0]["value"]


def a(data):
    f = panel("A", TITLES["A"])
    for j, (test, col) in enumerate([("A", "#28766c"), ("B", "#67799e")]):
        (num,) = data.loc[data.test_set == test, "n"].unique()
        ax = f.add_axes([0.211 + j * 0.386, 0.232, 0.331, 0.573])
        methods(ax)
        if j:
            ax.tick_params(axis="y", labelleft=False)
        ax.set(xlim=(-0.65, 9.15), xticks=[0, 2, 4, 6, 8])
        ax.tick_params(axis="x", length=2.4, pad=4)
        ax.set_axisbelow(True)
        ax.grid(axis="x", lw=0.6, color=GRID)
        for i, m in enumerate(METHODS.values()):
            val = int(get(data, m, test_set=test))
            ax.plot([0, val], [i, i], color=col, lw=1.8, zorder=2)
            ax.plot(
                val,
                i,
                marker="o" if j == 0 else "D",
                ms=5.6 if j == 0 else 5,
                color=col,
                markeredgewidth=0,
                ls="none",
                zorder=3,
            )
            ax.text(val + 0.52, i, str(val), va="center", color=INK, fontsize=8.5)
        ax.set_title(f"Test set {test} · {num} sites", fontsize=9.2, color=INK, pad=7.0)
    f.text(0.568, 0.071, "Sites outside 1.1–1.9 Å (count)", ha="center", fontsize=10)
    return f


def b(data):
    f = panel("B", TITLES["B"])
    ax = f.add_axes([0.211, 0.284, 0.748, 0.525])
    methods(ax)
    left = np.zeros(6)
    for key, (label, col) in OUTCOMES.items():
        val = np.array([get(data, m, metric=key) for m in METHODS.values()])
        ax.barh(range(6), val, height=0.70, left=left, color=col, edgecolor="white", lw=0.45)
        left += val
    assert np.allclose(left, 100)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.tick_params(axis="x", length=2.2, pad=2.5)
    ax.set_xlabel("Sites (%)", labelpad=3)
    handles = [Patch(facecolor=col, edgecolor="none", label=lab) for lab, col in OUTCOMES.values()]
    f.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(0.545, 0.072),
        ncol=2,
        fontsize=8.6,
        handlelength=1.36,
        handleheight=0.74,
        handletextpad=0.42,
        columnspacing=1.35,
        labelspacing=0.34,
        borderaxespad=0,
    )
    return f


def c(data):
    f = panel("C", TITLES["C"])
    ax = f.add_axes([0.211, 0.145, 0.751, 0.652])
    methods(ax)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlim(-0.55, 4.55)
    ax.tick_params(axis="x", length=0, pad=5, labelsize=8.5)
    ax.set_xticks(range(5), METRICS.values())
    for i, m in enumerate(METHODS.values()):
        for j, k in enumerate(METRICS):
            value = get(data, m, metric=k)
            fc = to_hex(HEAT(min(value / 75, 1)))
            ax.add_patch(
                Rectangle((j - 0.487, i - 0.459), 0.974, 0.918, facecolor=fc, edgecolor="none")
            )
            label = "0" if value == 0 else "<1" if value < 1 else str(int(np.floor(value + 0.5)))
            ax.text(
                j,
                i,
                label,
                va="center",
                ha="center",
                fontsize=9.6,
                color="white" if value >= 46 else INK,
                fontweight="bold" if value >= 41 else "normal",
            )
    return f
