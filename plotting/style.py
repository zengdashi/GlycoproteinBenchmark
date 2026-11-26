"""Shared typography and colours used in the submitted statistical panels."""

import matplotlib.pyplot as plt


METHODS = {
    "af3": "AlphaFold 3",
    "boltz2": "Boltz-2",
    "chai1": "Chai-1",
    "protenix": "Protenix",
    "rfaa": "RFAA",
    "sweetfold": "SweetFold",
}
COLORS = dict(zip(METHODS, ["#237c75", "#4f76a6", "#a66639", "#89649a", "#bf5351", "#6f8050"]))
INK = "#24363b"
GRID = "#dde5e4"
SIZE = (5.1, 3.8)


def panel(letter, title, size=SIZE):
    fig = plt.figure(figsize=size, facecolor="white")
    fig.text(0.025, 0.97, letter, weight="bold", fontsize=16, va="top")
    fig.text(0.104, 0.957, title, weight="bold", fontsize=12.3, va="top")
    return fig


def style(ax, grid=None):
    ax.set_axisbelow(True)
    ax.tick_params(length=2.6, width=0.65, pad=2.4)
    if grid:
        ax.grid(axis=grid, lw=0.55, color=GRID)


def configure_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelcolor": INK,
            "text.color": INK,
            "axes.edgecolor": "#91a2a1",
            "axes.linewidth": 0.7,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 9,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "pdf.fonttype": 42,
        }
    )
