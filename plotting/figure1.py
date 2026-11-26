"""Figure 1B. Composition of the two benchmark test sets."""
import json
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt

DATA = Path(__file__).resolve().parent.parent / "data" / "glycoproteins.json"


def topology(site):
    if len(site["sugars"]) == 1:
        return "Monosaccharide"
    degree = Counter()
    for edge in site["bonds"]:
        degree[edge["sugar_1"]] += 1
        degree[edge["sugar_2"]] += 1
    if degree[site["root_sugar"]] > 1 or max(degree.values()) > 2:
        return "Branched"
    return "Linear"


def test_set_panel():
    values = {"A": Counter(), "B": Counter()}
    for chain in json.loads(DATA.read_text(encoding="utf-8")):
        composition = Counter(topology(site) for site in chain["glycosylation_sites"])
        for name in chain["test_sets"]:
            values[name]["Inputs"] += 1
            values[name]["Sites"] += len(chain["glycosylation_sites"])
            values[name].update(composition)
    labels = ["Inputs", "Sites", "Monosaccharide", "Linear", "Branched"]
    data = [[name, values["A"][name], values["B"][name]] for name in labels]
    fig, axis = plt.subplots(figsize=(3.25, 2.6), layout="constrained")
    axis.axis("off")
    axis.set_title("b   Test sets", loc="left", fontweight="bold", fontsize=12)
    table = axis.table(cellText=data, colLabels=["", "Test set A", "Test set B"],
                       colWidths=[0.48, 0.26, 0.26], cellLoc="center", bbox=(0, 0.04, 1, 0.94))
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_facecolor("#d9e9e8" if row == 0 else "#f1f5f5" if row % 2 == 0 else "white")
        if col == 0:
            cell.get_text().set_ha("left")
    return fig
