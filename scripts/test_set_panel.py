"""Build Figure 1b from the six-method case lists and rooted glycan graphs."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def topology(glycan):
    if len(glycan["nodes"]) == 1:
        return "Monosaccharide"
    adjacency = defaultdict(set)
    for bond in glycan["bonds"]:
        a, b = bond["node_1"], bond["node_2"]
        adjacency[a].add(b)
        adjacency[b].add(a)
    root = glycan["root_node_id"]
    if len(adjacency[root]) > 1 or any(len(nodes) > 2 for nodes in adjacency.values()):
        return "Branched"
    return "Linear"


def counts():
    with (ROOT / "data/cases.tsv").open() as source:
        cases = {r["case_id"]: r for r in csv.DictReader(source, delimiter="\t") if r["common_six"] == "1"}
    graphs = json.loads((ROOT / "data/glycan_connections.json").read_text())
    by_set = {name: Counter() for name in ("A", "B")}
    for case in graphs["cases"]:
        info = cases.get(case["case_id"])
        if info is None:
            continue
        for name, stats in by_set.items():
            if info["in_" + name] == "1":
                stats["Inputs"] += 1
                stats["Sites"] += len(case["glycans"])
                stats.update(topology(glycan) for glycan in case["glycans"])
    return by_set


def render(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True, exist_ok=True)
    values = counts()
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
    for ext in ("pdf", "png"):
        fig.savefig(out / f"Fig1B.{ext}", dpi=400)
    plt.close(fig)
