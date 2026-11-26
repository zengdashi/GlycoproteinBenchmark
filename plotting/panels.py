"""Load the final plot values and render a statistical panel from Figures 2–5 or S1."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from . import figure2, figure3, figure4, figure5
from .style import configure_style


DATA = Path(__file__).resolve().parent / "data"
BENCHMARK = DATA.parents[1] / "data" / "glycoproteins.json"


def table(stem):
    return pd.read_csv(DATA / f"{stem}.tsv", sep="\t", keep_default_na=False,
                       dtype={"pdb_id": str, "protein_chain": str})


def with_sequences(observations):
    proteins = pd.DataFrame(
        [(r["pdb_id"], r["protein_chain"], r["sequence"], "B" in r["test_sets"])
         for r in json.loads(BENCHMARK.read_text())],
        columns=["pdb_id", "protein_chain", "sequence", "in_set_b"],
    )
    joined = observations.merge(proteins, on=["pdb_id", "protein_chain"], how="left",
                                validate="many_to_one", sort=False)
    if joined.sequence.isna().any():
        raise ValueError("A plotted protein chain is missing from the benchmark.")
    return joined


# The shared structure table contains one top-ranked prediction per method and chain.
RENDERERS = {
    "2A": (figure2.backbone, "fig2_structure"),
    "2B": (figure2.frames, "fig2_structure"),
    "2C": (figure2.assoc_panel, "fig2_structure"),
    "3A": (figure3.panel_a, "fig3A_bond_lengths"),
    "3B": (figure3.panel_b, "fig3B_stereochemistry"),
    "3C": (figure3.panel_c, "fig3C_attachment_angles"),
    "3D": (figure3.panel_d, "fig3D_disaccharides"),
    "4A": (figure4.a, "fig4A_root_bonds"),
    "4B": (figure4.b, "fig4B_joint_outcomes"),
    "4C": (figure4.c, "fig4C_defects"),
    "5A": (figure5.a, "fig5A_glycan_length"),
    "5B": (figure5.b, "fig5B_score_quartiles"),
    "5C": (figure5.c, "fig5C_candidate_choice"),
    "5D": (figure5.d, "fig5D_correlations"),
}
PANELS = (*RENDERERS, "S1A", "S1B")


def render_panel(key):
    plt.rcdefaults()
    if key.startswith("2"):
        figure2.configure_style()
    else:
        configure_style()
    if key in ("S1A", "S1B"):
        return figure3.torsions(key[-1], table("figS1_torsions"))
    render, source = RENDERERS[key]
    data = table(source)
    if key.startswith("2") or key == "3D":
        data = with_sequences(data)
    if key == "5D":
        return render(data, table("fig5D_AlphaFold3"))
    return render(data)
