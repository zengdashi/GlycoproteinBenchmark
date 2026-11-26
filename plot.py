"""Draw one statistical panel, or all panels in the glycoprotein benchmark."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plotting.figure1 import test_set_panel
from plotting.panels import PANELS, render_panel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    labels = ("1B", *PANELS)
    parser.add_argument("panel", type=str.upper, choices=(*labels, "ALL"))
    parser.add_argument("--out-dir", type=Path, default=Path("build/figures"))
    parser.add_argument("--dpi", type=int, help="PNG resolution (default 400 for 1B, 500 otherwise)")
    args = parser.parse_args()
    if args.dpi is not None and args.dpi <= 0:
        parser.error("--dpi must be positive")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for label in labels if args.panel == "ALL" else (args.panel,):
        if label == "1B":
            plt.rcdefaults()
            figure = test_set_panel()
        else:
            figure = render_panel(label)
        try:
            stem = args.out_dir / f"Fig{label}"
            figure.savefig(stem.with_suffix(".pdf"))
            figure.savefig(stem.with_suffix(".png"), dpi=args.dpi or (400 if label == "1B" else 500))
            print(stem)
        finally:
            plt.close(figure)


if __name__ == "__main__":
    main()
