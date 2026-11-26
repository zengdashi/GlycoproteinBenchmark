"""Render one statistical panel or every statistical panel in the paper."""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
from scripts.test_set_panel import render

PANELS = ("1B", "2A", "2B", "2C", "3A", "3B", "3C", "3D", "4A", "4B", "4C",
          "5A", "5B", "5C", "5D", "S1A", "S1B")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel", type=str.upper, choices=(*PANELS, "ALL"))
    parser.add_argument("--out-dir", type=Path, default=ROOT / "build/figures")
    args = parser.parse_args()
    output = args.out_dir.resolve()
    if args.panel in ("1B", "ALL"):
        render(output)
    if args.panel != "1B":
        choice = "all" if args.panel == "ALL" else args.panel
        subprocess.run([sys.executable, "-B", str(ROOT / "plotting/plot.py"), choice,
                        "--out-dir", str(output)], check=True)


if __name__ == "__main__":
    main()
