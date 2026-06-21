#!/usr/bin/env python3
"""Mean alignment time vs read length for each aligner × mutation grid parameter set.

x-axis: mutated read length (bp)
y-axis: mean alignment time (ns, log scale)

One panel per grid parameter combination, one line per aligner.

Usage (from eval-align-libs/):
    python scripts/plot_timing_by_len.py [--benchmarks benchmarks]
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plot_utils import (
    ALIGNER_STYLE, REDUCED_ALIGNERS, setup_style, grid_label,
    make_grid_fig, hide_unused, legend_below,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--benchmarks", default="benchmarks",
                    help="path to benchmarks directory (default: benchmarks)")
    args = ap.parse_args()

    bench   = Path(args.benchmarks)
    figures = bench.parent / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(bench / "summary_by_len.csv")
    df = df[df["read_len"].between(30, 150)]

    def make_plot(style_dict, out_path):
        setup_style()
        grid_ids = sorted(df["grid_idx"].unique())
        fig, axes, nrows, ncols = make_grid_fig(len(grid_ids))

        handles, labels = [], []

        for idx, gid in enumerate(grid_ids):
            ax       = axes[idx // ncols, idx % ncols]
            sub      = df[df["grid_idx"] == gid]
            grid_row = sub.iloc[0]

            for aligner, style in style_dict.items():
                adf = sub[sub["aligner"] == aligner].sort_values("read_len")
                if adf.empty or adf["mean_ns"].sum() == 0:
                    continue

                line, = ax.plot(
                    adf["read_len"], adf["mean_ns"],
                    color=style["color"], ls=style["ls"], lw=1.2,
                    marker="o", markersize=2, markeredgewidth=0,
                    label=style["label"],
                )
                if idx == 0:
                    handles.append(line)
                    labels.append(style["label"])

            ax.set_yscale("log")
            ax.set_title(grid_label(grid_row), fontsize=8)
            ax.set_xlabel("read length (bp)")
            ax.set_ylabel("mean time (ns)")

        hide_unused(axes, len(grid_ids), nrows, ncols)
        fig.suptitle("Alignment time vs read length", y=1.01, fontsize=11)
        legend_below(fig, handles, labels)
        fig.tight_layout()
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        print(f"Written: {out_path}")

    make_plot(ALIGNER_STYLE, figures / "timing_by_len.png")

    reduced_style = {k: ALIGNER_STYLE[k] for k in REDUCED_ALIGNERS}
    make_plot(reduced_style, figures / "timing_by_len_reduced.png")


if __name__ == "__main__":
    main()
