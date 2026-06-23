#!/usr/bin/env python3
"""Timing distribution histograms for each aligner × mutation grid parameter set.

Usage (from eval-align-libs/):
    python scripts/plot_timing.py [--benchmarks benchmarks]
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plot_utils import (
    ALIGNER_STYLE, setup_style, grid_label,
    make_grid_fig, hide_unused, legend_below, make_reduced_style, savefig,
)


def bucket_midpoint_us(lo: float, hi: float) -> float:
    """Geometric midpoint of a timing bucket, converted to µs."""
    if lo == 0:
        return hi / 2 / 1000       # [0, 1 µs) → midpoint at 0.5 µs
    return float(np.sqrt(lo * hi)) / 1000


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--benchmarks", default="benchmarks",
                    help="path to benchmarks directory (default: benchmarks)")
    args = ap.parse_args()

    bench   = Path(args.benchmarks)
    figures = bench.parent / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    df_t = pd.read_csv(bench / "timing.csv")
    df_s = pd.read_csv(bench / "summary.csv")

    df_t["midpoint_us"] = df_t.apply(
        lambda r: bucket_midpoint_us(r.bucket_low_ns, r.bucket_high_ns), axis=1
    )

    # (grid_idx, aligner) → mean_ns for the vertical mean marker
    mean_ns = df_s.set_index(["grid_idx", "aligner"])["mean_ns"]

    def make_plot(style_dict, out_path, outer_labels=False):
        setup_style()
        grid_ids = sorted(df_t["grid_idx"].unique())
        fig, axes, nrows, ncols = make_grid_fig(len(grid_ids))

        handles, labels = [], []

        for idx, gid in enumerate(grid_ids):
            ax       = axes[idx // ncols, idx % ncols]
            sub      = df_t[df_t["grid_idx"] == gid]
            grid_row = df_s[df_s["grid_idx"] == gid].iloc[0]

            for aligner, style in style_dict.items():
                adf = sub[sub["aligner"] == aligner].sort_values("midpoint_us")
                if adf.empty:
                    continue

                line, = ax.plot(
                    adf["midpoint_us"], adf["count"],
                    color=style["color"], ls=style["ls"], lw=1.5,
                    label=style["label"],
                )
                if idx == 0:
                    handles.append(line)
                    labels.append(style["label"])

                # Dotted vertical line at the mean
                key = (gid, aligner)
                if key in mean_ns.index:
                    ax.axvline(mean_ns[key] / 1000, color=style["color"],
                               ls=":", lw=0.8, alpha=0.5)

            ax.set_xscale("log")
            ax.set_title(grid_label(grid_row), fontsize=8)
            if not outer_labels or idx // ncols == nrows - 1:
                ax.set_xlabel("time (µs)")
            if not outer_labels or idx % ncols == 0:
                ax.set_ylabel("count")

        hide_unused(axes, len(grid_ids), nrows, ncols)
        fig.suptitle("Alignment timing distributions", y=1.01, fontsize=11)
        legend_below(fig, handles, labels)
        savefig(fig, out_path)

    make_plot(ALIGNER_STYLE, figures / "timing.png")

    make_plot(make_reduced_style(), figures / "timing_reduced.png", outer_labels=True)


if __name__ == "__main__":
    main()
