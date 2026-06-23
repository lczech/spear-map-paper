#!/usr/bin/env python3
"""Per-base score vs read length for each library × mutation grid parameter set.

x-axis: mutated read length (bp)
y-axis: mean score / read_len  (per-base score, linear scale)

Scores are only comparable within the same library (different scoring schemes),
so one PNG per library is produced matching plot_scores.py's LIBRARY_GROUPS.
Within each panel: one line per aligner variant.

For edlib (edit distance): per-base score = edit_distance / read_len (lower is better).
For all others (alignment score): per-base score = score / read_len (higher is better).

Usage (from eval-align-libs/):
    python scripts/plot_scores_by_len.py [--benchmarks benchmarks]
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plot_utils import (
    ALIGNER_STYLE, setup_style, grid_label,
    make_grid_fig, hide_unused, legend_below, savefig,
)
from plot_scores import LIBRARY_GROUPS


def compute_mean_score_by_len(df):
    """Count-weighted mean score per (grid_idx, aligner, read_len)."""
    def _wmean(g):
        total = g["count"].sum()
        if total == 0:
            return float("nan")
        return (g["score"] * g["count"]).sum() / total

    return (
        df.groupby(["grid_idx", "aligner", "read_len"])
        .apply(_wmean, include_groups=False)
        .rename("mean_score")
        .reset_index()
    )


def plot_library(mean_df, df_s, grid_ids, lib_name, group):
    setup_style()
    fig, axes, nrows, ncols = make_grid_fig(len(grid_ids))

    handles, labels = [], []

    for idx, gid in enumerate(grid_ids):
        ax       = axes[idx // ncols, idx % ncols]
        sub      = mean_df[mean_df["grid_idx"] == gid]
        grid_row = df_s[df_s["grid_idx"] == gid].iloc[0]

        for aligner in group["aligners"]:
            style = ALIGNER_STYLE[aligner]
            adf   = sub[sub["aligner"] == aligner].sort_values("read_len").dropna(subset=["mean_score"])
            if adf.empty:
                continue

            per_base = adf["mean_score"] / adf["read_len"]

            line, = ax.plot(
                adf["read_len"], per_base,
                color=style["color"], ls=style["ls"], lw=1.2,
                marker="o", markersize=2, markeredgewidth=0,
                label=style["label"],
            )
            if idx == 0:
                handles.append(line)
                labels.append(style["label"])

        ax.set_title(grid_label(grid_row), fontsize=8)
        ax.set_xlabel("read length (bp)")
        ax.set_ylabel(f"score / read_len  ({group['direction']})")

    hide_unused(axes, len(grid_ids), nrows, ncols)
    fig.suptitle(f"Per-base score vs read length — {lib_name}", y=1.01, fontsize=11)
    legend_below(fig, handles, labels)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--benchmarks", default="benchmarks",
                    help="path to benchmarks directory (default: benchmarks)")
    args = ap.parse_args()

    bench   = Path(args.benchmarks)
    figures = bench.parent / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    df   = pd.read_csv(bench / "scores_by_len.csv")
    df_s = pd.read_csv(bench / "summary_by_len.csv")
    df   = df[df["read_len"].between(30, 150)]
    grid_ids = sorted(df["grid_idx"].unique())

    mean_df = compute_mean_score_by_len(df)

    for lib_name, group in LIBRARY_GROUPS.items():
        fig = plot_library(mean_df, df_s, grid_ids, lib_name, group)
        out = figures / f"scores_by_len_{lib_name}.png"
        savefig(fig, out)


if __name__ == "__main__":
    main()
