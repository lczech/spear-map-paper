#!/usr/bin/env python3
"""Accuracy vs. speed scatter plot for each mutation grid parameter set.

x-axis: mean alignment time (ns, log scale) — lower is faster
y-axis: mean absolute offset (bp, averaged over start and end) — lower is more accurate

One panel per grid parameter combination (2×3 layout), one point per aligner.

Usage (from eval-align-libs/):
    python scripts/plot_accuracy_vs_speed.py [--benchmarks benchmarks]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plot_utils import (
    ALIGNER_STYLE, setup_style, grid_label,
    make_grid_fig, hide_unused, legend_below,
)


def weighted_mae(group, count_col):
    """Mean absolute offset weighted by count, or NaN if all counts are zero."""
    total = group[count_col].sum()
    if total == 0:
        return float("nan")
    return (group["offset_bp"].abs() * group[count_col]).sum() / total


def compute_mae(df_off):
    """Return DataFrame with columns [grid_idx, aligner, mae] (mean of start+end MAE)."""
    start_mae = (
        df_off.groupby(["grid_idx", "aligner"])
        .apply(lambda g: weighted_mae(g, "start_count"), include_groups=False)
        .rename("start_mae")
        .reset_index()
    )
    end_mae = (
        df_off.groupby(["grid_idx", "aligner"])
        .apply(lambda g: weighted_mae(g, "end_count"), include_groups=False)
        .rename("end_mae")
        .reset_index()
    )
    merged = start_mae.merge(end_mae, on=["grid_idx", "aligner"])
    merged["mae"] = (merged["start_mae"] + merged["end_mae"]) / 2
    return merged[["grid_idx", "aligner", "mae"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--benchmarks", default="benchmarks",
                    help="path to benchmarks directory (default: benchmarks)")
    args = ap.parse_args()

    bench   = Path(args.benchmarks)
    figures = bench.parent / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    df_s   = pd.read_csv(bench / "summary.csv")
    df_off = pd.read_csv(bench / "offsets.csv")

    mae_df = compute_mae(df_off)
    merged = df_s[["grid_idx", "aligner", "mean_ns"]].merge(
        mae_df, on=["grid_idx", "aligner"]
    )

    setup_style()
    grid_ids = sorted(merged["grid_idx"].unique())
    fig, axes, nrows, ncols = make_grid_fig(len(grid_ids))

    handles, labels = [], []

    for idx, gid in enumerate(grid_ids):
        ax       = axes[idx // ncols, idx % ncols]
        sub      = merged[merged["grid_idx"] == gid]
        grid_row = df_s[df_s["grid_idx"] == gid].iloc[0]

        for aligner, style in ALIGNER_STYLE.items():
            row = sub[sub["aligner"] == aligner]
            if row.empty:
                continue

            x = row["mean_ns"].values[0]
            y = row["mae"].values[0]

            marker = "o" if style["ls"] == "-" else "^"
            sc = ax.scatter(
                x, y,
                color=style["color"],
                marker=marker,
                s=60,
                zorder=3,
                label=style["label"],
            )
            if idx == 0:
                handles.append(sc)
                labels.append(style["label"])

            ax.annotate(
                style["label"].replace("ps-", ""),
                (x, y), fontsize=5, ha="left", va="bottom",
                xytext=(3, 3), textcoords="offset points",
            )

        ax.set_xscale("log")
        ax.set_title(grid_label(grid_row), fontsize=8)
        ax.set_xlabel("mean time (ns)")
        ax.set_ylabel("mean absolute offset (bp)")

    hide_unused(axes, len(grid_ids), nrows, ncols)
    fig.suptitle("Accuracy vs. speed  (lower-left = better)", y=1.01, fontsize=11)
    legend_below(fig, handles, labels)
    fig.tight_layout()

    out = figures / "accuracy_vs_speed.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
