#!/usr/bin/env python3
"""End-offset mean absolute error (MAE) vs read length for each aligner × mutation grid.

x-axis: mutated read length (bp)
y-axis: weighted MAE of end offset (bp, linear scale) — lower is more accurate

One panel per grid parameter combination, one line per aligner.
Uses end offset only (score-only aligners do not report a start position).

Usage (from eval-align-libs/):
    python scripts/plot_accuracy_by_len.py [--benchmarks benchmarks]
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plot_utils import (
    ALIGNER_STYLE, setup_style, grid_label,
    make_grid_fig, hide_unused, legend_below,
)


def compute_mae_by_len(df_off):
    """Weighted MAE of end offset grouped by (grid_idx, aligner, read_len)."""
    def _mae(g):
        total = g["end_count"].sum()
        if total == 0:
            return float("nan")
        return (g["offset_bp"].abs() * g["end_count"]).sum() / total

    return (
        df_off.groupby(["grid_idx", "aligner", "read_len"])
        .apply(_mae, include_groups=False)
        .rename("mae")
        .reset_index()
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--benchmarks", default="benchmarks",
                    help="path to benchmarks directory (default: benchmarks)")
    args = ap.parse_args()

    bench   = Path(args.benchmarks)
    figures = bench.parent / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    df_off = pd.read_csv(bench / "offsets_by_len.csv")
    df_s   = pd.read_csv(bench / "summary_by_len.csv")
    df_off = df_off[df_off["read_len"].between(30, 150)]

    mae_df = compute_mae_by_len(df_off)

    setup_style()
    grid_ids = sorted(mae_df["grid_idx"].unique())
    fig, axes, nrows, ncols = make_grid_fig(len(grid_ids))

    handles, labels = [], []

    for idx, gid in enumerate(grid_ids):
        ax       = axes[idx // ncols, idx % ncols]
        sub      = mae_df[mae_df["grid_idx"] == gid]
        grid_row = df_s[df_s["grid_idx"] == gid].iloc[0]

        for aligner, style in ALIGNER_STYLE.items():
            adf = sub[sub["aligner"] == aligner].sort_values("read_len").dropna(subset=["mae"])
            if adf.empty:
                continue

            line, = ax.plot(
                adf["read_len"], adf["mae"],
                color=style["color"], ls=style["ls"], lw=1.2,
                marker="o", markersize=2, markeredgewidth=0,
                label=style["label"],
            )
            if idx == 0:
                handles.append(line)
                labels.append(style["label"])

        ax.set_title(grid_label(grid_row), fontsize=8)
        ax.set_xlabel("read length (bp)")
        ax.set_ylabel("mean absolute end offset (bp)")

    hide_unused(axes, len(grid_ids), nrows, ncols)
    fig.suptitle("Alignment accuracy vs read length  (lower = better)", y=1.01, fontsize=11)
    legend_below(fig, handles, labels)
    fig.tight_layout()

    out = figures / "accuracy_by_len.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
