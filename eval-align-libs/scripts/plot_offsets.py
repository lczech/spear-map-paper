#!/usr/bin/env python3
"""Start and end offset distributions for each aligner × mutation grid parameter set.

Produces two combined PNGs (offsets_start.png, offsets_end.png) plus one pair of
per-aligner PNGs (offsets_start_<aligner>.png, offsets_end_<aligner>.png) for
detailed inspection without the clutter of overlapping lines.

Usage (from eval-align-libs/):
    python scripts/plot_offsets.py [--benchmarks benchmarks]
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


def plot_one_offset(df, count_col, title, figures, grid_ids, df_s, style_dict=None, outer_labels=False):
    if style_dict is None:
        style_dict = ALIGNER_STYLE
    setup_style()
    fig, axes, nrows, ncols = make_grid_fig(len(grid_ids))

    handles, labels = [], []

    for idx, gid in enumerate(grid_ids):
        ax       = axes[idx // ncols, idx % ncols]
        sub      = df[df["grid_idx"] == gid]
        grid_row = df_s[df_s["grid_idx"] == gid].iloc[0]

        for aligner, style in style_dict.items():
            adf = sub[sub["aligner"] == aligner].sort_values("offset_bp")
            if adf.empty or adf[count_col].sum() == 0:
                continue

            # Fill missing integer offsets with NaN so the line breaks at gaps
            # rather than drawing a false slope across them. NaN is required
            # instead of 0 because log scale silently skips zero values.
            counts = adf.set_index("offset_bp")[count_col].astype(float)
            full_range = range(int(counts.index.min()), int(counts.index.max()) + 1)
            counts = counts.reindex(full_range, fill_value=0.5)

            line, = ax.plot(
                counts.index, counts.values,
                color=style["color"], ls=style["ls"], lw=1.2,
                marker="o", markersize=2, markeredgewidth=0,
                label=style["label"],
            )
            if idx == 0:
                handles.append(line)
                labels.append(style["label"])

        ax.axvline(0, color="black", lw=0.6, ls="-", alpha=0.4)
        ax.set_yscale("log")
        ax.set_title(grid_label(grid_row), fontsize=8)
        if not outer_labels or idx // ncols == nrows - 1:
            ax.set_xlabel("offset (bp)")
        if not outer_labels or idx % ncols == 0:
            ax.set_ylabel("count")

    hide_unused(axes, len(grid_ids), nrows, ncols)
    fig.suptitle(title, y=1.01, fontsize=11)
    legend_below(fig, handles, labels)
    return fig


def plot_per_aligner_offset(df, count_col, title_prefix, grid_ids, df_s, aligner, style):
    """One figure for a single aligner: 2×3 grid of mutation-param panels, auto y-scale."""
    setup_style()
    fig, axes, nrows, ncols = make_grid_fig(len(grid_ids))

    for idx, gid in enumerate(grid_ids):
        ax       = axes[idx // ncols, idx % ncols]
        adf      = df[(df["grid_idx"] == gid) & (df["aligner"] == aligner)].sort_values("offset_bp")
        grid_row = df_s[df_s["grid_idx"] == gid].iloc[0]

        if not adf.empty and adf[count_col].sum() > 0:
            counts = adf.set_index("offset_bp")[count_col].astype(float)
            full_range = range(int(counts.index.min()), int(counts.index.max()) + 1)
            counts = counts.reindex(full_range, fill_value=0.5)
            ax.plot(
                counts.index, counts.values,
                color=style["color"], ls=style["ls"], lw=1.2,
                marker="o", markersize=2, markeredgewidth=0,
            )
            ax.axvline(0, color="black", lw=0.6, ls="-", alpha=0.4)
            ax.set_yscale("log")

        ax.set_title(grid_label(grid_row), fontsize=8)
        ax.set_xlabel("offset (bp)")
        ax.set_ylabel("count")

    hide_unused(axes, len(grid_ids), nrows, ncols)
    fig.suptitle(f"{title_prefix} — {style['label']}", y=1.01, fontsize=11)
    fig.tight_layout()
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--benchmarks", default="benchmarks",
                    help="path to benchmarks directory (default: benchmarks)")
    args = ap.parse_args()

    bench   = Path(args.benchmarks)
    figures = bench.parent / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    df     = pd.read_csv(bench / "offsets.csv")
    df_s   = pd.read_csv(bench / "summary.csv")
    grid_ids = sorted(df["grid_idx"].unique())

    for count_col, label, title_prefix, fname_prefix in [
        ("start_count", "Start offset distributions (aligner start − true start)",
         "Start offsets", "offsets_start"),
        ("end_count",   "End offset distributions (aligner end − true end)",
         "End offsets",   "offsets_end"),
    ]:
        # Combined plot — all aligners overlaid
        fig = plot_one_offset(df, count_col, label, figures, grid_ids, df_s)
        savefig(fig, figures / f"{fname_prefix}.png")

        # Combined plot — reduced aligner set
        fig = plot_one_offset(df, count_col, label, figures, grid_ids, df_s, make_reduced_style(), outer_labels=True)
        savefig(fig, figures / f"{fname_prefix}_reduced.png")

        # Per-aligner plots — one file per aligner
        for aligner, style in ALIGNER_STYLE.items():
            if df[df["aligner"] == aligner][count_col].sum() == 0:
                continue
            fig = plot_per_aligner_offset(
                df, count_col, title_prefix, grid_ids, df_s, aligner, style
            )
            safe_name = aligner.replace("/", "_")
            savefig(fig, figures / f"{fname_prefix}_{safe_name}.png")


if __name__ == "__main__":
    main()
