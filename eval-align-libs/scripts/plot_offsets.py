#!/usr/bin/env python3
"""Start and end offset distributions for each aligner × mutation grid parameter set.

Produces two PDFs: offsets_start.pdf and offsets_end.pdf.
Each has a 2×3 subplot grid (one panel per grid parameter combination).

Usage (from eval-align-libs/):
    python scripts/plot_offsets.py [--benchmarks benchmarks]
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


def plot_one_offset(df, count_col, title, figures, grid_ids, df_s):
    setup_style()
    fig, axes, nrows, ncols = make_grid_fig(len(grid_ids))

    handles, labels = [], []

    for idx, gid in enumerate(grid_ids):
        ax       = axes[idx // ncols, idx % ncols]
        sub      = df[df["grid_idx"] == gid]
        grid_row = df_s[df_s["grid_idx"] == gid].iloc[0]

        for aligner, style in ALIGNER_STYLE.items():
            adf = sub[sub["aligner"] == aligner].sort_values("offset_bp")
            if adf.empty or adf[count_col].sum() == 0:
                continue

            line, = ax.plot(
                adf["offset_bp"], adf[count_col],
                color=style["color"], ls=style["ls"], lw=1.5,
                label=style["label"],
            )
            if idx == 0:
                handles.append(line)
                labels.append(style["label"])

        ax.axvline(0, color="black", lw=0.6, ls="-", alpha=0.4)
        ax.set_title(grid_label(grid_row), fontsize=8)
        ax.set_xlabel("offset (bp)")
        ax.set_ylabel("count")

    hide_unused(axes, len(grid_ids), nrows, ncols)
    fig.suptitle(title, y=1.01, fontsize=11)
    legend_below(fig, handles, labels)
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

    for count_col, label, fname in [
        ("start_count", "Start offset distributions (aligner start − true start)", "offsets_start.pdf"),
        ("end_count",   "End offset distributions (aligner end − true end)",        "offsets_end.pdf"),
    ]:
        fig = plot_one_offset(df, count_col, label, figures, grid_ids, df_s)
        out = figures / fname
        fig.savefig(out, bbox_inches="tight")
        print(f"Written: {out}")


if __name__ == "__main__":
    main()
