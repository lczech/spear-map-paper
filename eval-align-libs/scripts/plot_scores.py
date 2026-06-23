#!/usr/bin/env python3
"""Score distributions for each library × mutation grid parameter set.

Produces one PNG per library (scores_edlib.png, scores_ksw2.png,
scores_parasail.png, scores_wfa2.png), each with a 2×3 panel grid.
Within each panel: one line per aligner variant, log y-scale, independent x-axis.

Usage (from eval-align-libs/):
    python scripts/plot_scores.py [--benchmarks benchmarks]
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

LIBRARY_GROUPS = {
    "edlib": {
        "aligners":  ["edlib"],
        "direction": "lower is better",
    },
    "ksw2": {
        "aligners":  ["ksw2-score", "ksw2-cigar"],
        "direction": "higher is better",
    },
    "parasail": {
        "aligners":  [
            "parasail-score-custom-hot",
            "parasail-score-custom-cold",
            "parasail-score-dnafull-hot",
            "parasail-score-dnafull-cold",
            "parasail-cigar-custom-hot",
            "parasail-cigar-custom-cold",
            "parasail-cigar-dnafull-hot",
            "parasail-cigar-dnafull-cold",
            "parasail-cigar-damage-hot",
        ],
        "direction": "higher is better",
    },
    "wfa2": {
        "aligners":  [
            "wfa2-score-exact",
            "wfa2-score-heuristic",
            "wfa2-cigar-exact",
            "wfa2-cigar-heuristic",
        ],
        "direction": "higher is better",
    },
}


def plot_library_scores(df, df_s, grid_ids, lib_name, group):
    setup_style()
    fig, axes, nrows, ncols = make_grid_fig(len(grid_ids))

    handles, labels = [], []

    for idx, gid in enumerate(grid_ids):
        ax       = axes[idx // ncols, idx % ncols]
        sub      = df[df["grid_idx"] == gid]
        grid_row = df_s[df_s["grid_idx"] == gid].iloc[0]

        for aligner in group["aligners"]:
            style = ALIGNER_STYLE[aligner]
            adf   = sub[sub["aligner"] == aligner].sort_values("score")
            if adf.empty or adf["count"].sum() == 0:
                continue

            line, = ax.plot(
                adf["score"], adf["count"],
                color=style["color"], ls=style["ls"], lw=1.2,
                marker="o", markersize=2, markeredgewidth=0,
                label=style["label"],
            )
            if idx == 0:
                handles.append(line)
                labels.append(style["label"])

        ax.set_yscale("log")
        ax.set_title(grid_label(grid_row), fontsize=8)
        ax.set_xlabel(f"score ({group['direction']})")
        ax.set_ylabel("count")

    hide_unused(axes, len(grid_ids), nrows, ncols)
    fig.suptitle(f"Score distributions — {lib_name}", y=1.01, fontsize=11)
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

    df   = pd.read_csv(bench / "scores.csv")
    df_s = pd.read_csv(bench / "summary.csv")
    grid_ids = sorted(df["grid_idx"].unique())

    for lib_name, group in LIBRARY_GROUPS.items():
        fig = plot_library_scores(df, df_s, grid_ids, lib_name, group)
        out = figures / f"scores_{lib_name}.png"
        savefig(fig, out)


if __name__ == "__main__":
    main()
