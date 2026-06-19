"""Shared constants and helpers for eval-align-libs plotting scripts."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Color groups: edlib=red, ksw2=orange, ps-loc-custom=blue, ps-loc-dnafull=navy,
#               ps-cigar-custom=green, ps-cigar-dnafull=darkgreen.
# Linestyle: solid=hot (or no hot/cold split), dashed=cold.
ALIGNER_STYLE = {
    "edlib": {
        "color": "#c0392b", "ls": "-", "label": "edlib",
    },
    "ksw2": {
        "color": "#e67e22", "ls": "-", "label": "ksw2",
    },
    "parasail-loc-custom-hot": {
        "color": "#2980b9", "ls": "-", "label": "ps-loc-custom (hot)",
    },
    "parasail-loc-custom-cold": {
        "color": "#2980b9", "ls": "--", "label": "ps-loc-custom (cold)",
    },
    "parasail-loc-dnafull-hot": {
        "color": "#1a5276", "ls": "-", "label": "ps-loc-dnafull (hot)",
    },
    "parasail-loc-dnafull-cold": {
        "color": "#1a5276", "ls": "--", "label": "ps-loc-dnafull (cold)",
    },
    "parasail-cigar-custom-hot": {
        "color": "#27ae60", "ls": "-", "label": "ps-cigar-custom (hot)",
    },
    "parasail-cigar-custom-cold": {
        "color": "#27ae60", "ls": "--", "label": "ps-cigar-custom (cold)",
    },
    "parasail-cigar-dnafull-hot": {
        "color": "#1e8449", "ls": "-", "label": "ps-cigar-dnafull (hot)",
    },
    "parasail-cigar-dnafull-cold": {
        "color": "#1e8449", "ls": "--", "label": "ps-cigar-dnafull (cold)",
    },
}


def setup_style():
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("seaborn-whitegrid")


def grid_label(row):
    return (
        f"sub={row.sub_rate:.3f}  indel={row.indel_rate:.3f}\n"
        f"dmg={row.damage_rate:.3f}  λ={row.decay_lambda:.3f}"
    )


def make_grid_fig(n_panels, ncols=3, panel_w=4.0, panel_h=3.5):
    """Create a figure with a 2×ncols subplot grid for n_panels panels."""
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(panel_w * ncols, panel_h * nrows),
        squeeze=False,
    )
    return fig, axes, nrows, ncols


def hide_unused(axes, n_used, nrows, ncols):
    for idx in range(n_used, nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)


def legend_below(fig, handles, labels, ncol=5):
    fig.legend(
        handles, labels,
        loc="lower center",
        ncol=ncol,
        bbox_to_anchor=(0.5, -0.04),
        fontsize=7,
        frameon=True,
    )
