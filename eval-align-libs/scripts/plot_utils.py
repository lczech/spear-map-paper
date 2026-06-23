"""Shared constants and helpers for eval-align-libs plotting scripts."""

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Output format for all figures. Change to "pdf" or "svg" for vector output.
OUTPUT_FORMAT = "png"


# Color groups: edlib=red, ksw2=orange, ps-score-custom=blue, ps-score-dnafull=navy,
#               ps-cigar-custom=green, ps-cigar-dnafull=darkgreen, ps-cigar-damage=teal,
#               wfa2-score=dark purple, wfa2-cigar=medium purple.
# Linestyle: solid=exact, dashed=heuristic (for wfa2); solid=hot, dashed=cold (for parasail);
#            solid=cigar, dashed=score (for ksw2).
REDUCED_ALIGNERS = [
    "edlib",
    "parasail-cigar-custom-hot",
    "parasail-cigar-dnafull-hot",
    "parasail-cigar-damage-hot",
    "wfa2-cigar-exact",
    "wfa2-cigar-exact-rescaled",
    "wfa2-cigar-damage",
    "wfa2-cigar-damage-rescaled",
]

ALIGNER_STYLE = {
    "edlib": {
        "color": "#c0392b", "ls": "-", "label": "edlib",
    },
    "ksw2-score": {
        "color": "#e67e22", "ls": "--", "label": "ksw2-score",
    },
    "ksw2-cigar": {
        "color": "#e67e22", "ls": "-", "label": "ksw2-cigar",
    },
    "parasail-score-custom-hot": {
        "color": "#2980b9", "ls": "-", "label": "ps-score-custom (hot)",
    },
    "parasail-score-custom-cold": {
        "color": "#2980b9", "ls": "--", "label": "ps-score-custom (cold)",
    },
    "parasail-score-dnafull-hot": {
        "color": "#1a5276", "ls": "-", "label": "ps-score-dnafull (hot)",
    },
    "parasail-score-dnafull-cold": {
        "color": "#1a5276", "ls": "--", "label": "ps-score-dnafull (cold)",
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
    "parasail-cigar-damage-hot": {
        "color": "#16a085", "ls": "-",  "label": "ps-cigar-damage (hot)",
    },
    "wfa2-cigar-damage": {
        "color": "#884ea0", "ls": "-",  "label": "wfa2-cigar-damage",
    },
    "wfa2-cigar-damage-rescaled": {
        "color": "#6c3483", "ls": "-",  "label": "wfa2-cigar-damage-rescaled",
    },
    "wfa2-cigar-exact-rescaled": {
        "color": "#c39bd3", "ls": "-",  "label": "wfa2-cigar-exact-rescaled",
    },
    "wfa2-score-exact": {
        "color": "#7d3c98", "ls": "-",  "label": "wfa2-score-exact",
    },
    "wfa2-score-heuristic": {
        "color": "#7d3c98", "ls": "--", "label": "wfa2-score-heuristic",
    },
    "wfa2-cigar-exact": {
        "color": "#a569bd", "ls": "-",  "label": "wfa2-cigar-exact",
    },
    "wfa2-cigar-heuristic": {
        "color": "#a569bd", "ls": "--", "label": "wfa2-cigar-heuristic",
    },
}


def setup_style():
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("seaborn-whitegrid")


def grid_label(row):
    return (
        f"sub={row.sub_rate:.3f}  indel={row.indel_rate:.3f}  dmg={row.damage_rate:.3f}"
    )


def make_grid_fig(n_panels, ncols=3, panel_w=4.0, panel_h=3.5):
    """Create a figure with ceil(n_panels/ncols) rows × ncols subplot grid."""
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


def savefig(fig, path):
    """Save fig with OUTPUT_FORMAT, close it, and print the path."""
    out = Path(path).with_suffix(f".{OUTPUT_FORMAT}")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Written: {out}")


def make_reduced_style():
    """Build the REDUCED_ALIGNERS style dict with '(hot)' stripped from labels."""
    style = {k: dict(ALIGNER_STYLE[k]) for k in REDUCED_ALIGNERS}
    for v in style.values():
        v["label"] = v["label"].replace(" (hot)", "").replace(" (cold)", "")
    return style


def legend_below(fig, handles, labels, ncol=None):
    # Auto-pick ncol: at most 2 rows, capped at 5 columns.
    if ncol is None:
        ncol = min(5, math.ceil(len(labels) / 2))
    # Reserve bottom space for the legend. Compute as absolute inches so the
    # margin stays consistent regardless of figure height.
    n_rows = math.ceil(len(labels) / ncol)
    legend_height_in = n_rows * 9 * 1.8 / 72 + 0.25  # font rows + padding
    bottom_margin = legend_height_in / fig.get_figheight()
    fig.legend(
        handles, labels,
        loc="lower center",
        ncol=ncol,
        bbox_to_anchor=(0.5, 0),
        fontsize=9,
        frameon=True,
    )
    fig.tight_layout(rect=[0, bottom_margin, 1, 1])
