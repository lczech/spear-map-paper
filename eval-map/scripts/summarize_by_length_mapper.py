# Build long-format TSV per mapper and metric ("mapped" or "correct"):
# cols: mapper  mm  length  total  <metric>  pct_<metric>

from pathlib import Path
import os
import re
import pandas as pd

len_min = int(snakemake.params.len_min)
len_max = int(snakemake.params.len_max)
metric  = snakemake.params.metric
pct_col = f"pct_{metric}"
lengths = list(range(len_min, len_max + 1))

def mm_from_path(p: str) -> int:
    m = re.search(r"mm(\d+)", os.path.basename(p))
    return int(m.group(1)) if m else None

# Load totals per mm
totals = {}
for p in snakemake.input.totals:
    mm = mm_from_path(p)
    df = pd.read_csv(p, sep="\t", header=None, names=["length", "total"])
    totals[mm] = df

# Load metric counts per mm for this mapper
metric_counts = {}
for p in snakemake.input.metric:
    mm = mm_from_path(p)
    if os.path.exists(p):
        df = pd.read_csv(p, sep="\t", header=None, names=["length", metric])
    else:
        df = pd.DataFrame({"length": lengths, metric: [0] * len(lengths)})
    metric_counts[mm] = df

mapper = snakemake.wildcards.mapper

rows = []
all_mm = sorted(set(totals.keys()) | set(metric_counts.keys()))
for mm in all_mm:
    df_t = totals.get(mm, pd.DataFrame({"length": lengths, "total": [0]*len(lengths)}))
    df_m = metric_counts.get(mm, pd.DataFrame({"length": lengths, metric: [0]*len(lengths)}))

    df_t = df_t.set_index("length").reindex(lengths, fill_value=0).astype(int)
    df_m = df_m.set_index("length").reindex(lengths, fill_value=0).astype(int)

    merged = df_t.join(df_m, how="outer").fillna(0)
    merged[pct_col] = merged.apply(
        lambda r: (100.0 * r[metric] / r["total"]) if r["total"] > 0 else 0.0, axis=1
    )

    for L, r in merged.iterrows():
        rows.append([mapper, mm, int(L), int(r["total"]), int(r[metric]), float(r[pct_col])])

out = pd.DataFrame(rows, columns=["mapper", "mm", "length", "total", metric, pct_col])
out.sort_values(["mm", "length"], inplace=True)
out.to_csv(snakemake.output[0], sep="\t", index=False)
