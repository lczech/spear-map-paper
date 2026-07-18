# Pivot long → wide heat matrix:
# rows = mm, cols = length, values = pct_<metric>

import pandas as pd

len_min = int(snakemake.params.len_min)
len_max = int(snakemake.params.len_max)
pct_col = f"pct_{snakemake.params.metric}"

df = pd.read_csv(snakemake.input[0], sep="\t")
df["mm"] = df["mm"].astype(int)
df["length"] = df["length"].astype(int)
df[pct_col] = df[pct_col].astype(float)

rows = sorted(df["mm"].unique())
cols = list(range(len_min, len_max + 1))

mat = df.pivot_table(index="mm", columns="length", values=pct_col, fill_value=0.0)
mat = mat.reindex(index=rows, columns=cols, fill_value=0.0)

mat.to_csv(snakemake.output[0], sep="\t", float_format="%.6f")

