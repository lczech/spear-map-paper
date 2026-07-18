from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

matrix_path = snakemake.input[0]
out_path    = snakemake.output[0]
title       = snakemake.params.title
cbar_label  = snakemake.params.cbar_label

df = pd.read_csv(matrix_path, sep="\t", index_col=0)
df.index = df.index.astype(int)
df = df.sort_index()
df.columns = df.columns.astype(int)
df = df[sorted(df.columns)]

arr = df.to_numpy(dtype=float)

fig = plt.figure(figsize=(10, 5.5))
ax = plt.gca()
im = ax.imshow(arr, aspect="auto", origin="lower", vmin=0, vmax=100)

ax.set_yticks(range(len(df.index)))
ax.set_yticklabels(df.index)
ax.set_ylabel("Mismatch rate (%)")

cols = df.columns.to_list()
# step = max(1, len(cols) // 15)
step = 10
xticks_idx = list(range(0, len(cols), step))
ax.set_xticks(xticks_idx)
ax.set_xticklabels([str(cols[i]) for i in xticks_idx])
ax.set_xlabel("Fragment length (bp)")

ax.set_title(title)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label(cbar_label)

fig.tight_layout()
Path(out_path).parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out_path, dpi=200)
plt.close(fig)
