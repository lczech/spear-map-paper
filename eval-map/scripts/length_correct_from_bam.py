# Count correctly-mapped reads per true length, by joining mapped BAM hits (including
# secondary/supplementary alignments, to catch multimapping) against the ground-truth origin.
# A read counts as correct if ANY of its reported hits matches its true chrom + start exactly
# (BAM POS is 1-based; truth start is 0-based).

import pandas as pd

len_min = int(snakemake.params.len_min)
len_max = int(snakemake.params.len_max)

hits = pd.read_csv(
    snakemake.input.hits, sep="\t", header=None,
    names=["read_id", "flag", "rname", "pos"],
)
truth = pd.read_csv(snakemake.input.truth, sep="\t")
truth["length"] = truth["end"] - truth["start"]

if hits.empty:
    correct_ids = set()
else:
    merged = hits.merge(truth[["read_id", "chrom", "start"]], on="read_id", how="inner")
    is_correct = (merged["rname"] == merged["chrom"]) & (merged["pos"] - 1 == merged["start"])
    correct_ids = set(merged.loc[is_correct, "read_id"].unique())

truth["correct"] = truth["read_id"].isin(correct_ids)
counts = (
    truth.loc[truth["correct"], "length"]
    .value_counts()
    .rename_axis("length")
    .reset_index(name="count")
)
counts = counts[(counts["length"] >= len_min) & (counts["length"] <= len_max)]
counts = counts.sort_values("length")

counts.to_csv(snakemake.output[0], sep="\t", header=False, index=False)
