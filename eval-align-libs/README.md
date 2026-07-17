# Evaluation of alignment libraries

This is a selfcontained benchmark to test which of the many existing C/C++ libraries for sequence alignment work best for our use case: aligning short reads with high error rates (in particular ancient DNA damage at the read ends) to genome subsequences that span the read length plus some extra intervals outside, of some hundred base pairs max in each direction of the read ends.

Compile:

```
cd spear-map-paper/eval-align-libs
cmake -B build && cmake --build build 8
```

Run:

```
./bin/eval-align-libs --max-chromosomes 1 --max-reads 25000 ../data/GCA_000001405.15_GRCh38_genomic.fna.gz
```

Plot:

```
python scripts/plot_accuracy_vs_speed.py
python scripts/plot_accuracy_by_len.py
python scripts/plot_offsets.py
python scripts/plot_scores.py
python scripts/plot_scores_by_len.py
python scripts/plot_timing.py
python scripts/plot_timing_by_len.py
```

```
python scripts/plot_accuracy_vs_speed.py ; python scripts/plot_accuracy_by_len.py ; python scripts/plot_offsets.py ; python scripts/plot_scores.py ; python scripts/plot_scores_by_len.py ; python scripts/plot_timing.py ; python scripts/plot_timing_by_len.py
```


## Summary of findings

We tested edlib, ksw2, parasail, and WFA2 for fitting (semi-global) alignment of short reads (30–150 bp, gamma-distributed) against reference windows of roughly 160–290 bp. Reads were drawn from a real genome, mutated with a grid of substitution rates, indel rates, and aDNA damage rates, and then re-aligned to the window they originated from. Alignment accuracy is measured as the offset of the reported start/end position relative to the true read origin.

### ksw2

Both `ksw2-score` and `ksw2-cigar` fail for approximately 36% of reads even on the clean (zero-mutation) baseline. The root cause is that `KSW_EZ_EXTZ_ONLY` performs **extension alignment** anchored at position (0, 0) of both query and target — it is not fitting alignment. To reach the read's true start position the DP must pay a deletion penalty to skip the preceding window context. When that penalty exceeds the read's match score (`4 + 2*(true_start − 1) > 2 * read_length`, i.e. roughly when `true_start ≥ read_length`), the extension finds a different, wrong path. For `ksw2-score` this produces negative end offsets (the extension terminates short of the true end); for `ksw2-cigar` the reverse-trick pass has the same flaw in the reverse direction, producing symmetric positive offsets. ksw2 would require a proper semi-global DP with a free-start row to handle this use case correctly.

### parasail

All parasail `sg_dx` variants (score and cigar, both matrices, hot and cold) perform true fitting alignment and get the correct position for >99.6% of reads on the clean baseline. The small fraction of non-zero offsets are genuine ties in repetitive reference regions. Under mutation, both negative (edlib-style leftmost) and positive offsets appear, but the absolute counts remain modest.

### edlib

Performs true fitting alignment (`EDLIB_MODE_HW`) and gives correct positions for >99.6% of reads on the clean baseline (the same tie cases as parasail). Under mutation, edlib shows a leftward bias: it returns the leftmost minimum-edit-distance alignment, so non-zero offsets are predominantly negative. This makes it look very accurate relative to the true origin, though this is partly a tie-breaking artifact rather than strictly better alignment quality. Edlib has no precomputation step and no hot/cold distinction, making per-call cost accounting straightforward.

### WFA2

Finding: wfadaptive is harmful for this use case. The adaptive heuristic prunes wavefront diagonals that are more than `max_distance_threshold` (50) steps from the current leading diagonal. For long-read near-diagonal alignment that is a useful approximation, but for fitting alignment of short reads in a window that is significantly wider than the read, the true alignment diagonal can be far from the main diagonal (offset up to `true_start ≈ 143`). The heuristic was pruning those diagonals, forcing WFA2 to find costly alternative paths with many extra mismatches and gaps — producing large positive end offsets and inflated penalties for ~29% of reads under mutation. Counterintuitively, the heuristic was also *slower* than the exact algorithm, because the wrong paths required more wavefront steps to complete.

### WFA2 scoring vs parasail

WFA2 operates in penalty space (match=0, lower=better), while parasail operates in score space (match=+2, higher=better). Although both use the same numerical values for gap_open=4 and gap_extend=2, they are not equivalent: the parasail match reward creates an asymmetry between insertions (query bases with no ref) and deletions (ref bases with no query) that WFA2's symmetric affine gap model cannot replicate with any single parameter choice. Concretely, parasail's effective gap/mismatch cost ratio is 8/6 ≈ 1.33, while WFA2's default is 6/4 = 1.5 — WFA2 is slightly more gap-averse. Raising WFA2's gap_open to 6 (mismatch=4, gap_open=6, gap_extend=2) increases the ratio to 2.0 and matches the BWA-MEM penalty scheme (Marco-Sola et al. 2021, Table 3). In practice this slightly improves offset accuracy at high divergence by discouraging spurious gap-opening near alignment boundaries.
