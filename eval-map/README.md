# Evaluation of mappers (short read aligners)

This is a first simple benchmark of short read mappers, comparing them to `spear map align`.

To this end, we use the human reference genome, build indices, simulate reads of different
lengths and divergence, and map them again, measuring how many are missed and how many landed
at their true origin.

## Setup

1. Build the local C++ tools (`generate_fragments`, the read simulator):

       cd ..
       ./build_all.sh
       cd eval-map

2. Build `spear` itself, in its own sibling repo, then point `spear_binary` in the `Snakefile`'s
   `config` dict at the resulting binary (default assumes `spear` and `spear-map-paper` are
   checked out side by side; adjust for your setup, e.g. to an absolute path on the cluster).
   Also set `module_load` there if your system needs an environment-module command to make the
   runtime libraries for `generate_fragments`/`spear` available (e.g. on the cluster); leave it
   as `""` on a plain local machine.

3. Create and activate the environment used to run Snakemake itself:

       conda env create -f envs/snakemake.yaml
       conda activate spear-eval-map

   The per-rule tool environments (`envs/mappers.yml`, `envs/python.yml`, `envs/utils.yml`) are
   created automatically by Snakemake via `--use-conda` (already enabled in the profiles below).

## Running

Do a dry run first to check the plan:

    snakemake --profile profiles/local -n

Then run for real, locally:

    snakemake --profile profiles/local

Or on the cluster (SLURM; adjust `profiles/config.yaml`'s account/partition first):

    snakemake --profile profiles
