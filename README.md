# canopy2

`canopy2` (`cc.bin`) clusters abundance profiles into co-abundant gene groups
using canopy clustering. It supports de novo clustering, correlation against
predefined guide profiles, MetaBAT2-derived guides, Pearson or Spearman
distance, and optional removal of autocorrelated samples.

This is a faster and more memory-efficient implementation of the canopy
clustering approach described by Nielsen et al. in
[Nature Biotechnology (2014)](https://doi.org/10.1038/nbt.2939).

Current version: **0.26**.

## Requirements and build

The supplied Makefile expects:

- a C++11 compiler with OpenMP support (normally GCC/G++);
- GNU Make;
- pthreads and the zlib development library;
- Bash, `awk`, and `cmp` to run the test suite.

On Debian or Ubuntu, the required build packages can be installed with:

```bash
sudo apt-get install build-essential zlib1g-dev
```

Build and verify the program:

```bash
git clone https://github.com/hildebra/canopy2.git
cd canopy2
make -j
./cc.bin --version
make test
```

Linux is the primary target. Windows users can build in WSL; other platforms
need a compiler/toolchain that provides OpenMP and zlib.

## Input format

The primary input is a tab-separated abundance matrix. The first line is a
header. Every other line contains a unique, non-empty profile ID followed by
the same two or more sample columns:

```text
profile_id	sample_1	sample_2	sample_3	sample_4
gene_1	1	2	3	4
gene_2	1.5	3	4.5	6
gene_3	4	3	2	1
```

Abundances must be finite, numeric, and non-negative. Empty fields, inconsistent
row lengths, duplicate IDs, and matrices with fewer than two sample columns are
rejected.

A file supplied with `--guide_matrix` has no header: each line is a guide ID
followed by numeric sample values. It must contain exactly the same number and
order of samples as the primary matrix.

`--referenceMB2` expects tab-separated `bin_id` and `profile_id` fields. All
profiles must occur in the primary matrix, and records for each bin must be
contiguous.

## Basic use

De novo clustering:

```bash
./cc.bin \
  -i profiles.tsv \
  -o clusters.tsv \
  -c cluster_profiles.tsv \
  -n 8 \
  --seed 42
```

The input and first two output paths may also be positional:

```bash
./cc.bin profiles.tsv clusters.tsv cluster_profiles.tsv -n 8 --seed 42
```

Guided clustering:

```bash
./cc.bin \
  -i profiles.tsv \
  --guide_matrix guides.tsv \
  -o guide_members.tsv \
  --seed 42
```

Sample-autocorrelation filtering (not used in guided mode):

```bash
./cc.bin \
  -i profiles.tsv \
  -o clusters.tsv \
  --sampleMinDist 0.15 \
  --sampleDistMatFile sample_distances.tsv \
  --sampleDistLog removed_samples.tsv \
  --seed 42
```

Run `./cc.bin --help` for the complete option reference, including all default
values and filter controls.

## Important settings

Distances are expressed as `1 - correlation`; smaller values therefore mean
more similar profiles. `--max_canopy_dist` controls membership and
`--max_merge_dist` controls merging. Pearson correlation is the default;
`--use_spearman` ranks each profile before correlation.

The default cluster profile is the sample-wise 75th percentile. Select
`median`, `mean`, or another supported percentile with `--profile_measure`.

Input profiles are filtered by minimum sample observations and top-three sample
contribution. Similar filters are applied to final de novo cluster profiles.
Use `--help` to review their defaults before changing or disabling them.

## Outputs

Output files are tab-separated and do not include headers:

- De novo membership output: `cluster_id`, `profile_id`.
- De novo profile output: `cluster_id`, followed by one value per original
  sample.
- Guided membership output: `guide_id`, `profile_id`, `1-correlation` distance.
- Partial-correlation output: the same three columns, using the partial distance.

In de novo mode, clusters with fewer than two members are omitted. The optional
input-filter report does include a header and records the profile ID and filter
reason.

## Reproducibility, memory, and interruption

Pass an explicit `--seed` for reproducible results. With the same input and
options, seeded clustering and sample-autocorrelation filtering are deterministic
across supported thread counts and between buffered and streaming input modes.
Without `--seed`, the current time is used.

By default the input file is buffered and profiles use a sparse representation.
`--dont_use_mmap` streams input lines to reduce peak memory. `--high_mem` selects
a dense profile representation and generally uses more memory.

The program writes per-seed progress statistics to `canopy_progress.out` unless
`--dont_create_progress_stat_file` is supplied. On `Ctrl+C`, it normally stops
clustering cleanly, merges the canopies already found, and writes partial
results. `--die_on_kill` requests immediate termination instead.

## Testing

Run the regression suite after changing the code:

```bash
make test
```

The suite covers de novo and guided clustering, MetaBAT2 guides, input modes,
thread-count determinism, autocorrelation, Pearson/Spearman behavior, filtering,
and invalid command-line options.

## License and attribution

This project is distributed under GPL-2.0-or-later; see [LICENSE](LICENSE).

The reimplementation is based on the original code by Piotr Dworzynski and has
since been improved and maintained by Falk Hildebrand.

Copyright (C) 2019-2024 Falk Hildebrand
(`Falk.Hildebrand@gmail.com`), Quadram Institute.

Copyright (C) 2013-2014 Piotr Dworzynski
(`piotr@cbs.dtu.dk`), Technical University of Denmark
