#!/usr/bin/env bash
set -euo pipefail

bin="${1:-./cc.bin}"
root_dir="$(cd "$(dirname "$0")" && pwd)"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

common_args=(
  -i "$root_dir/test/data/profiles.tsv"
  --filter_min_obs 0
  --filter_max_top3_sample_contribution 1
  --cag_filter_min_sample_obs 0
  --cag_filter_max_top3_sample_contribution 1
  --dont_create_progress_stat_file
  --stop_criteria 20
  --seed 7
)

"$bin" --version | grep -q "cc.bin v 0.26"

# Help is the public option reference: keep every accepted long option visible.
help_text="$("$bin" --help)"
for option in \
  --help --version --input_file_path --output_clusters_file_path \
  --output_cluster_profiles_file --output_clusters_partial_file_path \
  --cluster_name_prefix --input_filter_file --not_processed_profiles_file \
  --guide_matrix --referenceMB2 --maxMB2genes --redundant_guides \
  --max_canopy_dist --max_canopy_dist_part --max_merge_dist --max_close_dist \
  --min_step_dist --max_num_canopy_walks --profile_measure --use_spearman \
  --filter_min_obs --filter_max_top3_sample_contribution \
  --cag_filter_min_sample_obs --cag_filter_max_top3_sample_contribution \
  --sampleMinDist --sampleDistMatFile --sampleDistLog --num_threads --seed \
  --priority_reads_file_path --stop_criteria --progress_stat_file \
  --dont_create_progress_stat_file --show_progress_bar \
  --print_time_statistics --die_on_kill --dont_use_mmap --high_mem --verbosity
do
  grep -q -- "$option" <<< "$help_text"
done
grep -q "Input and output:" <<< "$("$bin")"

"$bin" "${common_args[@]}" --profile_measure mean \
  -o "$tmp_dir/mean.clusters" -c "$tmp_dir/mean.profiles"
test -s "$tmp_dir/mean.clusters"
test -s "$tmp_dir/mean.profiles"

# In-memory and streaming input, as well as different worker counts, must
# produce identical seeded clustering results.
"$bin" "${common_args[@]}" --profile_measure mean --dont_use_mmap -n 4 \
  -o "$tmp_dir/mean-streamed.clusters" -c "$tmp_dir/mean-streamed.profiles"
cmp "$tmp_dir/mean.clusters" "$tmp_dir/mean-streamed.clusters"
cmp "$tmp_dir/mean.profiles" "$tmp_dir/mean-streamed.profiles"

# Preserve the legacy short aliases accepted by the former Boost parser.
"$bin" "${common_args[@]}" -b -t -p MGS --profile_measure 75Q -n 2 \
  -o "$tmp_dir/short-options.clusters"
test -s "$tmp_dir/short-options.clusters"

"$bin" "${common_args[@]}" --profile_measure mean \
  -o "$tmp_dir/mean-repeat.clusters" -c "$tmp_dir/mean-repeat.profiles"
cmp "$tmp_dir/mean.clusters" "$tmp_dir/mean-repeat.clusters"
cmp "$tmp_dir/mean.profiles" "$tmp_dir/mean-repeat.profiles"

"$bin" "${common_args[@]}" --use_spearman \
  -o "$tmp_dir/spearman.clusters"
test -s "$tmp_dir/spearman.clusters"

"$bin" "${common_args[@]}" --sampleMinDist 0 \
  --sampleDistMatFile "$tmp_dir/sample-distances.tsv" \
  --sampleDistLog "$tmp_dir/removed-samples.txt" \
  -o "$tmp_dir/autocorrelation.clusters"
test -s "$tmp_dir/sample-distances.tsv"
awk -F '\t' 'NR == 1 { upper = $2 } NR == 2 && $1 != upper { exit 1 }' \
  "$tmp_dir/sample-distances.tsv"

# A-B and B-C are within the threshold, while A-C is not. The middle-depth
# sample must be removed without letting that removed sample eliminate C.
autocorrelation_args=(
  -i "$root_dir/test/data/autocorrelation-chain.tsv"
  --filter_min_obs 1
  --filter_max_top3_sample_contribution 1
  --cag_filter_min_sample_obs 0
  --cag_filter_max_top3_sample_contribution 1
  --dont_create_progress_stat_file
  --sampleMinDist 0.1
  --stop_criteria 20
  --seed 7
)
for threads in 1 4; do
  "$bin" "${autocorrelation_args[@]}" -n "$threads" \
    --sampleDistMatFile "$tmp_dir/autocorrelation-$threads.tsv" \
    --sampleDistLog "$tmp_dir/autocorrelation-$threads.log" \
    -o "$tmp_dir/autocorrelation-$threads.clusters" \
    -c "$tmp_dir/autocorrelation-$threads.profiles"
done
cmp "$tmp_dir/autocorrelation-1.tsv" "$tmp_dir/autocorrelation-4.tsv"
cmp "$tmp_dir/autocorrelation-1.log" "$tmp_dir/autocorrelation-4.log"
cmp "$tmp_dir/autocorrelation-1.clusters" "$tmp_dir/autocorrelation-4.clusters"
cmp "$tmp_dir/autocorrelation-1.profiles" "$tmp_dir/autocorrelation-4.profiles"
test "$(awk 'NR == 2 { print $1 }' "$tmp_dir/autocorrelation-1.log")" = "1"
awk -F '\t' 'NF != 4 { exit 1 }' "$tmp_dir/autocorrelation-1.profiles"

"$bin" -i "$root_dir/test/data/top3-boundary.tsv" \
  --filter_min_obs 1 --filter_max_top3_sample_contribution 0.75 \
  --cag_filter_min_sample_obs 0 --cag_filter_max_top3_sample_contribution 1 \
  --dont_create_progress_stat_file --input_filter_file "$tmp_dir/top3-filtered.tsv" \
  --stop_criteria 20 --seed 7 -o "$tmp_dir/top3.clusters"
test "$(wc -l < "$tmp_dir/top3-filtered.tsv")" -eq 2
grep -q $'^concentrated\tmax_top3_sample_contribution_filter$' \
  "$tmp_dir/top3-filtered.tsv"

"$bin" "${common_args[@]}" -n 1 \
  -g "$root_dir/test/data/guides.tsv" --redundant_guides \
  -o "$tmp_dir/guided.clusters" -r "$tmp_dir/partial.clusters"
grep -q '^guideA' "$tmp_dir/guided.clusters"
grep -q '^guideB' "$tmp_dir/guided.clusters"
test -s "$tmp_dir/partial.clusters"

"$bin" "${common_args[@]}" \
  --referenceMB2 "$root_dir/test/data/metabat-reference.tsv" \
  --redundant_guides -o "$tmp_dir/metabat-guided.clusters"
grep -q '^binA' "$tmp_dir/metabat-guided.clusters"
grep -q '^binB' "$tmp_dir/metabat-guided.clusters"

if "$bin" --num_threads invalid >/dev/null 2>&1; then
  echo "invalid numeric option was accepted" >&2
  exit 1
fi
if "$bin" --definitely-unknown >/dev/null 2>&1; then
  echo "unknown option was accepted" >&2
  exit 1
fi

echo "All regression checks passed."
