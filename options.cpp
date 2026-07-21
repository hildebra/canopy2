#include "options.h"

namespace {

void print_help(std::ostream& out)
{
	out << R"HELP(cc.bin - metagenomic co-abundance canopy clustering

Usage:
  cc.bin [OPTIONS] INPUT [CLUSTERS [CLUSTER_PROFILES]]
  cc.bin -i INPUT [OPTIONS]

INPUT is a tab-separated abundance matrix. Its first row is a header; each
remaining row contains a unique profile ID followed by two or more finite,
non-negative sample abundances. All rows must contain the same samples.

Input and output:
  -i, --input_file_path FILE       Input abundance matrix (required).
  -o, --output_clusters_file_path FILE
                                    Cluster-membership output [clusters_out].
  -c, --output_cluster_profiles_file FILE
                                    Optional cluster-profile output.
  -r, --output_clusters_partial_file_path FILE
                                    Optional partial-correlation output; requires
                                    --guide_matrix or --referenceMB2.
  -p, --cluster_name_prefix TEXT   Prefix for de novo cluster IDs [empty].
      --input_filter_file FILE     Write IDs and reasons for profiles removed by
                                    the input filters.
      --not_processed_profiles_file FILE
                                    Write profiles left when interrupted.

Guided clustering:
  -g, --guide_matrix FILE          Headerless guide profiles (ID plus values),
                                    with the same sample columns as INPUT.
      --referenceMB2 FILE          Tab-separated MetaBAT2 bin/profile assignments.
      --maxMB2genes N              Maximum profiles used per MetaBAT2 bin [1000].
      --redundant_guides           Permit a profile to remain assigned to more
                                    than one guide (default keeps its best guide).

Clustering and correlation:
  -d, --max_canopy_dist FLOAT      Maximum 1-correlation distance from a canopy
                                    center [0.1].
      --max_canopy_dist_part FLOAT Maximum partial-correlation distance [same as
                                    --max_canopy_dist].
      --max_merge_dist FLOAT       Merge canopies whose centers are closer than
                                    this distance [0.1].
      --max_close_dist FLOAT       Candidate-neighbour distance [0.6].
      --min_step_dist FLOAT        Minimum centroid movement during a walk [0.001].
      --max_num_canopy_walks N     Maximum walk iterations per seed [6].
      --profile_measure METHOD     Cluster profile: median, mean, 75Q, 80Q, 85Q,
                                    90Q, or 95Q [75Q].
      --use_spearman               Use Spearman instead of Pearson correlation.

Input filters:
      --filter_min_obs N           Require at least N non-zero samples; 0 disables
                                    the filter [3].
      --filter_max_top3_sample_contribution FLOAT
                                    Remove a profile if its three largest samples
                                    contribute more than this fraction; 1 disables
                                    the filter [0.9].

Output filters (de novo clustering):
      --cag_filter_min_sample_obs N
                                    Require at least N non-zero cluster-profile
                                    samples; 0 disables the filter [3].
      --cag_filter_max_top3_sample_contribution FLOAT
                                    Remove clusters exceeding the top-three sample
                                    contribution; 1 disables the filter [0.9].
  Clusters with fewer than two members are always omitted.

Autocorrelated-sample filtering:
      --sampleMinDist FLOAT        Remove redundant samples below this
                                    1-correlation distance; 2 disables [2].
      --sampleDistMatFile FILE     Write the deterministic sample-distance matrix.
      --sampleDistLog FILE         Write the indices of removed samples.
  This filter is disabled when guide profiles are used.

Execution, progress, and memory:
  -n, --num_threads N              Worker threads [1].
  -s, --seed N                     Random seed [current time]. A fixed seed makes
                                    results reproducible across thread counts and
                                    input modes.
      --priority_reads_file_path FILE
                                    Ordered profile IDs to consider before the
                                    remaining seeded/shuffled profiles.
      --stop_criteria N            Stop after processing N seeds; 0 disables
                                    early stopping [50000].
      --progress_stat_file FILE    Per-seed progress statistics
                                    [canopy_progress.out].
      --dont_create_progress_stat_file
                                    Do not create a progress-statistics file.
  -b, --show_progress_bar          Show a terminal progress bar.
  -t, --print_time_statistics      Print phase timing statistics.
      --die_on_kill                Exit immediately on interrupt instead of writing
                                    merged partial results.
      --dont_use_mmap              Stream input lines instead of buffering the
                                    complete input file (lower peak memory).
      --high_mem                   Store dense rather than sparse profiles.
  -v, --verbosity LEVEL            error, progress, warn, info, debug, debug1,
                                    debug2, or debug3 [info].

Other:
  -h, --help                       Show this help and exit.
      --version                    Show the version and exit (-v also does this
                                    when it is the only argument).

Examples:
  cc.bin -i profiles.tsv -o clusters.tsv -c cluster_profiles.tsv -n 8 --seed 42
  cc.bin -i profiles.tsv -g guides.tsv -o guide_members.tsv --seed 42
  cc.bin profiles.tsv clusters.tsv --dont_use_mmap --sampleMinDist 0.15
)HELP";
}

} // namespace



options::options(int argc, char** argv):
	input_file_path(""), input_filter_file(""),
	priority_reads_file_path(""), output_clusters_file_path("clusters_out"), output_clusters_partial_file_path(""),
	output_cluster_profiles_file(""), output_cluster_prefix(""), sampleDistMatFile(""),
	sampleDistLog(""),
	profile_measure_str("75Q"), guide_matrix_file(""), refMB2(""),refMB2_maxGenes(1000),
	num_threads(1), max_canopy_dist(0.1), max_canopy_dist_part(-1), max_close_dist(0.6),
	max_merge_dist(0.1), min_step_dist(0.001), sampleMinDist(2), verbosity_option("info"),
	filter_min_obs(3), filter_max_top3_sample_contribution(0.9), cag_filter_min_sample_obs(3),
	cag_filter_max_top3_sample_contribution(0.9),
	stop_after_num_seeds_processed(50000),
	dont_create_progress_stat_file(false), progress_stat_file("canopy_progress.out"), not_processed_profiles_file(""), show_progress_bar(false), print_time_statistics(false),
	die_on_kill(false), sparseMat(true), dont_use_mmap(false), use_spearman(false), max_num_canopy_walks(6), filter_redundant(true),
	RNG_Seed(-1)
{
	if (argc <= 1) {
		print_help(cout);
		exit(0);
	}

	int positional_arg = 0;
	for (int i = 1; i < argc; i++)
	{
		auto next_value = [&](const char* option_name) -> char* {
			if (i + 1 >= argc) {
				cerr << "Missing value for option: " << option_name << "\n";
				exit(2);
			}
			return argv[++i];
		};
		auto parse_int = [&](const char* option_name) -> int {
			const char* value = next_value(option_name);
			char* end = NULL;
			errno = 0;
			long parsed = strtol(value, &end, 10);
			if (errno == ERANGE || end == value || *end != '\0' || parsed < INT_MIN || parsed > INT_MAX) {
				cerr << "Invalid integer for option " << option_name << ": " << value << "\n";
				exit(2);
			}
			return static_cast<int>(parsed);
		};
		auto parse_double = [&](const char* option_name) -> double {
			const char* value = next_value(option_name);
			char* end = NULL;
			errno = 0;
			double parsed = strtod(value, &end);
			if (errno == ERANGE || end == value || *end != '\0' || !std::isfinite(parsed)) {
				cerr << "Invalid number for option " << option_name << ": " << value << "\n";
				exit(2);
			}
			return parsed;
		};

		if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
			print_help(cout);
			exit(0);
		}

		else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--input_file_path"))
			input_file_path = next_value(argv[i]);
		else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--output_clusters_file_path"))
			output_clusters_file_path = next_value(argv[i]);
		else if (!strcmp(argv[i], "--priority_reads_file_path"))
			priority_reads_file_path = next_value(argv[i]);
		else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--guide_matrix"))
			guide_matrix_file = next_value(argv[i]);
		else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--output_cluster_profiles_file"))
			output_cluster_profiles_file = next_value(argv[i]);
		else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--output_clusters_partial_file_path"))
			output_clusters_partial_file_path = next_value(argv[i]);
		else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--cluster_name_prefix"))
			output_cluster_prefix = next_value(argv[i]);
		else if (!strcmp(argv[i], "--profile_measure"))
			profile_measure_str = next_value(argv[i]);
		else if (!strcmp(argv[i], "--input_filter_file"))
			input_filter_file = next_value(argv[i]);
		else if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--num_threads"))
			num_threads = parse_int(argv[i]);
		else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--seed"))
			RNG_Seed = parse_int(argv[i]);

		else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbosity"))
			verbosity_option = next_value(argv[i]);
		else if (!strcmp(argv[i], "--sampleDistMatFile"))
			sampleDistMatFile = next_value(argv[i]);
		else if (!strcmp(argv[i], "--referenceMB2"))
			refMB2 = next_value(argv[i]);
		else if (!strcmp(argv[i], "--maxMB2genes"))
			refMB2_maxGenes = parse_int(argv[i]);
		else if (!strcmp(argv[i], "--sampleDistLog"))
			sampleDistLog = next_value(argv[i]);
		else if (!strcmp(argv[i], "--sampleMinDist"))
			sampleMinDist = parse_double(argv[i]);
		
		else if (!strcmp(argv[i], "--max_canopy_dist") || !strcmp(argv[i], "-d"))
			max_canopy_dist = parse_double(argv[i]);
		else if (!strcmp(argv[i], "--max_canopy_dist_part"))
			max_canopy_dist_part = parse_double(argv[i]);
		else if (!strcmp(argv[i], "--max_close_dist"))
			max_close_dist = parse_double(argv[i]);//this is fixed normally
		else if (!strcmp(argv[i], "--max_merge_dist"))
			max_merge_dist = parse_double(argv[i]);
		else if (!strcmp(argv[i], "--min_step_dist"))
			min_step_dist = parse_double(argv[i]); //also fixed
		else if (!strcmp(argv[i], "--filter_min_obs"))
			filter_min_obs = parse_int(argv[i]);
		else if (!strcmp(argv[i], "--filter_max_top3_sample_contribution"))
			filter_max_top3_sample_contribution = parse_double(argv[i]);


		else if (!strcmp(argv[i], "--cag_filter_min_sample_obs"))
			cag_filter_min_sample_obs = parse_int(argv[i]);
		else if (!strcmp(argv[i], "--cag_filter_max_top3_sample_contribution"))
			cag_filter_max_top3_sample_contribution = parse_double(argv[i]);
		else if (!strcmp(argv[i], "--stop_criteria"))
			stop_after_num_seeds_processed = parse_int(argv[i]);
		else if (!strcmp(argv[i], "--not_processed_profiles_file"))
			not_processed_profiles_file = next_value(argv[i]);
		else if (!strcmp(argv[i], "--progress_stat_file"))
			progress_stat_file = next_value(argv[i]);
		else if (!strcmp(argv[i], "--max_num_canopy_walks"))
			max_num_canopy_walks = parse_int(argv[i]);
		else if (!strcmp(argv[i], "--dont_create_progress_stat_file"))
			dont_create_progress_stat_file = true;
		else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--show_progress_bar"))
			show_progress_bar = true;
		else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--print_time_statistics"))
			print_time_statistics = true;
		else if (!strcmp(argv[i], "--die_on_kill"))
			die_on_kill = true;
		else if (!strcmp(argv[i], "--high_mem"))
			sparseMat = false;
		else if (!strcmp(argv[i], "--dont_use_mmap"))
			dont_use_mmap = true;
		else if (!strcmp(argv[i], "--redundant_guides"))
			filter_redundant = false;
		else if (!strcmp(argv[i], "--use_spearman"))
			use_spearman = true;
		else if (argv[i][0] != '-') {
			if (positional_arg == 0) input_file_path = argv[i];
			else if (positional_arg == 1) output_clusters_file_path = argv[i];
			else if (positional_arg == 2) output_cluster_profiles_file = argv[i];
			else {
				cerr << "Unexpected positional argument: " << argv[i] << "\n";
				exit(2);
			}
			positional_arg++;
		}
		else {
			cerr << "Unknown option: " << argv[i] << "\n";
			exit(2);
		}
	}
	if (max_canopy_dist_part == -1){//not initialized yet
		max_canopy_dist_part = max_canopy_dist;
	}
}
