#include "options.h"



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
	dont_create_progress_stat_file(false), progress_stat_file("canopy_progress.out"), not_processed_profiles_file(""), show_progress_bar(false), print_time_statistics(true),
	die_on_kill(true), sparseMat(true), use_spearman(false), max_num_canopy_walks(6), filter_redundant(true),
	RNG_Seed(-1)
{

	bool hasErr = false;

	if (argc <= 1) { cout << "No input args given, returning\n"; exit(0); }


	for (int i = 0; i < argc; i++)
	{
		if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--input_file_path"))
			input_file_path = argv[++i];
		else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--output_clusters_file_path"))
			output_clusters_file_path = argv[++i];
		else if (!strcmp(argv[i], "--priority_reads_file_path"))
			priority_reads_file_path = argv[++i];
		else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--guide_matrix"))
			guide_matrix_file = argv[++i];
		else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--output_cluster_profiles_file"))
			output_cluster_profiles_file = argv[++i];
		else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--output_clusters_partial_file_path"))
			output_clusters_partial_file_path = argv[++i];
		else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--cluster_name_prefix"))
			output_cluster_prefix = argv[++i];
		else if (!strcmp(argv[i], "--profile_measure"))
			profile_measure_str = argv[++i];
		else if (!strcmp(argv[i], "--input_filter_file"))
			input_filter_file = argv[++i];
		else if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--num_threads"))
			num_threads = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--seed"))
			RNG_Seed = atoi(argv[++i]);

		else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbosity"))
			verbosity_option = argv[++i];
		else if (!strcmp(argv[i], "--sampleDistMatFile"))
			sampleDistMatFile = argv[++i];
		else if (!strcmp(argv[i], "--referenceMB2"))
			refMB2 = argv[++i];
		else if (!strcmp(argv[i], "--maxMB2genes"))
			refMB2_maxGenes = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--sampleDistLog"))
			sampleDistLog = argv[++i];
		else if (!strcmp(argv[i], "--sampleMinDist"))
			sampleMinDist = atof(argv[++i]);
		
		else if (!strcmp(argv[i], "--max_canopy_dist") || !strcmp(argv[i], "-d"))
			max_canopy_dist = atof(argv[++i]);
		else if (!strcmp(argv[i], "--max_canopy_dist_part"))
			max_canopy_dist_part = atof(argv[++i]);
		else if (!strcmp(argv[i], "--max_close_dist"))
			max_close_dist = atof(argv[++i]);//this is fixed normally
		else if (!strcmp(argv[i], "--max_merge_dist"))
			max_merge_dist = atof(argv[++i]);
		else if (!strcmp(argv[i], "--min_step_dist"))
			min_step_dist = atof(argv[++i]); //also fixed
		else if (!strcmp(argv[i], "--filter_min_obs"))
			filter_min_obs = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--filter_max_top3_sample_contribution"))
			filter_max_top3_sample_contribution = atof(argv[++i]);


		else if (!strcmp(argv[i], "--cag_filter_min_sample_obs"))
			cag_filter_min_sample_obs = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--cag_filter_max_top3_sample_contribution"))
			cag_filter_max_top3_sample_contribution = atof(argv[++i]);
		else if (!strcmp(argv[i], "--stop_criteria"))
			stop_after_num_seeds_processed = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--not_processed_profiles_file"))
			not_processed_profiles_file = argv[++i];
		else if (!strcmp(argv[i], "--progress_stat_file"))
			progress_stat_file = argv[++i];
		else if (!strcmp(argv[i], "--max_num_canopy_walks"))
			max_num_canopy_walks = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--dont_create_progress_stat_file"))
			dont_create_progress_stat_file = !dont_create_progress_stat_file;
		else if (!strcmp(argv[i], "--show_progress_bar"))
			show_progress_bar = !show_progress_bar;
		else if (!strcmp(argv[i], "--print_time_statistics"))
			print_time_statistics = !print_time_statistics;
		else if (!strcmp(argv[i], "--die_on_kill"))
			die_on_kill = !die_on_kill;
		else if (!strcmp(argv[i], "--high_mem"))
			sparseMat = false;
		else if (!strcmp(argv[i], "--redundant_guides"))
			filter_redundant = false;
		else if (!strcmp(argv[i], "--use_spearman"))
			use_spearman = !use_spearman;


		//


	}
	if (max_canopy_dist_part == -1){//not initialized yet
		max_canopy_dist_part = max_canopy_dist;
	}
}