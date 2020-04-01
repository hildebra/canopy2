


#pragma once



#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "program_options_misc.hpp"





struct options
{
public:
	options(int argc, char** argv);

	//Prepare variables for command line input
	string input_file_path;
	string input_filter_file;
	string priority_reads_file_path; //This optional file contains names of the reads that will be used as canopy centers first, then clustering proceeds as before (points chosen by random)
	string output_clusters_file_path;
	string output_clusters_partial_file_path;
	string output_cluster_profiles_file;
	string output_cluster_prefix;
	string sampleDistMatFile;//write out matrix of sample distances
	string sampleDistLog;
	string profile_measure_str;
	string guide_matrix_file;
	string refMB2;
	int refMB2_maxGenes;
	int num_threads;
	PRECISIONT max_canopy_dist;
	PRECISIONT max_canopy_dist_part;
	PRECISIONT max_close_dist; //The value is hardcoded and the option to change it removed from CLI to not confuse users
	PRECISIONT max_merge_dist;
	PRECISIONT min_step_dist; //The value is hardcoded and the option to change it removed from CLI to not confuse users 
	PRECISIONT sampleMinDist;
	string verbosity_option;
	int filter_min_obs;
	double filter_max_top3_sample_contribution;
	int cag_filter_min_sample_obs;
	double cag_filter_max_top3_sample_contribution;
	double stop_after_num_seeds_processed;
	bool dont_create_progress_stat_file;
	string progress_stat_file;
	string not_processed_profiles_file;
	bool show_progress_bar;
	bool print_time_statistics;
	bool die_on_kill;
	bool sparseMat;
	bool dont_use_mmap;
	bool use_spearman; // simply converts input points to ranks
	int max_num_canopy_walks;
};
