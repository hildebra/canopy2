/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
 * Copyright (C) 2019, 2020, 2021 Falk Hildebrand (Falk.Hildebrand@gmail.com)
 *
 * This file is part of Metagenomics Canopy Clustering Implementation.
 *
 * Metagenomics Canopy Clustering Implementation is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Metagenomics Canopy Clustering Implementation is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

//#include <boost/program_options.hpp>
//#include <boost/type_index.hpp>

//#include <boost/iostreams/stream.hpp>
//#include <boost/iostreams/device/mapped_file.hpp>

//#include <boost/assign/std/vector.hpp>

#include "Canopy.hpp"
#include "CanopyClustering.hpp"

using namespace std;
//using namespace boost::program_options;
//using namespace boost::assign;



ProfileMeasureType profile_measure;

int main(int argc, char* argv[])
{
    //
    //Initialization
    //
    
    //Set initial logging level
    log_level = logINFO;


    //.22: fixed bug when more cores than samples were used
    //.24: added option to seed RNG
	//.25: conda release
    string ccbinVer = "0.25";

    //Preapre Time Profile
    TimeProfile time_profile;
    time_profile.start_timer("Total");

	//test_annoy();

    if (argc == 2 && (!strcmp(argv[1], "-v") || !strcmp(argv[1], "--version")) ) {
        cout << "cc.bin v " << ccbinVer << endl;
        exit(0);

    }

	options* opt = new options(argc, argv);
    cout << "cc.bin v " << ccbinVer << endl;

    int seed = opt->RNG_Seed;

    if (seed == -1) { seed = unsigned(std::time(NULL)); }
    std::srand(seed);


    //Prepare variables for command line input
	string input_file_path = opt->input_file_path;
    string input_filter_file = opt->input_filter_file;
    string priority_reads_file_path = opt->priority_reads_file_path; //This optional file contains names of the reads that will be used as canopy centers first, then clustering proceeds as before (points chosen by random)
	string output_clusters_file_path = opt->output_clusters_file_path;
	string output_clusters_partial_file_path = opt->output_clusters_partial_file_path;
	string output_cluster_profiles_file = opt->output_cluster_profiles_file;
	string profile_measure_str = opt->profile_measure_str;
	string guide_matrix_file = opt->guide_matrix_file;
	int num_threads = opt->num_threads;
	PRECISIONT max_canopy_dist = opt->max_canopy_dist;
	PRECISIONT max_canopy_dist_part = opt->max_canopy_dist_part;
	const PRECISIONT max_close_dist = opt->max_close_dist;  //The value is hardcoded and the option to change it removed from CLI to not confuse users
	PRECISIONT max_merge_dist = opt->max_merge_dist;
	const PRECISIONT min_step_dist = opt->min_step_dist;  //The value is hardcoded and the option to change it removed from CLI to not confuse users 
	string verbosity_option = opt->verbosity_option;
	
	int cag_filter_min_sample_obs = opt->cag_filter_min_sample_obs;
	double cag_filter_max_top3_sample_contribution = opt->cag_filter_max_top3_sample_contribution;
	double stop_after_num_seeds_processed = opt->stop_after_num_seeds_processed;
	string progress_stat_file = opt->progress_stat_file;
	string not_processed_profiles_file = opt->not_processed_profiles_file;
	bool show_progress_bar = opt->show_progress_bar;
	bool print_time_statistics = opt->print_time_statistics;
	bool die_on_kill = opt->die_on_kill;
	bool dont_use_mmap = opt->dont_use_mmap;
	bool use_spearman = opt->use_spearman;  // simply converts input points to ranks
	int max_num_canopy_walks = opt->max_num_canopy_walks;
	const string arr[] = { "median", "mean", "75Q", "80Q", "85Q", "90Q", "95Q" };
	const vector<string> valid_profile_measure_values (arr, arr + sizeof(arr) / sizeof(arr[0]) );

    bool filter_redundant = opt->filter_redundant;// true;//filter redundant genes assigned to multiple guides?


    //changed to simpler option parser to remove external lib dependency
    /*
    //Define and read command line options
    options_description all_options_desc("Allowed options");
    options_description options_shown_in_help_desc("Allowed options");
    options_description options_not_shown_in_help_desc("Options not shown in help");
    options_description general_options_desc("General");
    options_description algorithm_param_options_desc("Algorithm Parameters");
    options_description filter_in_options_desc("Input filter parameters");
    options_description filter_out_options_desc("Output filter parameters");
    options_description early_stop_options_desc("Early stopping");
    options_description misc_options_desc("Miscellaneous");


    general_options_desc.add_options()
        ("input_file_path,i", value<string>(&input_file_path), "Path to the input file")
        ("priority_reads_file_path", value<string>(&priority_reads_file_path)->default_value(""), "Path to (optional) file containing an ordered list (line by line) of read names according to which they should considered as cluster seeds or members. Reads from the input file not present in this file will be shuffled and considered second.")
		("guide_matrix,g", value<string>(&guide_matrix_file)->default_value(""), "Guide matrix that will be correlated against new matrix")
		("output_clusters_file_path,o", value<string>(&output_clusters_file_path)->default_value("clusters_out"), "Path to file to which clusters will be written")
		("output_clusters_partial_file_path,r", value<string>(&output_clusters_partial_file_path)->default_value(""), "Path to file to which partially correlating genes will be written")
		
		("output_cluster_profiles_file,c", value<string>(&output_cluster_profiles_file)->default_value(""), "Path to file to which cluster profiles will be written")
        ("cluster_name_prefix,p", value<string>(&output_cluster_prefix)->default_value("CAG"), "Prefix prepended to output cluster names")
        ("num_threads,n", value<int>(&num_threads)->default_value(4), "Number of cpu threads to use.")
        ("verbosity,v", value<string>(&verbosity_option)->default_value("info"), "Control how much information should be printed to the screen. Available levels according to their verbosity: error, progress, warn, info, debug, debug1.");

    algorithm_param_options_desc.add_options()
		
		("use_spearman", bool_switch(&use_spearman), "If set, calculates Spearman correlation coefficients instead of (default) Pearson correlation")
		("max_canopy_dist", value<PRECISIONT>(&max_canopy_dist)->default_value(0.1), "Max pearson correlation difference between a canopy center and a point included to the canopy")
        //This option is removed from CLI to avoid user confusion. The default value is hardcoded above
        //("max_close_dist", value<double>(&max_close_dist)->default_value(0.6), "Max pearson correlation difference between a canopy center and a point in which the point will be considered close to the canopy. As a heuristc, only points within this distance will be considered as potential neighbours during the canopy walk.")
        ("max_merge_dist", value<PRECISIONT>(&max_merge_dist)->default_value(0.1), "Max pearson correlation difference between two canopy centers in which the canopies should be merged. Please note, that the final canopy profiles are calculated after the merge step and consequently some final canopies might have profiles that are closer then max_merge_dist specifies.")
        //This option is removed from CLI to avoid user confusion. The default value is hardcoded above
        //("min_step_dist", value<double>(&min_step_dist)->default_value(0.001), "Min pearson correlation difference between canopy center and canopy centroid in which the centroid will be used as an origin for a new canpy (canopy walk). This is a stop criterion for canopy walk.")
        ("profile_measure", value<string>(&profile_measure_str)->default_value("75Q"), "Speicfies gene abundance measure should the algorithm use. Valid options are: \"median\", \"mean\", \"75Q\", \"80Q\", \"85Q\", \"90Q\", \"95Q\" where \"XXQ\" stands for XXth quantile measure.");

    filter_in_options_desc.add_options()
        ("filter_min_obs", value<int>(&filter_min_obs)->default_value(3), "Discard those profiles which have fewer than N non-zero samples. Setting it to 0 will disable the filter.")
        ("filter_max_top3_sample_contribution", value<double>(&filter_max_top3_sample_contribution)->default_value(0.9), "Discard those profiles for which top 3 samples constitute more than X fraction of the total signal. Setting it to 1 will disable the filter")
        ("input_filter_file", value<string>(&input_filter_file)->default_value(""), "The file to which profiles filtered out by either of the input filters will be written");

    filter_out_options_desc.add_options()
        ("cag_filter_min_sample_obs", value<int>(&cag_filter_min_sample_obs)->default_value(3), "Return only those canopies that have at least N non-zero cluster profile observations. Setting it to 0 will disable the filter.")
        ("cag_filter_max_top3_sample_contribution", value<double>(&cag_filter_max_top3_sample_contribution)->default_value(0.9), "Don't return canopies where top three(or less) samples constitute more than X fraction of the total profile signal. Setting it to 1 disables the filter.");

    early_stop_options_desc.add_options()
        ("stop_criteria", value<double>(&stop_after_num_seeds_processed)->default_value(50000), "Stop clustering after X number of seeds have been processed. Setting it to 0 will disable this stop criterion.");

    misc_options_desc.add_options()
        ("die_on_kill", bool_switch(&die_on_kill), "If set, after receiving a KILL signal, the program will die and no results will be produced. By default clustering will stop but clusters will be merged and partial results will be printed as usual.")
        ("dont_use_mmap", bool_switch(&dont_use_mmap), "If set, the program will not attempt to read in the entire file into memory but read it line by line. It will be slower but will potentially save lot of RAM.")
        ("not_processed_profiles_file", value<string>(&not_processed_profiles_file)->default_value(""), "Path to file to which unprocessed profiles will be dumped at KILL signal")
        ("print_time_statistics,t", bool_switch(&print_time_statistics), "Print wall clock time profiles of various analysis parts. This is not aggressive and won't increase compuatation time.")
        ("show_progress_bar,b", bool_switch(&show_progress_bar), "Show progress bar, nice if output is printed to console, don't use if you are redirecting to a file. Verbosity must be set to at least PROGRESS for it to have an effect.") 
        ("dont_create_progress_stat_file", bool_switch(&dont_create_progress_stat_file), "If set, the canopy progress file will not be created.")
        ("progress_stat_file", value<string>(&progress_stat_file)->default_value("canopy_progress.out"), "Name of the canopy size statistics file. To this file current progress after each processed seed profile will be dumped in format <index> <num_profiles_left> <this_canopy_size> <total_num_thread_collisions>")
        ("help", "write help message");

    options_not_shown_in_help_desc.add_options()
        //This option used to be removed from CLI to avoid user confusion and hardcoded above
        ("max_num_canopy_walks", value<int>(&max_num_canopy_walks)->default_value(6), "Max number of times the canopy will walk. This is a stop criterion for canopy walk.");

    all_options_desc.add(general_options_desc).add(algorithm_param_options_desc).add(filter_in_options_desc).add(filter_out_options_desc).add(early_stop_options_desc).add(misc_options_desc).add(options_not_shown_in_help_desc);
    options_shown_in_help_desc.add(general_options_desc).add(algorithm_param_options_desc).add(filter_in_options_desc).add(filter_out_options_desc).add(early_stop_options_desc).add(misc_options_desc);

    positional_options_description command_line_positional_desc;
    command_line_positional_desc.add("input_file_path",1);
    command_line_positional_desc.add("output_clusters_file_path",1);
    command_line_positional_desc.add("output_cluster_profiles_file",1);

    variables_map command_line_variable_map;
    store(command_line_parser(argc,argv).options(all_options_desc).positional(command_line_positional_desc).run(), command_line_variable_map);
    notify(command_line_variable_map);

    //
    //Verify command line input parameters
    //
    //verify_input_correctness(all_options_desc, command_line_variable_map);
    if (command_line_variable_map.count("help") || argc < 3) {
        cout << "Usage: cc.bin [options] PROFILES_INPUT_FILE CLUSTERS_OUTPUT_FILE" << endl << endl;;
        cout << options_shown_in_help_desc << "\n";
        exit(0);
    }
	*/


    check_if_file_is_readable("input_file_path",input_file_path);
    if(priority_reads_file_path != "")
        check_if_file_is_readable("priority_reads_file_path", priority_reads_file_path);
	if (guide_matrix_file != "")
		check_if_file_is_readable("guide_matrix_file", guide_matrix_file);
    //check_if_file_is_writable("output_clusters_file_path",output_clusters_file_path);
    //check_if_file_is_writable("output_cluster_profiles_file",output_cluster_profiles_file);
    check_if_file_is_writable("input_filter_file",input_filter_file);
	const string arr2[] = { "error", "progress", "warn", "info", "debug", "debug1", "debug2", "debug3" };
	const vector<string> valid_verbosities(arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]));

    check_if_one_of("verbosity_option",verbosity_option, valid_verbosities);
    check_if_within_bounds("num_threads",num_threads,1,999);//Not exactly future proof, but let's put foolproofness first
    check_if_within_bounds("max_canopy_dist",max_canopy_dist,0.0,1.0);
    check_if_within_bounds("max_close_dist",max_close_dist,0.0,1.0);
    check_if_within_bounds("max_merge_dist",max_merge_dist,0.0,1.0);
    check_if_within_bounds("min_step_dist",min_step_dist,0.0,1.0);
    check_if_within_bounds("max_num_canopy_walks",max_num_canopy_walks,0,100);
    check_if_one_of("profile_measure", profile_measure_str, valid_profile_measure_values);

    check_if_within_bounds("filter_min_obs",opt->filter_min_obs,0,10000);
    check_if_within_bounds("filter_max_top3_sample_contribution", opt->filter_max_top3_sample_contribution,0.0,1.0);
    check_if_within_bounds("cag_filter_min_sample_obs",cag_filter_min_sample_obs,0,10000);
    check_if_within_bounds("cag_filter_max_top3_sample_contribution",cag_filter_max_top3_sample_contribution,0.0,1.0);
	//bool dont_create_progress_stat_file = opt->dont_create_progress_stat_file;
	bool create_progress_stat_file = !opt->dont_create_progress_stat_file;
    if(create_progress_stat_file)
        check_if_file_is_writable("progress_stat_file",progress_stat_file);
    if(not_processed_profiles_file!= "")
        check_if_file_is_writable("not_processed_profiles_file",not_processed_profiles_file);

    //
    //Set appropriate profile measure method to the global var (ugh..)
    //
    if(profile_measure_str == "median"){
        profile_measure = MEDIAN;
    } else if(profile_measure_str == "mean"){
        profile_measure = MEAN;
    } else if(profile_measure_str == "75Q"){
        profile_measure = PERCENTILE_75;
    } else if(profile_measure_str == "80Q"){
        profile_measure = PERCENTILE_80;
    } else if(profile_measure_str == "85Q"){
        profile_measure = PERCENTILE_85;
    } else if(profile_measure_str == "90Q"){
        profile_measure = PERCENTILE_90;
    } else if(profile_measure_str == "95Q"){
        profile_measure = PERCENTILE_95;
    } else {
        cout << "Unknown type of profile measure method: \"" << profile_measure_str << "\"" << endl;
        cout << "This is most likely a programming error, please report this bug" << endl;
        exit(1);
    }

    //
    //Set user chosen logging level
    //
    if(verbosity_option == "error"){
        log_level = logERR;
    }else if(verbosity_option == "progress"){
        log_level = logPROGRESS;
    }else if(verbosity_option == "warn"){
        log_level = logWARN;
    }else if(verbosity_option == "info"){
        log_level = logINFO;
    }else if(verbosity_option == "debug"){
        log_level = logDEBUG;
    }else if(verbosity_option == "debug1"){
        log_level = logDEBUG1;
    }else if(verbosity_option == "debug2"){
        log_level = logDEBUG2;
    }else if(verbosity_option == "debug3"){
        log_level = logDEBUG3;
    }

    _log(logINFO) << "";
    _log(logINFO) << "cc.bin version " << ccbinVer;  
    _log(logINFO) << "Files:";
    _log(logINFO) << "input_file_path:\t " << input_file_path;
    _log(logINFO) << "priority_reads_file_path:\t" << priority_reads_file_path;
    _log(logINFO) << "output_cluster_profiles_file:\t " << output_cluster_profiles_file;
    _log(logINFO) << "progress_stat_file:\t " << progress_stat_file;
    _log(logINFO) << "not_processed_profiles_file:\t " << not_processed_profiles_file;
    _log(logINFO) << "input_filter_file:\t " << input_filter_file;
    _log(logINFO) << "";
    

    //Set signal handler
    if(die_on_kill) 
        signal(SIGINT, signal_callback_die_handler);
    else    
        signal(SIGINT, signal_callback_gentle_handler);

    //Set number of threads
    _log(logINFO) << "";
    _log(logINFO) << "General:";
    _log(logINFO) << "num_threads:\t " << num_threads;
   // _log(logINFO) << "precision_type:\t " << boost::typeindex::type_id<PRECISIONT>().pretty_name();
    _log(logINFO) << "";

    omp_set_num_threads(num_threads);
	bool sparseMat = opt->sparseMat;

    //
    //Parse priority point name file
    //
    vector<string> priority_read_names;
    if(priority_reads_file_path != ""){
		time_profile.start_timer("Loading priority reads");
		ifstream priority_reads_file (priority_reads_file_path);
        if(priority_reads_file.is_open())
        {
            string line;
            while(getline(priority_reads_file, line))
            {
                if(line.length()) {
                    priority_read_names.push_back(line);
                }
            }
            priority_reads_file.close();
        }
        _log(logINFO) << "";
        _log(logINFO) << "numbers of read priority reads:\t " << priority_read_names.size();
        _log(logINFO) << "";
        time_profile.stop_timer("Loading priority reads");
    }
	vector<Point*> guidePoints;
	if (guide_matrix_file != "") {
		time_profile.start_timer("Guide matrix");
		_log(logINFO) << "Reading guide matrix line by line";

		std::ifstream point_file(guide_matrix_file);
		std::string line;

		while (std::getline(point_file, line)) {
			if (line.length() < 2)
				break;
			Point * pp = new Point(line.c_str(), sparseMat);
			if (use_spearman) {
				pp->convert_to_rank();
			}
			pp->seal();
			guidePoints.push_back(pp);
			die_if_true(terminate_called);
		}
		point_file.close();


		_log(logINFO) << "";
		_log(logINFO) << "numbers entries in guide matrix:\t " << guidePoints.size();
		_log(logINFO) << "";
		time_profile.stop_timer("Guide matrix");
	}

    //
    //Parse point description file
    //

    vector<Point*> filtered_points;

//read input files..
    _log(logINFO) << "Reading file line by line";
    time_profile.start_timer("Loading file and reading profiles");

	//heavy IO routine
	vector<Point*> points(0);
	vector<PRECISIONT> sampleSums(0);
	readMatrix(points, sampleSums, input_file_path, sparseMat, use_spearman, num_threads);
	if (points.size() <= 1) {
		cerr << "Matrix is empty, aborting";
		exit(0);
	}

    time_profile.stop_timer("Loading file and reading profiles"); _log(logINFO) << "";
	readMB2preSet(opt, guidePoints, points);
	
	//rm autocorrelated samples..
	bool autocorr_sample_filter = true;
	if (guidePoints.size() > 0 || opt->sampleMinDist >= 2){
		autocorr_sample_filter = false;
	}
	vector<bool> rmSmpls; int sumRm(0);
	if (autocorr_sample_filter) {
		rmSmpls = autocorr_filter(opt, time_profile,points, sampleSums);
		sumRm = handleRms(points, rmSmpls);
	}

	filter(opt, time_profile, points, guidePoints, filtered_points);
    //Do not use "points" or points_filtered_out_due_to_three_point_proportion_filter or points_filtered_out_due_to_num_non_zero_samples_filter at this point
    
    
    die_if_true(terminate_called);
    die_if_true(filtered_points.size() < 1);

    //
    //This will precompute values for quicker pearson correlation calculation
    //This is a bit clumsy but helps prevent huge memory spikes
    //
	if (use_spearman)
		_log(logINFO) << "Precomputing Spearman correlation data to speed up distance calculations";
	else
		_log(logINFO) << "Precomputing Pearson correlation data to speed up distance calculations";

	time_profile.start_timer("Precomputing correlation data");

#pragma omp parallel for shared(filtered_points)
	for (size_t i = 0; i < filtered_points.size(); i++) {
		filtered_points[i]->allocate_and_precompute_pearson_data();
	}

#pragma omp parallel for shared(guidePoints)
	for (size_t i = 0; i < guidePoints.size(); i++) {
		guidePoints[i]->allocate_and_precompute_pearson_data();
	}


	time_profile.stop_timer("Precomputing correlation data");
    
    die_if_true(terminate_called);


    //
    //Run Canopy Clustering
    //
    std::vector<shared_ptr<Canopy>> canopies(guidePoints.size(), nullptr);
    bool guided = guidePoints.size() > 0;

	if (guided) {
		_log(logINFO) << "";
		_log(logINFO) << "Calculating genes correlationg to guide profiles";

		 multi_core_run_correlations(filtered_points, guidePoints, canopies,
			num_threads, max_canopy_dist,
			show_progress_bar, time_profile,false);
		_log(logINFO) << "Finished deep correlations";

        if (filter_redundant) {
            _log(logINFO) << "Filtering redundantly assigned genes";
            filter_redundant_genes(canopies, guidePoints, time_profile);
            _log(logINFO) << "Finished filtering redundantly assigned genes";
        }

	}
	else {
		_log(logINFO) << "";
		_log(logINFO) << "Calculating new canopies on gene matrix";
		canopies = multi_core_run_clustering_on(filtered_points,
			priority_read_names, opt, rmSmpls, sumRm,
			time_profile);
		_log(logINFO) << "Finished clustering";
	}


    //
    //Filter out canopies
    //

	if (!guided) {
		if (cag_filter_min_sample_obs) {
			time_profile.start_timer("Filtering canopies by minimum number of sample detections");
			filter_clusters_by_zero_medians(cag_filter_min_sample_obs, canopies);
			_log(logINFO) << "Finished filtering for minimum number of sample detections, number of canopies:" << canopies.size();
			time_profile.stop_timer("Filtering canopies by minimum number of sample detections");
		}


		if (cag_filter_max_top3_sample_contribution < 0.99999) { //It's due to a double comparison
			time_profile.start_timer("Filtering canopies by three sample signal contribution proportion");
			cag_filter_max_top3_sample_contributions(cag_filter_max_top3_sample_contribution, canopies);
			_log(logINFO) << "Finished filtering by three sample signal contribution proportion, number of canopies:" << canopies.size();
			time_profile.stop_timer("Filtering canopies by three sample signal contribution proportion");
		}

		{
			time_profile.start_timer("Filtering canopies by size");
			filter_clusters_by_size(canopies);
			_log(logINFO) << "Finished filtering by size(number of neighbours must be bigger than 1), number of canopies:" << canopies.size();
			time_profile.stop_timer("Filtering canopies by size");
		}
	}

	//from here on everything done wrt to clustering, restore the original size (samples)
	if (sumRm > 0) {
		_log(logPROGRESS) << "Expanding canopies to original sample size";
		for (uint i = 0; i < canopies.size(); i++) {
			canopies[i]->restore_rm(sumRm);
		}
	}

    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "#################### Writing Results ####################" ;
	ofstream* out_file(NULL); ofstream* out_file2(NULL);

    int num_digits = (int)ceil(log10(canopies.size()));
    //cout << std::setfill('0');


    out_file = new ofstream(output_clusters_file_path.c_str(), ios::out | ios::trunc);
	if (output_cluster_profiles_file != "") {
		out_file2 = new ofstream(output_cluster_profiles_file.c_str(), ios::out | ios::trunc);
	}
	if (!guided) {
		sort(canopies.begin(), canopies.end(), compare_canopy_ptrs_by_canopy_size);
	}
	for (int i = 0; i < (int)canopies.size(); i++) {
		//shared_ptr<Canopy> c = canopies[i];
		//cerr << "at cano " << i << "of" << canopies.size() << endl;
		if (canopies[i] == nullptr) { cerr << "Detected null pointer!! at " << i << " !!" << endl; continue; }
		canopies[i]->print2file(out_file, out_file2,opt,i, num_digits, guided);
	}
	(*out_file).close();
	_log(logPROGRESS) << "#################### Finished writing canopies ####################";

	if (output_cluster_profiles_file != "") {
		(*out_file2).close();
	}


	//partial correlations
	if (output_clusters_partial_file_path != "") {
        
		_log(logINFO) << "";
		_log(logINFO) << "Calculating genes PARTIALLY correlationg to guide profiles";
        std::vector<shared_ptr<Canopy>> canopies_par(guidePoints.size());
		 multi_core_run_correlations(filtered_points, guidePoints, canopies_par,
			num_threads, max_canopy_dist_part, show_progress_bar, time_profile,true);
	
		ofstream *OF;
		OF = new ofstream(output_clusters_partial_file_path.c_str(), ios::out | ios::trunc);
		_log(logPROGRESS) << "";
		_log(logPROGRESS) << "#################### Writing Results of partial correlations ####################";

		for (int i = 0; i < (int) canopies_par.size(); i++) {
			if (canopies_par[i] == nullptr) { cerr << "Detected null pointer par !! at " << i << " !!" << endl; continue; }
			canopies_par[i]->print2file(OF, NULL, opt, i, num_digits, true);
		}
		(*OF).close();
			//
 //Clean up
 //

	}

 

    time_profile.stop_timer("Total");
    //Write output statistics
    if(print_time_statistics){
        cout << time_profile << endl;
    }


    return 0;
}
