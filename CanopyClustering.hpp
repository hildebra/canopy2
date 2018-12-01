/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
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
#ifndef CANOPY_CLUSTERING
#define CANOPY_CLUSTERING


#include "Canopy.hpp"

/*Added to just run clustering on guide matrix, which are filtered canopies*/		
std::vector<shared_ptr<Canopy>> multi_core_run_correlations(vector< Point*>& points,
	vector< Point*>& guides,
		int num_threads, PRECISIONT max_canopy_dist,
		bool show_progress_bar, TimeProfile& time_profile,bool partial=false);

/**
    * Run the canopy clustering algorithm
    *
    * Parameters:
    * points - list of references to points to be clustered
	* guidePoints - if given only cluster to predefined set
    * num_threads - number of threads to be used in the calculation
    * max_canopy_dist, max_close_dist, max_merge_dist, max_num_canopy_walks, stop_after_num_seeds_processed, show_progress_bar - see program parameters description 
    * create_canopy_size_stats - boolean defining whether the progress file should be created
    * canopy_size_stats_fp - absolute file path to the canopy size statistics file
    * not_procesed_point_fp - absolute file path to the file containing not processed points if early stopped
    * time_profile - TimeProfile object instance for gathering statistics on time it took for each of the analysis steps
	*/
std::vector<shared_ptr<Canopy>> multi_core_run_clustering_on(vector< Point*>& points,
	vector<string>& priority_read_names, options* opt, const vector<bool> rm, int, TimeProfile& time_profile);

/**
    * Create canopy given an origin point
    *
    * Parameters:
    * origin - point being the canopy origin
    * points - list of all points
    * close_points - list of close_points
    * min_neighbour_correlation - minimum distance in correlation space between origin and a tested point for the point to be considered inside the canopy
    * min_close_correlation - minimum distance in correlation space between origin and a tested point for the point to be considered close to the canopy
    * sets_close_points - flag describing if the current execution of this function should set the close_points
    */
shared_ptr<Canopy> create_canopy(Point* origin, vector< Point*>& points, 
	vector< Point*>& close_points,
	PRECISIONT min_neighbour_correlation, PRECISIONT min_close_correlation, 
	bool sets_close_points,int);

shared_ptr<Canopy> create_canopy_singl(Point* origin, const vector< Point*>& points,
	PRECISIONT min_neighbour_correlation,bool partial=false, int origin_i=0);
		
/**
    * Execute the create_canopy function iteratively until a stable canopy is reached. 
    *
    * Parameters:
    * origin - point which is the canopy origin
    * points - list of all points
    * close_points - list of points which are within "close" distance to canopy
    * max_canopy_dist, max_close_dist, min_step_dist, max_num_canopy_walks  - see program parameters
    * num_canopy_jumps - number of times the create_canopy function was executed
    */
shared_ptr<Canopy> canopy_walk(Point* origin, vector< Point*>& points,
	vector< Point*>& close_points,
	PRECISIONT max_canopy_dist, PRECISIONT max_close_dist, PRECISIONT min_step_dist, 
	PRECISIONT max_num_canopy_walks, int& num_canopy_jumps,
	int deletedSmpls);

void shuffle_points(vector< Point*>& points, vector<string>& priority_read_names);
void filter_clusters_by_zero_medians(int min_num_non_zero_medians, vector<shared_ptr<Canopy>>& canopies_to_filter);
void cag_filter_max_top3_sample_contributions(PRECISIONT max_single_data_point_proportion, 
	vector<shared_ptr<Canopy>>& canopies_to_filter);
void filter_clusters_by_size(std::vector<shared_ptr<Canopy>>& canopies_to_filter);

void filter(options * opt, TimeProfile time_profile, vector<Point*>& points,
	vector<Point*>& guidePoints, vector<Point*>& filtered_points);

void writeMatrix(vector<vector<PRECISIONT>>& mat, string of);

vector<bool> autocorr_filter(options * opt, TimeProfile time_profile, 
	const vector<Point*>& points, const vector<PRECISIONT>& sampleSums);
int handleRms(vector<Point*>& points, const vector<bool>& rm);



#endif 
