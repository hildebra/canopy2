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
#include <iostream>
#include <fstream>
#include <ctime>

#include <algorithm>
#include <unordered_set>
#include <map>


#include <omp.h>

#include "Log.hpp"

#include "signal_handlers.hpp"
#include "prog_bar_misc.hpp"
#include "Canopy.hpp"
#include "Point.hpp"
#include "CanopyClustering.hpp"

shared_ptr<Canopy> create_canopy(Point* origin, vector< Point*>& points,
	vector< Point*>& close_points, PRECISIONT max_neighbour_dist, PRECISIONT max_close_dist,
	bool set_close_points, int deletedSmpls) {

	std::vector<Point*> neighbours;

	if (set_close_points) {
		Point* potential_neighbour;

		//Go through all points and set the close points to contain the ones that are "close"
		close_points.clear();//Will not reallocate
		for (size_t i = 0; i < points.size(); i++) {

			potential_neighbour = points[i];

			PRECISIONT dist = get_distance_between_points(origin, potential_neighbour);

			if (dist < max_close_dist) {

				close_points.push_back(potential_neighbour);

				if (dist < max_neighbour_dist) {

					neighbours.push_back(potential_neighbour);

				}
			}
		}

	}
	else {

		Point* potential_neighbour;

		for (size_t i = 0; i < close_points.size(); i++) {

			potential_neighbour = close_points[i];

			PRECISIONT dist = get_distance_between_points(origin, potential_neighbour);

			if (dist < max_neighbour_dist) {
				neighbours.push_back(potential_neighbour);
			}
		}

	}
	if (neighbours.size()) {
		return make_shared<Canopy>(neighbours, deletedSmpls);
	}
	else {
		return make_shared<Canopy>(origin, deletedSmpls);
	}
}

shared_ptr<Canopy> create_canopy_singl(Point* origin, const vector< Point*>& points,
	PRECISIONT max_neighbour_dist, bool partial, int origin_i) {

	std::vector< Point*> neighbours;
	list<PRECISIONT> corrs;

	Point* potential_neighbour;
	for (size_t  i = 0; i < points.size(); i++) {
		potential_neighbour = points[i];
		//cout << points[i]->id << " ";
		PRECISIONT dist;
		if (partial) {
			dist = get_partial_distance_between_points(origin, potential_neighbour);
		} else {
			dist = get_distance_between_points(origin, potential_neighbour);
		}

		if (dist < max_neighbour_dist) {
			neighbours.push_back(potential_neighbour);
			corrs.push_back(dist);
		}
	}

	
	shared_ptr<Canopy> cc = make_shared<Canopy>(origin,false); // first copy ori!
	cc->set_ori(origin_i);
	if (neighbours.size()>0) {
		//now set neighbor list
		cc->neighbours = neighbours;
		cc->corrs = corrs;
	}
	return cc;
}

shared_ptr<Canopy> canopy_walk(Point* origin, vector< Point*>& points,
	vector< Point*>& close_points, PRECISIONT max_canopy_dist, PRECISIONT max_close_dist,
	PRECISIONT min_step_dist, PRECISIONT max_num_canopy_walks, int& num_canopy_jumps,
	int deletedSmpls){

	shared_ptr<Canopy> c1;
	shared_ptr<Canopy> c2;


    c1 = create_canopy(origin, points, close_points, max_canopy_dist, max_close_dist, true,
		deletedSmpls);

    //special case for which there is no walking, return the canopy immediatelly
    if(max_num_canopy_walks == 0){
        return c1;                                                               
    }

    c2 = create_canopy(c1->center, points, close_points, max_canopy_dist, max_close_dist, false,
		deletedSmpls);
    
    PRECISIONT dist = get_distance_between_points(c1->center, c2->center);

    {
        _log(logDEBUG2) << "Canopy walking, first step";
        _log(logDEBUG2) << *c1;
        _log(logDEBUG2) << *c2;
        _log(logDEBUG2) << "First potential jump correlation dist: " << dist;
    }

    int num_canopy_jumps_local = 0;
    while((dist > min_step_dist) && (num_canopy_jumps_local <= max_num_canopy_walks )){
        c1=c2;

        num_canopy_jumps_local++; 

#pragma omp atomic
        num_canopy_jumps++;

        c2=create_canopy(c1->center, points, close_points, max_canopy_dist, max_close_dist, false,
			deletedSmpls);
        dist = get_distance_between_points(c1->center, c2->center); 
        {
            _log(logDEBUG2) << "Canopy walking, step: " << num_canopy_jumps_local;
            _log(logDEBUG2) << *c1;
            _log(logDEBUG2) << *c2;
            _log(logDEBUG2) << "distance: " << dist;
        }
    }

    //Now we know that c1 and c2 are close enough and we should choose the one that has more neighbours
	shared_ptr<Canopy> final_canopy;
    if(c1->neighbours.size() > c2->neighbours.size()){
        final_canopy = c1;
    } else {
        final_canopy = c2;
    }
    return final_canopy;
}

void filter_clusters_by_size(std::vector<shared_ptr<Canopy>>& canopies_to_filter){

    vector<int> canopy_indexes_to_remove;

    for(size_t  i=0; i < canopies_to_filter.size(); i++){
		shared_ptr<Canopy> canopy = canopies_to_filter[i];
        if(canopy->neighbours.size() < 2){
			delete canopy->center;
			canopy->center = NULL;
			canopy_indexes_to_remove.push_back(i);
        }
    }

    std::sort(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());
    std::reverse(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());

    for(size_t  i=0; i < canopy_indexes_to_remove.size(); i++)
        canopies_to_filter.erase(canopies_to_filter.begin() + canopy_indexes_to_remove[i]);

}

void cag_filter_max_top3_sample_contributions(PRECISIONT max_single_data_point_proportion,
	vector<shared_ptr<Canopy>>& canopies_to_filter){

    vector<int> canopy_indexes_to_remove;

    for(size_t  i=0; i < canopies_to_filter.size(); i++){
        Point* ccenter = canopies_to_filter[i]->center;
        //if(! ccenter->check_if_single_point_proportion_is_smaller_than(max_single_data_point_proportion) ){
        if(! ccenter->check_if_top_three_point_proportion_is_smaller_than(max_single_data_point_proportion) ){
			delete ccenter;
			canopies_to_filter[i]->center = NULL;
			canopy_indexes_to_remove.push_back(i);
        }
    }

    std::sort(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());
    std::reverse(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());

    for(size_t  i=0; i < canopy_indexes_to_remove.size(); i++)
        canopies_to_filter.erase(canopies_to_filter.begin() + canopy_indexes_to_remove[i]);

}

void filter_clusters_by_zero_medians(int min_num_non_zero_medians, 
	std::vector<shared_ptr<Canopy>>& canopies_to_filter){

    vector<int> canopy_indexes_to_remove;

    for(size_t  i=0; i < canopies_to_filter.size(); i++){
		Point* ccenter = canopies_to_filter[i]->center;
        if(! ccenter->check_if_num_non_zero_samples_is_greater_than_x(min_num_non_zero_medians) ){
            canopy_indexes_to_remove.push_back(i);
        }
    }

    std::sort(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());
    std::reverse(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());

    for(size_t  i=0; i < canopy_indexes_to_remove.size(); i++)
        canopies_to_filter.erase(canopies_to_filter.begin() + canopy_indexes_to_remove[i]);

}



void shuffle_points(vector< Point*>& points, vector<string>& priority_read_names){

    //Copy the read names from the vector into a map, we will be checking if the reads from the input file are in it and which position they should have
    std::map<string, int> priority_read_name__to_position;
    for(size_t  i=0; i<priority_read_names.size(); i++)
        priority_read_name__to_position[priority_read_names[i]] = i;

    //Sort the points vector so that those reads that are in priority_read_names come first and in the order of the priority_read_names
    sort(points.begin(), points.end(), [&priority_read_name__to_position](const Point* p1, const Point* p2) -> bool {
        //integers are compared through "<" to get ascending sort
        auto p1_map_it = priority_read_name__to_position.find(p1->id);
        auto p2_map_it = priority_read_name__to_position.find(p2->id);

        //If both are in the priority read map - then compare integers directly
        if((p1_map_it != priority_read_name__to_position.end()) && (p2_map_it != priority_read_name__to_position.end()))
            return p1_map_it->second < p2_map_it->second;
        //If first is in priority read map (we know the second must not be there) then the first one goes before the other
        else if(p1_map_it != priority_read_name__to_position.end())
            return true;
        else
            return false; //That is the second read was in the priority read map and first wasn't, then second one goes before the first one
    });

    //Go through all points and find the first point that is NOT in the priority read names (using the map)
    //This is a somewhat clever approach: find_if returns "An iterator to the first element in the range for which pred does not return false."
    auto first_non_priority_pint_it = find_if(points.begin(), points.end(), [&priority_read_name__to_position](const Point* p) -> bool {
        return priority_read_name__to_position.find(p->id) == priority_read_name__to_position.end();
    });

    //Now shuffle the non prioritized pionts
    std::random_shuffle(first_non_priority_pint_it, points.end());

}

std::vector<shared_ptr<Canopy>> multi_core_run_clustering_on(vector< Point*>& points,
	vector<string>& priority_read_names, options* opt, 
	const vector<bool> rm, int deletedSmpls,
	TimeProfile& time_profile){

	//establish parameters
	int num_threads = opt->num_threads;
	PRECISIONT max_canopy_dist = opt->max_canopy_dist;
	PRECISIONT max_close_dist = opt->max_close_dist;
	PRECISIONT max_merge_dist = opt->max_merge_dist;
	PRECISIONT min_step_dist = opt->min_step_dist;
	int max_num_canopy_walks = opt->max_num_canopy_walks;
	PRECISIONT stop_after_num_seeds_processed = opt->stop_after_num_seeds_processed;
	bool create_canopy_size_stats  = !opt->dont_create_progress_stat_file;
	bool show_progress_bar = opt->show_progress_bar;
	string canopy_size_stats_fp = opt->progress_stat_file;
	string not_processed_points_fp = opt->not_processed_profiles_file;

    _log(logINFO) << "";
    _log(logINFO) << "Algorithm Parameters:";
    _log(logINFO) << "max_canopy_dist:\t " << max_canopy_dist;
    _log(logINFO) << "max_close_dist:\t " << max_close_dist;
    _log(logINFO) << "max_merge_dist:\t " << max_merge_dist;
    _log(logINFO) << "min_step_dist:\t " << min_step_dist;
    _log(logINFO) << "max_num_canopy_walks:\t " << max_num_canopy_walks;
    _log(logINFO) << "";
    _log(logINFO) << "Early stopping:";
    _log(logINFO) << "stop_after_num_seeds_processed:\t " << stop_after_num_seeds_processed;
    _log(logINFO) << "";
    _log(logINFO) << "Priority reads";
    _log(logINFO) << "Number of reads with first priority:\t" << priority_read_names.size();

    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "############ Shuffling ############";
    time_profile.start_timer("Shuffling");
    shuffle_points(points, priority_read_names);
    time_profile.stop_timer("Shuffling");

    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "############ Creating Canopies ############";
    _log(logPROGRESS) << "To make the program stop and generate output send an interrupt signal by either:" ;
    _log(logPROGRESS) << "\t * running command \"kill -INT [ cc.bin PID ]\"";
    _log(logPROGRESS) << "\t * pressing \"CTRL + C\" in this terminal";
    _log(logPROGRESS) << "";

#pragma clang diagnostic push
#pragma ide diagnostic ignored "TemplateArgumentsIssues" //Clion is messing up, the set declaration is fine
    std::unordered_set< Point*> marked_points;//Points that should not be investigated as origins
#pragma clang diagnostic pop
    vector<unsigned int> canopy_size_per_origin_num;//Contains size of the canopy created from origin by it's number, so first origin gave canopy of size 5, second origin gave canopy of size 8 and so on
    int last_progress_displayed_at_num_points = 0;

    std::vector<shared_ptr<Canopy>> canopy_vector;

    vector< Point*> close_points;
    //close_points.reserve(points.size());

    //
    //Create canopies
    //
    time_profile.start_timer("Clustering");
        
    ofstream canopy_size_stats_file;

    if(create_canopy_size_stats)
        canopy_size_stats_file.open(canopy_size_stats_fp.c_str(), ios::out | ios::trunc);

    int canopy_stats_row_num = 0;

    int num_canopy_jumps = 0;
    int num_collisions = 0;
    int num_seeds_processed = 0;

    size_t first_non_processed_origin_due_interruption = points.size();

    //Disable stop criterion if set to zero
    if(stop_after_num_seeds_processed == 0){
        stop_after_num_seeds_processed = points.size();
    }

#pragma omp parallel for shared(marked_points, canopy_vector, num_canopy_jumps, canopy_size_per_origin_num, num_collisions, num_seeds_processed, canopy_stats_row_num, terminate_called, first_non_processed_origin_due_interruption) firstprivate(close_points, max_canopy_dist, max_close_dist, max_merge_dist, min_step_dist, last_progress_displayed_at_num_points) schedule(dynamic)
    for(size_t  origin_i = 0; origin_i < points.size(); origin_i++){
        //Early stopping after num of points
        if(num_seeds_processed > stop_after_num_seeds_processed){
            continue;
        }

        //Stop if exit signal received
        if(terminate_called){
            if(first_non_processed_origin_due_interruption > origin_i){
#pragma omp critical
                {
                    first_non_processed_origin_due_interruption = origin_i;
                }
            }
            continue;
        }

        Point* origin = points[origin_i]; 
		if (!origin->precomputed) { cerr << "Precomputation missing!\n"; exit(723); }

        if(marked_points.find(origin) != marked_points.end())
            continue;

        //Show progress bar
        {
            //Only master thread executes this
            if(omp_get_thread_num() == 0){
                if(log_level >= logPROGRESS && show_progress_bar){
                    if(marked_points.size() > last_progress_displayed_at_num_points + stop_after_num_seeds_processed/100){
                        printProgBar(marked_points.size(),(int)stop_after_num_seeds_processed * points.size());
                        last_progress_displayed_at_num_points = (int)stop_after_num_seeds_processed;
                    }
                }
            }
        }




        {
            _log(logDEBUG) << "Unmarked points count: " << points.size() - marked_points.size() << " Marked points count: " << marked_points.size();
            _log(logDEBUG) << "points.size: " << points.size() << " origin_i: " << origin_i << " origin->id: " << origin->id ;

            _log(logDEBUG1) << "Current canopy origin: " << origin->id;
        }
		//keepProfileStable
		shared_ptr<Canopy> final_canopy = canopy_walk(origin, points, close_points, max_canopy_dist,
			max_close_dist, min_step_dist, max_num_canopy_walks, num_canopy_jumps, deletedSmpls);

#pragma omp critical
        {
            //Do not commit anything if by chance another thread marked the current origin
            if(marked_points.find(origin) == marked_points.end()){

                //Add canopy
                marked_points.insert(origin);

                canopy_vector.push_back(final_canopy);

                for(Point* n : final_canopy->neighbours){
                    marked_points.insert(n);
                }

                //Statistics showing size of canopies per analyzed origin
                if(canopy_size_stats_file.is_open()){
                    canopy_size_stats_file << canopy_stats_row_num++  << "\t" << points.size() - marked_points.size() << "\t" << final_canopy->neighbours.size() << "\t" << num_collisions << endl;
                }

            } else {
                num_collisions++;
            }

            num_seeds_processed += 1;

        }

    }
    if(canopy_size_stats_file.is_open())
        canopy_size_stats_file.close();

    time_profile.stop_timer("Clustering");
    
    if(terminate_called && (not_processed_points_fp != "")){
        time_profile.start_timer("Saving unprocessed points file");
        _log(logERR) << "Received signal, clustering was stopped early, saving non processed points in file: " << not_processed_points_fp; 
        cout << "first_non_processed_origin_due_interruption:" << first_non_processed_origin_due_interruption << endl; 

        ofstream not_processed_points_file;
        not_processed_points_file.open(not_processed_points_fp.c_str(), ios::out | ios::trunc);
        for(size_t  i = first_non_processed_origin_due_interruption; i < points.size(); i++){
			Point* point = points[i];
            if(marked_points.find(point) == marked_points.end()){
                not_processed_points_file << point->id;
                for(int j = 0; j < point->num_data_samples; j++){
                    not_processed_points_file << "\t" << point->sample_data[j];
                }
                not_processed_points_file << "\n";
            }
        }
        not_processed_points_file.close();
        _log(logERR) << "Unprocessed points saved"; 
        time_profile.stop_timer("Saving unprocessed points file");
    }


    _log(logINFO) << "";
    _log(logINFO) << "Avg. number of canopy walks: " << num_canopy_jumps/((PRECISIONT)canopy_vector.size());
    _log(logINFO) << "Number of all canopies before merging: " << canopy_vector.size();

    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "############Removing canopies of size 1 to speed-up merging#############";
    int num_single_sample_canopies = 0;
    for(size_t  i=0; i < canopy_vector.size(); ){
        if(canopy_vector[i]->neighbours.size() == 1){
            canopy_vector.erase(canopy_vector.begin() + i);
            num_single_sample_canopies++;
        } else {
            i++;
        }
    }

    _log(logINFO) << "";
    _log(logINFO) << "Number of canopies which are removed due to having only 1 sample: " << num_single_sample_canopies;
    _log(logINFO) << "Number of canopies left after removal of single sample canpies: " << canopy_vector.size();

    int original_number_of_canopies = canopy_vector.size();

    //
    // Merge Canopies
    //
    std::vector<shared_ptr<Canopy>> merged_canopy_vector;

    time_profile.start_timer("Merging");
    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "############Merging canopies#############";


    //Actual merge 
    while(canopy_vector.size()){

        std::vector<shared_ptr<Canopy>> canopies_to_merge;

        //This is the canopy we will look for partners for
		shared_ptr<Canopy> c = *canopy_vector.rbegin();
        canopy_vector.pop_back();

        canopies_to_merge.push_back(c);

        //Get indexes of those canopies that are nearby
//#pragma omp parallel for shared(canopies_to_merge) 
        for(size_t  i = 0; i < canopy_vector.size(); i++){

			shared_ptr<Canopy> c2 = canopy_vector[i];

            PRECISIONT dist = get_distance_between_points(c->center, c2->center);

            if(dist < max_merge_dist){
//#pragma omp critical
                {
                    canopies_to_merge.push_back(c2);
                }
            }

        }

        //There are several canopies to merge, do it
        if( canopies_to_merge.size() > 1 ){

            vector< Point*> all_points_from_merged_canopies;
            
            for(shared_ptr<Canopy> canopy : canopies_to_merge){
                for(Point* n : canopy->neighbours){
                    if(std::find(all_points_from_merged_canopies.begin(), all_points_from_merged_canopies.end(), n) == all_points_from_merged_canopies.end()){ //If the element hasn't been added already
                        all_points_from_merged_canopies.push_back(n);
                    }
                }
            }

            //Create new canopy
			Point* temp_merged_canopy_origin = get_centroid_of_points(all_points_from_merged_canopies,
				deletedSmpls);
			shared_ptr<Canopy> merged_canopy = canopy_walk(temp_merged_canopy_origin, all_points_from_merged_canopies,
				close_points, max_canopy_dist, max_close_dist, min_step_dist, max_num_canopy_walks, 
				num_canopy_jumps, deletedSmpls);


            canopy_vector.push_back(merged_canopy);

            
            //Removed merged canopies //TODO might be slow
            for(shared_ptr<Canopy> canopy : canopies_to_merge){
                canopy_vector.erase(remove(canopy_vector.begin(), canopy_vector.end(), canopy), canopy_vector.end());
            }


        //If no canopies were merged remove the canopy we compared against the others
        } else {

            merged_canopy_vector.push_back(c);

            //Show progress bar
            {
                if(log_level >= logPROGRESS && show_progress_bar){
                    if(original_number_of_canopies - canopy_vector.size() % 1000)
                        printProgBar(original_number_of_canopies - canopy_vector.size(), original_number_of_canopies );
                }
            }
        }
    }
    time_profile.stop_timer("Merging");

    _log(logINFO) << "";
    _log(logINFO) << "Number of canopies after merging: " << merged_canopy_vector.size();

    return merged_canopy_vector;

}

std::vector<shared_ptr<Canopy>> multi_core_run_correlations(vector< Point*>& points,
	vector< Point*>& guides,
	int num_threads, PRECISIONT max_canopy_dist,
	bool show_progress_bar, TimeProfile& time_profile, bool partial) {

	_log(logINFO) << "";
	_log(logINFO) << "Algorithm Parameters:";
	_log(logINFO) << "max_canopy_dist:\t " << max_canopy_dist;
	_log(logINFO) << "Guide profile:\t" << guides.size();

	_log(logPROGRESS) << "";
	_log(logPROGRESS) << "############ Deep Canopies ############";
	_log(logPROGRESS) << "";

	vector<unsigned int> canopy_size_per_origin_num;//Contains size of the canopy created from origin by it's number, so first origin gave canopy of size 5, second origin gave canopy of size 8 and so on
	int last_progress_displayed_at_num_points = 0;

	std::vector<shared_ptr<Canopy>> canopy_vector(guides.size());

	//
	//Create canopies
	//
	time_profile.start_timer("ReClustering");

	ofstream canopy_size_stats_file;


	int totalSteps = guides.size();
	//thread pool
	vector < job > slots(num_threads);
	int j = 0; int origin_i = 0;
	while (true) {
		//show progress
		if (log_level >= logPROGRESS && show_progress_bar) {
			if (origin_i > last_progress_displayed_at_num_points + totalSteps / 100) {
				printProgBar(origin_i, totalSteps);
				last_progress_displayed_at_num_points = origin_i;
			}
		}
		//actual correlations..
		if (j >= num_threads) { j = 0; }
		if (origin_i >= totalSteps) { break; }
		if (slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
			slots[j].inUse = false;
			shared_ptr<Canopy> cc = slots[j].fut.get();
			canopy_vector[cc->get_ori()] = cc;
		}
		if (slots[j].inUse == false) {
			//create_canopy_singl(origin, points, max_canopy_dist, partial)
			slots[j].fut = async(std::launch::async, create_canopy_singl, guides[origin_i], ref(points), max_canopy_dist, partial, origin_i);
			origin_i++;
			slots[j].inUse = true;
		}
		j++;
	}
	for (int j = 0; j < num_threads; j++) {
		if (slots[j].inUse == true ) {
			shared_ptr<Canopy> cc = slots[j].fut.get();
			canopy_vector[cc->get_ori()] = cc;
			slots[j].inUse = false;
		}

	}
	//get remaining jobs


	/*
#pragma omp parallel for shared( canopy_vector) firstprivate( totalSteps, max_canopy_dist, max_close_dist, max_merge_dist, min_step_dist, last_progress_displayed_at_num_points) schedule(dynamic)
	for (int origin_i = 0; origin_i < totalSteps; origin_i++) {


		Point* origin = guides[origin_i];

		//Show progress bar
		{
			//Only master thread executes this
			if (omp_get_thread_num() == 0) {
				if (log_level >= logPROGRESS && show_progress_bar) {
					if (origin_i > last_progress_displayed_at_num_points + totalSteps / 100) {
						printProgBar(origin_i, totalSteps);
						last_progress_displayed_at_num_points = origin_i;
					}
				}
			}
		}




		{
			_log(logDEBUG) << "points.size: " << points.size() << " origin_i: " << origin_i << " origin->id: " << origin->id;
			_log(logDEBUG1) << "Current canopy origin: " << origin->id;
		}
		//keepProfileStable
		Canopy* final_canopy = create_canopy_singl(origin, points, max_canopy_dist, partial);

//		Canopy* final_canopy = canopy_walk(origin, points, close_points, max_canopy_dist,
//			max_close_dist, min_step_dist, max_num_canopy_walks, num_canopy_jumps);

		canopy_vector[origin_i] = final_canopy;

	}
	*/
	time_profile.stop_timer("ReClustering");

	
	_log(logINFO) << "Number of all canopies: " << canopy_vector.size();



	return canopy_vector;

}

void filter(options * opt, TimeProfile time_profile, vector<Point*>& points,
	vector<Point*>& guidePoints, vector<Point*>& filtered_points) {

	int filter_min_obs = opt->filter_min_obs;
	double filter_max_top3_sample_contribution = opt->filter_max_top3_sample_contribution;
	bool use_spearman = opt->use_spearman;
	string input_filter_file = opt->input_filter_file;

	_log(logINFO) << "Running basic validation of profiles";
	_log(logINFO) << "filter_max_top3_sample_contribution:\t " << filter_max_top3_sample_contribution;
	_log(logINFO) << "filter_min_obs:\t " << filter_min_obs;
	_log(logINFO) << "";

	time_profile.start_timer("Profiles validation");
	verify_proper_point_input_or_die(points, guidePoints);
	time_profile.stop_timer("Profiles validation");

	vector<Point*> points_filtered_out_due_to_three_point_proportion_filter;
	vector<Point*> points_filtered_out_due_to_num_non_zero_samples_filter;
	set<Point*> filtered_out_points;

	time_profile.start_timer("Input profiles filtering");
	if ((filter_min_obs == 0) && (filter_max_top3_sample_contribution > 0.9999)) {//all filters deactivated
		if (use_spearman) {
#pragma omp parallel for shared(points)
			for (size_t i = 0; i < points.size(); i++) {
				points[i]->convert_to_rank();
			}
		}
		filtered_points = points;

	}
	else {
#pragma omp parallel for shared(points_filtered_out_due_to_three_point_proportion_filter, points_filtered_out_due_to_num_non_zero_samples_filter, filtered_points, filtered_out_points)
		for (int i = 0; i < (int)points.size(); i++) {
			//Both filters are set
			if ((filter_min_obs > 0) && (filter_max_top3_sample_contribution < 0.9999)) {
				bool point_is_valid = true;

				if (!points[i]->check_if_num_non_zero_samples_is_greater_than_x(filter_min_obs))
				{
#pragma omp critical
					{
						points_filtered_out_due_to_num_non_zero_samples_filter.push_back(points[i]);
						filtered_out_points.insert(points[i]);
						point_is_valid = false;
					}
				}

				if (!points[i]->check_if_top_three_point_proportion_is_smaller_than(filter_max_top3_sample_contribution))
				{
#pragma omp critical
					{
						points_filtered_out_due_to_three_point_proportion_filter.push_back(points[i]);
						filtered_out_points.insert(points[i]);
						point_is_valid = false;
					}
				}

				if (point_is_valid) {
#pragma omp critical
					filtered_points.push_back(points[i]);
				}
			}
			else if (filter_min_obs > 0) {
				if (points[i]->check_if_num_non_zero_samples_is_greater_than_x(filter_min_obs)) {
#pragma omp critical
					filtered_points.push_back(points[i]);
				}
				else
				{
#pragma omp critical
					{
						points_filtered_out_due_to_num_non_zero_samples_filter.push_back(points[i]);
						filtered_out_points.insert(points[i]);
					}
				}
			}
			else if (filter_max_top3_sample_contribution < 0.9999) {
				if (points[i]->check_if_top_three_point_proportion_is_smaller_than(filter_max_top3_sample_contribution)) {
#pragma omp critical
					filtered_points.push_back(points[i]);
				}
				else
				{
#pragma omp critical
					{
						points_filtered_out_due_to_three_point_proportion_filter.push_back(points[i]);
						filtered_out_points.insert(points[i]);
					}
				}
			}
			if (use_spearman) {
				points[i]->convert_to_rank();
			}

		}
	}
	if (input_filter_file != "") {
		ofstream filtered_point_file;
		filtered_point_file.open(input_filter_file.c_str(), ios::out | ios::trunc);
		filtered_point_file << "#filtered_profile_id\tinput_filter_name\n";
		for (size_t i = 0; i < points_filtered_out_due_to_num_non_zero_samples_filter.size(); i++) {
			filtered_point_file << points[i]->id << "\t" << "min_observations_filter" << "\n";
		}
		for (size_t i = 0; i < points_filtered_out_due_to_three_point_proportion_filter.size(); i++) {
			filtered_point_file << points[i]->id << "\t" << "max_top3_sample_contribution_filter" << "\n";
		}
		filtered_point_file.close();
	}

	time_profile.stop_timer("Input profiles filtering");
	_log(logINFO) << "Number of profiles filtered out due to three sample signal contribution filter: " << points_filtered_out_due_to_three_point_proportion_filter.size();
	_log(logINFO) << "Number of profiles filtered out due to non zero samples number filter: " << points_filtered_out_due_to_num_non_zero_samples_filter.size();
	_log(logINFO) << "Number of profiles filtered out: " << points.size() - filtered_points.size();

	_log(logINFO) << "Finished input profiles processing";

	_log(logINFO) << "Number of profiles after filtering: " << filtered_points.size();

	//Sometimes these filters will remove over 50% of the dataset, let's free filtered out points now
	_log(logINFO) << "Relseasing filtered out points";
	for (auto p : filtered_out_points) {
		delete p;
	}

	_log(logINFO) << "Relseasing filtered out points: Done";

}



int handleRms(vector<Point*>& points, const vector<bool>& rm) {
	//action needed?
	int sumRm(false);

	for (uint i = 0; i < rm.size(); i++) {
		if (rm[i]) {
			sumRm ++;
			
		}
	}
	if (sumRm==0) { return sumRm; }
	//task here is to remove entries in each point, that should not be taken into account for creating MGS
	for (uint i = 0; i < points.size(); i++) {
		points[i]->pseudoRmSamples(rm, sumRm);
	}


	return sumRm;
}

//this routine scans for autocorrelated samples (and removes these)
vector<bool> autocorr_filter(options * opt, TimeProfile time_profile, const vector<Point*>& points,
	const vector<PRECISIONT>& sampleSums) {
	//look which samples are too close and mark down for filtering
	PRECISIONT bound = opt->sampleMinDist;
	if (bound >= 2) { vector<bool> rmT(0);  return rmT; }
	time_profile.start_timer("Input profiles filtering");
	_log(logINFO) << "Filtering autocorrelated samples";

	if (!points[0]->sparse) { cerr << "Autocorrelation only implemented for sparse data\n"; exit(232); }
	int nmSmpls = points[0]->num_data_samples;
	//reduce to most abundant points
	vector<unordered_map<int, PRECISIONT>> abPoints(nmSmpls);
	for (uint i = 0; i < points.size();i++) {
		for (auto x : points[i]->sp_data) {
			abPoints[x.first][i] = x.second; //basic transpose operation
		}
	}
	//srt descending
	typedef std::function<bool(std::pair<int, PRECISIONT>, std::pair<int, PRECISIONT>)> Comparator;
	Comparator compFunctor =
		[](std::pair<int, PRECISIONT> elem1, std::pair<int, PRECISIONT> elem2)
	{
		return elem1.second > elem2.second;
	};
	_log(logINFO) << "Transposed matrix";

	//filter out the X highest points..
	int MaxEntries = 50000;
	vector<unordered_map<int, PRECISIONT>> abPointsF(nmSmpls);
#pragma omp parallel for shared(abPoints, abPointsF)
	for (int j = 0; j < nmSmpls; j++) {
		if ((int)abPoints[j].size() <  MaxEntries) { 
			abPointsF[j] = abPoints[j];
			for (auto x : abPoints[j]) {
			}
			continue; 
		}
		set<pair<int, PRECISIONT>, Comparator> setOfWords(
			abPoints[j].begin(), abPoints[j].end(), compFunctor);
		int cnt = 0;
		for (pair<int, PRECISIONT> ele : setOfWords) {
			abPointsF[j][ele.first] = ele.second;
			cnt++;
			if (cnt > MaxEntries) { break; }
		}
	}

	//deleting big vector
	vector<unordered_map<int, PRECISIONT>>().swap(abPoints);

	
	_log(logINFO) << "Filtered "<< MaxEntries<<" largest entries";
	//now calculate correlations among samples..
	//abPointsF = abPoints;
	int num_threads = opt->num_threads;
	double nmSmpls_d = (double)nmSmpls;
//	double nrows_d = double(nrows);
	vector<bool> rm(abPointsF.size(), false);
	int cntRm(0);
	vector<vector<PRECISIONT>> corrs(abPointsF.size(), vector<PRECISIONT>(abPointsF.size(),0));
	vector < jobCor > slots(num_threads);
	int j = 0; 
	for (int i=0;i< nmSmpls;i++){
		if (rm[i]) { continue; }
		//get_distance_between_umaps_v (abPointsF, i, nmSmpls); continue;
		while (true) {
			if (!slots[j].inUse) {
				break;
			}
			else if (slots[j].fut.wait_for(std::chrono::milliseconds(5)) == std::future_status::ready) {
				slots[j].inUse = false;
				smplCor tmp = slots[j].fut.get();
				for (uint f = 0; f < tmp.dist.size(); f++) {
					corrs[tmp.i[f]][tmp.k[f]] = tmp.dist[f];
					if (tmp.dist[f] < bound) {//needs to be filtered
						//decide which sample has more depth
						if (sampleSums[tmp.i[f]] > sampleSums[tmp.k[f]]) {
							rm[tmp.k[f]] = true;
						}
						else {
							rm[tmp.i[f]] = true;
						}
						cntRm++;
					}
				}
			}
			else {
				j++;
				if (j >= num_threads) { j = 0; }
			}
		}
		if (slots[j].inUse == false) {
			//create_canopy_singl(origin, points, max_canopy_dist, partial)
			slots[j].fut = async(std::launch::async, get_distance_between_umaps_v, ref(abPointsF), i, nmSmpls);
			
			//get_distance_between_umaps(abPointsF[i], abPointsF[k], nmSmpls_d, i, k);
			slots[j].inUse = true;
			j++;
			if (j >= num_threads) { j = 0; }
		}
	}
	//write out matrix of sample correations (if requested)
	writeMatrix(corrs, opt->sampleDistMatFile);

	time_profile.stop_timer("Input profiles filtering");
	_log(logINFO) << "Removed " << cntRm << " sample due to correlation distance < " << bound;
	_log(logINFO) << "Done Filtering autocorrelated samples";
	vector<unordered_map<int, PRECISIONT>>().swap(abPointsF); // delete vector


	return(rm);
}

void writeMatrix(vector<vector<PRECISIONT>>& mat, string of){
	if (of == "") { return; }
	ofstream out(of.c_str(), ios::out);
	if (!out) {
		cerr << "Can't write to " << of << endl;
	}
	_log(logINFO) << "Writing sample correlation matrix";
	for (uint i = 0; i < mat.size(); i++) {
		out << mat[i][0];
		for (uint j = 1; j < mat.size(); j++) {
			out << "\t"<<mat[i][j];
		}
		out << endl;
	}
	out.close();

}