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
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <functional>
#include <numeric>
#include <limits>

//#include <boost/bind.hpp>
//#include <boost/functional/hash.hpp>
#//include <boost/algorithm/string.hpp>

#include "Point.hpp"
#include "Log.hpp"
#include "Stats.hpp"

using namespace std;



extern ProfileMeasureType profile_measure;


Point::Point( Point* p, int deletedSmpls):sample_data(NULL),
#ifdef PRECARRAY
sample_data_pearson_precomputed(NULL),
#else
SumD(0),StdDev(0),
#endif
			precomputed(false), sparse(p->sparse){
	
	id = p->id;
    num_data_samples = p->num_data_samples + deletedSmpls;

    sample_data = new PRECISIONT[num_data_samples];
	if (sparse) {
		for (int i = 0; i < num_data_samples; i++) {
			sample_data[i] = p->getDataSparse(i);
		}
	}
	else {
		for (int i = 0; i < num_data_samples; i++) {
			sample_data[i] = p->sample_data[i];
		}
	}
	for (auto x : sp_data_rm) {
		sp_data_rm[x.first] = x.second;
	}

    if(p->precomputed ){
        //The above if only checks if the point being copied has had its sample pearson data precomputed
        //In fact it should never copy a non precomputed point
        //Precomputing was added as a clumsy way to lower memory usage
#ifdef PRECARRAY
		if (sparse) {
			sp_data_precomp = p->sp_data_precomp;
		}
		else {
			sample_data_pearson_precomputed = new PRECISIONT[num_data_samples];
			for (int i = 0; i < num_data_samples; i++) {
				sample_data_pearson_precomputed[i] = p->sample_data_pearson_precomputed[i];
			}
		}
#else
		SumD = p->SumD; StdDev = p->StdDev;
#endif
		precomputed = true;
    } else {
#ifdef PRECARRAY
		sample_data_pearson_precomputed = NULL;
#endif
		precomputed = false;
    }
}


Point::~Point() {
	if (sample_data != NULL) {
		delete[] sample_data;
	}

#ifdef PRECARRAY
	if (sample_data_pearson_precomputed != NULL) {
		delete[] sample_data_pearson_precomputed; //Note that for some points this might be a null pointer, see constructor and delay_precomputing_pearson_data flag
	}
#endif
}
Point::Point(const char* line,bool sp):sample_data(NULL), 
#ifdef PRECARRAY
sample_data_pearson_precomputed(NULL),
#else
SumD(0), StdDev(0),
#endif
	precomputed(false),sparse(sp){
	char *next_token = NULL;
	const char* delim = "\t";

    //Copy line to private buffer - strtok will modify it
    char* private_line = new char[strlen(line) + 1];
	//strcpy_s(private_line, sizeof private_line, line);
	strcpy(private_line, line);
	_log(logDEBUG3)<< "Point constructor, got: \"" << line << "\"";

    //Read gene id - first word in the line
	//char* word = strtok_s(private_line, delim, &next_token);
	char* word = strtok(private_line, delim);
	id = string(word);
    _log(logDEBUG3)<< "Point constructor, point id: \"" << id << "\""; 

    //Fill vector with data samples
    std::vector<PRECISIONT> sample_data_vector;
    sample_data_vector.reserve(700);

    word = strtok(NULL, delim);
    while( word != NULL ){
        sample_data_vector.push_back((PRECISIONT)atof(word));
        word = strtok(NULL, delim);
    }

    //Get number of samples for this point
    num_data_samples = sample_data_vector.size();
    _log(logDEBUG3)<< "Point constructor, num data samples: \"" << num_data_samples << "\""; 


    //Allocate memory for sample_data (but not pearson precomputed)
    sample_data = new PRECISIONT[num_data_samples];
    for(size_t i = 0; i < sample_data_vector.size(); i++){
        sample_data[i] = sample_data_vector[i];
    }

    //Precomputing of point's pearson data should happen here, but it creates a huge memory spike when creating points of which many will be filtered out
    //Now it is up to user to precompute these using allocate_and_precompute_pearson_data
    //It is a bit clumsy but it does help with memory spikes
#ifdef PRECARRAY
	sample_data_pearson_precomputed = NULL;
#endif

    delete[] private_line;
}

vector<PRECISIONT> Point::rankSort(const PRECISIONT* v_temp, const size_t size) {
	vector<pair<PRECISIONT, size_t> > v_sort(size);

	for (size_t i = 0U; i < size; ++i) {
		v_sort[i] = make_pair(v_temp[i], i);
	}

	sort(v_sort.begin(), v_sort.end());

	pair<double, size_t> rank;
	vector<PRECISIONT> result(size);

	for (size_t i = 0U; i < size; ++i) {
		if (v_sort[i].first != rank.first) {
			rank = make_pair(v_sort[i].first, i);
		}
		result[v_sort[i].second] = rank.second;
	}
	return result;
}

void Point::seal() {
	if (sparse) {
		for (int i = 0; i < num_data_samples; i++) {
//			if (sample_data[i] != 0) {
			if (sample_data[i] > 1e-15) { //take compputational inaccuracy into account..
				sp_data[i] = sample_data[i];
			}
		}
		delete[] sample_data;
		sample_data = NULL;
	}

}

void Point::convert_to_rank(){
	vector<PRECISIONT> tmp = rankSort(sample_data,num_data_samples);
	for (int i = 0; i < num_data_samples; i++) {
		sample_data[i] = tmp[i];
	}

}


void  Point::precompute_pearson_data() {
	//calc both for debugging..
#ifdef PRECARRAY
	if (sparse) {
		precompute_pearson_data_sparse();
	}
	else {
		precompute_pearson_data_array();
	}
#else
	if (sparse) {
		precompute_pearson_data_sparse_2();
	}
	else {
		precompute_pearson_data_array_2();
	}
#endif
	precomputed = true;

}

void  Point::restore_rm(int sr) {
	for (auto x : sp_data_rm) {
		if (sparse) {
			sp_data[x.first] = x.second;
		}
		else {
			sample_data[x.first] = x.second;
		}
	}
}
void Point::pseudoRmSamples(const vector<bool> & rm, int sumRm) {
	if (sparse) {
		unordered_map<int, PRECISIONT>::iterator  x = sp_data.begin();
		while(x != sp_data.end() ) {
			int idx = x->first;
			if (rm[x->first]) {
				sp_data_rm[x->first] = x->second;
				x= sp_data.erase(x);
			}
			else {
				x++;
			}
		}
	}
	else {
		for (int i = 0; i < num_data_samples; i++) {
			if (rm[i] && sample_data[i] > 0) {
				sp_data_rm[i] = sample_data[i];
				sample_data[i] = -1;
			}
		}
	}

	num_data_samples -= sumRm;
}

void Point::addToVec(vector<PRECISIONT>& sms) {
	if ((int)sms.size() < num_data_samples) {
		sms.resize(num_data_samples, 0);
	}
	if (sparse) {
		for (auto x : sp_data) {
			sms[x.first] += x.second;
		}
	}
	else {
		for (int i = 0; i < num_data_samples; i++) {
			sms[i] += sample_data[i];
		}
	}

}


#ifndef PRECARRAY
PRECISIONT Point::getDist_precomp(Point* oth) {
	double sum_XY(0);
	if (sparse) {
		const unordered_map<int, PRECISIONT>& v2 = oth->sp_data;
		for (auto x : sp_data) {
			auto fnd = v2.find(x.first);
			if (fnd == v2.end()) {//v2===Yi == 0
				continue;
			}
			sum_XY += x.second * fnd->second;
		}
	}
	else {
		for (int i = 0; i < num_data_samples; i++)
		{
			// sum of elements of array X/ Y
			sum_XY += sample_data[i] * oth->sample_data[i];
		}

	}
	PRECISIONT dist = 1 - ((PRECISIONT)((double)num_data_samples * sum_XY - SumD * oth->SumD)
		/ sqrt(StdDev * oth->StdDev));
			//* (n * squareSum_Y - sum_Y * sum_Y)));
	return dist;
}


void  Point::precompute_pearson_data_array_2() {

	//Calculate sum and average of data samples
	SumD = 0; StdDev = 0;
	for (int i = 0; i < num_data_samples; i++) {
		if (sample_data[i] < 0) {continue;}
		SumD += sample_data[i];
	}
	//Calculate standard deviation of data samples
	double factor_sum = 0;
	for (int i = 0; i < num_data_samples; i++) {
		if (sample_data[i] < 0) { continue; }
		factor_sum += pow((sample_data[i]), 2);
	}
	
	StdDev = ((double)num_data_samples * factor_sum) - 
		(SumD * SumD);

}
void  Point::precompute_pearson_data_sparse_2() {

	//Calculate sum and average of data samples
	SumD = 0; StdDev = 0;
	for (auto x : sp_data)
		SumD += x.second;

	//Calculate standard deviation of data samples
	double factor_sum = 0;
	for (auto x : sp_data)
		factor_sum += pow((x.second ), 2);
	StdDev = ((double)num_data_samples * factor_sum) -
		(SumD * SumD);

}


#else
void  Point::precompute_pearson_data_array() {

	//Calculate sum and average of data samples
	double sum = 0, avg = 0, num_data_samples_d = num_data_samples;
	for (int i = 0; i < num_data_samples; i++)
		sum += sample_data[i];

	avg = sum / num_data_samples_d;

	//Calculate standard deviation of data samples
	double factor_sum = 0;
	for (int i = 0; i < num_data_samples; i++)
		factor_sum += pow((sample_data[i] - avg), 2);

	double stddev = 0;
	stddev = sqrt(factor_sum / num_data_samples_d);

	//Precompute pearson data
	for (int i = 0; i < num_data_samples; i++) {
		if (fabs(stddev) < 2 * std::numeric_limits<PRECISIONT>::min())
			sample_data_pearson_precomputed[i] = 0;
		else
			double tmp = (sample_data[i] - avg) / (stddev * num_data_samples_d);
			sample_data_pearson_precomputed[i] = (sample_data[i] - avg) / (stddev * num_data_samples_d);
	}
}
void  Point::precompute_pearson_data_sparse() {

	//Calculate sum and average of data samples
	double sum = 0, avg = 0, num_data_samples_d = num_data_samples;
	for (auto x: sp_data)
		sum += x.second;

	avg = sum / num_data_samples_d;

	//Calculate standard deviation of data samples
	double factor_sum = 0;
	for (auto x : sp_data)
		factor_sum += pow((x.second - avg), 2);

	double stddev = 0;
	stddev = sqrt(factor_sum / num_data_samples_d);

	//Precompute pearson data
	for (auto x : sp_data){
		if (fabs(stddev) > 2 * std::numeric_limits<PRECISIONT>::min()) {
			sp_data_precomp[x.first] = (x.second - avg) / (stddev * num_data_samples_d);
		}
	}
}
#endif

void Point::allocate_and_precompute_pearson_data(){
#ifdef PRECARRAY
    if(sample_data_pearson_precomputed != NULL){
        //When this is thrown it means that the point on which this function is executing had already had it's sample_data allocated
        //this function should only be called on Points with points not precomputed
        throw "Precomputing already existing pearson data, this is a memory leak";
    }

    //Allocate and copy samples into array
    sample_data_pearson_precomputed = new PRECISIONT[num_data_samples]; 
#endif
    precompute_pearson_data();  

}



bool Point::check_if_num_non_zero_samples_is_greater_than_x(int x){
    int num_non_zero_samples = 0;
	if (sparse) {
		for (auto y : sp_data) {
			if (y.second > std::numeric_limits<PRECISIONT>::min()) {
				num_non_zero_samples++;
				if (num_non_zero_samples >= x)
					return true;
			}
		}
	}
	else {
		for (int i = 0; i < num_data_samples; i++) {
			if (sample_data[i] > std::numeric_limits<PRECISIONT>::min()) {
				num_non_zero_samples++;
				if (num_non_zero_samples >= x)
					return true;
			}
		}
	}
    return false;

}

bool Point::check_if_top_three_point_proportion_is_smaller_than(PRECISIONT x){

    vector<PRECISIONT> temp_data_samples;
    temp_data_samples.resize(num_data_samples, 0.0);//pseudo vector, exat filling doesn't matter
	int cnt(0);
	if (sparse) {
		for (auto y : sp_data){
			temp_data_samples[cnt] = y.second;
			cnt++;
		}
	}
	else {
		for (int i = 0; i < num_data_samples; i++) {
			temp_data_samples[i] = sample_data[i];
		}
	}

    std::sort(temp_data_samples.begin(), temp_data_samples.end(), std::greater<PRECISIONT>());
    //std::reverse(temp_data_samples.begin(), temp_data_samples.end());

    PRECISIONT sum_data_samples = std::accumulate(temp_data_samples.begin(), temp_data_samples.end(), 0.0 );
    PRECISIONT sum_top_three = temp_data_samples[0] + temp_data_samples[1] + temp_data_samples[2]; 

    if(sum_data_samples > std::numeric_limits<PRECISIONT>::min()){
        return (sum_top_three / sum_data_samples) < x - std::numeric_limits<PRECISIONT>::min();
    } else {
        //All samples have 0 value - can't divide by 0
        return false;
    }

}

void verify_proper_point_input_or_die(const std::vector< Point*>& points,
	const std::vector< Point*>& gp){
    
    //Verify all points have the same number of samples
    int num_samples = points[0]->num_data_samples;
	for (const Point* point : points) {
		assert(point->num_data_samples == num_samples);
	}
	for (const Point* point : gp) {
		assert(point->num_data_samples == num_samples);
	}

    _log(logINFO) << "Finished reading profiles input file";
    _log(logINFO) << "Observed number of samples per profile: " << num_samples;
    _log(logINFO) << "Number of profiles read: " << points.size();

}


PRECISIONT get_partial_distance_between_points(const Point* p1, const Point* p2) {
	// function that returns correlation coefficient. 
	double sum_X(0); double sum_Y(0); double sum_XY(0), n(0);
	double squareSum_X(0); double squareSum_Y(0);

	bool sparse = p1->sparse;

	if (!sparse) {
		PRECISIONT* X = p1->sample_data;
		PRECISIONT* Y = p2->sample_data;
		for (int i = 0; i < p1->num_data_samples; i++)
		{
			//partial part: skip entries in Y that are absent
			if (Y[i] == 0 && X[i] != 0) {
				continue;
			}
			PRECISIONT Xi = X[i];
			PRECISIONT Yi = Y[i];
			// sum of elements of array X/ Y
			sum_X += Xi;
			sum_Y += Yi;
			sum_XY += Xi * Yi;
			// sum of square of array elements. 
			squareSum_X += Xi * Xi;
			squareSum_Y += Yi * Yi;
			n += 1.0;
		}
	} else {
		const unordered_map<int, PRECISIONT>& v1 = p1->sp_data;
		const unordered_map<int, PRECISIONT>& v2 = p2->sp_data;
		for (auto x : v1) {
			auto fnd = v2.find(x.first);
			if (fnd == v2.end()) {//v2===Yi == 0
				continue;
			}
			PRECISIONT Xi = x.second;
			PRECISIONT Yi = fnd->second;

			sum_X += Xi;
			sum_Y += Yi;
			sum_XY += Xi * Yi;
			// sum of square of array elements. 
			squareSum_X += Xi * Xi;
			squareSum_Y += Yi * Yi;
			n += 1.0;
		}
	}
	// use formula for calculating correlation coefficient. 
	//"1-" makes a distance from corr
	PRECISIONT dist = 1 - ((PRECISIONT)(n * sum_XY - sum_X * sum_Y)
		/ sqrt((n * squareSum_X - sum_X * sum_X)
			* (n * squareSum_Y - sum_Y * sum_Y)));
	//cout << corr << " "<< n<< " "<< n * squareSum_Y - sum_Y * sum_Y<<"X ";

	if (dist < 0) {
		int x = 0;
	}

	return dist;
}



smplCor get_distance_between_umaps_v(const vector<unordered_map<int, PRECISIONT>>& vs,
	uint i, int nmSmpls) {
	smplCor ret;
	for (int k = i + 1; k < nmSmpls; k++) {
		ret.i.push_back(i);
		ret.k.push_back(k);
		PRECISIONT dd = get_distance_between_umaps(vs[i], vs[k]);
		ret.dist.push_back(dd);
	}
	return ret;
}

PRECISIONT get_distance_between_umaps(const unordered_map<int, PRECISIONT>& v1,
	const unordered_map<int, PRECISIONT>& v2) {
	// function that returns correlation coefficient. 
	double sum_X(0); double sum_Y(0); double sum_XY(0);
	double squareSum_X(0); double squareSum_Y(0);
	double n = v1.size() + v2.size();
	auto v2end = v2.end();
	for (auto x : v1) {
		PRECISIONT Xi = x.second;
		sum_X += Xi;
		squareSum_X += Xi * Xi;

		auto fnd = v2.find(x.first);
		if (fnd == v2end) {//v2===Yi == 0
			continue;
		}
		n--;
		PRECISIONT Yi = fnd->second;

		sum_Y += Yi;
		sum_XY += Xi * Yi;
		// sum of square of array elements. 
		squareSum_Y += Yi * Yi;
	}
	auto v1end = v1.end();
	for (auto y : v2) {
		auto fnd = v1.find(y.first);
		if (fnd != v1end) {//already counted in loop above
			continue;
		}
		PRECISIONT Yi = y.second;
		sum_Y += Yi;
		squareSum_Y += Yi * Yi;
	}
	
	// use formula for calculating correlation coefficient. 
	//"1-" makes a distance from corr
	PRECISIONT dist = ((PRECISIONT)(n * sum_XY - sum_X * sum_Y  )
		/ sqrt((n * squareSum_X - sum_X * sum_X)
			* (n * squareSum_Y - sum_Y * sum_Y)));
	//cout << corr << " "<< n<< " "<< n * squareSum_Y - sum_Y * sum_Y<<"X ";
	dist = 1 - dist;
	return dist;
}


PRECISIONT get_distance_between_points(Point* p1, Point* p2) {
#ifndef PRECARRAY
	PRECISIONT dist = p1->getDist_precomp(p2);
#else
	int len = p1->num_data_samples;
	PRECISIONT dist = 1 - pearsoncorr_from_precomputed(len, p1->sample_data_pearson_precomputed,
		p2->sample_data_pearson_precomputed);
#endif // !PRECARRAY

    //if(log_level >= logDEBUG3){
    //    _log(logDEBUG3) << "<<<<<<DISTANCE<<<<<<";
    //    _log(logDEBUG3) << "point: " << p1->id;
    //    for(int i=0; i < p1->num_data_samples; i++){
    //        _log(logDEBUG3) << "\t"<<p1->sample_data[i];
    //    }
    //    _log(logDEBUG3) << "point: " << p2->id;
    //    for(int i=0; i < p2->num_data_samples; i++){
    //        _log(logDEBUG3) << "\t"<<p2->sample_data[i];
    //    }
    //    _log(logDEBUG3) << "distance: " << dist;
    //}

    return dist; 
}

Point* get_centroid_of_points(const std::vector< Point*>& points,int deletedSmpls){

    assert(points.size());

	Point* centroid = new Point(points[0], deletedSmpls);
    centroid->id = "!GENERATED!";

    const int num_samples = points[0]->num_data_samples;
    const int num_points = points.size();
	double num_points_d = (PRECISIONT)num_points;
	
	bool sparse = points[0]->sparse;


	//Number which multiplied with length of the vector
 //will give us the element corresponding to the percentile
	PRECISIONT percentile_multiplier = -1;

	//The reason for the ignore here is that when profile_measure is MEAN then code above is ececuted (which is much faster)
	//The value below will thus never be equal to MEAN
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wswitch"
	switch (profile_measure) {
	case MEDIAN:
		percentile_multiplier = 0.5;
		break;
	case PERCENTILE_75:
		percentile_multiplier = 0.75;
		break;
	case PERCENTILE_80:
		percentile_multiplier = 0.80;
		break;
	case PERCENTILE_85:
		percentile_multiplier = 0.85;
		break;
	case PERCENTILE_90:
		percentile_multiplier = 0.90;
		break;
	case PERCENTILE_95:
		percentile_multiplier = 0.95;
		break;
	}
#pragma clang diagnostic pop
	//make correction for counting from 0
	 //so median in vector of length 5 would be (5 - 1)*0.5 = 2
	PRECISIONT target_element_i = (num_points - 1)*percentile_multiplier; //We cannot use it since it might (and usually will be) a float like 2.25
	int lower_element_i = (int) floor(target_element_i);
	int upper_element_i = (int) ceil(target_element_i);

	//now we want to take value which is proportional to the percentile
	PRECISIONT lower_to_upper_proportion = target_element_i - lower_element_i;

	assert(percentile_multiplier != -1);

    //_log(logDEBUG4) << "num samples: " << num_samples;

	std::vector<PRECISIONT> point_samples(num_points, 0);
	for (int i = 0; i < num_samples + deletedSmpls; i++) {
		PRECISIONT percentile = getMedian(points, point_samples, lower_element_i,
			upper_element_i, num_points, i, lower_to_upper_proportion, num_points_d, true);
		centroid->sample_data[i] = percentile;
	}
	centroid->seal();
//also take care of eventually deleted points
	for (int i = 0; i < num_samples + deletedSmpls; i++) {
		PRECISIONT percentile = getMedian(points, point_samples, lower_element_i,
			upper_element_i, num_points, i, lower_to_upper_proportion, num_points_d, false);
		if (percentile > 0) {
			centroid->sp_data_rm[i] = percentile;
		}
	}


	centroid->precompute_pearson_data();
    
    return centroid;
}

/*
std::size_t hash_value(const Point& p){
    boost::hash<std::string> hasher;
    return hasher(p.id);
}
*/


PRECISIONT getMedian(const vector<Point*>& points, vector<PRECISIONT>& point_samples, 
	int lower_element_i, int upper_element_i, const int num_points, int i,
	PRECISIONT lower_to_upper_proportion, PRECISIONT num_points_d,
	bool real)
{
	// for(const Point* p : points){
	if (!real) {
		for (int j = 0; j < num_points; j++) {
			point_samples[j] = points[j]->getDataRm(i);
		}
	}
	else if (points[0]->sparse) {
		for (int j = 0; j < num_points; j++) {
			point_samples[j] = points[j]->getDataSparse(i);
		}
	}
	else {
		for (int j = 0; j < num_points; j++) {
			//TODO: this is slow 
			//point_samples.push_back(p->sample_data[i]);
			point_samples[j] = points[j]->sample_data[i];

		}
	}
	PRECISIONT percentile(0);
	if (profile_measure == MEAN) {
		PRECISIONT sum = std::accumulate(point_samples.begin(), point_samples.end(), 0);
		percentile = sum / num_points_d;

	}
	else {
		std::sort(point_samples.begin(), point_samples.end());

		PRECISIONT lower_element_val = point_samples[lower_element_i];
		PRECISIONT upper_element_val = point_samples[upper_element_i];

		percentile = lower_element_val + lower_to_upper_proportion * (upper_element_val - lower_element_val);

	}
	return percentile;
}

std::ostream& operator<<(std::ostream& ost, const Point& p)
{
        ost << "============================" << std::endl;
        ost << "Point: " << p.id << std::endl;
        for(int i=0; i < p.num_data_samples; i++){
            ost << p.sample_data[i] << "\t" ;
        }
        ost << std::endl;
        ost << "============================" << std::endl;
        
        return ost;
}

/*
Point* create_Point(string l) {
	Point * pp = new Point(l.c_str());
	return pp;
}
*/
