/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
 * Copyright (C) 2018 Falk Hildebrand (falk.hildebrand@gmail.com), EMBL
 * Copyright (C) 2020 Falk Hildebrand (falk.hildebrand@gmail.com), QIB/EI
 * This file is part of Metagenomics Canopy Clustering Implementation.
 *
 * Metagenomics Canopy Clustering Implementation is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
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
//#include <boost/algorithm/string.hpp>

#include "Point.hpp"
#include "Log.hpp"
#include "Stats.hpp"

using namespace std;



extern ProfileMeasureType profile_measure;
struct job2 {
	std::future<std::unique_ptr<Point>> fut;
	bool inUse = false;
};

std::unique_ptr<Point> line2point(string line, bool sparseMat, int lcnt) {
	std::unique_ptr<Point> pp(new Point(line, sparseMat, lcnt));
	pp->seal();
	return pp;

}

void readMatrix(vector<Point*>& points, vector<std::unique_ptr<Point>>& point_owners,
	vector<PRECISIONT>& sampleSums, string input_file_path, bool sparseMat,
	int num_threads, bool dont_use_mmap) {

	std::unique_ptr<std::istream> point_file;

	if (isGZfile(input_file_path)) {
#ifdef _gzipread
		point_file.reset(new igzstream(input_file_path.c_str(), ios::in));
		cout << "Reading gzip input\n";
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
}
	else {
		point_file.reset(new ifstream(input_file_path.c_str()));
	}
	if (!point_file || !(*point_file)) {
		cerr << "Could not open input matrix: " << input_file_path << "\n";
		exit(1);
	}

	std::string line;
	vector<string> buffered_lines;
	size_t buffered_line_i = 0;
	if (dont_use_mmap) {
		_log(logINFO) << "Reading input matrix line by line";
	}
	else {
		_log(logINFO) << "Reading input matrix into memory";
		while (getline(*point_file, line)) {
			buffered_lines.push_back(line);
		}
	}
	auto get_next_line = [&](string& next_line) -> bool {
		if (dont_use_mmap) {
			return static_cast<bool>(getline(*point_file, next_line));
		}
		if (buffered_line_i >= buffered_lines.size()) {
			return false;
		}
		next_line = buffered_lines[buffered_line_i++];
		return true;
	};

	// The matrix format requires a header row.
	if (!get_next_line(line)) {
		cerr << "Input matrix is empty\n";
		return;
	}

	vector<job2> slots(num_threads);
	int j(0); int lcnt(0);
	auto get_completed_point = [](
		std::future<std::unique_ptr<Point>>& future) -> std::unique_ptr<Point> {
		try {
			return future.get();
		}
		catch (const std::exception& error) {
			cerr << error.what() << "\n";
			exit(1);
		}
	};
	//vector<job2> fut(num_threads); int ji = 0;
	while (true) {
		
		if (j >= num_threads) { j = 0; }
		if (slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(1)) == std::future_status::ready) {
			slots[j].inUse = false;
			std::unique_ptr<Point> pp = get_completed_point(slots[j].fut);
			//cout << pp->lineCnt << " ";
			if (pp->lineCnt >= static_cast<int>(points.size())) {
				points.resize(pp->lineCnt+1);
				point_owners.resize(pp->lineCnt+1);
			}
			points[pp->lineCnt] = pp.get();
			point_owners[pp->lineCnt] = std::move(pp);
			//points.push_back(pp);
		}
		if (slots[j].inUse == false) {
			bool have_data_line = false;
			while (get_next_line(line)) {
				if (line.length() >= 2) {
					have_data_line = true;
					break;
				}
			}
			if (!have_data_line) { break; }
			//line2point(line,sparseMat, use_spearman);
			//cout << line<<endl;
			string lineC = line;
			slots[j].fut = async(std::launch::async, line2point, lineC, sparseMat, lcnt);
			slots[j].inUse = true;
			lcnt++;
		}
		j++;

		//points.push_back(new Point(line.c_str(), sparseMat));		
	}

	for (j = 0; j < num_threads; j++) {
		if (slots[j].inUse == true) {
			slots[j].inUse = false;
			std::unique_ptr<Point> pp = get_completed_point(slots[j].fut);
			if (pp->lineCnt >= static_cast<int>(points.size())) {
				points.resize(pp->lineCnt + 1, nullptr);
				point_owners.resize(pp->lineCnt + 1);
			}
			points[pp->lineCnt] = pp.get();
			point_owners[pp->lineCnt] = std::move(pp);
		}
	}

	// Aggregate in input order rather than future-completion order. Besides
	// producing repeatable floating-point sums, this makes autocorrelation
	// sample priorities identical between input modes and thread counts.
	sampleSums.clear();
	for (Point* point : points) {
		if (point != NULL) {
			point->addToVec(sampleSums);
		}
	}

	cerr << "Read " << lcnt << " rows\n";
}

Point::Point( Point* p, int deletedSmpls):sample_data(NULL),
#ifdef PRECARRAY
sample_data_pearson_precomputed(NULL),
#else
SumD(0),StdDev(0),
#endif
			num_data_samples(p->num_data_samples), lineCnt(p->lineCnt),
			precomputed(false), sparse(p->sparse){
	(void)deletedSmpls;

	id = p->id;

	if (sparse) {
		// Sparse keys retain their original sample indexes while
		// num_data_samples records only the observations used for
		// correlations. Copy the maps directly so those two concepts do not
		// become conflated after autocorrelation filtering.
		sp_data = p->sp_data;
	}
	else {
		sample_data = new PRECISIONT[num_data_samples];
		for (int i = 0; i < num_data_samples; i++) {
			sample_data[i] = p->sample_data[i];
		}
	}
	sp_data_rm = p->sp_data_rm;

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
Point::Point(string line,bool sp, int lc):sample_data(NULL),
#ifdef PRECARRAY
sample_data_pearson_precomputed(NULL),
#else
SumD(0), StdDev(0),
#endif
num_data_samples(0), lineCnt(lc), precomputed(false), sparse(sp)
{
	auto parse_error = [&](const string& message) {
		throw runtime_error("Invalid profile on input row " + to_string(lc + 2) + ": " + message);
	};
	if (line.empty()) {
		parse_error("empty row");
	}
	if (line[line.size() - 1] == '\t') {
		parse_error("trailing empty sample field");
	}
	stringstream ss;
	ss << line;
	int cnt2(-2);
	const char sep = '\t';
	string segments;
	std::vector<PRECISIONT> sample_data_vector;

	while (getline(ss, segments, sep)) {
		cnt2++;
		if (cnt2 == -1) {//this is the row ID
			id = segments;
			if (id.empty()) {
				parse_error("empty profile ID");
			}
			continue;
		}
		char* end = NULL;
		errno = 0;
		const char* value = segments.c_str();
		double parsed = strtod(value, &end);
		if (segments.empty() || end == value || *end != '\0' || errno == ERANGE ||
			!std::isfinite(parsed) || parsed < 0) {
			parse_error("invalid non-negative numeric value '" + segments + "'");
		}
		const PRECISIONT stored_value = static_cast<PRECISIONT>(parsed);
		if (!std::isfinite(stored_value)) {
			parse_error("numeric value is outside the configured precision range");
		}
		sample_data_vector.push_back(stored_value);
	}
	if (sample_data_vector.empty()) {
		parse_error("profile has no sample values");
	}
	num_data_samples = sample_data_vector.size();
	sample_data = new PRECISIONT[num_data_samples];

	num_data_samples = sample_data_vector.size();
	for (size_t i = 0; i < sample_data_vector.size(); i++) {
		sample_data[i] = sample_data_vector[i];
	}

#ifdef PRECARRAY
	sample_data_pearson_precomputed = NULL;
#endif
/*
	char *next_token = NULL;
	const char* delim = "\t";

    //Copy line to private buffer - strtok will modify it
	char* private_line = &line[0];// new char[strlen(line) + 1];
	//strcpy_s(private_line, sizeof private_line, line);
	//strcpy(private_line, line);
	//_log(logDEBUG3)<< "Point constructor, got: \"" << line << "\"";

    //Read gene id - first word in the line
	//char* word = strtok_s(private_line, delim, &next_token);
	char* word = strtok(private_line, delim);
	id = string(word);
    //_log(logDEBUG3)<< "Point constructor, point id: \"" << id << "\""; 

    //Fill vector with data samples
    std::vector<PRECISIONT> sample_data_vector;
    //sample_data_vector.reserve(700);

    word = strtok(private_line, delim);
    while( word != NULL ){
        sample_data_vector.push_back((PRECISIONT)atof(word));
        word = strtok(private_line, delim);
    }

    //Get number of samples for this point
    num_data_samples = sample_data_vector.size();
    //_log(logDEBUG3)<< "Point constructor, num data samples: \"" << num_data_samples << "\""; 
	//return;


    //Allocate memory for sample_data (but not pearson precomputed)
	//assert(num_data_samples > 0);
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

    //delete[] private_line;
	*/
}

vector<PRECISIONT> Point::rankSort(const PRECISIONT* v_temp, const size_t size) {
	vector<pair<PRECISIONT, size_t> > v_sort(size);

	for (size_t i = 0U; i < size; ++i) {
		v_sort[i] = make_pair(v_temp[i], i);
	}

	sort(v_sort.begin(), v_sort.end());

	vector<PRECISIONT> result(size);

	for (size_t group_begin = 0; group_begin < size;) {
		size_t group_end = group_begin + 1;
		while (group_end < size && v_sort[group_end].first == v_sort[group_begin].first) {
			group_end++;
		}
		const PRECISIONT average_rank = static_cast<PRECISIONT>(group_begin + group_end - 1) /
			static_cast<PRECISIONT>(2);
		for (size_t i = group_begin; i < group_end; ++i) {
			result[v_sort[i].second] = average_rank;
		}
		group_begin = group_end;
	}
	return result;
}

void Point::seal() {
	if (sparse && sample_data != NULL) {
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
	vector<bool> removed_samples(num_data_samples, false);
	convert_to_rank(removed_samples);
}

void Point::convert_to_rank(const vector<bool>& removed_samples){
	if (removed_samples.size() != static_cast<size_t>(num_data_samples)) {
		throw runtime_error("Rank conversion received an invalid sample-removal mask");
	}
	vector<PRECISIONT> original(num_data_samples, 0);
	if (sample_data != NULL) {
		for (int i = 0; i < num_data_samples; i++) {
			original[i] = sample_data[i];
		}
	}
	else if (sparse) {
		for (auto value : sp_data) {
			if (value.first >= 0 && value.first < num_data_samples) {
				original[value.first] = value.second;
			}
		}
	}
	else {
		throw runtime_error("Cannot rank a point without sample data");
	}

	vector<PRECISIONT> retained_values;
	vector<int> retained_indexes;
	retained_values.reserve(num_data_samples);
	retained_indexes.reserve(num_data_samples);
	for (int i = 0; i < num_data_samples; ++i) {
		if (!removed_samples[i]) {
			retained_values.push_back(original[i]);
			retained_indexes.push_back(i);
		}
	}
	if (retained_values.empty()) {
		throw runtime_error("Cannot rank a profile after removing every sample");
	}
	vector<PRECISIONT> retained_ranks = rankSort(retained_values.data(), retained_values.size());
	// Preserve full-profile ranks for removed samples so they can be restored
	// for output, while clustering uses ranks recomputed over retained samples.
	vector<PRECISIONT> ranked = rankSort(original.data(), original.size());
	for (size_t i = 0; i < retained_indexes.size(); ++i) {
		ranked[retained_indexes[i]] = retained_ranks[i];
	}
	if (sample_data != NULL) {
		for (int i = 0; i < num_data_samples; i++) {
			sample_data[i] = ranked[i];
		}
	}
	else {
		sp_data.clear();
		for (int i = 0; i < num_data_samples; i++) {
			if (ranked[i] > static_cast<PRECISIONT>(1e-15)) {
				sp_data[i] = ranked[i];
			}
		}
	}
	precomputed = false;
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
	num_data_samples += sr;
}
void Point::pseudoRmSamples(const vector<bool> & rm, int sumRm) {
	if (sparse) {
		mvec::iterator  x = sp_data.begin();
		while(x != sp_data.end() ) {
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
		const mvec& v2 = oth->sp_data;
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
	const double denominator_squared = StdDev * oth->StdDev;
	if (!(denominator_squared > 0.0) ||
		!std::isfinite(denominator_squared)) {
		return std::numeric_limits<PRECISIONT>::infinity();
	}
	double correlation = ((double)num_data_samples * sum_XY - SumD * oth->SumD) /
		sqrt(denominator_squared);
	if (!std::isfinite(correlation)) {
		return std::numeric_limits<PRECISIONT>::infinity();
	}
	correlation = std::max(-1.0, std::min(1.0, correlation));
	return static_cast<PRECISIONT>(1.0 - correlation);
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
	if (temp_data_samples.empty()) {
		return false;
	}
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
	const size_t top_count = std::min<size_t>(3, temp_data_samples.size());
	PRECISIONT sum_top_three = std::accumulate(temp_data_samples.begin(),
		temp_data_samples.begin() + top_count, static_cast<PRECISIONT>(0));

    if(sum_data_samples > std::numeric_limits<PRECISIONT>::min()){
		const PRECISIONT proportion = sum_top_three / sum_data_samples;
		const PRECISIONT tolerance = std::numeric_limits<PRECISIONT>::epsilon() *
			std::max(static_cast<PRECISIONT>(1), std::fabs(x));
		// The option rejects profiles whose contribution is greater than the
		// threshold, so equality (within floating-point precision) must pass.
		return proportion <= x + tolerance;
    } else {
        //All samples have 0 value - can't divide by 0
        return false;
    }

}

void verify_proper_point_input_or_die(const std::vector< Point*>& points,
	const std::vector< Point*>& gp){
	if (points.empty() || points[0] == NULL) {
		_log(logERR) << "The input matrix contains no valid profiles.";
		exit(1);
	}
	const int num_samples = points[0]->num_data_samples;
	if (num_samples < 2) {
		_log(logERR) << "At least two sample columns are required.";
		exit(1);
	}
	unordered_set<string> profile_ids;
	for (const Point* point : points) {
		if (point == NULL || point->num_data_samples != num_samples || point->id.empty()) {
			_log(logERR) << "Input profiles must be non-null, named, and have identical sample counts.";
			exit(1);
		}
		if (!profile_ids.insert(point->id).second) {
			_log(logERR) << "Duplicate profile ID in input matrix: " << point->id;
			exit(1);
		}
	}
	unordered_set<string> guide_ids;
	for (const Point* point : gp) {
		if (point == NULL || point->num_data_samples != num_samples || point->id.empty()) {
			_log(logERR) << "Guide profiles must be non-null, named, and match the input sample count.";
			exit(1);
		}
		if (!guide_ids.insert(point->id).second) {
			_log(logERR) << "Duplicate guide profile ID: " << point->id;
			exit(1);
		}
	}

    _log(logINFO) << "Finished reading profiles input file";
    _log(logINFO) << "Observed number of samples per profile: " << num_samples;
    _log(logINFO) << "Number of profiles read: " << points.size();

}


PRECISIONT get_partial_distance_between_points(const Point* p1, const Point* p2) {
	double sum_X(0), sum_Y(0), sum_XY(0), n(0);
	double squareSum_X(0), squareSum_Y(0);
	const int sample_count = std::min(p1->num_data_samples, p2->num_data_samples);

	// A partial correlation is calculated only where the target profile (Y)
	// is observed. Using getData keeps sparse and dense modes identical.
	for (int i = 0; i < sample_count; i++) {
		const PRECISIONT Yi = p2->getData(i);
		if (Yi == 0) {
			continue;
		}
		const PRECISIONT Xi = p1->getData(i);
		sum_X += Xi;
		sum_Y += Yi;
		sum_XY += Xi * Yi;
		squareSum_X += Xi * Xi;
		squareSum_Y += Yi * Yi;
		n += 1.0;
	}

	if (n < 2.0) {
		return std::numeric_limits<PRECISIONT>::infinity();
	}
	const double denominator_squared = (n * squareSum_X - sum_X * sum_X) *
		(n * squareSum_Y - sum_Y * sum_Y);
	if (!(denominator_squared > 0.0) ||
		!std::isfinite(denominator_squared)) {
		return std::numeric_limits<PRECISIONT>::infinity();
	}
	double correlation = (n * sum_XY - sum_X * sum_Y) / sqrt(denominator_squared);
	if (!std::isfinite(correlation)) {
		return std::numeric_limits<PRECISIONT>::infinity();
	}
	correlation = std::max(-1.0, std::min(1.0, correlation));
	return static_cast<PRECISIONT>(1.0 - correlation);
}
smplCor get_distance_between_umaps_v(vector<mvec2>& vs,
	uint i, int nmSmpls, int num_observations) {
	smplCor ret;
	for (int k = i + 1; k < nmSmpls; k++) {
		ret.i.push_back(i);
		ret.k.push_back(k);
		PRECISIONT dd = get_distance_between_umaps(vs[i], vs[k], num_observations);
		ret.dist.push_back(dd);
	}
	return ret;
}

PRECISIONT get_distance_between_umaps(const mvec2& v1,
	const mvec2& v2, int num_observations) {
	// function that returns correlation coefficient. 
	double sum_X(0); double sum_Y(0); double sum_XY(0);
	double squareSum_X(0); double squareSum_Y(0);
	const double n = static_cast<double>(num_observations);
	auto v2end = v2.end();
	for (auto x : v1) {
		PRECISIONT Xi = x.second;
		sum_X += Xi;
		squareSum_X += Xi * Xi;

		auto fnd = v2.find(x.first);
		if (fnd == v2end) {//v2===Yi == 0
			continue;
		}
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
	
	const double denominator_squared = (n * squareSum_X - sum_X * sum_X) *
		(n * squareSum_Y - sum_Y * sum_Y);
	if (n < 2.0 || !(denominator_squared > 0.0) ||
		!std::isfinite(denominator_squared)) {
		return std::numeric_limits<PRECISIONT>::infinity();
	}
	double correlation = (n * sum_XY - sum_X * sum_Y) / sqrt(denominator_squared);
	if (!std::isfinite(correlation)) {
		return std::numeric_limits<PRECISIONT>::infinity();
	}
	correlation = std::max(-1.0, std::min(1.0, correlation));
	return static_cast<PRECISIONT>(1.0 - correlation);
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

std::unique_ptr<Point> get_centroid_of_points(
	const std::vector<Point*>& points, int deletedSmpls) {

    assert(points.size());

	std::unique_ptr<Point> centroid(new Point(points[0], 0));
    centroid->id = "!GENERATED!";
	centroid->precomputed = false;
	centroid->sp_data.clear();
	centroid->sp_data_rm.clear();

    const int num_samples = points[0]->num_data_samples;
	const int original_sample_count = num_samples + deletedSmpls;
	if (!points[0]->sparse && deletedSmpls != 0) {
		throw runtime_error("Removed-sample restoration is only supported for sparse profiles");
	}
    const int num_points = points.size();
	double num_points_d = (PRECISIONT)num_points;
	
	//Number which multiplied with length of the vector
 //will give us the element corresponding to the percentile
	PRECISIONT percentile_multiplier = -1;

	//The reason for the ignore here is that when profile_measure is MEAN then code above is ececuted (which is much faster)
	//The value below will thus never be equal to MEAN
	switch (profile_measure) {
	case MEAN:
		// The indices are unused by getMedian for means, but a valid value keeps
		// assertion-enabled and release builds on the same control path.
		percentile_multiplier = 0;
		break;
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
	for (int i = 0; i < original_sample_count; i++) {
		PRECISIONT percentile = getMedian(points, point_samples, lower_element_i,
			upper_element_i, num_points, i, lower_to_upper_proportion, num_points_d, true);
		if (centroid->sparse) {
			if (percentile > static_cast<PRECISIONT>(1e-15)) {
				centroid->sp_data[i] = percentile;
			}
		}
		else {
			centroid->sample_data[i] = percentile;
		}
	}
//also take care of eventually deleted points
	for (int i = 0; i < original_sample_count; i++) {
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
		PRECISIONT sum = std::accumulate(point_samples.begin(), point_samples.end(),
			static_cast<PRECISIONT>(0));
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
			ost << p.getData(i) << "\t" ;
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
