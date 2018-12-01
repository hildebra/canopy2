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

#include "options.h"
#include "signal_handlers.hpp"

//#define PRECARRAY

#ifndef POINT
#define POINT




enum ProfileMeasureType { MEDIAN, MEAN, PERCENTILE_75, PERCENTILE_80, PERCENTILE_85, PERCENTILE_90, PERCENTILE_95 };




using namespace std;
class Point {
    public:
        Point( Point* p, int);
        Point(const char* line,bool =false);
        virtual ~Point();

        PRECISIONT* sample_data;
		//added for sparse representation
		unordered_map<int, PRECISIONT> sp_data;
#ifdef PRECARRAY
        PRECISIONT* sample_data_pearson_precomputed;
		unordered_map<int, PRECISIONT> sp_data_precomp;
#else
		//added summary stats for corr
		double SumD; //just sum of samples
		double StdDev; // n*SumD2 - Sum2D
#endif
		int num_data_samples;

		bool precomputed;
		bool sparse;

		//pseud removals of samples
		void pseudoRmSamples(const vector<bool> & rm, int sumRm);
		void restore_rm(int);
		unordered_map<int, PRECISIONT> sp_data_rm;

		//finish precomputations, convert to sparse after ranks etc are calculated..
		void seal();
		void addToVec(vector<PRECISIONT>& sms);
		void convert_to_rank();
		vector<PRECISIONT> rankSort(const PRECISIONT* v_temp, const size_t size);
        
        std::string id;

        void allocate_and_precompute_pearson_data();
		void precompute_pearson_data();
		PRECISIONT getDist_precomp( Point* oth);//pearson distance, new implementation..

		PRECISIONT getDataSparse(int x) {
			auto fnd = sp_data.find(x);
			if (fnd != sp_data.end()) {
				return fnd->second;
			}
			else {
				return 0;
			}
		}
		PRECISIONT getDataRm(int x) { // only check for already deleted data
			auto fnd = sp_data_rm.find(x);
			if (fnd != sp_data_rm.end()) {
				return fnd->second;
			}
			else {
				return 0;
			}
		}
		PRECISIONT getData(int x) {
			if (sparse) {
				return getDataSparse(x);
			}
			else {
				return sample_data[x];
			}
		}

        bool check_if_num_non_zero_samples_is_greater_than_x(int x);
        bool check_if_top_three_point_proportion_is_smaller_than(PRECISIONT x);

//        friend PRECISIONT* precompute_pearson_data(PRECISIONT* sample_data);
       // friend std::size_t hash_value(const Point &p);
        friend std::ostream& operator<<(std::ostream& ost, const Point& ls);

		friend PRECISIONT get_distance_between_points(Point* p1, Point* p2);
		friend PRECISIONT get_partial_distance_between_points(const Point* p1, const Point* p2);
		friend Point* get_centroid_of_points(const std::vector< Point*>& points,int);
		friend PRECISIONT getMedian(const vector<Point*>& points, vector<PRECISIONT>& point_samples,
			int lower_element_i, int upper_element_i, const int, int, PRECISIONT, PRECISIONT,
			bool);
        friend void verify_proper_point_input_or_die(const std::vector<Point*>& points, 
			const std::vector< Point*>& gp);
	private:
		//fills "sample_data_pearson_precomputed
#ifdef PRECARRAY
		void precompute_pearson_data_array();
		void precompute_pearson_data_sparse();
#endif
		//new versions that don't rely on array "sample_data_pearson_precomputed" any longer..
		void precompute_pearson_data_array_2();
		void precompute_pearson_data_sparse_2();
};

//Point* create_Point(string l);
struct smplCor {
	vector<PRECISIONT> dist = vector<PRECISIONT>(0);
	vector<uint> i = vector<uint>(0);
	vector<uint> k = vector<uint>(0);
};
struct jobCor {
	std::future <smplCor> fut;
	bool inUse = false;
};


//note that "n" is derrived from entries, this is too imprecise for some cases
PRECISIONT get_distance_between_umaps(const unordered_map<int, PRECISIONT>& v1,
	const unordered_map<int, PRECISIONT>& v2);
smplCor get_distance_between_umaps_v(const vector<unordered_map<int, PRECISIONT>>& vs,
	uint i, int);

#endif
