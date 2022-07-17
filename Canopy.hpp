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
#ifndef CANOPY 
#define CANOPY


#include "Point.hpp"


using namespace std;




/**
 * Represents a canopy
 */
class Canopy {
    public:
        //Constructor - copies the given point as center, neighbours are empty (used only in rare cases) 
        Canopy(Point* center_to_copy,  int deletedSmpls, bool = true);

        //Constructor - assigns the neighour points and creates new point for center (representing canopy profile) 
        //only used temporarily
		Canopy(vector< Point*>& neighbours, int deletedSmpls);

        //Destructor - deletes only the center point - not the neighours
        virtual ~Canopy();
		
		//handle output to file
		void print2file(ofstream* out_file_memb, ofstream* out_file_pro,
			options*, int, int, bool );

        //Set's the profile for the center point
        void find_and_set_center(int deletedSmpls);

		//origin number of submission, used for correlating index
		void set_ori(int x) { ori_index = x; }
		int get_ori(void){ return(ori_index); }
		int ori_index;
		
		void restore_rm(int sr) { center->restore_rm(sr); }
        
        //Center point representing the canopy profile
		Point* center;

        //List of points belonging to the canopy
        std::vector< Point*> neighbours;
		//list of correlations corresponding to neighbor entries..
		vector<PRECISIONT> corrs;

        //Debugging printout of the canopy
        friend std::ostream& operator<<(std::ostream& ost, const Canopy& c);
        
        //Comparison operator of Canopy objects by the number of their neighbours
        //friend bool compare_canopy_ptrs_by_canopy_size(const Canopy* a, const Canopy* b);
};


/*
class CanopyCC : public Canopy {
	CanopyCC(Point *p):Canopy(p){}
};
*/
bool compare_canopy_ptrs_by_canopy_size(const shared_ptr<Canopy> a, const shared_ptr<Canopy> b);

//internal structures
struct job {
	std::future <shared_ptr<Canopy>> fut;
	bool inUse = false;
};
struct job2 {
	std::future <Point*> fut;
	bool inUse = false;
};


inline vector<string> splitByComma(const string& fileS, bool requireTwo, char SrchStr = ',') {
	string::size_type pos = fileS.find(SrchStr);
	if (pos == string::npos) {
		if (requireTwo) {
			cerr << fileS << endl;
			cerr << "Could not find \"" << SrchStr << "\" in input file (required for paired end sequences, as these are two files separated by \",\")\n";
			exit(14);
		}
		else {
			return vector<string>(1, fileS);
		}
	}
	vector<string> tfas(2, "");
	tfas[0] = fileS.substr(0, pos);
	tfas[1] = fileS.substr(pos + 1);
	return tfas;
}

inline vector<string> splitByCommas(const string& fileS, char SrchStr = ',') {
	if (fileS.find(SrchStr) == string::npos) { return vector<string>(1, fileS); }
	vector<string> res = splitByComma(fileS, true, SrchStr);
	vector<string> ret(0); ret.push_back(res[0]);
	while (res[1].find(SrchStr) != string::npos) {
		res = splitByComma(res[1], true, SrchStr);
		ret.push_back(res[0]);
	}
	ret.push_back(res[1]);
	return ret;
}


#endif
