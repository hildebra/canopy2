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
 /*
 Metagenomics Canopy Clustering Implementation v2
 all copyrights lie with Falk Hildebrand (falk.hildebrand [ta] gmail dot com)
 */


#include "Canopy.hpp"


Canopy::Canopy(Point* center_to_copy, int deletedSmpls, bool replID):neighbours(0),corrs(0){
    center = new Point(center_to_copy, deletedSmpls);
	if (replID) {
		center->id = "!GENERATED!";
	}
}
Canopy::~Canopy() {
	delete center;
}

Canopy::Canopy(std::vector< Point*>& neighbours, int deletedSmpls): 
	neighbours(neighbours), corrs(neighbours.size()){
    find_and_set_center(deletedSmpls);
}

void Canopy::print2file(ofstream* out_file_memb, ofstream* out_file_pro,
	options* opt, int kk, int num_digits,bool guided) {
	string output_cluster_prefix = opt->output_cluster_prefix;
	if (!guided) {
		//members
		for (Point* p : this->neighbours) {
			*out_file_memb << output_cluster_prefix << std::setw(num_digits) << std::setfill('0') << kk << "\t";
			*out_file_memb << p->id;
			*out_file_memb << "\n";
		}

		if (out_file_pro != NULL) { //classical canopy write out
			//profile
			*out_file_pro << output_cluster_prefix << std::setw(num_digits) << std::setfill('0') << kk;

			for (int j = 0; j < this->center->num_data_samples; j++) {
				*out_file_pro << "\t" << this->center->getData(j);
			}

			*out_file_pro << "\n";

		}
	}
	else {//this is a post correlation analysis..
			string curID = this->center->id;
			//cerr << "Writing canopy " << curID << ", size = "<< neighbours.size() <<endl;
			//vector<Point*> nei = this->neighbours;
			//vector<PRECISIONT> corres = this->corrs.begin();
			if (neighbours.size() != corrs.size()) {
				cerr << curID << "\n";
				cerr << "Found unequal neighbours(" << neighbours.size() << ") and corrs(" << corrs.size() << ") size\n";
				exit(932);
			}
			size_t j = 0;
			while ( j < neighbours.size() && j < corrs.size() ) {
				Point* p = neighbours[j];
				//PRECISIONT dist = get_distance_between_points(c->center, p);
				*out_file_memb << curID << "\t";
				*out_file_memb << p->id;
				*out_file_memb << "\t" << corrs[j] << "\n";
				j++;
			}
			//cerr << "  wrote " << j << " genes out\n";
	}

}

void Canopy::find_and_set_center(int deletedSmpls){

    center = get_centroid_of_points(neighbours, deletedSmpls);

}

uint Canopy::cleanUp() {
	bool hasnull = false;
	uint clndGs = 0;
	uint osize = neighbours.size();
	for (size_t i = 0; i < neighbours.size(); i++) {
		if (neighbours[i] == nullptr){
			//hasnull = true;// break;
			clndGs++;
		}
	}
	if (!clndGs) { return clndGs; }

	//clean up//
	std::vector< Point*> neighbours2(osize- clndGs);
	vector<PRECISIONT> corrs2(osize - clndGs);
	uint i2 = 0;
	for (size_t i = 0; i < neighbours.size(); i++) {
		if (neighbours[i] != nullptr) {
			//neighbours2.push_back(neighbours[i]);
			//clndGs++;
			neighbours2[i2] = neighbours[i];
			corrs2[i2] = corrs[i];
			i2++;
		}
	}
	neighbours = neighbours2;
	corrs = corrs2;
	return clndGs;
}

std::ostream& operator<<(std::ostream& ost, const Canopy& c)
{
    ost << ">>>>>>>>>>Canopy>>>>>>>>" << std::endl;
    ost << "Center:" << std::endl;
    if(c.center != NULL)
        ost << *c.center;
    else
        ost << "===NONE===" << endl;
    ost << "Neighbours: " << c.neighbours.size() << std::endl;
    //BOOST_FOREACH(const Point* p, c.neighbours)
    //    ost << p->id << "\t";
    //ost << std::endl;
    ost << ">>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    return ost;
}

bool compare_canopy_ptrs_by_canopy_size(const shared_ptr<Canopy> a, const shared_ptr<Canopy> b){
    return (a->neighbours.size() > b->neighbours.size());
}


