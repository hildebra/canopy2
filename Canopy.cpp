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


#include "Canopy.hpp"

Canopy::Canopy(Point* center_to_copy, int deletedSmpls, bool replID):neighbours(0){
    center = new Point(center_to_copy, deletedSmpls);
	if (replID) {
		center->id = "!GENERATED!";
	}
}
Canopy::~Canopy() {
	delete center;
}

Canopy::Canopy(std::vector< Point*> neighbours, int deletedSmpls): 
	neighbours(neighbours){
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
			vector<Point*> nei = this->neighbours;
			list <PRECISIONT> ::iterator it = this->corrs.begin();
			for (size_t j = 0; j < nei.size(); j++) {
				Point* p = nei[j];
				//PRECISIONT dist = get_distance_between_points(c->center, p);
				*out_file_memb << curID << "\t";
				*out_file_memb << p->id << "\t" << *it << "\n";
				it++;
			}
		

	}

}

void Canopy::find_and_set_center(int deletedSmpls){

    center = get_centroid_of_points(neighbours, deletedSmpls);

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


