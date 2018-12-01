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

#include "Stats.hpp"


PRECISIONT pearsoncorr_from_precomputed(int n, const PRECISIONT*   v1, const PRECISIONT*   v2) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += v1[i] * v2[i];
	}
	return (PRECISIONT)sum * n;
}

PRECISIONT pearsoncorr_from_precomputed(int n, unordered_map<int,PRECISIONT> &  v1, unordered_map<int,PRECISIONT> &  v2) {
	double sum = 0;
	for (auto x : v1) {
		auto fnd = v2.find(x.first);
		if (fnd != v2.end()) {
			sum += x.second * fnd->second;
		}
		//else {
			//do nothing
			//sum += x.second * 0;
		//}
		//v2 only entries don't need to be checked, since these are 0...
	}
	return (PRECISIONT)sum * n;
}


