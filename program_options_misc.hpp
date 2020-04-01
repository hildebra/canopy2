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

#pragma once

//centralized library load


#include <cstring>
//#include <stdio.h>
#include <string>
//#include <string.h>
#include <fstream>
#include <iomanip>

//#include <vector>
#include <list>
#include <set>
#include <unordered_map> 
#include <unordered_set>
#include <map>

#include <omp.h>
#include <future>


//#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include <math.h>
#include <fcntl.h>
//#include <algorithm> 

#include "TimeProfile.hpp"

#include <iostream>
//#include <stdio.h>
//#include <math.h>
//#include <limits>





#ifdef SINGLEPRECISION
typedef float PRECISIONT;
#else
typedef double PRECISIONT;
#endif

typedef unsigned int uint;
using namespace std;


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#define _gziprea//d
#pragma warning(disable:4996)
#else
#define _gzipread
#endif


#ifdef _gzipread
#include "gzstream.h"
#endif


using namespace std;

static bool check_if_within_bounds(string option_name, double value, double lower, double higher){
    if( value >= lower - numeric_limits<double>::epsilon() && value <= higher + numeric_limits<double>::epsilon()){
        return true;
    }else{ 
        _log(logERR) << "Option: \"" << option_name << "\" must be a value within range: <" << lower << ";" << higher << ">";
        exit(1);
    }
}

static bool check_if_within_bounds(string option_name, int value, int lower, int higher){
    if( value >= lower && value <= higher ){
        return true;
    }else{ 
        _log(logERR) << "Option: \"" << option_name << "\" must be a value within range: <" << lower << ";" << higher << ">";
        exit(1);
    }
}
static bool isGZfile(const std::string fi) {
	std::string subst = fi.substr(fi.length() - 3);
	if (subst == ".gz") {
		return true;
	}
	return false;
}

static bool check_if_file_is_readable(string option_name, string path) {
    ifstream test_file(path);
    if (test_file.good()){
        return true;
    }
    else {
        _log(logERR) << "Option: \"" << option_name << "\" must be accessible and readable.";
        exit(1);
    }
}

static bool check_if_file_is_writable(string option_name, string path){
    ofstream file;
    try{
        file.open(path.c_str(), ios::out | ios::trunc);
        file.close();
        return true;
    } catch (ios_base::failure){
        _log(logERR) << "Option: \"" << option_name << "\" must be accessible and writable.";
        exit(1);
    }
}



static bool check_if_one_of(string option_name, string value, vector<string> valid_options){
    for(string valid_opt: valid_options){
        if( value == valid_opt){
            return true;
        }
    }
	_log(logERR) << "Option: \"" << option_name << "\" must be one of:<defunct>\n";
	//<< boost::algorithm::join(valid_options, ", ");
    exit(1);
}
 




