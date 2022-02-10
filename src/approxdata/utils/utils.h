/*
 * ApproxData Library
 * Copyright (c) 2018 Michael Shekelyan <michael.shekelyan@gmail.com>
 */

/*
Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions 
of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef HEADER_UTILS_H
#define HEADER_UTILS_H

#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif

#include <memory>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <time.h>

#include <iostream>
#include <algorithm>
#include <limits>


#include <assert.h>
#include <vector>
#include <random>
#include <fstream>
#include <string>

using std::ifstream;
using std::ostream;
using std::ofstream;
using std::array;
using std::string;
using std::stringstream;
using std::istringstream;

using std::shared_ptr;
using std::unique_ptr;

using std::ostringstream;

#include <bitset>

//using namespace std;

using std::vector;
using std::array;
using std::unordered_map;
using std::map;
using std::unordered_set;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::numeric_limits;
using std::to_string;


#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <initializer_list>

#include <iostream>     // std::cout
#include <algorithm>    // std::shuffle
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono> 

struct Exception : public std::exception
{
   std::string s;
   Exception(std::string ss) : s(ss) {}
   ~Exception() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
   
   inline void print(){
   
   		cout << "exception " << what() << endl;
   }
};

// classes:
// namespaces: 
#include <approxdata/utils/UtilCollection.h>

// deprecated classes: Stats, Interval, BinarySearcher, NormGenerator, RealFunction
// classes: BinarySearch, RandomGenerator, OneDimEstimator, Histogram1d, Equiwidth1d, BasicAggregator
// namespaces: UtilMath, UtilHist, UtilKS
// enums: QueryMode

#include <approxdata/utils/UtilMath.h>


#include <approxdata/utils/UtilString.h>

#include <approxdata/utils/UtilTime.h>
#include <approxdata/utils/UtilHash.h>
#include <approxdata/utils/UtilPoint.h>




//#include <approxdata/utils/Timer.h>

#include <approxdata/utils/UtilBits.h>
#include <approxdata/utils/UtilGrid.h>


#include <approxdata/utils/ZCurveCoords.h>

#include <approxdata/utils/MultiSet.h>


#include <approxdata/utils/UtilJson.h>




#endif
