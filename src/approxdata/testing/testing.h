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

#ifndef SHEKELYAN_TESTING_H
#define SHEKELYAN_TESTING_H


//#include <approxdata/utils/utils.h>
//#include <approxdata/data/data.h>
#include <approxdata/datasummaries/datasummaries.h>

//#include <approxdata/visualization/DataImage.h>

#include <approxdata/testing/DataImport.h>
#include <approxdata/testing/BoxIndex.h>
#include <approxdata/testing/QueryBuilder.h>
#include <approxdata/testing/Experiments.h>


class Tester{


};

inline void evaluate(const DataSet& data, const DataSummary& summary, const RangeQueries& queries, int num=1000000000, double scale=1.0){
	
	Timer timer("runQueries");
	timer.start(0, queries.ranges.size()/10, queries.ranges.size() );
	
	Stats ests("est");
	
	Stats counts("tru");
	
	Stats absErr("aerr", "%");
	Stats absWidth("awid", "%");
	
	Stats relErr("err", "%");
	Stats relWidth("wid", "%");
	
	int c = 0;
	for (auto it = queries.ranges.begin(); it != queries.ranges.end(); it++){
	
		const double t = it->count;
	
		const double lb = summary.boxCount( it->min, it->max, QueryMode::LB )*scale;
		const double est = summary.boxCount( it->min, it->max, QueryMode::EST )*scale;
		const double ub = summary.boxCount( it->min, it->max, QueryMode::UB )*scale;
		
		if ( timer.tick() || (ub < t) || (lb > t))
			cout << "LB " << lb << " EST " << est << " UB " << ub << " T " << t << " qmin " << it->min << " qmax " << it->max << endl;
		
		//if ( (ub < t) || (lb > t) ){
		//	cout << data.boxCount(it->min, it->max) << endl;
		//}
		
		assert (lb <= t);
		assert (ub >= t);
		
		ests.add(est);
		counts.add(t);
		
		absErr.add( 100*UtilMath::absVal(est-t)/data.getSize() );
		relErr.add( 100*UtilMath::absVal(est-t)/t );
		relWidth.add( 100*UtilMath::absVal(ub-lb)/t );
		
		absWidth.add( 100*UtilMath::absVal(ub-lb)/data.getSize() );
		
		c++;
		
		if (c >= num)
			break;
	}
	
	timer.end();
	cout << "/runQueries" << endl;
	cout << absErr << endl;
	cout << absWidth << endl;
	cout << relErr << endl;
	cout << relWidth << endl;
	cout << ests << endl;
	cout << counts << endl;

}

inline void evaluate(shared_ptr<DataSet> ds, const DataSummary& summary, const RangeQueries& queries, int num=1000000000, double scale=1.0){

	evaluate( (*ds), summary, queries, num, scale);
}

#endif