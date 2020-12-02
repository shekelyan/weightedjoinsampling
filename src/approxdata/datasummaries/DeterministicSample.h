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

#ifndef SHEKELYAN_RANDOM_SAMPLE_H
#define SHEKELYAN_RANDOM_SAMPLE_H

#include <approxdata/datasummaries/datasummaries.h>


#include <list>
#include <queue>

using namespace std; 

class DeterministicSampleSummary : public DataSummary{
	
	
public: 
	
	const int DIMS;
	
private:

	const Params params;
	shared_ptr<DataSet> sample;
	const long datasize;
	
	double samplePointWeight;
	
	double maxAbsErr = -1;
	
public:
	
	inline double getSamplePointWeight() const{
	
		return samplePointWeight;
	}
	
	inline double getSampleSize() const{
	
		return sample->getSize();
	}
	
	inline double getBytes() const override{
	
		return getSampleSize()*8*DIMS;
	
	}
	
	inline double getDataSize() const{
		
		return datasize;
	}
	
	
	shared_ptr<vector<Point>> getBetween( const DataSet& data, const Point& a , const Point& b) const{
	
		shared_ptr<vector<Point>> ret(new vector<Point>() );
		
		return ret;
	}
	
	inline void sort(vector<Point>& v, const PointComparator& pc, long first = 0, long last=-1){
	
		if (last == -1)
			last = v.size();
	
		sort(v.begin()+first, v.begin()+last, pc);
		
	}
	
	
	
	inline DeterministicSampleSummary(const DataSet& data, string paramStr) : DIMS(data.getDims()), params(paramStr), datasize(data.getSize() ){

		const double kbPerSample = DIMS*8.0/(1000.0);
		const double sampleSize = UtilMath::minVal<long>(datasize, params.get("kb")/kbPerSample );
		
		assert (sampleSize > 0);
		
		sample = data.getRandomSubset(sampleSize, to_string(params.get("seed")) );
		
		sample->setVerbose(false);
		samplePointWeight = 1.0/getSampleSize()*getDataSize();
		
		
		long quantileInserts = 0;
		long quantileLookups = 0;
		
		if (params.get("verifykb") > 0){
			
			for (int k = 1; k <= 2; k++){
			
				if (DIMS == 2 && k == 2)
					continue;
				
				SliceHistSummary sh(sample, "slicehist -k "+to_string(k)+" -kb "+to_string(params.get("verifykb") ));
		
				const double epsn = sh.getEpsilonDistance(data);
				
				if (k == 1 || epsn < maxAbsErr){
					
					maxAbsErr = epsn;
					
					observer["eps"] = maxAbsErr/getDataSize();
					observer["k"]  = k;
					observer["verifybytes"] = sh.getBytes();
					
				}
				
				quantileInserts += sh.getQuantileInserts();
				quantileLookups += sh.getQuantileLookups();
			}
		}
		
		observer["quantileInserts"] = quantileInserts;
		observer["quantileLookups"] = quantileLookups;
	}
	
	inline RandomSampleSummary(shared_ptr<DataSet> data, string paramStr) : RandomSampleSummary( (*data), paramStr){
		
	}
	
	inline RandomSampleSummary(unique_ptr<DataSet> data, string paramStr) : RandomSampleSummary( (*data), paramStr){
		
	}
	
	inline void drawRandomPoint(RandomGenerator& rg, Point& p) const override{
	
		p = (*sample)[rg.randomInteger(0, sample->getSize() )];
	}
	
	~RandomSampleSummary() override {};
	
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode = QueryMode::EST) const override{
		
		if (getDataSize() == 0)
			return 0;
			
		double ret = 0;
		
		if (mode == QueryMode::LB){
		
			if (maxAbsErr >= 0 && maxAbsErr <= getDataSize() )
				ret = sample->boxCount(qmin, qmax)*samplePointWeight-maxAbsErr;
			else 
				ret = sample->boxCount(qmin, qmax);
				
		} else if (mode == QueryMode::UB){
		
			if (maxAbsErr >= 0 && maxAbsErr <= getDataSize() )
				ret = sample->boxCount(qmin, qmax)*samplePointWeight+maxAbsErr;
			else
				ret = datasize-(getSampleSize()-sample->boxCount(qmin, qmax));
				
		} else {
			
			ret = sample->boxCount(qmin, qmax)*samplePointWeight;
		}
		
		return UtilMath::makeBetween<double>(0, getDataSize(), ret);
	}
	
	inline const string toString() const override{
		
		if (params.get("verifykb") > 0)
			return "verifiedsampling "+params.toString();
		
		if (params.get("vkb") > 0)
			return "vsampling "+params.toString();
		
		
		return "randomsampling "+params.toString();
	}
};

#endif