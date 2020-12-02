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
#ifndef SHEKELYAN_RANKTRANSFORM_H
#define SHEKELYAN_RANKTRANSFORM_H

#include <approxdata/external/gk.h>
#include <approxdata/external/gk2.h>
#include <approxdata/datasummaries/datasummaries.h>

class RankTransform{
public:
	vector<unique_ptr<RankSummary>> pointRankSummaries;
	
	//vector<unique_ptr<ValueRankSummary<double> >> valueRankSummaries;
	
	bool finalised = false;
	
	const int DIMS;
	
	
	inline void init(const RankSummaryParameters& pars){
	
		for (int j = 0; j < pointRankSummaries.size(); j++){
		
			if (pointRankSummaries.at(j) ){
			
				pointRankSummaries.at(j)->init(pars);
			}
		}
	}
	
	/*
	inline long getMaxCountOverestimation() const{
	
		long max = 0;
		for (int j = 0; j < pointRankSummaries.size(); j++){
		
			if (pointRankSummaries.at(j) )
				max = UtilMath::maxVal( max, pointRankSummaries.at(j)->getMaxCountOverestimation() );
		}
		
		return max;
	}*/
	
	inline void debug(){ //const DataSet& data){
	
		/*
		for (int j = 0; j < pointRankSummaries.size(); j++){
		
			if (pointRankSummaries.at(j) )
				pointRankSummaries.at(j)->debug(100);
		}*/
		
	}
	
	inline RankTransform(int dims) : DIMS(dims){
	
		for (int i = 0; i < DIMS; i++){
			
			Point dmin(DIMS);
			Point dmax(DIMS);
			
			dmin.id = UtilMath::lowestDoubleInteger();
			dmax.id = UtilMath::largestDoubleInteger();
			
			for (int ii = 0; ii < DIMS; ii++){
			
				dmin[ii] = 0.0;
				dmax[ii] = 1.0;
			}
			
			//cout << "DOMAIN MIN " << dmin << endl;
			//cout << "DOMAIN MAX " << dmax << endl;
			
			
			unique_ptr< RankSummary> q(new RankSummary(i,dmin, dmax, PointComparator(i, true)));
			pointRankSummaries.push_back( std::move(q) );
		}
	}
	
	inline void add(const Point& p){
	
		assert (!finalised);
		
		for (int i = 0; i < DIMS; i++)
			pointRankSummaries.at(i)->add(p);
			
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		const int size = pointRankSummaries.size();
	
		for (int i = 0; i < size; i++){
			
			s << "quantile summary in i = " << to_string(i) << endl;
			
			s << (*pointRankSummaries[i]);
			
		}
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const RankTransform& m) { 
		
		return os << m.toString();
	}
	
	
	inline void finalise(){
		
		assert (!finalised);
	
		for (int i = 0; i < DIMS; i++){
		
			pointRankSummaries[i]->finalise();
		}
		
		finalised = true;
	}
	
	inline double getRank(const Point& p, int i, QueryMode mode=QueryMode::EST) const{
	
		return pointRankSummaries[i]->getRank(p, mode);
		
	}
	
	inline double getNormalizedRank(const Point& p, int i, QueryMode mode=QueryMode::EST) const{
	
		
		
		return pointRankSummaries[i]->getNormalizedRank(p, mode);
	}
	
	inline void rankTransform(const Point& in, Point& out) const{
	
		assert (finalised);
	
		for (int i = 0; i < DIMS; i++)
			out[i] = getRank(in, i);
	}
	
	inline double getBytes() const{
	
		double ret = 0;
		
		for (int i = 0; i < DIMS; i++){
		
			if (pointRankSummaries.size() > 0)
				ret += pointRankSummaries[i]->getBytes();
		}
		return ret;
	}
	
	inline bool less(int i, const Point& a, const Point& b) const{
	
		if (pointRankSummaries.size() > 0)
			return pointRankSummaries[i]->compare(a, b);
		
		return a[i] < b[i];
	}
	
	inline void normRankTransform(const Point& in, Point& out, QueryMode mode=QueryMode::EST) const{
	
		assert (finalised);
	
		for (int i = 0; i < DIMS; i++)
			out[i] = getNormalizedRank(in, i, mode);
			
		out.id = in.id;
	}
	
	
	inline double getQuantileNum(int i) const{
	
		return pointRankSummaries[i]->getQuantileNum();
		
	}
	
	inline double getQuantileInfimum(const long q, int i) const{
	
		return pointRankSummaries[i]->getQuantileInfimum(q)[i];
		
	}
	
	inline double getQuantileMax(const long q, int i) const{
	
		return pointRankSummaries[i]->getQuantileMax(q)[i];
		
	}
	
	inline void normRankTransformInv(const Point& in, Point& out) const{
	
		assert (finalised);
	
		for (int i = 0; i < DIMS; i++){
		
			const double q = in[i]*getQuantileNum(i);
			
			const long qlong = q;
			
			const double w = q-qlong;
			
			out[i] = getQuantileInfimum(q,i)*(1-w)+getQuantileMax(q,i)*w;
		}
	}
	
	
	
	inline void normRankTransform(const QueryBox& in, QueryBox& out) const{
	
		Point tempPoint(DIMS);
	
		assert (finalised);
		
		tempPoint = in.min;
		
		for (int i = 0; i < DIMS; i++){
		
			if (in.isMinUnbounded(i)){
				tempPoint[i] = 0;
				tempPoint.id = UtilMath::lowestDoubleInteger();
			}
		}
		
		normRankTransform(tempPoint, out.min);
		
		tempPoint = in.max;
		
		for (int i = 0; i < DIMS; i++){
		
			if (in.isMaxUnbounded(i)){
				tempPoint[i] = 1;
				tempPoint.id = UtilMath::largestDoubleInteger();
			}
		}
		
		normRankTransform(tempPoint, out.max);
		
		out.mode = in.mode;
		out.flags = in.flags;
		
		for (int i = 0; i < DIMS; i++){
		
			if (in.isMinUnbounded(i))
				out.min[i] = 0;
			
			if (in.isMaxUnbounded(i))
				out.max[i] = 1;
		}
	}
	
	inline void normRankTransform(const Point& in1, const Point& in2, Point& out1, Point& out2, QueryMode mode=QueryMode::EST) const{
	
		assert (finalised);
		
		Point p1 = in1;
		Point p2 = in2;
		
		p1.id = UtilMath::lowestDoubleInteger();
		p2.id = UtilMath::largestDoubleInteger();
		
		if (mode == QueryMode::LB){
		
			normRankTransform(p1, out1, QueryMode::UB);
			normRankTransform(p2, out2, QueryMode::LB);
			
		} else if (mode == QueryMode::EST){
		
			normRankTransform(p1, out1);
			normRankTransform(p2, out2);
				
		} else if (mode == QueryMode::UB){
		
			normRankTransform(p1, out1, QueryMode::LB);
			normRankTransform(p2, out2, QueryMode::UB);	
		}
		
		
		
	}
	inline void normRankTransform(const Box& in, Box& out,QueryMode mode=QueryMode::EST) const{
	
		normRankTransform(in.min, in.max, out.min, out.max,mode);
	}
	
	
};

#endif