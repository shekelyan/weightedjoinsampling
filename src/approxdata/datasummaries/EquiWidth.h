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

#ifndef SHEKELYAN_EQUIWIDTH_H
#define SHEKELYAN_EQUIWIDTH_H

#include <approxdata/datasummaries/datasummaries.h>

class EquiwidthHist{

public:
	DenseGrid g;
	const int DIMS;
	long state = -1;
	
	inline EquiwidthHist(int dims, long n) : DIMS(dims), g(dims, n, false){
	
		
	}
	
	inline double getBytes(){
	
		return g.getBytes();
	}

	inline void setState(int s){
	
		assert (s == state+1);
	
		state = s;
		
		if (state == 0)
			g.allocate();
			
		if (state == 1)
			g.finalise();
	}


	inline long getGridCoord(const Point& p) const{
		
		long ret = 0;
		
		for (int i = 0; i < DIMS; i++)
			ret = g.lc.setDimCoord(ret, i, UtilHist::discretize(p[i], g.getCellsPerDim() ) );
		
		return ret;
	}

	inline void add(const Point& p){
	
		assert (state == 0);
	
		g.addCount( getGridCoord(p), 1 );	
	}
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const{
		
		assert (state == 1);
		
		Box b(qmin, qmax);
		
		const long ubMin = getGridCoord(qmin);
		const long ubMax = getGridCoord(qmax);
		
		const double ubCount = UtilMath::isIn(-1L, ubMin, ubMax) ? 0 : g.sum(ubMin, ubMax);
		
		if (ubCount == 0)
			return 0;
		
		if (mode == QueryMode::UB)
			return ubCount;
		
		const long lbMin = g.lc.incAllCoords(ubMin);
		const long lbMax = g.lc.decAllCoords(ubMax);
		
		const long lbCells = g.lc.getArea(lbMin, lbMax);
		
		const double lbCount = lbCells > 0 ? ( UtilMath::isIn(-1L, lbMin, lbMax) ? 0 : g.sum(lbMin, lbMax) ) : 0;
		
		if (mode == QueryMode::LB)
			return lbCount;
		
		const double qVol = b.getVolume();
		
		if (qVol == 0)
			return 0;
			
		const long ubCells = g.lc.getArea(ubMin, ubMax);
		
		const double lbVol = UtilMath::isIn(-1L, lbMin, lbMax) ? 0 : lbCells*1.0/pow(g.getCellsPerDim(), DIMS);
		const double ubVol = UtilMath::isIn(-1L, ubMin, ubMax) ? 0 : ubCells*1.0/pow(g.getCellsPerDim(), DIMS);
		
		return lbCount+(ubCount-lbCount)*(qVol-lbVol)/(ubVol-lbVol);
	}
};

class EquiwidthSummary : public DataSummary{

public:

	const int DIMS;
	Params params;
	shared_ptr<EquiwidthHist> h;
	
	inline double getBytes() const override{
	
		return h->getBytes();
	}
	
	inline EquiwidthSummary(const DataSet& data, const string s) : DIMS(data.getDims() ), params(s) {
		
		
		double bytes = params.get("kb")*1024;
		double bytesPerBucket = 8;
		
		const double buckets = bytes/bytesPerBucket;
		
		const double bucketsPerDim = pow(buckets, 1.0/DIMS);
		
		h.reset( new EquiwidthHist(DIMS, bucketsPerDim));
		
		h->setState(0);
			
		for (auto it = data.begin(); it != data.end(); it++)
			h->add( (*it) );
		
		h->setState(1);
	}
	
	inline EquiwidthSummary(shared_ptr<DataSet>& data, const string s) : EquiwidthSummary( (*data), s) {
	
	}
	
	inline virtual double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override{
	
		return h->boxCount(qmin, qmax, mode);
	}
	
	inline const string toString() const override{
		
		return "equiwidth "+params.toString();
	}
	
};


#endif