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

#ifndef SHEKELYAN_DYADICSKETCH_H
#define SHEKELYAN_DYADICSKETCH_H

#include <approxdata/datasummaries/datasummaries.h>






class DyadicSketchSummary : public DataSummary{

public:
	
	unique_ptr<DyadicCMSketch> sketches;
	Params params;
	
	double eps1;
	double eps;
	
	const int DIMS;
	
	
	
	
	inline long getSize() const{
	
		if (sketches)
			return sketches->size;
	
		return 0;
	}
	
	inline double getBytes() const override{
	
		return getSize()*8;
	}
	
	const long datasize;
	
	inline void drawRandomPoint(RandomGenerator& rg, Point& p) const override{
	
		
	}
	
	inline void convert(const Point& p, vector<long>& v) const{
		
		
		for (int i = 0; i < DIMS; ++i)
			v[i] = sketches->toInt(p[i]);
			
			/*
			const float f = static_cast<float>(p[i]);
			
			const int& f_int = reinterpret_cast<const int&>(f);
			
			v[i] = f_int;*/
		
	}
	
	inline DyadicSketchSummary( const DataSet& data, const string& paramStr) : datasize(data.getSize() ), DIMS(data.getDims() ),  params(paramStr){
	
		eps = params.get("eps")/100.0;
		const double kb = params.get("kb");
		
		
		const double nofast = params.get("nofast");
		
		
		const double m = params.get("m");
		
		
		sketches.reset( new DyadicCMSketch(DIMS, UtilMath::makeMultipleOfUB(32, m), kb, eps, 0.01, nofast == 0, m == 0 ? 1 : m) );
		
		
		observer["eps"] = eps;
		
		cout << toString() << endl;
		
		const int dataScans = 1;
		
		Timer t("dyadicsketch");
		
		t.start(10, 1000000, data.getSize()*dataScans);
		
		vector<long> p(DIMS);
		
		for (int i = 0; i < dataScans; i++){
			
			const string key = "datascan"+to_string(i+1);
			
			observer.begin(key);
		
			for (auto it = data.begin(); it != data.end(); it++){
			
				t.tick();
			
				convert(*it, p);
				sketches->add(p);
			}
			
			observer.end(key);
		}
		
		t.end();
		
		
		cout << "getSize() " << getSize() << endl;
		
		observer["bytes"] = getBytes();

	}
	
	inline DyadicSketchSummary( const shared_ptr<DataSet>& data, const string& paramStr) : DyadicSketchSummary( (*data), paramStr){
	
	}
	
	inline DyadicSketchSummary( const unique_ptr<DataSet>& data, const string& paramStr) : DyadicSketchSummary( (*data), paramStr){
	
	}
	
	inline double boxCount(const Point& bmin, const Point& bmax, QueryMode mode=QueryMode::EST) const override{
		
		if (sketches){
		
			vector<long> qmin(DIMS);
			vector<long> qmax(DIMS);
		
			convert(bmin, qmin);
			convert(bmax, qmax);
			
			return mode == QueryMode::LB ? sketches->lowerbound(qmin, qmax) : mode == QueryMode::UB ? sketches->upperbound(qmin, qmax) : sketches->estimate(qmin, qmax);
			
		} else {
		
			return 0;
		}
	}
	
	inline void print() const{
		
	}
	
	inline void test() const{
	
	
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const DyadicSketchSummary& m) { 
    	
    	return os << m.toString();
	}
	
	inline const string toString() const {
		
		return "dyadicsketch "+params.toString();
	}
};

#endif