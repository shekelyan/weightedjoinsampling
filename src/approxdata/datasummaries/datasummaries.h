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

#ifndef SHEKELYAN_DATASUMMARIES_H
#define SHEKELYAN_DATASUMMARIES_H

#include <approxdata/utils/utils.h>


class TheoreticalSummaryModel{

public:
	
	inline virtual const string toString() const { return ""; }

		
	inline virtual long double getBytes() const = 0;

	inline long double getMegabytes(){
	
		const long double kb = getBytes()/1000.0L;
		const long double mb = kb/1000.0L;
	
		return mb;
	}
	
	inline virtual void setEps(long double eps) = 0;
	
	inline virtual void reset() {
	
	
	}
	
	inline virtual long double getAsymptotic(long double eps) const{
	
		return 0;
	}
	
	inline virtual long double getAsymptotic2(long double eps) const{
	
		return 0;
	}
	
	inline virtual int getDims() const = 0;
	
	inline virtual long double getEps() const = 0;
	
	inline virtual long getQuantileOps() const {
	
		return 0;
	}
	
	inline virtual bool setQueryDimensionality(int d){ return false;}
	
	inline virtual void setBytes(long double bytes) = 0;
	
	inline virtual void setDataSize(long n) = 0;
	
	inline virtual long getDataSize() const = 0;
	
	inline virtual const string getName() const = 0;
	
	
	inline virtual const string getParameterLatex() const{
	
		return "";
	};
	
	inline virtual void setDataBytes(long double bytes){
	
		long double bytesPerPoint = 8*getDims();
	
		setDataSize( ceil(bytes/bytesPerPoint) );
	}
};



#include <approxdata/data/data.h>


#include <approxdata/datasummaries/DataSummary.h>

#include <approxdata/datasummaries/UtilDyadic.h>
#include <approxdata/datasummaries/DyadicSketch.h>
#include <approxdata/datasummaries/EquiWidth.h>

#include <approxdata/datasummaries/RankSummary.h>
#include <approxdata/datasummaries/RankTransform.h>
#include <approxdata/datasummaries/GridHist.h>
#include <approxdata/datasummaries/DigitHist.h>
#include <approxdata/datasummaries/EquiDepth.h>
#include <approxdata/datasummaries/DyadicHist.h>
#include <approxdata/datasummaries/SliceHist.h>
#include <approxdata/datasummaries/RandomSample.h>
#include <approxdata/datasummaries/VaryWidth.h>


class TheoreticalHybrid{

public:

	int parentSummary;
	long p;
	int DIMS;
	
	long double bytes;
	long double eps;
	
	
	const int EQUIDEPTH = 1;
	const int DYADIC = 2;
	const int SLICEHIST = 3;
	const int HYBRID = 4;
	
	inline TheoreticalHybrid(int dims) : DIMS(dims){
	
		
	
	}
	
	inline const string toString() const{
 	
 		stringstream s;
 	
 		s << "hybrid";
 		
 		s << "\t -size " << UtilString::removeAll(UtilString::bytesToString(bytes), ' ');
 		
 		s << "\t -eps " << UtilString::doubleToString(eps*100, 10) << "%";
 		//s << "\t -eps1 " << UtilString::doubleToString(eps1*100, 10) << "%";;
 		
 		s << "\t -p " << p;
 		
 		
 		return s.str();
 	}

	
	
	inline long double getBytes(int mode, long double eps, int dims) const{
	
		if (dims == 1){
		
			return (ceil(2.0/eps)+1)*(DIMS+1)*8;
		}
	
		if (mode == 1){
		
			TheoreticalEquidepth t(dims);
			t.setEps(eps);
			return t.getBytes();
		}
		
		if (mode == 2){
		
			TheoreticalDyadicHist t(dims);
			t.setEps(eps);
			return t.getBytes();
		}
		
		if (mode == 3){
		
			TheoreticalSliceHist t(dims);
			t.setEps(eps);
			return t.getBytes();
		}
		
		if (mode == 4){
		
			TheoreticalHybrid t(dims);
			t.setEps(eps);
			
			return t.getBytes();
		}
		
		assert(false);
	}
	
	
	inline int getBestMode(double eps, int dims) const{
	
		if (dims == 1)
			return EQUIDEPTH;
	
		
		int bestMode = 0;
		long double bestBytes = std::numeric_limits<double>::max();
	
		for (int i = 1; i <= 4; i++){
		
			const long double b = getBytes(i, eps, dims);
			
			if (b < bestBytes){
			
				bestMode = i;
				bestBytes = b;
			}
		}
		
		return bestMode;
	}
	
	
	inline long double getBytes() const{
	
		return bytes;
	}

	inline long double getMegabytes() const{
	
		return bytes/1000000;
	}
	
	inline void setEps(long double eps){
	
		
		long double bestBytes = std::numeric_limits<long double>::max();
		int bestMode = 0;
		int bestP = 0;
		
		if (false)
		for (int i = 1; i <= 3; i++){
		
			const long double b = getBytes(i, eps, DIMS);
			
			if (b < bestBytes){
			
				bestMode = i;
				bestBytes = b;
			}
		}
		
		this->eps = eps;
		
		for (int p = 4; p < 30; p++){
		
			const long normParts = 4;
			const long subParts = 2*(p-2);
			
			assert (subParts >= 0);
			
			const long double q = 1.0L/(1L << p);
			
			if (4*q >= eps)
				continue;
			
			long double b = ((1L << p)+1)*(DIMS+1)*8;
			
			// eps*N = 4*q*N + 2*sum[i=2, p-1] (q*N)*2^i *x
			
			// y = 2^i x
			
			// eps*N = 4*q*N + sum[i=2, p-1] (q*N)*y
			
			// eps*N = 4*q*N + 2*(p-2)* (q*N)*y
			
			// eps*N = (4+2*(p-2)y) q*N
			
			// eps/q = 4+2*(p-2)y
			
			// y = (eps/q-4 )/2*(p-2)
			
			// x= ((eps/q-4 )/2*(p-2))/2^i
			
			
			for (int i = 2; i <= (p-1); i++){
			
				const long len = 1L << i;
				const long double childEps = ((eps/q-4)/(2*p-4))/len;
				
				const long num = 1L << (p-i);
				
				assert (childEps > 0);
				
				if (childEps < 1){
					
					const int bestChildMode = getBestMode(childEps, DIMS-1);
					b += num * getBytes(bestChildMode, childEps, DIMS-1);
				}
			}
			
			if (b < bestBytes){
			
				bestMode = HYBRID;
				bestP = p;
				bestBytes = b;
			}
		}
		
		parentSummary = bestMode;
		p = bestMode == HYBRID ? bestP : -bestMode;
		bytes = bestBytes;	
	}


};

inline shared_ptr<DataSummary> newDataSummary(const DataSet& data, const string s){

	shared_ptr<DataSummary> ret;
	
	//try {
	
		const Params p(s);
		
		if (UtilString::equals(p.name, "equiwidth"))
			ret.reset(new EquiwidthSummary(data, s) );
			
		if (UtilString::equals(p.name, "equidepth"))
			ret.reset(new EquidepthSummary(data, s) );
		
		if (UtilString::equals(p.name, "randomsampling"))
			ret.reset(new RandomSampleSummary(data, s) );
		
		if (UtilString::equals(p.name, "randomsample"))
			ret.reset(new RandomSampleSummary(data, s) );
			
			
		if (UtilString::equals(p.name, "digithist"))
			ret.reset(new DigitHistSummary(data, s) );
			
		if (UtilString::equals(p.name, "varywidth"))
			ret.reset(new VarywidthSummary(data, s) );
		
		if (UtilString::equals(p.name, "dyadichist"))
			ret.reset(new DyadicHistSummary(data, s) );
			
		if (UtilString::equals(p.name, "dyadicsketch"))
			ret.reset(new DyadicSketchSummary(data, s) );
		
		if (UtilString::equals(p.name, "slicehist"))
			ret.reset(new SliceHistSummary(data, s) );
		/*
	} catch (...){
	
		cout << "exception during summary construction!!" << endl;
	
	}*/
	
	return ret;
}

class ElemHist{

	
	
	inline vector<Box> getBox(const Box& box, long xres){
	
	
		
	}


};

/*
inline void equidepth(vector<Point>& data, long xres, long yres){

	const int dims = 2;
	
	PointComparator xcomp(0, true);
	PointComparator ycomp(1, true);
	
	


}*/

shared_ptr<DataSet> getSample(const DataSet& data, long size){
	
	const long datasize = data.getSize();
	
	long L = floor( pow(size, 1.0/3) );
	
	long m = datasize / (L*L*L);
	
	
	shared_ptr<vector<Point>> sample(new vector<Point>() );
	
	long c = 0;
	
	while ( (m*(L)*(L)*(L+c+1)) <= datasize)
		c++;
		
	cout << "L " << L << " c " << c << " m " << m << " data " << datasize << endl;
	
	
	const long quantileNum = (L+c);
		
	const long s = (L)*(L)*quantileNum;
	
	cout << datasize-(m*s) << endl;
	
	
	RandomSkipper skipper(datasize, datasize-(m*s), "test");
	
	
	const long pointsPerBigBlock = m*(L)*quantileNum;
	const long pointsPerSmallBlock = m*quantileNum;
	
	const long pointsPerBucket = m;
	
	Point min(2, 0);
	Point max(2, 1);
	
	cout << min << endl;
	cout << max << endl;
	
	PointComparator xcomp(0, true);
	PointComparator ycomp(1, true);
	
	vector<Point> lastOfBigBlocks;
	
	const double eps1 = 0.01;
	
	/*
	if (true){
	
		
		RankSummary qs(0, min, max, xcomp);
		
		RankSummaryParameters params;
		
		params.eps = eps1;
		
		qs.init(params);

		qs.setState(1);
		
		long c = 0;
		
		
		skipper.reset();
		
		for (auto it = data.begin(); it != data.end(); it++){
		
			if (skipper.skip())
				continue;
		
			if (c%100000 == 0)
				cout << c << endl;
		
			qs.add( (*it) );
			c++;
		}
		
		qs.setState(2);
		
		for (long L1 = 1; L1 <= L; L1++){
		
			lastOfBigBlocks.push_back( qs.getQuantileForRank( L1*pointsPerBigBlock + 2*eps1*data.getSize() ) );
		}
	}*/
	
	vector<Point> q1;
	vector<Point> q2;
		
	vector<Point> ps(pointsPerBigBlock*2);
	
	Point next(2, -1);
	
	RandomGenerator g("test");
	
	BipartiteMatcher matcher(quantileNum, quantileNum);
	
	
	
	for (long L1 = 0; L1 < L; L1++){
	
		const Point& a = next;
		const Point& b = lastOfBigBlocks[L1];
		
		ps.clear();
		
		long pre = 0;
		
		skipper.reset();
		
		for (auto it = data.begin(); it != data.end(); it++){
		
			if (skipper.skip())
				continue;
		
			const Point& p = (*it);
			
			if (xcomp(p, a)){
			
				pre++;
				continue;
			}
			
			if (L1 < (L-1) && xcomp(b, p))
				continue;
				
				
			ps.push_back(p);
		}
		
		
		cout << "pre " << pre << " between " << ps.size() << " pointsPerBigBlock " << pointsPerBigBlock << endl;
		
		xcomp.sort(ps.begin(), ps.end() );
		
		if (L1 < (L-1)){
		
			next = ps[pointsPerBigBlock];
		}
		
		ps.resize(pointsPerBigBlock);
		
		ycomp.sort(ps.begin(), ps.end());
		
		for (long L2 = 0; L2 < L; L2++){
		
			const long s = L2 * pointsPerSmallBlock;
			
			q1.clear();
			q2.clear();
			
			auto it1 = ps.begin()+s;
			auto it2 = ps.begin()+(s+pointsPerSmallBlock);
			
			{
				xcomp.sort(it1, it2);
				
				for (long k = pointsPerBucket; k < pointsPerSmallBlock; k += pointsPerBucket)
					q1.push_back( ps[s+k] );
			}
			
			{
				ycomp.sort(it1, it2);
				
				for (long k = pointsPerBucket; k < pointsPerSmallBlock; k += pointsPerBucket)
					q2.push_back( ps[s+k] );
			}
			
			MultiArray<bool> qc(quantileNum, quantileNum);
			
			qc.fill(false);
			
			for (long k = 0; k < pointsPerSmallBlock; k++){
			
				const Point& p = ps[s+k];
				
				const long qx = xcomp.countLess(q1, p);
				const long qy = ycomp.countLess(q2, p);
				
				qc.set(qx, qy, true);
			}
			
			matcher.getMatching(qc);
			
			g.shuffle(it1, it2);
			
			for (long k = 0; k < pointsPerSmallBlock; k++){
			
				const Point& p = ps[s+k];
			
				const long qx = xcomp.countLess(q1, p);
				const long qy = ycomp.countLess(q2, p);
				
				if (qc.get(qx, qy)){
					
					sample->push_back(p);
					qc.set(qx, qy, false);
				}
			}
		}
	}
		
	return data.getCopy(sample);
}



#endif
