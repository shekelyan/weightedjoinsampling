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

#ifndef SHEKELYAN_GRIDHIST_H
#define SHEKELYAN_GRIDHIST_H

#include <approxdata/datasummaries/datasummaries.h>



template<bool CONTAIN>
class IntersectionFunction : public RealFunction{

public:
	Point& a;
	Point& b;
	
	inline IntersectionFunction(Box& box) : a(box.min), b(box.max){
	
		
	}
	
	inline ~IntersectionFunction() override {}
	
	inline double operator()(const double x) const override{
		
		const int DIMS = a.size();
		
		const double r = 0.5*pow(x, 1.0/DIMS);
		
		double y = 1.0/(1.0-r);
		
		if (CONTAIN){
		
			for (int i = 0; i < DIMS; i++)
				y *= UtilMath::minVal<double>(b[i]+r, 1-r)-UtilMath::maxVal<double>(a[i]-r, r);
		
		} else {
		
			for (int i = 0; i < DIMS; i++)
				y *= UtilMath::minVal<double>(a[i]+r, 1-r)-UtilMath::maxVal<double>(b[i]-r, r);
		}
		
		return y;
	}
};

class SparseGridHist : public MultiSet<long, long>{

private:
	
	UnifOneDimEstimator unif;
	vector<OneDimEstimator*> tempOneDimEsts;
	
	vector<long> ind;
	
public:

	ZCurveCoords z;
	const int DIMS;
	const long maxSparse;
	
	~SparseGridHist() override {};
	
	inline SparseGridHist(int dims, long sparse, long maxLod=62) : MultiSet(sparse), DIMS(dims), z(dims, maxLod), ind(dims), maxSparse(sparse){
		
		init(maxLod);
		
		for (int i = 0; i < DIMS; i++)
			tempOneDimEsts.push_back(&unif);
		
	}
	
	inline SparseGridHist(SparseGridHist& h) : SparseGridHist(h.DIMS, h.getLod(), h.maxSparse){
		
		elems = h.elems;
		counts = h.counts;
		total = h.total;
		
	}
	
	inline void write(IntegerOutputStream& out) const{
	
		//cout << "write" << endl;
		long last = -1;
		
		for (long j = 0; j < elems.size(); j++){
		
		
			const long v = elems[j];
			const long gap = v-last; 
			const long count = counts[j];
			
			out.write(gap);
			out.write(count);
			
			last = elems[j];
		}
		
		out.write(0);
		
		//cout << "/write" << endl;
	}
	
	inline double getBytes() const{
	
		double ret = 0;
	
		long last = -1;
		
		for (long j = 0; j < elems.size(); j++){
		
			const long v = elems[j];
			const long gap = v-last; 
			const long count = counts[j];
			
			ret += UtilBits::getNumOfVariableBytes(gap);
			ret += UtilBits::getNumOfVariableBytes(count);
			
			last = elems[j];
		}
		
		ret++;
		
		return ret;
	}
	
	
	inline void read(IntegerInputStream& in){
	
		//cout << "read" << endl;
		long last = -1;
		
		elems.clear();
		counts.clear();
		
		while (true){
		
			const long gap = in.read();
			
			if (gap == 0){
			
				//cout << "/read" << endl;
				return;
			}
			
			const long v = gap+last; // v-last == gap => gap+last == v
			const long count = in.read();
			
			elems.push_back(v);
			counts.push_back(count);
			
			last = v;
		}
	}
	
	inline long getCellNum() const{
	
		return 1L << z.getLod();
	}
	
	inline long getCellsPerDim(int i) const{
	
		return z.getCellsPerDim(i);
	}
	
	inline double getCellVol() const{
	
		return 1.0/getCellNum();
	}
	
	inline long getBucketNum() const{
	
		return elems.size();
	}
	
	inline long getBucketCount() const{
	
		return UtilMath::sum<long>(counts);
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		s << "gridhist ";
		s << "2^(" << z.getLod(0) << ")";
		
		for (int i = 1; i < DIMS; i++)
			s << " x " << "2^(" << z.getLod(i) << ")";
		
		s << " = " << "2^(" << z.getLod() << ")";
		s << " [non empty " << getBucketNum() << "]";
		
		s << " count " << getBucketCount();
		
		s << " elems " << UtilString::listToString(elems);
		s << " counts " << UtilString::listToString(counts);
		
		return s.str();
	}
	
	inline double volSelErr() const{
		
		double ret = 0.0;
		
		const double bucketVol = getCellVol();
		
		for (long j = 0; j < counts.size(); j++)
			ret += counts.at(j) * bucketVol;
		
		return ret;
	}
	
	
	inline double getUError(int steps=30) const{
	
		double ret = 0;
	
		//cout << "GG1" << endl;
	
		for (long j = 0; j < getBucketNum(); j++){
		
			Box b = getBox(j);
			
			double m = 0;
			
			//cout << "GG2" << endl;
			
			for (int i = 0; i < DIMS; i++){
				
				const double diff = b.max[i]-b.min[i];
				if (diff > m)
					m = diff;

			}
			
			//cout << "GG3" << endl;
			
			IntersectionFunction<true> f(b);
			
			//cout << "GG4" << endl;
			
			IntersectionFunction<false> g(b);
			
			//cout << "GG5" << endl;
			
			const double prob = UtilMath::integrate(f, 0, 1, steps)- UtilMath::integrate(g, m, 1, steps);
			
			//cout << "GG6" << endl;
			
			ret += counts.at(j)*prob;
			
			//cout << "GG7" << endl;
		}
		
		return ret;
	
	}
	
	inline Box getBox(long index) const{
	
		Point pmin(DIMS);
		Point pmax(DIMS);
		
		for (int i = 0; i < DIMS; i++){
		
			pmin[i] = UtilHist::bucketMin( z.getDimCoord(index, i), z.getCellsPerDim(i) );
			pmax[i] = UtilHist::bucketMax( z.getDimCoord(index, i), z.getCellsPerDim(i) );
		}
		
		Box b;
		b.enclose(pmin);
		b.enclose(pmax);
		
		return b;
	}
	
	inline void init(int lod){
	
		z.init(lod);
	}
	
	inline int getLod() const{
	
		return z.getLod();
	}
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode, const vector<OneDimEstimator*>& oneDimEsts) const{
	
		assert (oneDimEsts.size() == DIMS);
	
		const long zmin = getZ(qmin);
		const long zmax = getZ(qmax);
		
		assert (zmin != -1);
		assert (zmax != -1);
		
		const long zmin2 = z.incAllCoords(zmin);
		const long zmax2 = z.decAllCoords(zmax);
		
		Point tempmin(DIMS);
		Point tempmax(DIMS);
		
		for (int i = 0; i < DIMS; i++){
				
			const double bmin = UtilHist::bucketMin( z.getDimCoord(zmin, i), z.getCellsPerDim(i) );
			const double bmax = UtilHist::bucketMax( z.getDimCoord(zmin, i), z.getCellsPerDim(i) );
			
			const double imin = UtilMath::maxVal<double>( qmin[i], bmin);
			const double imax = UtilMath::minVal<double>( qmax[i], bmax);
			
			if (imin > imax){
			
				tempmin[i] = 0;
				
			} else {
			
				const double a = oneDimEsts[i]->intervalCount(imin, imax);
				const double b = oneDimEsts[i]->intervalCount(bmin, bmax);
				
				tempmin[i] =  b > 0 ? UtilMath::makeBetween<double>(0.0,1.0, (a/b)) : 0;
			}
		}
		
		for (int i = 0; i < DIMS; i++){
				
			const double bmin = UtilHist::bucketMin( z.getDimCoord(zmax, i), z.getCellsPerDim(i) );
			const double bmax = UtilHist::bucketMax( z.getDimCoord(zmax, i), z.getCellsPerDim(i) );
					
			const double imin = UtilMath::maxVal<double>( qmin[i], bmin);
			const double imax = UtilMath::minVal<double>( qmax[i], bmax);
			
			if (imin > imax){
			
				tempmax[i] = 0;
				
			} else {
			
				const double a = oneDimEsts[i]->intervalCount(imin, imax);
				const double b = oneDimEsts[i]->intervalCount(bmin, bmax);
				
				tempmax[i] = b > 0 ? UtilMath::makeBetween<double>(0.0,1.0, (a/b)) : 0;
			}
		}
		
		double ret = 0;
		
		for (int k = 0; k < elems.size(); k++){
		
			const long e = elems[k];
			const long c = counts[k];
			
			if (z.contains(zmin2, zmax2, e)){
				ret += c;
				continue;
			}
			
			if (mode == QueryMode::LB)
				continue;
			
			if (!z.contains(zmin, zmax, e))
				continue; 
			
			if (mode == QueryMode::UB){
				ret += c;
				continue;	
			}
			
			if (mode == QueryMode::EST){
			
				double frac = 1.0;
			
				for (int i = 0; i < DIMS; i++){
				
					if (z.sameDimCoord(zmin, e, i) )
						frac *= tempmin[i];
					else if (z.sameDimCoord(zmax, e, i))
						frac *= tempmax[i];
				}
				
				assert (frac >= 0);
				assert (frac <= 1);
				
				if (frac > 0)
					ret += frac*c;
			}
		}
		
		if (mode == QueryMode::EST){
		
			assert (ret <= boxCount(qmin, qmax, QueryMode::UB, oneDimEsts) );
			assert (ret >= boxCount(qmin, qmax, QueryMode::LB, oneDimEsts) );
		}
		
		return ret;
		
	}
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode) const{
	
		return boxCount(qmin, qmax, mode, tempOneDimEsts);
	}
	
	
	
	inline void reduceLod(int red=1){
	
		processUnsorted();
		processSorted();
		
		/*
		long check = 0;
		
		for (long k = 0; k < elems.size(); k++)
			check += (elems[k] >> red)*counts[k];
		*/
		
		for (int k = 0; k < elems.size(); k++)
			elems[k] >>= red;
		
		mergeDuplicates(elems, counts, tempElems, tempCounts);
		
		clear(elems, counts);
		
		elems = tempElems;
		counts = tempCounts;
		
		elems.shrink_to_fit();
		counts.shrink_to_fit();
		
		//for (long k = 0; k < elems.size(); k++)
		//	check -= elems[k]*counts[k];
		
		//assert (check == 0);
		
		clear(tempElems, tempCounts);
		
		init(getLod()-red);
	}
	
	inline long getZ(const Point& p) const{
	
		long ret = 0;
	
		for (int i = 0; i < DIMS; i++)
			ret = z.setDimCoord(ret, UtilHist::discretizePowerOfTwo(p[i], z.getLod(i) ), i);
			
		return ret;
		
	}
	
	inline void addZ(long z, long c){
	
		this->add(z, c);
		
		while (elems.size() > maxSparse)
			reduceLod();
	}
	
	inline void addSortedZ(long z, long c){
	
		this->addSorted(z, c);
		
		while (elems.size() > maxSparse)
			reduceLod();
	}
	
	inline void addPoint(const Point& p){
	
		this->add(getZ(p) );
		
		while (elems.size() > maxSparse)
			reduceLod();
	}
	
	inline void finalize(){
	
		processUnsorted();
		processSorted();
		
		while (elems.size() > maxSparse)
			reduceLod();
	}

};

#endif