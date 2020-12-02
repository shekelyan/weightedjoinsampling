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

#ifndef SHEKELYAN_DIGITHIST_H
#define SHEKELYAN_DIGITHIST_H

#include <approxdata/utils/utils.h>
#include <approxdata/data/data.h>
#include <approxdata/datasummaries/datasummaries.h>

/*

Class for the digit histograms composing DigitHist summaries.

*/

class DigitHistogram : public DataSummary{

typedef unique_ptr<DigitHistogram> DH;
typedef unique_ptr<SparseGridHist> Hist;
typedef unique_ptr<EquiWidth1d> Marg;

protected:

	long mul = 1; // (Digit)-multiplier of all counts
	const int DIMS; // Number of dimensions
	bool AVI = false; // Intra-bucket distribution according to attribute value independence assumption?
	vector<OneDimEstimator*> oneDimEsts; // Intra-bucket distribution estimators (pointers to 1D histograms)
	
	inline DigitHistogram(int dims, bool avi=false) : DIMS(dims), AVI(avi){
		
		mul = 1;
	}
	
	inline DigitHistogram(int dims, int s, Hist& h, bool avi) : DigitHistogram(dims, avi){
		
		gridhist = std::move(h);
		mul = 1L << s;
	}
	
public:
	
	
	vector<Marg> marginals; // 1D histograms for intra-bucket distribution
	Hist gridhist; // multidimensional grid histogram
	
	// Constructor only used internally (by other class in the library)
	inline DigitHistogram(const DataSet& data, Params& params) : DigitHistogram(data.getDims() ){
		
		summarize(data, params);
	}
	
	
	inline const string toString() const override{
		
		return "gridhist";
	}
	
	inline double getBytes() const override{
	
	
		double ret = gridhist->getBytes();
		
		for (auto it = marginals.begin(); it != marginals.end(); it++)
			ret += (*it)->getBytes();
	
		return ret;
	}
	
	// Constructor receiving a dataset and string with parameters
	
	inline DigitHistogram(const DataSet& data, string paramStr="-kb 100") : DigitHistogram(data.getDims() ){
	
		Params params(paramStr);
		summarize(data, params);
	}
	
	// All counts of the grid histogram are multiples of {mul} (they are stored in gridhist divided by mul)
	inline long getMultiplier() const{
	
		return mul;
	}
	
	~DigitHistogram() override {};
	
	inline void checkConsistency(){
		
		assert (AVI);
		assert (marginals.size() > 0);
		
		for (int i = 0; i < DIMS; i++){
		
			marginals[i]->makeCumulative(false);
			
			const int margLod = gridhist->z.dimMaskBits[i];
			const int srcLod = UtilMath::getPowerOfTwo(marginals[i]->counts.size() );
			
			const long buckets = gridhist->getBucketNum();
			const long margBuckets = 1L << margLod;
			const long srcBuckets = 1L << srcLod;
			
			vector<long>& src = marginals[i]->counts;
			
			if (src.size() != srcBuckets)
				cout << "src.size" << src.size() << " != " << srcBuckets << " srcBuckets" << endl;
			
			
			assert (src.size() == srcBuckets);
			
			
			for (long jsrc = 0; jsrc < srcBuckets; jsrc++)
				assert(src[jsrc] >= 0);
			
			if (margLod >= srcLod){
			
				const int lodDiff = margLod-srcLod;
				
				for (long j = 0; j < buckets; j++){
			
					const long jsrc = gridhist->z.getDimCoord(gridhist->elems[j], i) >> lodDiff;
					const long c = gridhist->counts[j]*mul;
					
					assert (jsrc >= 0);
					assert (jsrc < src.size() );
					
					assert (src[jsrc] >= 0);
					
					assert (c >= 0);
					assert (c <= src[jsrc]);
					
					src[jsrc] -= c;
				}
				
			} else {
				
				vector<long> tmp(margBuckets, 0);
				
				for (long j = 0; j < buckets; j++){
			
					const long x = gridhist->z.getDimCoord(gridhist->elems[j], i);
					const long y = gridhist->counts[j] * mul;
					
					assert (x >= 0);
					assert (x < tmp.size() );
					
					assert (y > 0);
					
					tmp[x] += y;
				}
				
				const int lodDiff = srcLod-margLod;
				
				vector<double> v( 1L << lodDiff );
			
				for (long jmarg = 0; jmarg < margBuckets; jmarg++){
					
					long b = tmp[jmarg];
					
					if (b == 0)
						continue;
					
					double vsum = 0;
					
					const long jsrc = jmarg << lodDiff;
					
					for (long jv = 0; jv < v.size(); jv++){
						v[jv] = src[jsrc+jv];
						vsum += v[jv];
					}
					
					if (vsum == 0)
						continue;
					
					//if (b > vsum)
						//b = vsum;
					
					const double sumDiv = b*1.0/vsum;
					
					assert(sumDiv > 0);
					assert(sumDiv <= 1);
					
					for (long jv = 0; jv < v.size(); jv++)
						v[jv] *= sumDiv;
					
					UtilMath::optimalRounding(v);
					
					for (long jv = 0; jv < v.size(); jv++){
					
						const long c = v[jv];
						
						if (c == 0)
							continue;
						
						if (c > src[jsrc+jv]){
							cout << "c " << c << " > " << src[jsrc+jv] << endl;
							cout << endl;
						}
						
						assert(c <= src[jsrc+jv]);
						
						src[jsrc+jv] -= c;
					}
				}
			}
			
			for (long jsrc = 0; jsrc < srcBuckets; jsrc++)
				assert(src[jsrc] >= 0);
			
			assert (UtilMath::sum<long>(src) == 0 );
		}
		
		cout << "check passed!" << endl;
	}
	
	inline void take(vector<Marg>& initMarginals){
		
		assert (AVI);
		assert (initMarginals.size() > 0);
		
		for (int i = 0; i < DIMS; i++){
		
			initMarginals[i]->makeCumulative(false);
			
			const int margLod = gridhist->z.dimMaskBits[i];
			
			assert (UtilMath::isPowerOfTwo( initMarginals[i]->counts.size() ) );
			
			const int srcLod = UtilMath::getPowerOfTwo(initMarginals[i]->counts.size() );
			const int dstLod = srcLod;
			
			const long buckets = gridhist->getBucketNum();
			const long margBuckets = 1L << margLod;
			const long srcBuckets = 1L << srcLod;
			const long dstBuckets = 1L << dstLod;
			
			Marg m(new EquiWidth1d(dstBuckets) );
			
			vector<long>& src = initMarginals[i]->counts;
			vector<long>& dst = m->counts;
			
			if (src.size() != srcBuckets)
				cout << "src.size" << src.size() << " != " << srcBuckets << " srcBuckets" << endl;
			
			
			assert (src.size() == srcBuckets);
			assert (dst.size() == dstBuckets);
			
			for (long jsrc = 0; jsrc < srcBuckets; jsrc++)
				assert(src[jsrc] >= 0);
			
			if (margLod >= srcLod){
			
				const int lodDiff = margLod-srcLod;
				
				for (long j = 0; j < buckets; j++){
			
					const long jsrc = gridhist->z.getDimCoord(gridhist->elems[j], i) >> lodDiff;
					const long c = gridhist->counts[j]*mul;
					
					assert (jsrc >= 0);
					assert (jsrc < dst.size() );
					assert (jsrc < src.size() );
					
					assert (src[jsrc] >= 0);
					
					assert (c >= 0);
					assert (c <= src[jsrc]);
					
					dst[jsrc] += c;
					src[jsrc] -= c;
				}
				
			} else {
				
				vector<long> tmp(margBuckets, 0);
				
				for (long j = 0; j < buckets; j++){
			
					const long x = gridhist->z.getDimCoord(gridhist->elems[j], i);
					const long y = gridhist->counts[j] * mul;
					
					assert (x >= 0);
					assert (x < tmp.size() );
					
					assert (y > 0);
					
					tmp[x] += y;
				}
				
				const int lodDiff = srcLod-margLod;
				
				vector<double> v( 1L << lodDiff );
			
				for (long jmarg = 0; jmarg < margBuckets; jmarg++){
					
					long b = tmp[jmarg];
					
					if (b == 0)
						continue;
					
					double vsum = 0;
					
					const long jsrc = jmarg << lodDiff;
					
					for (long jv = 0; jv < v.size(); jv++){
						v[jv] = src[jsrc+jv];
						vsum += v[jv];
					}
					
					if (vsum == 0)
						continue;
					
					if (b > vsum)
						b = vsum;
					
					const double sumDiv = b*1.0/vsum;
					
					assert(sumDiv > 0);
					assert(sumDiv <= 1);
					
					for (long jv = 0; jv < v.size(); jv++)
						v[jv] *= sumDiv;
					
					UtilMath::optimalRounding(v);
					
					for (long jv = 0; jv < v.size(); jv++){
					
						const long c = v[jv];
						
						if (c == 0)
							continue;
						
						if (c > src[jsrc+jv]){
							cout << "c " << c << " > " << src[jsrc+jv] << endl;
							cout << endl;
						}
						
						assert(c <= src[jsrc+jv]);
						
						src[jsrc+jv] -= c;
						dst[jsrc+jv] += c;
					}
				}
			}
			
			for (long jsrc = 0; jsrc < srcBuckets; jsrc++)
				assert(src[jsrc] >= 0);
			
			m->makeCumulative(true);
			
			marginals.push_back(std::move(m) );
			oneDimEsts.push_back( marginals[i].get() );	
		}
	}
	
	inline void add(const Point& p){
	
		gridhist->addPoint(p);
			
		if (AVI){
			for (int i = 0; i < DIMS; i++)
				marginals[i]->add(p[i] );
		}
	}
	
	bool operator<(const DigitHistogram& other) const{
	
		if (mul == other.mul)
			return gridhist->total < other.gridhist->total;
		
		return mul < other.mul;
    }
	
	inline void reduceLod(){
	
		gridhist->reduceLod();
	}
	
	inline long getLargestCount() const{
	
		long ret = 0;
		
		for (long j = 0; j < gridhist->counts.size(); j++)
			if (gridhist->counts[j] > ret)
				ret = gridhist->counts[j];
		
		return ret;
	}
	
	
	inline int getLod() const{
	
		return gridhist->z.getLod();
	}
	
	inline double getUError(int steps=30) const{
	
		return gridhist->getUError()*mul;
	}
	
	/*
	inline const string toString() const{
	
		stringstream s;
		
		s << gridhist->toString();
		
		s << " " << mul;
		
		return s.str();
	}*/
	
	
	
	
	inline DH createDigitHistogram(int s, int w, long maxBuckets, int lod=62){
		
		Hist retgrid( new SparseGridHist(gridhist->DIMS, maxBuckets, UtilMath::minVal<int>(lod, gridhist->getLod())) );
		
		assert (gridhist->getBucketNum() >= 0);
		
		for (long j = 0; j < gridhist->getBucketNum(); j++){
			
			const long c = (gridhist->counts[j] >> s) & UtilBits::mask(w);
			
			if (c == 0)
				continue;
			
			const int lodDiff = (gridhist->getLod()-retgrid->getLod() );
			
			assert (lodDiff >= 0);
			
			const long v = gridhist->elems[j] >> lodDiff;
			retgrid->addSortedZ(v, c);
		}
		
		retgrid->finalize();
		
		DH ret(new DigitHistogram(DIMS, s, retgrid, AVI) );
		
		return ret;
	}
	
	inline long count(){
	
		return UtilMath::sum<long>(gridhist->counts);
	}
	
	inline void summarize(const DataSet& data, Params& params){
		
		const int initLod = params.get("lod") ? params.get("lod") : 62;
		
		const double maxSparse = params.get("kb")/( 8.0/(1024) );
		const double maxMarg = params.get("margkb")/( 8.0/(1024) );
		
		AVI = !params.get("nomarg") && maxMarg > 0;
		mul = 1;
		
		//cout << "create summary DigitHistogram " << params << endl;
		
		long marginalBuckets = UtilMath::powerOfTwoUB( maxMarg/DIMS );
		long gridBuckets =  maxSparse;
		
		/*
		if (params.get("check")){
		
			initLod = params.get("check")*DIMS;
			
			AVI = true;
			
			marginalBuckets = 1L << (10+ (int) params.get("check"));
			gridBuckets = 10000000;
		}*/
		
		assert (gridBuckets > 0);
		
		if (AVI)
			assert (marginalBuckets > 0);
		
		gridhist.reset(new SparseGridHist(DIMS, gridBuckets, initLod) );
		
		if (AVI){
			
			for (int i = 0; i < DIMS; i++){
			
				Marg m(new EquiWidth1d(marginalBuckets));
			
				marginals.push_back( std::move(m) );
				oneDimEsts.push_back( marginals[i].get() );
			}
		}
		
		for (auto it = data.begin(); it != data.end(); it++)
			add( (*it) );
		
		gridhist->finalize();
		
		cout << gridhist->toString() << endl;
		
		if (AVI)
		for (int i = 0; i < DIMS; i++){
		
			marginals[i]->makeCumulative(true);
			//cout << "i " << i << " " << marginals[i]->toString() << endl;
		}
	}
	
	inline unique_ptr<RangeQueries> getBuckets() const{
	
		unique_ptr<RangeQueries> ret(new RangeQueries() );
		
		for (auto it = gridhist->elems.begin(); it != gridhist->elems.end(); it++){
		
			Box b = gridhist->getBox( (*it) );
			ret->add( b );
		}
		
		return ret;
	}
	
	inline unique_ptr<RangeQueries> getBuckets(const Box& box) const{
	
		unique_ptr<RangeQueries> ret(new RangeQueries() );
		
		const long zmin = gridhist->getZ(box.min);
		const long zmax = gridhist->getZ(box.max);
		
		long k = 0;
		for (auto it = gridhist->elems.begin(); it != gridhist->elems.end(); it++, k++){
		
			if ( !gridhist->z.contains(zmin, zmax, (*it)) )
				continue;
			
			Box b = gridhist->getBox( (*it) );
			
			b.count = gridhist->counts[k];
			ret->add(b);
		}
		
		return ret;
	}
	
	inline double marginalBoxCountUpperBound(int i, double qmin, double qmax) const {
	
		if (AVI)
			return marginals[i]->intervalCount(qmin, qmax, QueryMode::UB);
		else
			return mul*gridhist->getTotalNum();
	}
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode = QueryMode::EST) const override{
		
		if (AVI){
		
			assert (oneDimEsts.size() == DIMS);
			return mul*gridhist->boxCount(qmin, qmax, mode, oneDimEsts);
			
		} else{
			
			return mul*gridhist->boxCount(qmin, qmax, mode);
		}
	}
	
	#if defined(SHEKELYAN_SELIMAGE)  
	
	inline void visualize(DataImage& img) const{
	
		array<float,3> black = {0,0,0};
		
		array<float,3> green = {0,100,0};
		array<float,3> violet = {255,0,255};
		array<float,3> red = {255,0,0};
		array<float,3> blue = {0,0,255};
		array<float,3> gray = {100,100,100};
		
		double sum = 1.0;
		
		shared_ptr<RangeQueries> buckets = getBuckets(img.getDataBox());
			
		for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
			sum += it->count;
		
		Box b;
		
		for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++){
			
			b = (*it);
			b.count /= sum;
				
			img.fillRect(b, black);	
		}
	}
	
	#endif
	
};


// DigitHist data summary as described in http://www.vldb.org/pvldb/vol10/p1514-shekelyan.pdf (except two-level representation)

// brief overview of terminology:
// lod = level of detail (binary logarithm of number of cells/buckets)

class DigitHistSummary : public DataSummary{

typedef unique_ptr<DigitHistogram> DH;
typedef unique_ptr<DigitHistogram> Hist;

public:
	
	~DigitHistSummary() override {};
	
	vector<DH> digitHistograms;
	
	int DIMS;
	
	Params params;
	
	inline double getBytes() const override{
	
		double ret = 0;
		
		for (auto it = digitHistograms.begin(); it != digitHistograms.end(); it++)
			ret += (*it)->getBytes();
			
		return ret;
	}
	
	inline DigitHistSummary(const DataSet& data, string paramStr) : params(paramStr), DIMS(data.getDims() ){
		
		summarize(data, params);
	}
	
	inline DigitHistSummary(shared_ptr<DataSet> data, string paramStr) : params(paramStr), DIMS(data->getDims() ){
		
		summarize( (*data), params);
	}
	
	inline DigitHistSummary(unique_ptr<DataSet> data, string paramStr) : params(paramStr), DIMS(data->getDims() ){
		
		summarize( (*data), params);
	}
	
	inline const string toString() const override{
		
		return "digithist "+params.toString();
	}
	
	inline void summarize(const DataSet& data, Params& params){
		
		Params paramsInit;
		
		if (params.get("tkb") == 0)
			params.set("tkb", 100000);
		
		paramsInit.set("kb", params.get("tkb") );
		
		const bool AVI = !params.get("nomarg");
		
		if (!AVI)
			paramsInit.set("nomarg", 1);
			
		observer["avi"] = AVI ? 1 : 0;
		
		const long bytes = params.get("kb")*1024;
		
		const long margBytes = AVI ? bytes*0.25 : 0;
		
		paramsInit.set("margkb", margBytes/1024.0);
		
		cout << "create summary DigitHist " << params << endl;
		
		const int DIMS = data.getDims();
		
		
		// 5.1.1 Initial Histograms
		
		// 1. Create an initial d-dimensional histogram H and, for each
		// dimension i, a one-dimensional marginal histogram.
		
		DH init( new DigitHistogram(data, paramsInit) );
		// both H and M1, ... Md are stored in "init"
		
		observer["initlod"] = init->getLod();
		
		cout << "init lod " << init->getLod() << endl;
		
		const int K = params.get("k") > 0 ? params.get("k") : 4;
		
		const long digitBytes = bytes-margBytes;
		
		observer["digitbytes"] = digitBytes;
		
		const double minBytesPerElem = 2;
		const double maxBytesPerElem = 16;
		
		const long minElems = floor(digitBytes/maxBytesPerElem);
		const long maxElems = ceil(digitBytes/minBytesPerElem);
		
		cout << "maxElems " << maxElems << endl;
		
		// 5.1.2 Digit Histogram Compression
		
		// 2. Decompose and compress initial histogram into K digit histograms of maximum S bytes
		
		vector<int> lods = {init->getLod() };
		
		//cout << "A" << endl;
		
		{
		
			if (init->gridhist->getBytes() <= digitBytes){
		
				digitHistograms.push_back( std::move(init) );
				return;
			}
			
			//cout << "B" << endl;
		
			double maxProfit = std::numeric_limits<double>::lowest();
			DH bestInit = std::move( init->createDigitHistogram(0, 62, maxElems, init->getLod() ) );
			
			Timer t("digithist");
		
			t.start(10, 1, init->getLod()*init->getLod()*K );
		
			// Algorithm 1: DigitHistCompress, Line 2
		
			//cout << "C" << endl;
		
			while (init->gridhist->getBytes() > digitBytes ){
				
				const long ub = UtilMath::powerOfTwoUB( init->getLargestCount() );
				const int countLod = UtilMath::getPowerOfTwo(ub);
				const int step = ceil(countLod*1.0/K);
				
				//cout << "D" << endl;
				
				if (step <= 0){
				
					cout << "step <= 0" << endl;
					init->reduceLod();
					continue;
				}
				
				// brute-force solver for multiple-choice knapsack problems
				MultipleChoiceKnapsack knap(K, init->getLod()+1);
			
				// Algorithm 1, Line 3: decompose H into (D[K1], . . . , D[0])
			
				for (int k = 0; k < K; k++){
					
					//cout << "E" << endl;
					
					DH dh = std::move( init->createDigitHistogram(k*step, step, maxElems) );
					
					//cout << "F" << endl;
					
					while (dh->gridhist->getBytes() > digitBytes && dh->gridhist->getLod() > 1)
						dh->gridhist->reduceLod();
						
					//cout << "G" << endl;
				
					while (dh->gridhist->getLod() > 1){
				
						t.tick();
						
						//cout << "G1" << endl;
				
						const int lod = dh->gridhist->getLod();
				
						//cout << "G2" << endl;
				
						const double profit = -dh->getUError();
						
						//cout << "G3" << endl;
						
						const double weight = dh->gridhist->getBytes();
				
						//cout << "H" << endl;
				
						knap.add(k, lod, profit, weight );
						dh->gridhist->reduceLod();
						
						//cout << "I" << endl;
				
					}
				}
			
				// Algorithm 1, Line 4 and 5
				
				//cout << "J" << endl;
				
				vector<int> candLods = knap.maximizeProfit(digitBytes);
			
				//cout << "K" << endl;
			
				const double profit = knap.getProfit(candLods);
				const double weight = knap.getWeight(candLods);
				
				assert (weight <= digitBytes);
			
				// Algorithm 1, Line 6
			
				if (profit >= maxProfit){
			
					// Algorithm 1, Line 7
					
					lods = candLods;
					maxProfit = profit;
					bestInit = std::move( init->createDigitHistogram(0, 62, maxElems, init->getLod() ) );	
				}
				
				//cout << "L" << endl;
				
				// Algorithm 1, Line 8
				
				init->reduceLod();
				
				//cout << "M" << endl;
			}
			
			//cout << "N" << endl;
		
			t.end();
		
			
			//long total = 0;
		
			for (int k = 0; k < lods.size(); k++){
		
				const long ub = UtilMath::powerOfTwoUB( bestInit->getLargestCount() );
				const int countLod = UtilMath::getPowerOfTwo(ub);
				const int step = ceil(countLod*1.0/K);
			
				DH dh = std::move( bestInit->createDigitHistogram(k*step,step, maxElems, lods.at(k)) );
			
				//for (long j = 0; j < dh->gridhist->counts.size(); j++)
				//	total += dh->gridhist->counts[j] * dh->getMultiplier();
				
				digitHistograms.push_back( std::move(dh) );
			}
		}
		
		// 5.1.3 Marginal Histograms
		
		// 3. Distribute each initial marginal histogram over the marginal histograms
		
		if (AVI){
		
			Sorter s;
		
			for (int k = 0; k < lods.size(); k++)
				s.add(k, -lods.at(k) );
			
			for (int kk = 0; kk < lods.size(); kk++)
				digitHistograms[s.getSortedIndex(kk)]->take( init->marginals );
		}
		
		observer["bytes"] = getBytes();
		observer["k"] = K;
		
		//cout << "total " << total << " data size " << data.getSize() << endl;
		
		/*
		assert (total <= data.getSize() );
		assert (total >= data.getSize() );
		
		if (true){
			
			const long siz = bestInit->gridhist->counts.size();
			
			map<long, long> c;
			
			for (int k = 0; k < K; k++){
			
				const long ub = UtilMath::powerOfTwoUB( bestInit->getLargestCount() );
				const int countLod = UtilMath::getPowerOfTwo(ub);
				const int step = ceil(countLod*1.0/K);
				
				DH dh = std::move( bestInit->createDigitHistogram(k*step, step, siz+1 ) );
				
				assert (dh->getLod() == bestInit->getLod() );
				
				for (long j = 0; j < dh->gridhist->counts.size(); j++){
			
					assert (dh->gridhist->counts[j] > 0);
 					
					c[dh->gridhist->elems[j]] += dh->gridhist->counts[j] * dh->getMultiplier();
				}
			}
			
			assert (c.size() <= bestInit->gridhist->counts.size() );
			assert (c.size() >= bestInit->gridhist->counts.size() );
			
			for (long j = 0; j < bestInit->gridhist->counts.size(); j++){
				
				if ( c[bestInit->gridhist->elems[j]] != bestInit->gridhist->counts[j])
					cout << c[bestInit->gridhist->elems[j]] << " != " << bestInit->gridhist->counts[j] << endl;
				
				assert (c[bestInit->gridhist->elems[j]] == bestInit->gridhist->counts[j]);
			}
		}
		
		
		
		if (AVI)
		for (int i = 0; i < DIMS; i++)
			cout << "i " << i << " strays " << UtilMath::sum(init->marginals[i]->counts) << endl;
		
		//std::sort(digitHistograms.begin(), digitHistograms.end() );
		
		cout << "finished" << endl;*/
	}
	
	// approximated count in box is equal to sum of component counts
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override {
		
		double ret = 0;
		
		for (int k = 0; k < digitHistograms.size(); k++)
			ret += digitHistograms[k]->boxCount(qmin, qmax, mode);
		
		for (int i = 0; i < DIMS; i++){
		
			double r = 0;
		
			for (int k = 0; k < digitHistograms.size(); k++)
				r += digitHistograms[k]->marginalBoxCountUpperBound(i, qmin[i], qmax[i]);
		
			if (r < ret)
				return r;
		}
		
		if (mode == QueryMode::EST){
		
			assert (ret <= boxCount(qmin, qmax, QueryMode::UB) );
			assert (ret >= boxCount(qmin, qmax, QueryMode::LB) );
		}
		
		return ret;
	}
	

	#if defined(SHEKELYAN_SELIMAGE) 
	
	inline void visualize1(DataImage& img) const{
	
		array<float,3> black = {0,0,0};
		array<float,3> green = {0,100,0};
		array<float,3> violet = {255,0,255};
		array<float,3> red = {255,0,0};
		array<float,3> blue = {0,0,255};
		array<float,3> gray = {100,100,100};
		
		const int K = digitHistograms.size();
		
		
		double sum = 1.0;
		
		
		for (int k = 0; k < K; k++){
			
			unique_ptr<RangeQueries> buckets = std::move( digitHistograms[k]->getBuckets(img.getDataBox()) );
			
			for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
				sum += it->count;
		}
		
		for (int k = 0; k < K; k++){
			
			unique_ptr<RangeQueries> buckets = std::move( digitHistograms[k]->getBuckets(img.getDataBox()) );
			
			Box b;
			
			for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++){
			
				b = (*it);
				b.count /= sum;
				
				img.fillRect(b, black);	
			}
		}
	}
	
	inline void visualize(DataImage& img) const{
	
		array<float,3> black = {0,0,0};
		array<float,3> green = {0,100,0};
		array<float,3> violet = {255,0,255};
		array<float,3> red = {255,0,0};
		array<float,3> blue = {0,0,255};
		array<float,3> gray = {100,100,100};
		
		const int K = digitHistograms.size();
		
		Sorter s;
		
		for (int k = 0; k < K; k++)
			s.add(k, -digitHistograms[k]->getMultiplier() );
		
		for (int kk = 0; kk < K; kk++){
		
			for (int k = 0; k < K; k++){
			
				if (k != s.getSortedIndex(kk))
					continue;
			
				shared_ptr<RangeQueries> buckets = digitHistograms[k]->getBuckets(img.getDataBox());
				
				if (k == 3){
				
					for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
						img.drawRect( (*it), green, 1.0);
				}
				
				if (k == 2){
			
					for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
						img.drawRect( (*it), violet, 1.0);
				}
			
				if (k == 1){
			
					for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
						img.drawRect( (*it), blue, 0.2);
				}
			
				if (k == 0){
			
					for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
						img.drawRect( (*it), red, 1.0);
				}
			}
		}
	}
	
	#endif
	
};
#endif