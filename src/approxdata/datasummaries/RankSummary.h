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
#ifndef SHEKELYAN_RANKSUMMARY_H
#define SHEKELYAN_RANKSUMMARY_H

#include <approxdata/external/gk.h>
#include <approxdata/external/gk2.h>
#include <approxdata/datasummaries/datasummaries.h>


/*
	[-inf] < x1 <= [b1] < x2 <= [b2]
	
	gk gives me for any z both count(z <= x) and x
	
	x	0	0.1		0.2	0.3 0.4	0.5	1
	r	0	1		2	3	4	5	6
	
	for any x, tell me a < x <= b s.t. count(z <= a) <= count(z <= x) <= count(z <= b)
	
	
	
*/


template<typename E = Point, typename F=PointComparator>
class PointEquiwidth{

public:
	
	const int dim;
	
	const E dmin;
	const E dmax;
	
	E pmin;
	E pmax;
	
	const F comp;

	const double vmin;
	const double vmax;
	
	const long buckets;
	
	const double vdiv;

	const long maxBucketCount;
	
	const double ldiv;
	
	vector<bool> counts0;
	vector<uint8_t> counts1;
	vector<uint16_t> counts2;
	vector<uint32_t> counts4;
	vector<uint64_t> counts8;
	
	int mode;
	
	
	long count = 0;
	
	int state = -1;
	
	inline long getBucketNum() const {
	
		return buckets;
	}
	
	inline void setState(int newstate){
	
		assert (state+1 == newstate);
		
		state = newstate;
	
	}
	
	inline E& getBucketInfimum(long q) const {
	
		assert (false);
	}
	
	/*
	inline E& getBucketMax(long q) const {
	
		assert (false);
	}*/
	
	int getMode(long s) const {
	
		if ( ((bool) s) == s)
			return 0;
		
		if ( ((uint8_t) s) == s)
			return 1;
			
		if ( ((uint16_t) s) == s)
			return 2;
			
		if ( ((uint32_t) s) == s)
			return 3;
			
		if ( ((uint64_t) s) == s)
			return 4;
			
		assert (false);
		
	}
	
	inline void incrementCount(long bucket){
	
		switch(mode){
		
			case 0:
				assert(counts0[bucket] == 0);
				counts0[bucket] = 1;
				break;
			
			case 1:
				++counts1[bucket];
				break;
		
			case 2:
				++counts2[bucket];
				break;
			
			case 3:
				++counts4[bucket];
				break;
			
			case 4:
				++counts8[bucket];
				break;
				
			default:
				assert(false);
		}
		
	}
	
	inline void addCount(long bucket, long c){
	
		if (mode == 0){
		
			assert ( counts0[bucket]+c <= 1);
		
			if (c == 1)
				counts0[bucket] = 1;
			
			return;
		}
	
		if (mode == 1)
			counts1[bucket] += (uint8_t) c;
		
		else if (mode == 2)
			counts2[bucket] += (uint16_t) c;
			
		else if (mode == 3)
			counts4[bucket] += (uint32_t) c;
			
		else if (mode == 4)
			counts8[bucket] += c;
		
		else	
			assert(false);
	}
	
	inline long getBucketCount(long bucket) const{
	
		if (bucket < 0)
			return 0;
	
		assert (bucket < getBucketNum());
	
		switch(mode){
		
			case 0:
				return counts0.at(bucket);
			
			case 1:
				return counts1.at(bucket);
		
			case 2:
				return counts2.at(bucket);
			
			case 3:
				return counts4.at(bucket);
			
			case 4:
				return counts8.at(bucket);
		}
		
		assert(false);
	}
	
	bool cumulativeCounts = false;
	long cacheBucket = -1;
	long cacheBucketCount = 0;
	
	inline long getCumulativeBucketCount(long bucket){
	
		assert (state == 2);
	
		if (!cumulativeCounts){
		
			if (bucket < 0)
				return 0;
		
			if (bucket < cacheBucket){
			
				if (bucket == cacheBucket-1)
					return cacheBucketCount-getBucketCount(cacheBucket);
				else if (bucket == cacheBucket-2)
					return cacheBucketCount-getBucketCount(cacheBucket)-getBucketCount(cacheBucket-1);
				else if (bucket == cacheBucket-3)
					return cacheBucketCount-getBucketCount(cacheBucket)-getBucketCount(cacheBucket-1)-getBucketCount(cacheBucket-2);
					
				cacheBucket = -1;
				cacheBucketCount = 0;
			}
			
			while (bucket > cacheBucket){
			
				++cacheBucket;
				cacheBucketCount += getBucketCount(cacheBucket);
			}
			
			assert (bucket == cacheBucket);
				
			return cacheBucketCount;
		}
	
		return getBucketCount(bucket);
	}
	
	inline long getBucketCount(long q1, long q2) const{
	
		assert (cumulativeCounts);
	
		return getBucketCount(q2)-getBucketCount(q1-1);
	}
	
	inline PointEquiwidth(int dim, long bs, const E& dmin, const E& dmax, const F& comp, long maxBucketCount = 1L << 60) : 
		dim(dim), dmin(dmin), dmax(dmax),
		pmin(dmin), pmax(dmax), 
		comp(comp),
		vmin(dmin[dim]), vmax(dmax[dim]), 
		buckets(vmin == vmax ? 4 : bs), 
		vdiv(1.0/(vmax-vmin)), 
		maxBucketCount(maxBucketCount), 
		ldiv(1.0/(buckets-2))
		{
		
		setState(0);
		
		assert (buckets >= 4);
		
		mode = getMode(maxBucketCount);
		
		//cout << "mode " << mode << " maxBucketCount " << maxBucketCount << " buckets " << buckets << endl;
	
		if (mode == 0){
		
			counts0.reserve(buckets);
			for (long j = 0; j < (buckets); ++j)
				counts0.push_back(0);
				
			assert (counts0.size() == getBucketNum() );
		
		} else if (mode == 1){	
		
			counts1.reserve(buckets);
			for (long j = 0; j < (buckets); ++j)
				counts1.push_back(0);
				
			assert (counts1.size() == getBucketNum() );
				
		} else if (mode == 2){	
		
			counts2.reserve(buckets);
			for (long j = 0; j < (buckets); ++j)
				counts2.push_back(0);
				
			assert (counts2.size() == getBucketNum() );
				
		} else if (mode == 3){	
		
			counts4.reserve(buckets);
			for (long j = 0; j < (buckets); ++j)
				counts4.push_back(0);
			
			assert (counts4.size() == getBucketNum() );
				
		} else if (mode == 4){	
		
			counts8.reserve(buckets);
			for (long j = 0; j < (buckets); ++j)
				counts8.push_back(0);
				
			assert (counts8.size() == getBucketNum() );
			
		} else {
		
			assert (false);
		}
		
		setState(1);
	}
	
	inline double getBytes() const{
	
		assert (state == 2);
	
		return (mode == 0 ? 0.125 : mode == 1 ? 1.0 : mode == 2 ? 2.0 : mode == 3 ? 4.0 : mode == 4 ? 8.0 : -1 )*(buckets)+dmin.size()*4*8;
	}
	
	inline void finalise(bool cumulative = true){
	
		if (!cumulative){
		
			cumulativeCounts = false;
			setState(2);
			return;
		}
		
		int oldmode = mode;
		
		long sum = 0;
	
		for (long j = 0; j < getBucketNum(); ++j)
			sum += (long) getBucketCount(j);
		
		int newmode = getMode(sum);
		
		if (newmode != oldmode){
		
			assert (newmode != 0);
			
			if (newmode == 1){
				
				counts1.reserve(getBucketNum() );
			
				for (long j = 0; j < getBucketNum(); ++j)
					counts1.push_back(getBucketCount(j));
			}
				
			if (newmode == 2){
			
				counts2.reserve(getBucketNum() );
			
				for (long j = 0; j < getBucketNum(); ++j)
					counts2.push_back(getBucketCount(j));
			}
		
			if (newmode == 3){
			
				counts4.reserve(getBucketNum() );
			
				for (long j = 0; j < getBucketNum(); ++j)
					counts4.push_back(getBucketCount(j));
			}
		
			if (newmode == 4){
			
				counts8.reserve(getBucketNum() );
			
				for (long j = 0; j < getBucketNum(); ++j)
					counts8.push_back(getBucketCount(j));
					
			}
			
			mode = newmode;
		
			if (oldmode == 0){
			
				for (long j = 0; j < counts0.size(); ++j)
					assert (counts0[j] == getBucketCount(j));
			
				counts0.resize(0);
			}
			
			if (oldmode == 1){
			
				for (long j = 0; j < counts1.size(); ++j)
					assert (counts1[j] == getBucketCount(j));
			
				counts1.resize(0);
			}
				
			if (oldmode == 2){
			
				for (long j = 0; j < counts2.size(); ++j)
					assert (counts2[j] == getBucketCount(j));
			
				counts2.resize(0);
			}
				
			if (oldmode == 3){
			
				for (long j = 0; j < counts4.size(); ++j)
					assert (counts4[j] == getBucketCount(j));
			
				counts4.resize(0);
				
			}
		}
	
		if (mode == 0){
		
			for (long j = 1; j < counts0.size(); ++j)
				counts0[j] = counts0[j] | counts0[j-1];
				
			assert (counts0[counts0.size()-1] == sum);
				
		} else if (mode == 1){
		
			for (long j = 1; j < counts1.size(); ++j)
				counts1[j] += counts1[j-1];
				
			assert (counts1.back() == sum);
				
		} else if (mode == 2){
			
			for (long j = 1; j < counts2.size(); ++j)
				counts2[j] += counts2[j-1];	
				
			assert (counts2.back() == sum);
				
		} else if (mode == 3){
		
			for (long j = 1; j < counts4.size(); ++j)
				counts4[j] += counts4[j-1];	
				
			assert (counts4.back() == sum);
				
		} else if (mode == 4){
		
			for (long j = 1; j < counts8.size(); ++j)
				counts8[j] += counts8[j-1];	
				
			assert (counts8.back() == sum);
		}
				
		for (long j = 1; j < getBucketNum(); ++j){
		
			if (! (getBucketCount(j) >= getBucketCount(j-1)))
				cout << getBucketCount(j) << " >= " << getBucketCount(j-1) << endl;
			assert ( getBucketCount(j) >= getBucketCount(j-1) );
		}
		
		cumulativeCounts = true;
		setState(2);
	}
	
	inline const long numOfEndpoints() const{
	
		return buckets-1;
	}
	
	inline const E getBucketMax(long j) const{
	
		if (j == 0){
					
			return dmin;
		
		} else if ( (j-1) < buckets ){
		
			E p = dmax;
			p[dim] = UtilHist::bucketMin(j-1, buckets-2);
			
			return p;
			
		} else {
		
			assert (false);
		}
	}
	
	
	inline void finalise(vector<E>& endpoints, vector<double>& ranks, long absEps = -1){
		
		E p = dmax;
		
		long res = 0;
		
		if (absEps == -1)
			absEps = maxBucketCount;
		
		assert (absEps >= maxBucketCount);
		
		
		for (int k = 0; k <= 1; ++k){
		
			long last = 0;
		
			long sum = 0;
		
			for (long j = 0; j < buckets-1; ++j){
		
				sum += getBucketCount(j);
		
				if (j < getBucketNum()-1 && ( (sum-last)+getBucketCount(j+1) <= absEps))
					continue;
					
				if (k == 0){
				
					++res;
					
				} else {
				
					endpoints.push_back( getBucketMax(j) );
					ranks.push_back(sum );
					
					assert (sum-last <= maxBucketCount);
			
				}
				
				last = sum;
				
			}
			
			if (k == 0){
			
				endpoints.clear();
				ranks.clear();
			
				endpoints.reserve(res);
				ranks.reserve(res+1);
			}
		}
		
		ranks.push_back(count);
	}
	
	
	inline double getRank(const E& p, QueryMode mode=QueryMode::EST) {
		
		
		if (count == 0)
			return count;
		
		const long q = getBucket(p);
		
		const double s1 = getCumulativeBucketCount(q-1);
		const double s2 = getCumulativeBucketCount(q);
		
		if (mode == QueryMode::LB)
			return s1;
		
		if (mode == QueryMode::UB)
			return s2;
		
		if (s1 == s2)
			return s2;
		
		assert (p[dim] >= UtilHist::bucketMin(q-1, buckets-2 ));
		
		
		double w = 0.5;
		
		if (q == 0){
		
			assert (!comp(dmin, p));
			w = comp.interpolation(pmin, dmin, p);
		
		} else if (q == getBucketNum()-1){
		
			assert (comp(dmax, p));
			w = comp.interpolation(dmax, pmax, p);
			
		} else {
		
			w = ((p[dim]-vmin)*vdiv)*(buckets-2)-q;
		}
		
		return w <= 0 ? std::nexttoward(s1, s2) : w >= 1 ? std::nexttoward(s2, s1) : s1 + w*(s2-s1);
	}
	
	inline long getBucket(const E& p) const{
	
		const double v = p[dim];
	
		if (v <= vmin && comp(p, dmin) )
			return 0;
		
		if (v >= vmax && comp(dmax, p))
			return getBucketNum()-1;
		
		return 1+UtilHist::discretize( (p[dim]-vmin)*vdiv, buckets-2);
	}
	
	
	inline bool add(const E& p){
	
		assert (state == 1);
	
		const long bucket = getBucket(p);
		
		if (getBucketCount(bucket) == maxBucketCount)
			return false;
			
		const double v = p[dim];
		
		if (v <= vmin && comp(p, pmin))
			pmin = p;
		
		if (v >= vmax && comp(pmax, p))
			pmax = p;
			
		++count;
		
		incrementCount(bucket);
		
		return true;
	}

};

enum RankSummaryMode{

	Q_EQUI, Q_SKETCH, Q_GK, Q_NAIVE, Q_MERGE
};

struct RankSummaryParameters{

	// USER INPUT
	double eps = 0;
	long datasize = -1;
	long datasizeLB = 0;
	long datasizeUB = 1L << 40;
	
	double equibudget = 1;
	double sketchbudget = 1;
	double gkbudget = 1;
	double naivebudget = 1;

	long sketchbase = 4;
	double sketchdelta = 0.05;
	
	bool noid = false;
	
	int merge = 0;
	int naive = 10;
	int gk = 5;
	int equi = 0;
	int sketch = 0;
	int index = 0;
	bool keepEqui = false;
	bool oneDim = false;
	
	inline void init(){
	
		if (datasize != -1){
		
			datasizeLB = datasize;
			datasizeUB = datasize;
		}
		
		if (datasizeUB > 0 && datasizeLB == 0)
			datasizeLB = 1;
			
	
		const double totalErr = datasizeLB == 0 ? 1 : floor(eps*datasizeLB);
		
		bufferSize = totalErr == 0 ? datasizeUB : naive * ceil(datasizeUB / totalErr);
		
		if (datasizeUB < 10000 || bufferSize > datasizeUB)
			bufferSize = datasizeUB;
		
		if ( bufferSize >= datasizeUB){
		
			merge = 0;
		
			naiveErr = totalErr == 0 ? 1 : totalErr;
			return;
		}
		
		
		if (merge > 0){
			
			// 
				
			mergeErr2 = floor(totalErr * (merge/(1.0+merge) ) );
			mergeErr1 = totalErr-mergeErr2;
			
			if (mergeErr2 * mergeErr1 == 0){
			
				merge = 0;
				mergeErr1 = 0;
				mergeErr2 = 0;
				
				naiveErr = totalErr;
				return;
			
			} else {
			
				if (equi > 0){
				
					equiErr = ceil( mergeErr1 * 0.5);
			
					equibuckets = equi * ceil(datasizeUB/equiErr);
				
					if (equibuckets < 100)
						equibuckets = 100;
				
					if (equibuckets >= 0.1*datasizeUB){
						equiErr = 0;
						equi = 0;
						equibuckets = 0;
					}
					
				}
			
				const double gkErr = mergeErr1-equiErr;
			
				if (gk > 0){
			
					gkErr1 = floor( (mergeErr1 - equiErr)*(gk/(1.0+gk)) );
				
					gkErr2 = gkErr-gkErr1;
				
					if ( (1.0*datasizeUB/gkErr1) >= 0.1*datasizeUB)
						gkErr1 = 0;
				}
			
				if (gkErr1 * gkErr2 == 0){
			
					gkErr1 = 0;
					gkErr2 = 0;
				
					naiveErr = gkErr;
				}
				
				return;
			}
			
			
		}
		
		
		
		
		
		
		for (int i = 0; i < 1000; ++i){
		
			double err = 0;
		
			
		
			double budgetsum = equibudget+sketchbudget+gkbudget+naivebudget;
		
		
			
			
			if (equibudget * equi > 0){
		
				const double equieps = eps * (equibudget / budgetsum);
		
			
				equiErr = floor( equieps * datasizeLB);
			
			
				equibuckets = equi * ceil(datasizeUB/equiErr);
			
				if (equibuckets < 10)
					equibuckets = 10;
			
				/*
				if (equiErr < 1000){
					if (equibuckets < 100000)
						equibuckets = 100000;
				}*/
			
				err += equiErr;
			}
		
			if (equiErr == 0){
		
				equi = 0;
				equibuckets = 0;
				
				if (equibudget != 0){
					
					equibudget = 0;
					continue;
				}
				budgetsum -= equibudget;
			}
		
			if (sketchbudget * sketch > 0){
		
				const double sketcheps = eps * (sketchbudget / budgetsum);
		
				sketchErr = floor(sketcheps * datasizeLB);
				
				err += sketchErr;
			}
			
			if (sketchErr == 0){
			
				if (sketchbudget != 0){
					
					sketchbudget = 0;
					continue;
				}
			}
		
			if (gkbudget * gk > 0){
		
				const double gkeps = eps * (gkbudget / budgetsum);
				const double gkErr = floor(gkeps * datasizeLB);
		
				gkErr1 = floor(gkErr / gk);
				gkErr2 = gkErr1 == 0 ? 0 : gkErr - gkErr1;
			
				err += gkErr1+gkErr2;
			}
			
			if (gkErr1 == 0){
			
				if (gkbudget != 0){
					
					gkbudget = 0;
					continue;
				}
			}
		
			naiveErr = totalErr-err;
			
		
			break;
		}
		
		print();
		assert (naiveErr > 0);
	}
	
	inline void print() const{
	
		cout << endl;
		cout << endl;
		cout << "eps " << (eps) << endl;
		cout << "datasize " << (datasize) << endl;
		cout << "datasizeLB " << (datasizeLB) << endl;
		cout << "datasizeUB " << (datasizeUB) << endl;
		cout << endl;
		cout << "bufferSize " << bufferSize << endl;
		cout << "totalErr " << (eps*datasizeLB) << endl;
		cout << endl;
		cout << "merge " << merge << endl;
		cout << "mergeErr1 " << mergeErr1 << endl;
		cout << "mergeErr2 " << mergeErr2 << endl;
		cout << endl;
		cout << "equiErr " << equiErr << endl;
		cout << "sketchErr " << sketchErr << endl;
		cout << "gkErr1 " << gkErr1 << endl;
		cout << "gkErr2 " << gkErr2 << endl;
		cout << "naiveErr " << naiveErr << endl;
		cout << endl;
		cout << endl;
	}
	
	// AUTOMATICALLY DETERMINED
	
	long equibuckets = 0;
	
	double sketchErr = 0;
	double equiErr = 0;
	double gkErr1 = 0;
	double gkErr2 = 0;
	double naiveErr = 0;
	
	double mergeErr1 = 0;
	double mergeErr2 = 0;
	
	double bufferSize = 0;
};


template<typename E = Point, typename F=PointComparator>
class QuantileSketch{

public:

	const int dim;
	const int dims;
	
	const E domainMin;
	const E domainMax;
	
	unique_ptr<DyadicCMSketch> sketches;
	
	const F comp;

	long count = 0;

	vector<long> temp1;
	vector<long> temp2;
	
	E temp;
	
	shared_ptr<RankSummaryParameters> params;
	
	long maxCard;
	
	int state = -1;
	
	inline QuantileSketch<E,F>(int dim, const E& domainMin, const E& domainMax, const F& comp) : temp1(1), temp2(1), dim(dim), domainMin(domainMin), domainMax(domainMax), dims(domainMin.size() ), comp(comp) {
		
		assert (comp(domainMin, domainMax));
		setState(0);
	}
	
	
	inline void setState(int newstate){
	
		assert (state+1 == newstate);
		state = newstate;
	}
	
	inline void init(shared_ptr<RankSummaryParameters> par ){
		
		params = par;
		
		const int dims0 = 1;
		const long logUniverse = 32;
		
		const double eps0 = params->sketchErr / params->datasizeUB;
		
		const double kbPerValue = 8 / (1024.0);
		
		const double kb = (params->sketch / eps0 ) * kbPerValue;
		
		const double delta0 = params->sketchdelta;
		const bool fastHashes = false;
		const int m = params->sketchbase;
	
		sketches.reset( new DyadicCMSketch(dims0, logUniverse, kb, eps0 / params->sketch, delta0, fastHashes, m) );
		
		maxCard = params->sketchErr;
		
		setState(1);
	}
	
	inline uint64_t getLong(const double d1, const long l1) const{
	
		assert (d1 >= 0);
		const float f2 = static_cast<float>(d1);
		const int* i2 = reinterpret_cast<const int*>(&f2);
		const uint64_t l2 = *i2;
		
		//return l2;
		
		const uint64_t ret = l2;
		
		//const uint64_t ret = (l2 << 32) | (l1 & ((1L << 32)-1L));
		
		assert (ret >= 0);
		
		return ret;
	}
	
	inline double getDouble(const long l1) const{
	
		assert (l1 >= 0);
		const int i2 = static_cast<int>(l1);
		const float* f2 = reinterpret_cast<const float*>(&i2);
		
		return static_cast<double>(*f2);
	}
	
	inline uint64_t getLong(const E& p){
	
		return getLong(p[0], p[p.size()-1]);
	}
	
	inline const vector<long>& assign(vector<long>& a, const long v) const{
	
		a[0] = v;
		
		return a;
	}
	
	inline const E& assign(E& p, const long v) const{
	
		p[0] = getDouble(v);
		p[p.size()-1] = 0;
	
		//p[0] = getDouble(v >> 32);
		//p[p.size()-1] = v & ((1L << 32)-1L);
		return p;
	}
	
	inline bool add(const E& p){
	
		assert (state == 1);
		
		
		//cout << "P ";
		
		/*
		for (int i = 0; i < p.size(); ++i)
			cout << p[i] << " ";
		cout << endl;*/
	
		assign(temp1, getLong(p));
		
		//cout << "temp1 " << temp1[0] << endl;
	
		if ( sketches->add(temp1 , maxCard)){
		
		
			++count;
			return true;
		}
		
		
		return false;
	}

	
	inline void finalise(vector<E>& endpoints, vector<double>& ranks, long absEps = -1){
	
	
		if (absEps == -1)
			absEps = maxCard;
	
		endpoints.reserve( 4.0 * params->datasizeUB / params->sketchErr);
		ranks.reserve( endpoints.capacity() );
		
		const long x1 = getLong(domainMin);
		const long x2 = getLong(domainMax);
		
		long x = x1;
		
		
		
		long sum = 0;
		
		temp1[0] = x1;
		
		for (int i = 0; i < temp.size(); ++i)
			temp[i] = 0;
		
		while ( (x <= x2) && (sum+absEps < count) ){
		
			const long x1_ = x;
			
			BinarySearch<long, long> s( x1_, x2, sum+absEps);
			
			while ( s.find( sketches->upperbound(temp1, assign(temp2, s.x)) ) ){
			
			
			}
			
			
			
			
			if (s.lb < x1){
			
				x = x1_;
				sum += absEps;
			} else {
				
				x = s.lb;
				sum = s.lb_y;
			
			}
			
			//cout << "x1 " << x1 << " x2 " << x2 << " absEps " << absEps << " count " << count << endl;
			//cout << "x " << x << " sum " << sum << " temp ";
			
			
			assert (x >= x1_);
			
			assign(temp, x);
			
			/*
			for (int i = 0; i < temp.size(); ++i)
				cout << temp[i] << " ";
			cout << endl;*/
			
			
			endpoints.push_back( temp );
			ranks.push_back(sum);
			
			++x;
		}
		
		//cout << "1/eps " << params->datasize / params->sketchErr << " s " << endpoints.size() << endl;
		
		ranks.push_back(count);
		
		endpoints.shrink_to_fit();
		ranks.shrink_to_fit();
	
		setState(2);
	
	}	
};

template<typename E, typename F=std::less<E>>
class BasicRankSummary {
	

public:
	unique_ptr< GK2<E,F> > gk;
	
	unique_ptr< PointEquiwidth<E,F> > equi;
	
	unique_ptr< QuantileSketch<E,F> > sketch;
	
	unique_ptr< PointEquiwidth<E,F> > index;
	
	
	shared_ptr<RankSummaryParameters> params;
	
	
	inline long numOfEndpoints() const{
	
		if (equi)
			return equi->numOfEndpoints();
	
		return endpoints.size();
	}
	
	inline const E& getEndpoint(long pos){
	
		if (equi){
			temp = equi->getBucketMax(pos);
			return temp;
		}
		
		return endpoints.at(pos);
	}
	
	E temp;
	
	vector<E> endpoints;
	
	vector<double> counts;
	vector<double> ranks;
	
	/*
	vector<double> ranksLB;
	vector<double> ranksUB;*/
	
	long count = 0;
	
	const F comp;
	
	const E domainMin;
	const E domainMax;
	
	E pmin;
	E pmax;
	const int dim;
	
	const int dims;
	
	double absRankError = 0;
	
	double gkeps = 0;

	const RankSummaryMode mode;
	
	inline BasicRankSummary<E,F>(RankSummaryMode mode, int dim, const E& domainMin, const E& domainMax, const F& comp) : dim(dim), domainMin(domainMin), domainMax(domainMax), dims(domainMin.size() ), mode(mode), comp(comp) {
		
		assert (comp(domainMin, domainMax));
		
		pmin = domainMax;
		pmax = domainMin;
		
		temp = domainMin;
	}
	
	inline bool merge(vector< BasicRankSummary<E,F>* >& summaries){
	
		class MergeIter{
		
		public:
		
			vector<E> points;
			
			const E dmax;
			
			E last;
			
			vector<long> positions;
			
			long lb;
			long ub;
			
			vector<long> lbs;
			vector<long> ubs;
			
			vector<long> sizes;
			
			vector< BasicRankSummary<E,F>* > s;
			
			int current;
			const F* ptrcomp;
			
			inline BasicRankSummary<E,F>& currentSummary(){
			
				return *s.at(current);
			}
			
			inline int getMin(){
			
				int ret = 0;
				
				const F& comp = *ptrcomp;
				
				for (int i = 1; i < points.size(); ++i)
					if (comp( points.at(i), points.at(ret)))
						ret = i;
				
				return ret;
			}
			
			inline long getTotalCount(){
			
				long ret = 0;
				
				for (int i = 0; i < s.size(); ++i)
					ret += s.at(i)->count;
					
				return ret;
			}
		
			inline MergeIter(vector< BasicRankSummary<E,F>* >& summaries, const E& domainMin, const E& domainMax, const F& comp) : last(domainMin), dmax(domainMax), ptrcomp(&comp){
			
			
				for (int i = 0; i < summaries.size(); ++i)
					if (summaries.at(i)->numOfEndpoints() > 0)
						s.push_back(summaries.at(i) );
			
				for (int i = 0; i < s.size(); ++i){
		
					sizes.push_back( s.at(i)->numOfEndPoints() );
					positions.push_back(0);
					points.push_back( s.at(i)->getEndpoint(0) );					
				}
				
				current = getMin();
				
				const E& p = getCurrentPoint();
				
				lb = 0;
				ub = 0;
				
				assert (ub >= lb);
				
				for (int i = 0; i < s.size(); ++i){
				
					lbs.push_back(0);
					ubs.push_back(0);
					
					long& lbi = lbs.at(i);
					long& ubi = ubs.at(i);
					
					BasicRankSummary<E,F>& sm = *s.at(i);
					
					lbi = 0;
					ubi = sm.mode == RankSummaryMode::Q_NAIVE ? 1 : ceil(sm.getCumulativeBucketCount(0));
					
					assert (sm.absRankError >= 0);
					assert (ubi >= lbi);
					
					lbi -= ceil(sm.absRankError);
					ubi += ceil(sm.absRankError);
					
					if (lbi <= 0)
						lbi = 0;
					
					if (ubi >= sm.count)
						ubi = sm.count;
						
					assert (ubi >= lbi);
					
					lb += lbi;
					ub += ubi;
					
					assert (ub >= lb);
				}
				
			}
			
			bool hasReachedEnd(){
			
				const F& comp = *ptrcomp;
			
				return !comp(points.at(current), dmax);
			}
			
			inline void advance(){
			
				if (positions.at(current) > sizes.at(current) )
					return;
					
				const int i = current;
				
				
				++positions.at(i);
				
				BasicRankSummary<E,F>& sm = *s.at(i);
				
				const long bucket = positions.at(i);
				
				if (bucket == sizes.at(i) ){
				
					points.at(i) = dmax;
					
					current = getMin();
					
					assert (ub >= lb);
					
					lb -= lbs.at(i);
					ub -= ubs.at(i);
					
					assert (ub >= lb);
				
					lbs.at(i) = sm.count;
					ubs.at(i) = sm.count;
					
					assert (ubs.at(i) >= lbs.at(i) );
					
					lb += lbs.at(i);
					ub += ubs.at(i);
					
					assert (ub >= lb);
					
				} else if ( bucket < sizes.at(i)){
					
					long& lbi = lbs.at(i);
					long& ubi = ubs.at(i);
					
					assert (ubi >= lbi);
					
					points.at(i) = sm.getEndpoint(bucket);
					
					current = getMin();
				
					const E& p = getCurrentPoint();
					
					assert (ub >= lb);
					
					lb -= lbi;
					ub -= ubi;
					
					if ( !(ub >= lb)){
					
						cout << "i " << i << " lb " << lb << " ub " << ub << " lbi " << lbi << " ubi " << ubi << endl;
					}
					
					assert (ub >= lb);
					
					lbi = sm.mode == RankSummaryMode::Q_NAIVE ? bucket : floor(sm.getCumulativeBucketCount(bucket-1));
					ubi = sm.mode == RankSummaryMode::Q_NAIVE ? bucket+1 : ceil(sm.getCumulativeBucketCount(bucket));
					
					assert (sm.absRankError >= 0);
					assert (ubi >= lbi);
					
					lbi -= sm.absRankError;
					ubi += sm.absRankError;
					
					if (lbi <= 0)
						lbi = 0;
					
					if (ubi >= sm.count)
						ubi = sm.count;
						
					assert (ubi >= lbi);
					
					/*if (s.at(old)->ranks.size() > 0){
					
						lbs.at(old) = bucket <= 0 ? 0 : s.at(old)->ranks.at(bucket-1);
						ubs.at(old) = bucket < 0 ? 0 : s.at(old)->ranks.at(bucket);
						
					} else {
					
						lbs.at(old) = s.at(old)->getRank(p, QueryMode::LB);
						ubs.at(old) = s.at(old)->getRank(p, QueryMode::UB);
					}*/
					
					lb += lbi;
					ub += ubi;
					
					assert (ub >= lb);
				}
				
			}
			
			inline const E& getCurrentPoint(){
			
				return points.at(current);
			}
			
			inline double getRank(const QueryMode& qmode){
			
				if (qmode == QueryMode::LB)
					return lb;
			
				if (qmode == QueryMode::UB)
					return ub;
			
				const E& p = getCurrentPoint();
				
				double ret = 0;
			
				for (int i = 0; i < s.size(); ++i)
					ret += s.at(i)->getRank(p, qmode);
				
				return ret;
			}		
		};
	
		
		for (int i = 0; i < summaries.size(); ++i)
			summaries.at(i)->premerge();
		
		count = 0;
		
		for (int i = 0; i < summaries.size(); ++i)
			count += summaries.at(i)->count;
			
		for (int k = 0; k <= 1; ++k){
		
			long reserv = 0;
			
			if (k == 1){
			
				ranks.reserve(reserv);
				endpoints.reserve(reserv);
			}
		
			
			double lb1 = 0;
			double ub1 = 0;
		
			double lb2 = 0;
			double ub2 = 0;
			E p2 = domainMin;
			
			
			const long totalErr = params->eps * count; //mergeErr1+params->mergeErr2;
		
			
		
			for (MergeIter m(summaries, domainMin, domainMax, comp); ; m.advance() ){
		
				const bool last = m.hasReachedEnd();
			
				const long lb3 = m.lb < 0 ? 0 : m.lb > count ? count : m.lb; //last ? 0 : m.getRank(QueryMode::LB);
				const long ub3 = m.ub < 0 ? 0 : m.ub > count ? count : m.ub; //last ? 0 : m.getRank(QueryMode::UB);
				
				assert (ub3 >= lb3);
			
				const E& p3 = m.getCurrentPoint();
				
				if (last ||(ub3-lb1 > totalErr) ){
					
					assert (ub2 >= lb2);
					
					const long rankErr = ceil( (ub2-lb2)*0.5 );
				
					if (k == 1 && rankErr > absRankError)
						absRankError = rankErr;
			
					if (rankErr > params->mergeErr1){
					
						cout << (ub2-lb2)*0.5 << " > " << params->mergeErr1 << endl;
						
						assert (false);
					
						return false;
					}
					
					++reserv;	
					if (k == 1){

						ranks.push_back( lb2*0.5+ub2*0.5 );
						endpoints.push_back( p2);
					}
				
					if (last){
				
						++reserv;
						if (k == 1){
							
							ranks.push_back(count);
							
							setState(2);
							
							return true;
							
						} else {
						
							break;
						}
						
					}
					
					lb1 = lb2;
					ub1 = ub2;
				
					if (k == 1 && (ub3-lb1 > totalErr) ){
					
						cout << m.currentSummary().mode << " NAIV " << RankSummaryMode::Q_NAIVE << " EQUI " << RankSummaryMode::Q_EQUI << endl;
						cout << "p2 " << UtilString::iterableToString(p2) << " lst " << UtilString::iterableToString(endpoints.back() ) << " p3 " << UtilString::iterableToString(p3) << " ub3 " << ub3 << "-" << lb1 << " lb1 > " << params->mergeErr2 << endl;
					
						assert (false);
					
						return false;
					}
				}
				
				if (!comp(p2,p3))
					continue;
			
				lb2 = lb3;
				ub2 = ub3;
				p2 = p3;
			
			}
		}
		
		assert (false);
	
		
	}
	
	inline void init(shared_ptr<RankSummaryParameters> par){
		
		params = par;
		
		//cout << "EPS " << params.eps << endl;
		
		/*
		if (params.modeGK && params.datasize >= 0){
		
		
			if (params.gk == 0 ||params.datasize == 0 ||params.eps/params.gk < (params.naive/params.datasize) ){
				
				params.modeGK = false;
				params.modeNaive = true;
			}
		}*/
		
		/*
		assert ( params.modeGK+params.modeNaive != 2);
		assert ( params.modeGK + (params.modeEqui && !params.keepEqui) != 2 );
		
		assert ( params.modeNaive + (params.modeEqui && !params.keepEqui) != 2);*/
		
		setState(0);
		
		
	
		E v1 = domainMin;
		E v2 = domainMax;

		for (int i = 0; i < dims; i++){
			v1[i] = 0;
			v2[i] = 1;
		}
			
		v1[dims-1] = UtilMath::lowestDoubleInteger();
		v2[dims-1] = UtilMath::largestDoubleInteger();
				
		if (mode == RankSummaryMode::Q_NAIVE){
			
			//endpoints.reserve(params.datasize);
			
		} else if (mode == RankSummaryMode::Q_GK){
		
			gkeps = params->gkErr1*0.5 / params->datasizeUB;
		
			gk.reset( new GK2<E, F>( gkeps, v1, v2, comp) );
			
		} else if (mode == RankSummaryMode::Q_EQUI){
		
			assert (params->datasizeUB > 0);
		
			equi.reset( new PointEquiwidth<E, F>(dim, params->equibuckets, v1, v2, comp, params->equiErr) );
		
		} else if (mode == RankSummaryMode::Q_SKETCH){
		
			assert (params->datasizeUB > 0);
		
			sketch.reset( new QuantileSketch<E, F>(dim, v1, v2, comp) );
			
			sketch->init(params);
		}
		
		setState(1);
	}
	
	inline long getBytes() const {
	
		long ret = (endpoints.size()+1)*(dims*8+8);
		
		if (gk)
			assert (false);
		
		
		if (sketch)
			assert (false);
	
		if (equi)
			ret += equi->getBytes();
			
		return ret;
	}
	
	inline void reserve(long n){
	
		endpoints.reserve(n);
	}
	
	
	inline void startCounting(){
	
		counts.clear();
	}
	
	inline void countPoint(const E& p){
		
		if (counts.size() == 0){
			
			counts.reserve(ranks.size() );
			
			for (long j = 0; j < ranks.size(); j++)
				counts.push_back(0);
		}
		
		counts.at(getBucket(p))++;
	}
	
	inline double getAbsRankError() const {
	
		return absRankError;
	}
	
	inline string toString(RankSummaryMode m){
	
		switch(m){
		
			case Q_NAIVE:
				return "NAIVE";
		
			case Q_EQUI:
				return "EQUI";
		
			case Q_SKETCH:
				return "SKETCH";
			
			case Q_GK:
				return "GK";
				
			case Q_MERGE:
				return "MERGE";
		
		
			default:
				return "UNKNOWN";
		};
	
	}
	
	inline void endCounting(){
		
		for (long j = 1; j < counts.size(); j++)
			counts[j] += counts[j-1];
		
		for (long j = 0; j < counts.size(); j++){
		
			if (ranks[j]-absRankError > counts[j])
				cout << toString(mode) << " ranks[j] " << ranks[j] << "-" << absRankError << " > " << counts[j] << " counts[j]" << endl;
			
			if (ranks[j]+absRankError < counts[j])
				cout << toString(mode) << " ranks[j] " << ranks[j] << "+" << absRankError << " < " << counts[j] << " counts[j]" << endl;
			
			//assert (ranks[j]-absRankError <= counts[j]);
			
			
			
			//assert (ranks[j]+absRankError >= counts[j]);
			
			ranks[j] = counts[j];
		}
		
		counts.resize(0);
		counts.shrink_to_fit();
	
		absRankError = 0;
	}
	
	template<typename T>
	inline string str(const T& arr) const{
	
		stringstream ss;
		for (int i = 0; i < arr.size(); ++i)
			ss << " "<< arr[i];
		
		return ss.str();
	}
	
	inline void premerge(){
	
		if (mode == RankSummaryMode::Q_GK)
			finalise();
		else
			setState(2);
			
		
		if (mode == RankSummaryMode::Q_NAIVE){
		
			/*
			ranks.reserve(endpoints.size()+1);
			
			ranks.clear();
	
			for (long j = 1; j <= endpoints.size(); ++j)
				ranks.push_back(j);
		
			ranks.push_back(count);*/
		
			std::sort(endpoints.begin(), endpoints.end(), comp );
		}
			
		if (mode == RankSummaryMode::Q_EQUI)
			equi->finalise(false);
	}
	
	inline void finalise(long absEps = -1){
	
		assert (count <= 1 || comp(pmin, pmax));
		
		assert (absEps == -1);
		
		if (mode == RankSummaryMode::Q_SKETCH){
		
			if (count == 0){
			
				ranks.push_back(0);
				
				sketch.reset();	
				
			} else {
	
				sketch->finalise(endpoints, ranks, params->sketchErr);
			
				sketch.reset();
			}
		}
		
		if (mode == RankSummaryMode::Q_EQUI){
		
			if (count == 0){
			
				ranks.push_back(0);
				
				equi.reset();	
			} else {
		
				
				if (params->keepEqui){
			
					equi->finalise();
					
				} else {
			
					equi->finalise(endpoints, ranks, params->equiErr);
			
					equi.reset();
				}
			}
		}
		
		if (mode == RankSummaryMode::Q_GK){
		
			if (count == 0){
			
				ranks.push_back(0);
				gk.reset();
				
			} else {
		
				gk->finalize();
		
				const double step = params->gkErr2/count;
			
				absRankError = ceil(gkeps*count);
				
				/*
				if (step != step ||step <= 0){
					
					cout << "step " << step << " params.eps " << params.eps << " params.gk " << " " << params.gk << " datasize " << params.datasize << " count " << count << endl;
				}*/
			
				assert (step > 0);
				
				long res = 0;
			
				for (double rank = step; rank <= 1.0-step; rank += step)
					++res;
			
				endpoints.reserve(res);
				ranks.reserve(res+1);
				
				for (double rank = step; rank <= 1.0-step; rank += step){
					endpoints.push_back( gk->query_for_value(rank) );
					ranks.push_back(rank * count);
				}
				
				ranks.push_back(count);
			
				gk.reset();
			}	
		}
	
		if (mode == RankSummaryMode::Q_NAIVE){
		
			const long absErr = params->naiveErr;

			std::sort( endpoints.begin(), endpoints.end(), comp);
	
			assert(endpoints.size() == count);
			
			if (ranks.size() == 0){
				
				ranks.reserve(endpoints.size()+1);
			
				ranks.clear();
			
				for (long j = 1; j <= endpoints.size(); ++j)
					ranks.push_back(j);
				
				ranks.push_back(count);
			}
			
			{
	
				long last = 0;
	
				long k = 0;
	
				for (long j = 0; j < endpoints.size(); ++j){
				
					if (j < (endpoints.size()-1)){
				
						/*if (!cmp(endpoints.at(j), endpoints.at(j+1)))
							continue;*/
						
						if (ranks.at(j+1)-last <= absErr)
							continue;
					}
			
					last = ranks.at(j);
					ranks.at(k) = last;
					endpoints.at(k) = endpoints.at(j);
			
					++k;
				}
			
				ranks.at(k) = count;
				++k;
		
				ranks.resize(k);
				endpoints.resize(k-1);
		
				ranks.shrink_to_fit();
				endpoints.shrink_to_fit();
			
			}
			
			assert (ranks.back() == count);
		}
		
		for (long j = 0; j < endpoints.size(); ++j){
		
			const E& e = j > 0 ? endpoints.at(j-1) : domainMin;
		
			if (comp(endpoints.at(j), e)){
			
				cout << str(e) << " " << str(endpoints.at(j)) << endl;
				assert (false);
			}
		}
		
		if (endpoints.size() > 1){
			assert (!comp(domainMax, endpoints.back()) );
			
			//assert (comp(pmin, pmax) );
		}
		
		if (params->index > 0){
		
			index.reset( new PointEquiwidth<E, F>(dim, params->index * endpoints.size(), pmin, pmax, comp) );
			
			for (long j = 0; j < endpoints.size(); ++j)
				index->add(endpoints[j]);
			
		}
		
		setState(2);
	}
	
	int state = -1;
	
	inline void setState(int j){
	
		assert (j == state+1);
		
		state = j;
	}
	
	inline bool add(vector<E>& vec, long& count){
	
		assert (state > 0);
		assert (state < 2);
	
		long k = 0;
		
		//cout << "A " << vec.size() << " MODE " << mode << " NAIVE/EQUI/SKETCH/GK " << Q_NAIVE << "/" << Q_EQUI << "/" << Q_SKETCH << "/" << Q_GK << endl;
		
		//cout << "Z " << vec.size() << " MODE " << mode << " NAIVE/EQUI/SKETCH/GK " << Q_NAIVE << "/" << Q_EQUI << "/" << Q_SKETCH << "/" << Q_GK << endl;
		
		for (long j = 0; j < vec.size(); ++j){
		
			//cout << "K " << k << " MODE " << mode << " NAIVE/EQUI/SKETCH/GK " << Q_NAIVE << "/" << Q_EQUI << "/" << Q_SKETCH << "/" << Q_GK << endl;
		
			if (add(vec.at(j)) ){
			
				--count;
			
			} else {
				
				vec.at(k) = vec.at(j);
				++k;
			}
			
			
		}
		
		//cout << "K " << k << " MODE " << mode << " NAIVE/EQUI/SKETCH/GK " << Q_NAIVE << "/" << Q_EQUI << "/" << Q_SKETCH << "/" << Q_GK << endl;
		
		vec.resize(k);
	
		return true;
	}
	
	inline bool add(const E& p){
		
		if (state == 1){
		
			if (comp(p,pmin)) // p < min
				pmin = p;
		
			if (comp(pmax,p)) // p > max
				pmax = p;
				
			//cout << "B MODE " << mode << " NAIVE/EQUI/SKETCH/GK " << Q_NAIVE << "/" << Q_EQUI << "/" << Q_SKETCH << "/" << Q_GK << endl;
			
			switch (mode){
			
				case Q_EQUI:
				
					if (equi->add(p)){
				
						++count;
						return true;
					}
					
					break;
					
				case Q_SKETCH:
				
					//cout << "C MODE " << mode << " NAIVE/EQUI/SKETCH/GK " << Q_NAIVE << "/" << Q_EQUI << "/" << Q_SKETCH << "/" << Q_GK << endl;
				
					if (sketch->add(p)){
					
						++count;
						return true;
					}
					
					break;
					
				case Q_GK:
				
					gk->feed(p);
					++count;
					return true;
					
				case Q_NAIVE:
				
					endpoints.push_back(p);
					++count;
					return true;
					
				default:
					break;
			}
			
			return false;
		}
		
		assert (state == 2);
		
		if (absRankError > 0){
		
		
			countPoint(p);
		}
		
		return true;	
	}
	
	inline long getBucket(const E& p) const{
	
		if (equi)
			return equi->getBucket(p);
		
		if (count == 0)
			return 0;
		
		if (endpoints.size() == 0)
			return 0;
		
		if (endpoints.size() == 1)
			return comp( endpoints.at(0), p) ? 1 : 0;
		
		assert (endpoints.size() >= 1);
		
		auto s = endpoints.begin();
		auto t = endpoints.end();
		
		if (index){
		
			const long b = index->getBucket(p);
			const long c1 = getCumulativeBucketCount(b-1);
			const long c2 = getCumulativeBucketCount(b);
		
			s = s+c1;
			t = s+(c2-c1);
		}
		
		auto it = std::lower_bound(s,t, p ,comp); // p <= it;
		
		if (it != endpoints.begin() ){
		
			assert( comp(*(it-1), p));
		}
		if (it == endpoints.end())
			return endpoints.size();
		
		return it-endpoints.begin();
		
	}
	
	inline double getCumulativeBucketCount(const long q) const{
		
		if (equi)
			return equi->getCumulativeBucketCount(q);
		
		const long r = ranks.size();
		
		if (r == 0 && mode == RankSummaryMode::Q_NAIVE && endpoints.size() > 0)
			return q+1;
		
		/*
		if (q == r)
			return ranks.back();*/
	
		if (! (q<r) )
			cout << "q " << q << " >= ranks " << r << endl;
		
		assert (q == q); // check for not a number ...
		assert (q < r);
		
		return q < 0 ? 0 : ranks.at(q);
	}
	
	inline double getBucketCount(const long q1, const long q2) const{
	
		if (equi)
			return equi->getBucketCount(q1, q2);
	
		return q1 <= 0 ? getCumulativeBucketCount(q2) : getCumulativeBucketCount(q2)-getCumulativeBucketCount(q1-1);
	}
	
	inline double getMaxRank() const{
		
		if (equi)
			return equi->count;
		
		return count; // ranks.size() == 0 ? 0 : ranks.back();
	}
	
	// getBucketMin( getBucket(x) ) < x <= getBucketMax( getBucket(x) )
	inline const E& getBucketMax(const long q){
		
		if (equi){
		
			temp = equi->getBucketMax(q);
			return temp;
		}
		
		if (q < 0)
			return domainMin;
		
		if (q >= endpoints.size() )
			return domainMax;
		
		return endpoints.at(q);
		
	}
	
	inline const long getBucketNum() const{
		
		if (equi)
			return equi->getBucketNum();
		
		return endpoints.size()+1;
		
	}
	
	inline const long numOfEndPoints() const{
		
		return getBucketNum()-1;
		
	}
	
	inline const E& getBucketInfimum(const long q){
		
		if (equi)
			return equi->getBucketInfimum(q);
		
		return getBucketMax(q-1);
	}
	
	inline double getRank(const E& p, QueryMode qmode=QueryMode::EST){
		
		assert (state == 2);
		
		if (equi)
			return equi->getRank(p, qmode);
		
		
		const long q = getBucket(p);
		
		/*
		if (mode == Q_MERGE && ranksLB.size() > 0){
		
			if (qmode == QueryMode::LB)
				return ranksLB.at(q);
			
			if (qmode == QueryMode::UB)
				return ranksUB.at(q);
		
			if (ranksLB.at(q) == ranksUB.at(q))
				return ranksLB.at(q);
		}*/
		
		
		const double s1 = getBucketCount(0, q-1);
		const double s2 = getBucketCount(0, q);
		
		assert (s1 <= s2);
		
		if (qmode == QueryMode::LB)
			return UtilMath::makeBetween<double>(0, s1, s1-absRankError);
		
		if (qmode == QueryMode::UB)
			return UtilMath::makeBetween<double>(s2, getMaxRank(), s2+absRankError);
		
		if (s1 == s2)
			return s2;
		
		const E& p1 = getBucketInfimum(q);
		const E& p2 = getBucketMax(q);
		
		const double w = comp.interpolation(p1, p2, p);
		
		return w <= 0 ? std::nexttoward((double) s1, (double) s2) : w >= 1 ? std::nexttoward((double) s2, (double) s1) : s1 + w*(s2-s1);
		
	}
};



class PointRankSummary{

	int dims;
	int dim;
	
	int DIMS = -1;

	unique_ptr< BasicRankSummary< Point, PointComparator > > qn;

	Point temp;
	PointComparator comp;
	
	Point domainMin;
	Point domainMax;
	
	typedef std::array<double, 1> arr1;
	typedef std::array<double, 2> arr2;
	typedef std::array<double, 3> arr3;
	typedef std::array<double, 4> arr4;
	typedef std::array<double, 5> arr5;

	arr1 temp1;
	arr2 temp2;
	arr3 temp3;
	arr4 temp4;
	arr5 temp5;
	
	arr1 ttemp1;
	arr2 ttemp2;
	arr3 ttemp3;
	arr4 ttemp4;
	arr5 ttemp5;
	
	typedef ArrayComparator<1> comp1;
	typedef ArrayComparator<2> comp2;
	typedef ArrayComparator<3> comp3;
	typedef ArrayComparator<4> comp4;
	typedef ArrayComparator<5> comp5;

	unique_ptr< BasicRankSummary<arr1, comp1> > q1;
	unique_ptr< BasicRankSummary<arr2, comp2> > q2;
	unique_ptr< BasicRankSummary<arr3, comp3> > q3;
	unique_ptr< BasicRankSummary<arr4, comp4> > q4;
	unique_ptr< BasicRankSummary<arr5, comp5> > q5;
	

public:
	inline PointRankSummary(int dim, const Point& domainMin, const Point& domainMax, const PointComparator& comp) : dim(dim), domainMin(domainMin), domainMax(domainMax), dims(domainMin.size() ), comp(comp) {
	
		
	}
	
	template<int DIMS>
	inline Point& fromArray(const std::array<double, DIMS>& a, Point& p) const{
		
		if (a.size() == p.size()+1){
		
			int k = 1;
		
			for (int i = 0; i < p.size(); ++i){
		
				if (i == dim)
					continue;
		
				p[i] = a[k];
				++k;
			}
			
		} else {
		
			p = domainMin;
		}
		
		p[dim] = a[0];
		p.id = a[a.size()-1];
		
		return p;
		
	}
	
	template<int DIMS>
	inline std::array<double, DIMS>& toArray(const Point& p, std::array<double, DIMS>& ret) const{
		
		switch (DIMS){
		
			case 1:
				ret[0] = p[dim];
				return ret;
				
		
			case 2:
				ret[0] = p[dim];
				ret[1] = p.id;
				return ret;
				
			case 3:
				ret[0] = p[dim];
				ret[1] = p[1-dim];
				ret[2] = p.id;
				return ret;
				
			case 4:
				
				switch (dim){
				
					case 0:
						ret[0] = p[dim];
						ret[1] = p[1];
						ret[2] = p[2];
						ret[3] = p.id;
						return ret;
					
					case 1:
						ret[0] = p[dim];
						ret[1] = p[0];
						ret[2] = p[2];
						ret[3] = p.id;
						return ret;
						
					case 2:
						ret[0] = p[dim];
						ret[1] = p[0];
						ret[2] = p[1];
						ret[3] = p.id;
						return ret;
						
					default:
						break;
				
				}
				
				break;
				
			case 5:
				
				switch (dim){
				
					case 0:
						ret[0] = p[dim];
						ret[1] = p[1];
						ret[2] = p[2];
						ret[3] = p[3];
						ret[4] = p.id;
						return ret;
					
					case 1:
						ret[0] = p[dim];
						ret[1] = p[0];
						ret[2] = p[2];
						ret[3] = p[3];
						ret[4] = p.id;
						return ret;
						
					case 2:
						ret[0] = p[dim];
						ret[1] = p[0];
						ret[2] = p[1];
						ret[3] = p[3];
						ret[4] = p.id;
						return ret;
						
					case 3:
						ret[0] = p[dim];
						ret[1] = p[0];
						ret[2] = p[1];
						ret[3] = p[2];
						ret[4] = p.id;
						return ret;
					
					default:
						break;
				}
				
				break;
				
			default:
				break;
		}
		
		
		ret[0] = p[dim];
		
		int k = 1;
		
		for (int i = 0; i < dim; ++i){
		
			ret[k] = p[i];
			++k;
		}
		
		for (int i = dim+1; i < (DIMS-1); ++i){
		
			ret[k] = p[i];
			++k;
		}
		
		ret[k] = p.id;
		
		return ret;
	}
	
	
	inline void init(const RankSummaryMode& mode, shared_ptr<RankSummaryParameters> params){
	
		DIMS = (params->noid ? 0 : 1)+(params->oneDim ? 1 : dims);
		
		if (DIMS == 1){
		
			q1.reset( new BasicRankSummary< arr1, comp1 >(mode, 0, toArray<1>(domainMin, temp1), toArray<1>(domainMax, ttemp1), comp1() ));
			
			q1->init(params);
			
		} else if (DIMS == 2){
		
			q2.reset( new BasicRankSummary< arr2, comp2 >(mode, 0, toArray<2>(domainMin, temp2), toArray<2>(domainMax, ttemp2), comp2() ));
			
			q2->init(params);
			
		} else if (DIMS == 3){
		
			q3.reset( new BasicRankSummary< arr3, comp3 >(mode, 0, toArray<3>(domainMin, temp3), toArray<3>(domainMax, ttemp3), comp3() ));
			
			q3->init(params);
		
		} else if (DIMS == 4){
		
			q4.reset( new BasicRankSummary< arr4, comp4 >(mode, 0, toArray<4>(domainMin, temp4), toArray<4>(domainMax, ttemp4), comp4() ));
			
			q4->init(params);
		
		} else if (DIMS == 5){
		
			q5.reset( new BasicRankSummary< arr5, comp5 >(mode, 0, toArray<5>(domainMin, temp5), toArray<5>(domainMax, ttemp5), comp5() ));
			
			q5->init(params);
		
		} else {
		
			qn.reset( new BasicRankSummary< Point, PointComparator >(mode, dim, domainMin, domainMax, comp ));
			
			qn->init(params);
		}
	
	}
	
	inline double getRank(const Point& p, QueryMode qmode=QueryMode::EST){
	
		switch(DIMS){
		
			case 1:
				return q1->getRank( toArray<1>(p, temp1), qmode );
		
			case 2:
				return q2->getRank( toArray<2>(p, temp2), qmode );
			
			case 3:
				return q3->getRank( toArray<3>(p, temp3), qmode );
				
			case 4:
				return q4->getRank( toArray<4>(p, temp4), qmode );
			
			case 5:
				return q5->getRank( toArray<5>(p, temp5), qmode );
				
			default:
				return qn->getRank(p, qmode);
		}
		
		return 0;
	}
	
	inline double getBytes() const{
	
		switch(DIMS){
		
			case 1:
				return q1->getBytes();
		
			case 2:
				return q2->getBytes();
			
			case 3:
				return q3->getBytes();
				
			case 4:
				return q4->getBytes();
			
			case 5:
				return q5->getBytes();
				
			default:
				return qn->getBytes();
		}
		
		return 0;
	}
	
	inline double getAbsRankError() const{
	
		switch(DIMS){
		
			case 1:
				return q1->getAbsRankError();
		
			case 2:
				return q2->getAbsRankError();
			
			case 3:
				return q3->getAbsRankError();
				
			case 4:
				return q4->getAbsRankError();
			
			case 5:
				return q5->getAbsRankError();
				
			default:
				return qn->getAbsRankError();
		}
	}
	
	inline bool add(const Point& p){
	
		switch(DIMS){
		
			case 1:
				return q1->add( toArray<1>(p, temp1) );
		
			case 2:
				return q2->add( toArray<2>(p, temp2) );
			
			case 3:
				return q3->add( toArray<3>(p, temp3) );
				
			case 4:
				return q4->add( toArray<4>(p, temp4) );
			
			case 5:
				return q5->add( toArray<5>(p, temp5) );
				
			default:
				return qn->add(p);
		}
		
		return false;
	}
	
	
	inline bool merge(vector<PointRankSummary*>& v){
		
		
		
		switch(DIMS){
		
			case 1:
				{
				vector<BasicRankSummary< arr1, comp1 >*> v1;
			
				for (long i = 0; i < v.size(); ++i)
					v1.push_back( v.at(i)->q1.get() );
				
				return q1->merge(v1);
				}
		
			case 2:
				{
				vector<BasicRankSummary< arr2, comp2 >*> v2;
			
				for (long i = 0; i < v.size(); ++i)
					v2.push_back( v.at(i)->q2.get() );
				
				return q2->merge(v2);
				}
			
			case 3:
				{
				vector<BasicRankSummary< arr3, comp3 >*> v3;
			
				for (long i = 0; i < v.size(); ++i)
					v3.push_back( v.at(i)->q3.get() );
				
				return q3->merge(v3);
				}
			
			case 4:
				{
				vector<BasicRankSummary< arr4, comp4 >*> v4;
			
				for (long i = 0; i < v.size(); ++i)
					v4.push_back( v.at(i)->q4.get() );
				
				return q4->merge(v4);
				}
				
			case 5:
				{
				vector<BasicRankSummary< arr5, comp5 >*> v5;
			
				for (long i = 0; i < v.size(); ++i)
					v5.push_back( v.at(i)->q5.get() );
				
				return q5->merge(v5);
				}
			default:
				{
				vector<BasicRankSummary< Point, PointComparator >*> vn;
			
				for (long i = 0; i < v.size(); ++i)
					vn.push_back( v.at(i)->qn.get() );
				
				return qn->merge(vn);
				}
		}
	}
	
	
	inline bool add(PointRankSummary& s){
		
		switch(DIMS){
		
			case 1:
				return q1->add( s.q1->endpoints, s.q1->count );
		
			case 2:
				return q2->add( s.q2->endpoints, s.q2->count );
			
			case 3:
				return q3->add( s.q3->endpoints, s.q3->count );
				
			case 4:
				return q4->add( s.q4->endpoints, s.q4->count );
			
			case 5:
				return q5->add( s.q5->endpoints, s.q5->count );
				
			default:
				return qn->add( s.qn->endpoints, s.qn->count );
		}
		
		return false;
	}
	
	inline void finalise(long x = -1){
	
		switch(DIMS){
		
			case 1:
				return q1->finalise(x);
		
			case 2:
				return q2->finalise(x);
			
			case 3:
				return q3->finalise(x);
				
			case 4:
				return q4->finalise(x);
			
			case 5:
				return q5->finalise(x);
				
			default:
				return qn->finalise(x);
		}
	}
	
	inline void startCounting(){
	
		switch(DIMS){
		
			case 1:
				return q1->startCounting();
		
			case 2:
				return q2->startCounting();
			
			case 3:
				return q3->startCounting();
				
			case 4:
				return q4->startCounting();
			
			case 5:
				return q5->startCounting();
				
			default:
				return qn->startCounting();
		}
	}
	
	inline void countPoint(const Point& p){
	
		switch(DIMS){
		
			case 1:
				return q1->countPoint(toArray<1>(p, temp1));
		
			case 2:
				return q2->countPoint(toArray<2>(p, temp2));
			
			case 3:
				return q3->countPoint(toArray<3>(p, temp3));
				
			case 4:
				return q4->countPoint(toArray<4>(p, temp4));
			
			case 5:
				return q5->countPoint(toArray<5>(p, temp5));
				
			default:
				return qn->countPoint(p);
		}
	}
	
	inline void reserve(long n){
	
		switch(DIMS){
		
			case 1:
				return q1->reserve(n);
		
			case 2:
				return q2->reserve(n);
			
			case 3:
				return q3->reserve(n);
				
			case 4:
				return q4->reserve(n);
			
			case 5:
				return q5->reserve(n);
				
			default:
				return qn->reserve(n);
		}
	}
	
	inline void endCounting(){
	
		switch(DIMS){
		
			case 1:
				return q1->endCounting();
		
			case 2:
				return q2->endCounting();
			
			case 3:
				return q3->endCounting();
				
			case 4:
				return q4->endCounting();
			
			case 5:
				return q5->endCounting();
				
			default:
				return qn->endCounting();
		}
	}
	/*
	inline void setState(int j){
	
		assert (DIMS != -1);
	
		switch(DIMS){
		
			case 2:
				return q2->setState(j);
			
			case 3:
				return q3->setState(j);
				
			case 4:
				return q4->setState(j);
			
			case 5:
				return q5->setState(j);
				
			default:
				return qn->setState(j);
		}
	}*/
	
	inline long getBucket(const Point& p){
	
		switch(DIMS){
		
			case 1:
				return q1->getBucket( toArray<1>(p, temp1) );
		
			case 2:
				return q2->getBucket( toArray<2>(p, temp2) );
			
			case 3:
				return q3->getBucket( toArray<3>(p, temp3) );
				
			case 4:
				return q4->getBucket( toArray<4>(p, temp4) );
			
			case 5:
				return q5->getBucket( toArray<5>(p, temp5) );
				
			default:
				return qn->getBucket(p);
		}
	}
	
	inline long getBucketNum() const{
	
		switch(DIMS){
		
			case 1:
				return q1->getBucketNum();
		
			case 2:
				return q2->getBucketNum();
		
			case 3:
				return q3->getBucketNum();
				
			case 4:
				return q4->getBucketNum();
			
			case 5:
				return q5->getBucketNum();
				
			default:
				return qn->getBucketNum();
		}
	}
	
	inline const Point& getBucketInfimum(const long q){
		
		switch(DIMS){
		
			case 1:
				return fromArray<1>(q1->getBucketInfimum(q), temp);
		
			case 2:
				return fromArray<2>(q2->getBucketInfimum(q), temp);
			
			case 3:
				return fromArray<3>(q3->getBucketInfimum(q), temp);
				
			case 4:
				return fromArray<4>(q4->getBucketInfimum(q), temp);
			
			case 5:
				return fromArray<5>(q5->getBucketInfimum(q), temp);
				
			default:
				return qn->getBucketInfimum(q);
		}
	}
	
	inline const Point& getBucketMax(const long q){
		
		switch(DIMS){
		
			case 1:
				return fromArray<1>(q1->getBucketMax(q), temp);
		
			case 2:
				return fromArray<2>(q2->getBucketMax(q), temp);
			
			case 3:
				return fromArray<3>(q3->getBucketMax(q), temp);
				
			case 4:
				return fromArray<4>(q4->getBucketMax(q), temp);
			
			case 5:
				return fromArray<5>(q5->getBucketMax(q), temp);
				
			default:
				return qn->getBucketMax(q);
		}
	}
	
	inline double getCumulativeBucketCount(long q){
		
		switch(DIMS){
		
			case 1:
				return q1->getCumulativeBucketCount(q);
		
			case 2:
				return q2->getCumulativeBucketCount(q);
			
			case 3:
				return q3->getCumulativeBucketCount(q);
				
			case 4:
				return q4->getCumulativeBucketCount(q);
			
			case 5:
				return q5->getCumulativeBucketCount(q);
				
			default:
				return qn->getCumulativeBucketCount(q);
		}
		
	}

	inline long getCount(){
		
		switch(DIMS){
		
			case 1:
				return q1->count;
				
			case 2:
				return q2->count;
			
			case 3:
				return q3->count;
				
			case 4:
				return q4->count;
			
			case 5:
				return q5->count;
				
			default:
				return qn->count;
		}
		
	}
	
	inline double getBucketCount(const long q1, const long q2){
	
		return q1 <= 0 ? getCumulativeBucketCount(q2) : getCumulativeBucketCount(q2)-getCumulativeBucketCount(q1-1);
	}
	
};



// rank summary that uses an equi-width histogram and a GK-summary for overflowing buckets
class RankSummary{

	typedef Point E;
	typedef PointComparator F;

	private:
	
		//unique_ptr< GK<vector<double>,VectorComparator> > gk;
		
		unique_ptr< PointRankSummary > merged;
		
		unique_ptr< PointRankSummary > buffer;
		
		unique_ptr< PointRankSummary > exact;
		
		unique_ptr< PointRankSummary > gk;
		
		unique_ptr< PointRankSummary > equi;
		
		unique_ptr< PointRankSummary > sketch;
		
		bool gkonly = true;
		
		shared_ptr<RankSummaryParameters> params;
		
		//vector<E> borders;
		
		long count = 0;
		
		E domainMin;
		E domainMax;
		
		F comp;
		
		const int dim;
		
		//const double eps; // maximal quantile selectivity
		//double rankError;
		
		bool finalised = false;
		
		/*
		double gkErr1; // quantile precision of GK summary
		double gkErr2; // precision of equidepth histogram based on GK summary
		double eqErr; // precision of equidepth histogram based on equiwidth histogram
		*/
		
		long bufferSize = 10000;
		
	public:
	
		inline void transferPoints(unique_ptr< PointRankSummary >& src, unique_ptr< PointRankSummary >& dst){
		
			dst->add(*src);
			
		}
		
		// initialise equi-width histogram, datasize needed to detect overflow, param
		inline void init(const RankSummaryParameters& pars){
		
			/*
			gkErr1 = rankError/5;
			eqErr = 0;
			gkErr2 = rankError-gkErr1-eqErr;*/
			
			params = std::make_shared<RankSummaryParameters>();
			
			(*params) = pars;
			
			params->init();
			
			//params->print();
			
			buffer.reset( new PointRankSummary(dim, domainMin, domainMax, comp));
			
			
			
			buffer->init(RankSummaryMode::Q_NAIVE, params);
			
			bufferSize = params->bufferSize;
			
			buffer->reserve(bufferSize);
			
			
			
			//bufferSize = //;params.eps == 0 ? 1L << 60 : ceil(params.naive*(1.0 / params.eps) );
			
			count = 0;
		}
	
		inline void shrink(){
		
			if (buffer->getCount() == 0)
				return;
			
			if (equi){
			
			} else if (params->equiErr > 0){
			
				equi.reset( new PointRankSummary(dim, domainMin, domainMax, comp) );
				equi->init( RankSummaryMode::Q_EQUI, params);
				
				gkonly = false;
			}
			
			if (equi)
				transferPoints(buffer, equi);
			
			if (buffer->getCount() == 0)
				return;		
				
			if (sketch){
			
			
			} else if (params->sketchErr > 0){
			
				sketch.reset( new PointRankSummary(dim, domainMin, domainMax, comp) );
				sketch->init( RankSummaryMode::Q_SKETCH, params);
				
				gkonly = false;	
			}
			
			if (sketch)
				transferPoints(buffer, sketch);
			
			if (buffer->getCount() == 0)
				return;		
				
			if (gk){
			
			} else if (params->gkErr1 > 0){
			
				gk.reset( new PointRankSummary(dim, domainMin, domainMax, comp));
			
				gk->init(RankSummaryMode::Q_GK, params);
			}
			
			if (gk)
				transferPoints(buffer, gk);
				
			if (buffer->getCount() == 0)
				return;		
			
			if (exact){
			
			
			} else {
			
				exact.reset( new PointRankSummary(dim, domainMin, domainMax, comp));
				exact->init(RankSummaryMode::Q_NAIVE, params);
			}
			
			exact->reserve(exact->getBucketNum()+1+buffer->getCount());
			transferPoints(buffer, exact);	
			
			if (buffer->getCount() == 0)
				return;		
			
			assert (false);
		}
	
		
		inline void finalise(){
		
			shrink();
			
			assert (buffer->getCount() == 0);
			buffer.reset();		
			
			if (params->merge * params->mergeErr1 * params->mergeErr2 > 0){
			
				merged.reset( new PointRankSummary(dim, domainMin, domainMax, comp) );
				merged->init( RankSummaryMode::Q_MERGE, params);
			
				vector<PointRankSummary*> summaries;
			
				if (exact)
					summaries.push_back( exact.get() );
					
				if (equi)
					summaries.push_back( equi.get() );
					
				if (sketch)
					summaries.push_back( sketch.get() );
					
				if (gk)
					summaries.push_back( gk.get() );
			
				//cout << "MERGE" << endl;
			
				const bool b = merged->merge(summaries);
				
				assert (b);
				
				//cout << "/MERGE" << endl;
				
				exact.reset();
				equi.reset();
				sketch.reset();
				gk.reset();
			
			} else {
		
				if (exact){
				
					exact->finalise();
					
					merged = std::move(exact);
			
				} else if (equi){
				
					equi->finalise();
			
				} else if (sketch){
				
					sketch->finalise();
				
				} else if (gk){
				
					gk->finalise();
					
					merged = std::move(gk);
				
				} else {
				
					assert (count == 0);
					merged.reset( new PointRankSummary(dim, domainMin, domainMax, comp));
					merged->init(RankSummaryMode::Q_NAIVE, params);
					
					merged->finalise();
				}
				
			}
			
				
				
			finalised = true;
		}
		
		
		inline RankSummary(int dim, const E& domainMin, const E& domainMax, const F& comp) : domainMin(domainMin), domainMax(domainMax), dim(dim), comp(comp) {
		
			
		}
		
		inline void startCounting(){
		
			if (merged)
				merged->startCounting();
		}
		
		inline void countPoint(const Point& p){
		
			if (merged)
				merged->countPoint(p);
		}
		
		inline void endCounting(bool b = true){
		
			if (!b)
				return;
		
			if (merged)
				merged->endCounting();
		}
		
		
		inline double getMaxQuantileCount() const{
		
			double max = 0;
			
			for (long q = 0; q < getQuantileNum(); q++){
			
				const double c = getQuantileCount(q,q, QueryMode::LB);
				
				if (c > max)
					max = c;
			}
			
			return max;
		}
		
		
		inline double getBytes() const{
	
			if (merged)
				return merged->getBytes();
	
			long ret = 0;
			
			if (buffer)
				assert(false);
			
			if (gk)
				ret += gk->getBytes();
				
			if (equi)
				ret += equi->getBytes();
			
			if (sketch)
				ret += sketch->getBytes();
				
			if (exact)
				ret += exact->getBytes();
			
	
			return ret;
		}
		
		
		inline double getEps() const{
		
			assert (false);
			//return params.eps;
		}
		
		inline long getQuantile(const E& p) const{
		
			assert (finalised);
			
			if (buffer)
				assert(false);
			
			if (merged)
				return merged->getBucket(p);
			
			assert(gkonly);
		
			return gk->getBucket(p);
		}
		
		inline const E& getQuantileInfimum(long q){
		
			assert (finalised);
			
			if (buffer)
				assert(false);
			
			if (merged)
				return merged->getBucketMax(q-1);
			
			assert(gkonly);
		
			return gk->getBucketMax(q-1);
		}
		
		inline const E& getQuantileMax(long q){
		
			assert (finalised);
			
			if (buffer)
				assert(false);
			
			if (merged)
				return merged->getBucketMax(q);
			
			assert(gkonly);
		
			return gk->getBucketMax(q);
		}
		
		inline long getQuantileNum() const{
		
			assert (finalised);
			
			if (buffer)
				assert(false);
		
			if (merged)
				return merged->getBucketNum();
			
			assert(gkonly);
			
		
			return gk->getBucketNum();
		}
		
		inline double getRank(const E& p, QueryMode qmode=QueryMode::EST) const{
		
			assert (finalised);
			
			if (buffer)
				assert(false);
		
			if (merged){
			
				/*if (mode == QueryMode::LB)
					return UtilMath::makeBetween<double>(0, count, merged->getRank(p, mode)-merged->getAbsRankError() );
				else if (mode == QueryMode::UB)
					return UtilMath::makeBetween<double>(0, count, merged->getRank(p, mode)+merged->getAbsRankError() );
				else*/
				
				return UtilMath::makeBetween<double>(0, count, merged->getRank(p, qmode) ); 
			}
			
			/*
			if (merged)
				return UtilMath::makeBetween<double>(0, count, merged->getRank(p, qmode));*/
		
			double r = 0;
			
			if (exact)
				r += exact->getRank(p, qmode);
			
			if (gk)
				r += gk->getRank(p, qmode);
			
			if (equi)
				r += equi->getRank(p, qmode);
			
			if (sketch)
				r += sketch->getRank(p, qmode);
			
			
			return UtilMath::makeBetween<double>(0, count, r);
		}
		
		inline double getNormalizedRank(const E& p, QueryMode qmode=QueryMode::EST) const{
		
		
			return count > 0 ? getRank(p, qmode)/count : 0;
		}
		
		
		
		inline double getQuantileCount(long q1, long q2, QueryMode mode) const{
		
			assert (finalised);
			
			if (merged){
			
				if (mode == QueryMode::LB)
					return UtilMath::makeBetween<double>(0, count, merged->getBucketCount(q1,q2)-2*merged->getAbsRankError() );
				else if (mode == QueryMode::UB)
					return UtilMath::makeBetween<double>(0, count, merged->getBucketCount(q1,q2)+2*merged->getAbsRankError() );
				else
					return UtilMath::makeBetween<double>(0, count, merged->getBucketCount(q1,q2) ); 
			}
				
			
			assert (gkonly);
			
			if (gk)
				return UtilMath::makeBetween<double>(0, count, gk->getBucketCount(q1,q2) );
			else
				return UtilMath::makeBetween<double>(0, count, exact->getBucketCount(q1,q2) );
				
		}
		
		
		inline double getQuantileCount(long q, QueryMode mode) const{
		
			return getQuantileCount(q, q, mode);
				
		}
		
		inline double getCount() const{
		
			return count;
		}
		
		
		inline double intervalCount(const E& qmin, const E& qmax, QueryMode mode) const{
		
			assert (finalised);
		
			double r = 0;
			
			if ( mode == QueryMode::UB)
				r = getRank(qmax, QueryMode::UB)-getRank(qmin, QueryMode::LB);
			else if ( mode == QueryMode::LB)
				r = getRank(qmax, QueryMode::LB)-getRank(qmin, QueryMode::UB);
			else
				r = getRank(qmax, QueryMode::EST)-getRank(qmin, QueryMode::EST);
			
			return UtilMath::makeBetween<double>(0, count, r);
			
		}
		
		inline bool add(const E& p){
		
			assert (!finalised);
		
			++count;
			
			if (buffer){
				
				if (buffer->getCount() >= params->bufferSize)
					shrink();
				
				if ( buffer->add(p) )
					return true;
				
			}
			
			/*
			if (equi)
				if (equi->add(p))
					return true;
			
			if (sketch)
				if (sketch->add(p))
					return true;
			
			if (gk)
				if (gk->add(p))
					return true;*/
			
			--count;
			
			return false;
			
		}
		
		inline const string toString() const{
		
			stringstream s;
			
			s <<  "min " << domainMin << endl;
			s <<  "max " << domainMax << endl;
		
			s << "count " << count << endl;
		
			//s << h.toString();
		
			return s.str();
		}
		
		friend inline std::ostream& operator<<(std::ostream &os, const RankSummary& m) { 
			
			return os << m.toString();
		}
		
		inline int compare(const Point& a, const Point& b) const{
		
			return comp(a,b);
		}
		
		/*
		inline void finalise(){
			
			assert (!finalised);
			
			
			if (equi)
				equi->finalise();
				
			if (gk)
				gk->finalise();
			
			finalised = true;
		}*/
		
		inline void print() const{
		
			assert (finalised);
			
			//cout << h << endl;
		}
};


#endif