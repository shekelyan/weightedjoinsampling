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

#ifndef SHEKELYAN_UTILDYADIC_H
#define SHEKELYAN_UTILDYADIC_H


#include <approxdata/datasummaries/datasummaries.h>

namespace UtilDyadic{

	
	inline long getFirst(long e, long p, int m = 1){
	
		return e << (p*m);
	}
	
	inline long getEnd(long e, long p, int m = 1){
	
		return (e+1) << (p*m);
	}
	
	inline long getLast(long e, long p, int m =1){
	
		return getEnd(e, p, m)-1;
	}
	
	inline long getBinaryPow(long x){
	
		assert (x > 0);
		
		return 63-__builtin_clzl( ((uint64_t) x) );
	}
	
	inline long getBinaryPowLB(long x){
	
		const int p =  getBinaryPow(x);
	
		return x == 1L << p ? p : p-1;
	}
	
	inline long getBinaryPowUB(long x){
	
		const int p =  getBinaryPow(x);
	
		return x == 1L << p ? p : p+1;
	}
	
	inline long getPowerOfTwoLB(long x){
	
		return 1L << getBinaryPowLB(x);
	}
	
	inline long getPowerOfTwoUB(long x){
	
		return 1L << getBinaryPowUB(x);
	}
	
	
	inline long getLevel(long x){
	
		return x == 0 ? 0 : __builtin_ctzl(x);
	}
	
	inline const string toString(long e, long p){
	
		stringstream ss;
		
		ss << "[" << getFirst(e, p) << "," << getLast(e,p) << "]";
	
		return ss.str();
	}
	
	inline bool decomposePre(long a, long b, long& e, long& p, int m = 1){
	
		if (a == 0){
		
			e = 0;
			p = b == 0 ? 0 : getBinaryPow(b+1);
			
			if (m > 1)
				p /= m;
			
			return true;
			
		} else if (a == b){
		
			e = a;
			p = 0;
		
			return true;
		}
	
		if (a > b)
			return false;
		
		p = getLevel(a);
		
		if (m > 1){
			p /= m;
			
			e = a >> (p*m);
			
		} else {
		
			e = a >> p;
		}
		
		return getLast(e,p,m) <= b;
	}

	inline bool decomposePost(long a, long b, long& e, long& p, int m = 1){
	
		assert (a != 0);
	
		if (a > b){
		
			return false;
			
		} else if (a == b){
		
			e = a;
			p = 0;
			
			return true;
		}
	
		p = getLevel(b+1)/m;
		e = b >> (p*m);
		
		return getFirst(e,p, m) >= a;
	}



};



class ComplementIterator{

private:
	int dim = 0;
	
	const int dims;
	
	int k = -1;
	
	bool ready = false;
	
	const vector<long> qmin;
	const vector<long> qmax;
	
	vector<long> bordermin;
	vector<long> bordermax;
	
	inline void advance(){
	
		++k;
	
		ready = false;
		
		if (k >= (1L << dims))
			return;
		
		if (k % 2 == 0){
		
			if (qmin.at(dim) == bordermin.at(dim) )
				return;
			
			min = bordermin;
			max = bordermax;
			
			max.at(dim) = qmin.at(dim)-1;
			bordermin.at(dim) = qmin.at(dim);
			ready = true;
			
		} else {
		
			if (qmax.at(dim) == bordermax.at(dim) )
				return;
			
			min = bordermin;
			max = bordermax;
			
			min.at(dim) = qmax.at(dim)+1;
			bordermax.at(dim) = qmax.at(dim);

			ready = true;
			++dim;
		}
		
	}
	
public:
	vector<long> min;
	vector<long> max;
	
	inline ComplementIterator(const vector<long>& amin, const vector<long>& amax, const vector<long>& qmin, const vector<long>& qmax) : dims(qmin.size() ), qmin(qmin), qmax(qmax){
	
		bordermin = amin;
		bordermax = amax;
		min = amin;
		max = amax;
		
		increment();
	}
	
	inline bool endReached() const{
	
		return !ready;
	}

	inline void increment(){
		
		ready = false;
		while ( (k < (1L << dims)) && !ready){
		
			advance();
		}
		
	}	

};

struct DyadicIterator{

	long a = 0;
	long b = 0;
	
	long e = 0; // bucket
	long p = 0; // resolution
	
	long count = 0;
	
	int m;
	
	bool ready = false;

	inline DyadicIterator(long a, long b, int m = 1) : a(a), b(b), m(m){
	
		ready = a <= b;
		
		if (ready)
			increment();
	}
	
	inline long getFirst() const{
	
		return UtilDyadic::getFirst(e,p, m);
	}
	
	inline long getLast() const{
	
		return UtilDyadic::getLast(e,p, m);
	}
	
	inline bool endReached() const{
	
		return !ready;
	}
	
	inline void increment(){
	
		assert (count < 256);
	
		if ( (count % 2 == 0) && UtilDyadic::decomposePre(a,b, e, p,m)){
			
			a = UtilDyadic::getLast(e, p,m)+1;
			
		} else {
		
			if (UtilDyadic::decomposePost(a,b, e, p,m)){
		
				b = UtilDyadic::getFirst(e, p,m)-1;
			
			} else if ( UtilDyadic::decomposePre(a,b, e, p,m) ){
		
				a = UtilDyadic::getLast(e, p,m)+1;
			
			} else {
				
				ready = false;
			}
		}
		
		++count;
	}

};



class CMSketch{
public:
	
	const long w; // buckets per histogram
	const long k; // histograms
	
	vector<long> data;
	
	long totalcount = 0;
	
	bool hash = true;
	
	vector<long> dummy;
	
	
	inline long size() const{
	
		return w*k;
	}
	
	inline long buckets() const{
	
		return w;
	}
	
	inline long independentCounts() const{
	
		return k;
	}
	
	inline void add(long x, const vector<long>& hashes){
	
		
	
		if (hash == false){
		
			++data[x];
		
		} else if (hashes.size() > 0){
		
			for (int kk = 0; kk < k; ++kk)
				++data[ getIndex(kk, hashes[kk] % w) ];
			
		} else {
			
			for (int kk = 0; kk < k; ++kk)
				++data.at( getIndex(kk, getHash(kk, x) ) );	
		}
		
		++totalcount;
	}
	
	
	inline void add(long x){
	
		add(x, dummy);
	}
	
	
	inline long count() const{
	
		return totalcount;
	}
	
	
	inline long getIndex(int kk, long x) const{
	
		return kk*w+x;
	}
	
	inline long getIndex(long j, long x, const vector<long>& hashes) const{
	
		if (hashes.size() > 0){
		
			return getIndex(j, hash ? hashes.at(j) % w : x );
		}
	
		return getIndex(j, hash ? getHash(j, x) : x);
	}
	
	inline long count(long x, const vector<long>& hashes) const{
	
		long ret = totalcount+1;
		
		for (int i = 0; i < k; ++i){
		
			const long c = data.at( getIndex(i,x, hashes) );	
			
			if (c < ret)
				ret = c;
		}
		
		assert (ret <= totalcount);
		
		return ret;
	}
	
	inline long count(long x) const{
	
		return count(x,dummy);
	}
	
	inline long getHash(long j, long x) const{
		
		return UtilHash::hash32_32(x, j) % w;
		
	}
	
	
	
	inline CMSketch(long k, long w) : w(w), k(k), data(w*k, 0){
	
		//cout << "k " << k << " w " << w << endl;	
	}
	
	
	inline CMSketch(long w) : CMSketch(1, w){
	
		hash = false;
	}

	inline CMSketch(double eps, double delta) : CMSketch((long) ceil(log(1.0/delta)), (long) ceil(exp(1.0)/eps)){
	
		
	}

};


class DyadicCMSketch{

public:
	
	const int m;
	const int dims;
	const long logUniverse;
	const int maxp;
	
	vector<CMSketch> sketches;
	
	vector<MultiDimIndex<true, true>> sketchInds;
	
	MultiDimIndex<> ind;
	
	MultiDimIndex<true, true> hashesInd;
	
	LongCoords lc;
	
	
	vector<long> resolutionCoords;
	vector<long> bucketCoords;
	
	vector<long> tempmin;
	vector<long> tempmax;
	
	vector<long> bordermin;
	vector<long> bordermax;
	
	vector<long> hashes;
	vector<long> dimHashes;
	
	const bool fastHashes;
	
	int maxc = 0;
	
	const int sh;
		
	const long mod;
	const long mask;
		
	vector<BasicAggregator> hashCheck;
	
	vector<std::set<long>> hashCheckSets;
	
	inline long getBucketNum(const vector<long>& resolutionCoords) const{
	
		long ret = 1;
		
		for (int i = 0; i < resolutionCoords.size(); ++i)
			ret <<= resolutionCoords.at(i);
		
		return ret;
	
	}
	
	long size = 0;
	
	vector<vector<long>> allResolutionCoords;
	
	//vector<int> allResolutionIndices;
	
	
	int queries = 0;
	
	inline DyadicCMSketch(int dims, int logUniverse, double kb, double eps0, double delta0, bool fastHashes = true, int m = 1): m(m),  fastHashes(fastHashes), dims(dims), sh(32/dims), mod(1L << sh), mask((1L << sh)-1), lc(dims), logUniverse(logUniverse),maxp(logUniverse/m), hashesInd(3), ind(dims), tempmin(dims), tempmax(dims){
	
		ind.setAllWidths( maxp+1);
		
		const long decompositionCardinality = 2*(1+((logUniverse-1)/m)*((1L << m)-1) );//pow( logUniverse-1, dims) * pow(2*m, dims);
		
		const long numberOfSketches = pow( maxp+1, dims);
		
		double delta = 1.0-pow(1-delta0, 1.0/decompositionCardinality); ///decompositionCardinality;
		double eps = eps0/( ( (1 << dims)+1)* decompositionCardinality );
		
		if (kb > 0){
		
			delta = 1.0-pow(1-delta0, 1.0/decompositionCardinality);
		
			const long w = (long) ceil(log(1.0/delta));
			
			const double kbPerSketch = kb/numberOfSketches;
			
			const long h = (kbPerSketch*1024/8)/(w*2);
		
			eps = 1.0/h;
			
		}
		
		cout << "eps " << eps << " delta " << delta << endl;
		
		for (int i = 0; i < dims; ++i){
		
			bordermin.push_back(0);
			bordermax.push_back( (1L << logUniverse)-1 );
		}
		
		{
			vector<long> resolutionCoords1;
		
			for (ind.begin(resolutionCoords1); !ind.endReached(resolutionCoords1); ind.increment(resolutionCoords1) ){
			
				assert (ind.getIndex(resolutionCoords1) == allResolutionCoords.size() );
		
				allResolutionCoords.push_back(resolutionCoords1);
			}
		}
		
		for (int j = 0; j < allResolutionCoords.size(); ++j){
		
			hashCheck.push_back(BasicAggregator() );
			hashCheckSets.push_back(std::set<long>() );
		
			const vector<long>& resolutionCoords = allResolutionCoords.at(j);
		
			assert (j == sketches.size() );
			assert (j == sketchInds.size() );
			
			vector<long> v;
			
			for (int i = 0; i < dims; ++i)
				v.push_back( 1L << (resolutionCoords.at(i)*m) );
			
			sketchInds.push_back( MultiDimIndex<true, true>(v) );
		
			if (eps*delta != 0 && (sketchInds.back().getLog2Length() > 10) )
				sketches.push_back( CMSketch( eps, delta) );
			else
				sketches.push_back( CMSketch( sketchInds.back().getLength() ) );
				
			size += sketches.back().size();
			
			if (sketches.back().independentCounts() > maxc)
				maxc = sketches.back().independentCounts();
		}
		
		hashesInd.setWidth(0, 1L << UtilDyadic::getBinaryPowUB(maxc)  );
		hashesInd.setWidth(1, 1L << UtilDyadic::getBinaryPowUB(dims) );
		hashesInd.setWidth(2, 1L << UtilDyadic::getBinaryPowUB( maxp+1) );
		
		cout << "maxc " << maxc << " dimhashes " << maxc*dims*( maxp+1) << " hind " << hashesInd.getLength() << endl;
		
		hashes.reserve(maxc);
		dimHashes.reserve(hashesInd.getLength() );
	}
	
	
	inline long translate(long x, long p) const{
	
		/*
		assert (x >= 0);
		assert (p >= 0);
		assert (p <= maxp);*/
	
		return x >> ((maxp-p)*m);
	}
	
	inline void translate(const vector<long>& v, const vector<long>& translator, vector<long>& translation){
	
		for (int i = 0; i < translator.size(); ++i)
			translation.at(i) = translate( v.at(i), translator.at(i) );
	}
	
	
	
	inline long getHash(long k, int i, int p, long bucketcoord){
	
		return UtilHash::hash32_32(bucketcoord, maxc * dims * p + k * dims + i);
	}
	
	
	inline long getHash(const vector<long>& resolutionCoords, const vector<long>& bucketCoords, long k){
	
		long h = getHash(k,0, resolutionCoords.at(0), bucketCoords.at(0) );		
		
		for (int i = 1; i < dims; i++){
			
			h = combineHashes(h, getHash(k,i, resolutionCoords.at(i), bucketCoords.at(i) ));
			
			/*
			h <<= sh;
			
			h |= mask & getHash(k,i, resolutionCoords.at(i), bucketCoords.at(i) );*/
		}
		
		//cout << UtilString::toString(bucketCoords) << " " << h << endl;
		
		return h;
	}
	
	inline void getHashesQuerying(const vector<long>& resolutionCoords, const vector<long>& bucketCoords, vector<long>& hashes){
	
		hashes.clear();
		
		for (int k = 0; k < maxc; ++k)
			hashes.push_back(getHash(resolutionCoords, bucketCoords, k));
		
	}
	
	size_t hash_combine( size_t lhs, size_t rhs ) const {
  		lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
  		return lhs;
	}
	
	inline size_t combineHashes(size_t h1, size_t h2) const{
	
		return hash_combine(h1, h2) % (1L << 32);
		//return h1 ^h2; // ^ (h2 + 0x9e3779b9 + (h1<<6) + (h1>>2));

	}
	
	inline void getHashesConstruction(const vector<long>& dimhashes, const vector<long>& resolutionCoords, vector<long>& hashes){
	
		hashes.clear();
		
		for (int k = 0; k < maxc; ++k){
		
			size_t h = dimhashes.at( hashesInd.getIndex(k, 0, resolutionCoords.at(0)  ) );
			
			for (int i = 1; i < dims; i++){
			
				h = combineHashes(h, dimhashes.at( hashesInd.getIndex(k, i, resolutionCoords.at(i)  ) ));
			
				//h <<= sh;
				
				//h |= mask & dimhashes.at( hashesInd.getIndex(k, i, resolutionCoords.at(i)  ) );
			}
			
			hashes.push_back(h);
		}
		
	}
	
	
	inline void getHashesConstruction(const vector<long>& c, vector<long>& dimhashes){
	
		
		if (dimhashes.size() == 0)
			for (int i = hashesInd.getLength(); i >= 1; --i)
				dimhashes.push_back(-1);
		
		for (int k = 0; k < maxc; ++k){
		
			for (int i = 0; i < dims; i++){
			
				const long ci = c[i];
				
				for (int p = 0; p <= maxp; ++p){
		
					dimhashes.at(hashesInd.getIndex(k,i,p)) = getHash(k,i,p, ci >> ((maxp-p)*m) );
				}
			}
		}
	}
	
	inline void checkHash(const vector<long>& c){
	
		if (fastHashes)	
			getHashesConstruction(c, dimHashes);
		
		for (int j = 0; j < allResolutionCoords.size(); ++j){
		
			if (!sketches.at(j).hash)
				continue;
		
			const vector<long>& resolutionCoords1 = allResolutionCoords.at(j);
			
			for (int i = 0; i < resolutionCoords1.size(); ++i)
				assert( resolutionCoords1.at(i) >= 0);
			
			translate(c, resolutionCoords1, tempmin);
			
			const long ind = sketchInds.at(j).getIndex(tempmin);
			
			if (hashCheckSets.at(j).count(ind) > 0)
				continue;
				
			hashCheckSets.at(j).insert(ind);
			
			CMSketch& s = sketches.at(j);
			
			hashes.clear();
			if (fastHashes && s.hash)
				getHashesConstruction(dimHashes,resolutionCoords1 , hashes);
			
			
			for (auto it = hashes.begin(); it != hashes.end(); ++it){
			
				
				hashCheck.at(j).add( (*it) );
			}
			
		}
	}
	
	inline void checkHash(){
	
		for (int j = 0; j < allResolutionCoords.size(); ++j){
		
			if (!sketches.at(j).hash)
				continue;
				
			auto it = &hashCheck.at(j);
			
			cout << "count " << ((long)(*it).count()) << "\t min " << ((long)(*it).min()) << "\t max " << ((long)(*it).max()) << endl;
		}
	
		
	}
	
	
	inline bool add(const vector<long>& c, long maxCard = -1){
	
		if (fastHashes)	
			getHashesConstruction(c, dimHashes);
		
		for (int j = allResolutionCoords.size()-1; j >= 0; --j){
		
			const vector<long>& resolutionCoords1 = allResolutionCoords.at(j);
			
			//assert( ind.getIndex(resolutionCoords) == j);
			
			translate(c, resolutionCoords1, tempmin);
			
			CMSketch& s = sketches.at(j);
			
			hashes.clear();
			if (fastHashes && s.hash)
				getHashesConstruction(dimHashes,resolutionCoords1 , hashes);
			
			const long bucket = sketchInds.at(j).getIndex(tempmin);
			
			//assert (bucket >= 0);
			
			//if (!s.hash)
			//	assert (bucket < sketchInds.at(j).getLength() );
			
			if (maxCard >= 0 && (j == allResolutionCoords.size()-1) && s.count(bucket, hashes) >= maxCard )
				return false;
			
			
			s.add( bucket , hashes);
		}
		
		return true;
	}
	
	/*
	inline void add(const vector<long>& c){
	
		if (fastHashes)	
			getHashesConstruction(c, dimHashes);
		
		
		
		for (int j = 0; j < allResolutionCoords.size(); ++j){
		
			const vector<long>& resolutionCoords1 = allResolutionCoords.at(j);
			
			//assert( ind.getIndex(resolutionCoords) == j);
			
			translate(c, resolutionCoords1, tempmin);
			
			CMSketch& s = sketches.at(j);
			
			hashes.clear();
			if (fastHashes && s.hash)
				getHashesConstruction(dimHashes,resolutionCoords1 , hashes);
			
			//assert (hashes.at(0) == getHash(resolutionCoords, tempmin, 0));
			
			//cout << "* " << hashes.at(0) << " * " << getHash(tempmin, 0) << endl;
			
			
			//cout << UtilString::toString(c)<< " _ " << UtilString::toString(resolutionCoords)<< " _ " << UtilString::toString(tempmin) << endl;
			
			//cout << UtilString::toString(sketchInds.at(j).getLengths() )<< endl;
			
			s.add( sketchInds.at(j).getIndex(tempmin) , hashes);
		}
	}*/
	
	
	
	inline long upperboundInternal(const vector<long>& resolutionCoords, const vector<long>& bucketCoords){
		
		++queries;
		//cout << UtilString::toString(resolutionCoords)<< " _ " << UtilString::toString(bucketCoords) << endl;
		
		
		const int j = ind.getIndex(resolutionCoords) ;
		
		const long bucket = sketchInds.at(j).getIndex(bucketCoords);
		
		//getHashes(resolutionCoords, bucketCoords, hashes);
		
		hashes.clear();
		if (fastHashes && sketches.at(j).hash)
			getHashesQuerying(resolutionCoords, bucketCoords, hashes);
		
		
		const long totalCount = sketches.at(0).count(0);
		
		const long ret = sketches.at(j).count( bucket, hashes );
		
		return ret > totalCount ? totalCount : ret;
	}
	
	
	inline long upperbound( const vector<long>& q){
	
		return upperboundInternal(allResolutionCoords.back(), q);
	}
	
	
	inline long lowerbound(const vector<long>& qmin, const vector<long>& qmax){
	
		long ret = sketches.at(0).count(0); 	
		
		for (ComplementIterator it(bordermin, bordermax, qmin, qmax); !it.endReached(); it.increment() ){
		
			ret -= upperbound(it.min, it.max);
			
			if (ret < 0)
				return 0;	
		}
		
		return ret;
	}
	
	inline double estimate(const vector<long>& qmin, const vector<long>& qmax){
	
		//const double w1 = 1.0/( (1L << dims)+1);
		//const double w2 = ((double)(1L << dims))/( (1L << dims)+1);
		
		queries = 0;
		
		const long lb = lowerbound(qmin, qmax);
		
		const double q1 = queries;
		
		queries = 0;
		
		const long ub = upperbound(qmin, qmax);
		
		const double q2 = queries;
		
		assert (q1 > 0);
		assert (q2 > 0);
		
		return (q2/(q1+q2))*lb + (q1/(q1+q2))*ub;
	}
	
	inline double toFloat(int ff) const{
		
		const int f = static_cast<int>(ff);
			
		const float& f_int = reinterpret_cast<const float&>(f);
		
		return f_int;
	}
	
	template<typename F>
	inline int toLong(const F& ff) const{
		
		const double f = static_cast<double>(ff);
			
		const long& f_int = reinterpret_cast<const long&>(f);
		
		return f_int;
	}
	
	template<typename F>
	inline int toInt(const F& ff) const{
		
		const float f = static_cast<float>(ff);
			
		const int& f_int = reinterpret_cast<const int&>(f);
		
		return f_int;
	}
	
	inline const vector<long>& assignDim(vector<long>& v, int dim, long x) const{
	
		v[dim] = x;
		return v;
	}
	
	/*
	inline double find(const vector<long>& qmin, int dim, long qmax, vector<long>& v, long m){
	
		v = qmin;
	
		for (BinarySearch<long,long> s(qmin.at(dim), qmax.at(dim), m, 32); s.hasConverged(); s.find( upperbound(qmin, assigned(v, dim, s.x() )) ) ){
		
			v.at(dim) = s.x();
		
		}
		
		v.at(dim) = s.x();
	}*/
	
	
	inline long upperbound(const vector<long>& qmin, const vector<long>& qmax){
	
		assert (dims <= 4);
		
		long ret = 0;
		
		ind.begin(resolutionCoords);
		ind.begin(bucketCoords);
		
		
		//cout << "qmin " << UtilString::toString(qmin)<< " _ qmax " << UtilString::toString(qmax) << endl;
		
		const long totalCount = sketches.at(0).count(0);
		
		for (DyadicIterator it0(qmin.at(0), qmax.at(0), m); !it0.endReached(); it0.increment() ){
		
			resolutionCoords.at(0) = maxp-it0.p;
			bucketCoords.at(0) = it0.e;
			
			if (dims == 1){
			
				
				ret += upperboundInternal(resolutionCoords, bucketCoords);
				
				if (ret >= totalCount)
					return totalCount;
				
				continue;
			}
		
			for (DyadicIterator it1(qmin.at(1), qmax.at(1), m); !it1.endReached(); it1.increment() ){
			
				resolutionCoords.at(1) = maxp-it1.p;
				bucketCoords.at(1) = it1.e;
				
				if (dims == 2){
			
					ret += upperboundInternal(resolutionCoords, bucketCoords);
					
					
					if (ret >= totalCount)
						return totalCount;
					
					continue;
				}
				
				for (DyadicIterator it2(qmin.at(2), qmax.at(2), m); !it2.endReached(); it2.increment() ){
		
					resolutionCoords.at(2) = maxp-it2.p;
					bucketCoords.at(2) = it2.e;
		
					if (dims == 3){
			
						ret += upperboundInternal(resolutionCoords, bucketCoords);
						
						if (ret >= totalCount)
							return totalCount;
						
						continue;
					}
					
					for (DyadicIterator it3(qmin.at(3), qmax.at(3)); !it3.endReached(); it3.increment() ){
		
						resolutionCoords.at(3) = maxp-it3.p;
						bucketCoords.at(3) = it3.e;
		
						if (dims == 4){
			
							ret += upperboundInternal(resolutionCoords, bucketCoords);
							
							if (ret >= totalCount)
								return totalCount;
							
							continue;
						}
						
						assert(false);
					}
				}
			}
		}
	
		return ret;
	}


};



#endif
