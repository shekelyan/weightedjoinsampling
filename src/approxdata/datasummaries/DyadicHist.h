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

#ifndef HEADER_DYADICHIST_H
#define HEADER_DYADICHIST_H

#include <approxdata/datasummaries/datasummaries.h>

class TheoreticalDyadicHist : public TheoreticalSummaryModel{


private:
	
	unsigned long bytes = 0;
	long double eps = -1;
	long p = -1;
	long datasize = 0;

	const int DIMS;
	
	const int BYTES_PER_VALUE = 8;
	
	const int VALUES_PER_BOUNDARY;
	

	inline long values(int d, long p) const{
		
		assert (p >= 1);
		
		if (p == 0)
			return 1;
		
		long ret = 0;
		
		ret += VALUES_PER_BOUNDARY*((1L << p)-1); // bucket boundaries
		
		if (d == 1)
			return ret;
		
		for (int q = (d == 2 ? 2 : 3); q <= p; q++) // children of bucket
			ret += (1L << (p-q)) * values(d-1, q);
		
		return ret;
	}
	
	inline long intersectedBuckets(int d, long p) const{
		
		assert (p >= 0);
		
		if (p == 0)
			return 1;
		
		if (d == 1)
			return 2;
		
		if (p <= 2)
			return 1L << p;
		
		long ret = 4;
		
		for (int pp = 1; pp <= (p-2); pp++)
			ret += 2*intersectedBuckets(d-1, pp);
		
		return ret;
	}
	
	
	
	inline long double getBytes(long p) const{
		
		return values(DIMS, p)* BYTES_PER_VALUE; // 
	}
	
	inline long double getEps(long p){
		
		return intersectedBuckets(DIMS, p)*(1.0/(1L << p));
	}

public:

	inline TheoreticalDyadicHist(int dims) : DIMS(dims), VALUES_PER_BOUNDARY(dims){
	
		
	}
	
	inline long double getAsymptotic(long double eps) const override{
	
		return pow( log(log(1/eps)/eps), 2*DIMS-2)/eps;
	}
	
	inline virtual const string getParameterLatex() const override{
	
		return "$\\begin{bmatrix} p = "+to_string(p)+"\\end{bmatrix}$";
	};
	
	inline const string toString() const override{
 	
 		stringstream s;
 	
 		s << "dyadichist";
 		
 		
 		
 		s << "\t -size " << UtilString::removeAll(UtilString::bytesToString(bytes), ' ');
 		
 		s << "\t -eps " << UtilString::doubleToString(eps*100, 10) << "%";
 		s << "\t -eps1 " << UtilString::doubleToString(getEps1()*100, 10) << "%";;
 		
 		s << "\t -p " << p;;
 		
 		
 		return s.str();
 	}
	
	inline void reset() override{
	
		bytes = 0;
		eps = -1;
		p = -1;
		datasize = 0;
	}
	
	inline long getP() const{
	
		return p;
	}
	
	/*
	inline void setEps1(long double e1){
		
		if (e1 <= 0){
	
			bytes = 0;
			eps = 1;
			
			return;
		}
		
		p = ceil(log2(1.0/e1));
		
		bytes = getBytes(p);
		eps = getEps(p);
	}*/
	
	inline long double getEps() const override{
	
		return eps > 1 ? 1 : eps;
	}
	
	
	inline long double getEps1() const{
	
		return 1.0/(1L << p);
	}
	
	inline long getQuantileOpsPerPoint(long pp, int d) const{
	
		long ret = 0;
		
		if (d == 1)
			return 1;
		
		for (long j = 1; j <= pp; j++){
			
			ret += 1;
			ret += getQuantileOpsPerPoint(j, d-1);
		}
		
		return ret;
	}
	
	inline int getDims() const override{
	
		return DIMS;
	}
	
	inline const string getName() const override{
	
		return "dyadichist";
	}
	
	inline long getQuantileOps() const override{
		
		long n = datasize;
		
		long ret = 0;
		
		for (int i = 1; i <= DIMS; i++)
			ret += getQuantileOpsPerPoint(p, i);
		
		return ret*n; //UtilMath::multichoose<long>(DIMS, p);
	}
	
	inline void setDataSize(long n) override{
	
		datasize = n;
	}
	
	inline void setP(long newP){
	
		p = newP;
		bytes = getBytes(p);
		eps = getEps(p);
		
		assert (bytes >= 0);
	}
	
	inline void setEps(long double EPS) override{
	
		// find smallest p, s.t. eps < EPS
		for (long j = 1; j < 60; j++){
			
			setP(j);
			
			if ( eps < EPS)
				return;
		}
		
		assert (false);
		
	}
	
	inline void setBytes(long double size) override{
		
		// find largest p, s.t. bytes <= size
		
		for (long j = 2; j < 60; j++){
		
			setP(j);
		
			if ( bytes > size){
			
				setP(j-1);
				return;
			}
		}
		
		assert (false);
	}
	
	inline long double getBytes() const override{
	
		return bytes;
	}
	
	
	inline long getDataSize() const override{
	
		return datasize;
	}

};

/* 
class TheoreticalDyadicHist{

private:
	
	double bytes = -1;
	double eps = -1;
	double eps1 = -1;

	const int DIMS;

	inline long values(long buckets, int d) const{
		
		assert (buckets >= 2);
		
		assert (UtilBits::countBits(buckets) == 1);
		
		long ret = 0;
		
		ret += (DIMS+1)*(buckets+1); // bucket boundaries
		
		assert (ret >= 0);
		
		if (d == 1)
			return ret;
		
		for (int p = 1; (1L << p) <= buckets; p++){ // children of bucket
			
			const long bucketLen = 1L << p; // 2^p
			const long bucketNum = buckets >> p; // buckets/(2^p);
		
			ret += bucketNum * values(bucketLen, d-1 );
		}
		
		return ret;
	}
	
	inline long intersectedBuckets(long buckets, int d) const{
		
		assert (buckets >= 2);
		
		assert (UtilBits::countBits(buckets) == 1);
		
		if (d == 1)
			return 2;
		
		if (buckets <= 4)
			return buckets;
		
		long ret = 4;
		
		const long L = buckets >> 1;
		
		for (int p = 1; (1L << p) < L; p++){
			
			ret += 2*intersectedBuckets(1 << p, d-1);
		}
		
		return ret;
	}
	
	
	inline double getBytes(double e1) const{
	
		if (e1 <= 0)
			return numeric_limits<double>::max();
			
		if (e1 > 0.5)
			return getBytes(0.5);
		
		const long buckets = 1L << ( (long) ceil(log2(1.0/e1)) );
		
		assert (UtilBits::countBits(buckets) == 1);
		
		return values(buckets, DIMS)*8;
	}
	
	inline double getEps(double e1){
		
		if (e1 <= 0)
			return 0;
			
		if (e1 > 0.5)
			return getEps(0.5);
		
		const long buckets = 1L << ( (long) ceil(log2(1.0/e1)) );
		
		assert (UtilBits::countBits(buckets) == 1);
		
		return intersectedBuckets(buckets, DIMS)*(1.0/buckets);
	}

public:

	inline TheoreticalDyadicHist(int dims) : DIMS(dims){
	
		
	}
	
	inline double getEps1() const{
	
		return eps1;
	}
	
	inline void setEps1(double e1){
		
		if (e1 <= 0){
		
			eps1 = 0;
			bytes = numeric_limits<double>::max();
			eps = 0;
			
			return;
		}
		
		const long buckets = 1L << ( (long) ceil(log2(1.0/e1)) );
		
		eps1 = 1.0/buckets;
		
		bytes = getBytes(eps1);
		eps = getEps(eps1);
	}
	
	const double getEps() const{
	
		return eps > 1 ? 1 : eps;
	}
	
	inline void setEps(double eps){
	
		{
			BinarySearcher s;
		
			s.setOutput(eps);
		
			s.inputGreaterThan(0);
			s.inputLessThan(eps);
			
			for (int i = 0; i < 100; i++){
		
				const double x = s.getInput();
				const double y = getEps(x);
				
				//cout << "x " << x << " y " << y << endl;
				
				s.feed(x,y);
			}
			
			setEps1(s.getInputLB());
		}
	}
	
	inline void setBytes(double size){
	
		
		BinarySearcher s;
	
		s.setOutput(size);
	
		s.inputGreaterThan(0);
		s.inputLessThan(1);
		
		for (int i = 0; i < 1000; i++){
	
			const double x = s.getInput();
			
			setEps1(x);
			
			const double y = bytes;
			
			s.feed(x,y);
		}
		
		setEps1(s.getInputLB());
	}
	
	inline double getBytes() const{
	
		return bytes;
	}
	
	inline double getMegabytes() const{
	
		return getBytes()/1000000.0;
	}

};
 */


/*
class DyadicHistCalc{

private:
	class DyadicStorageFunction : public RealFunction{
	public:
		const int DIMS ;
	
		inline DyadicStorageFunction(int dims): DIMS(dims){
		
		}
		
		inline long multichoose(long n, long k) const{
		
			return UtilMath::choose<long>( n+k-1, k);
		}
		
		inline double nextPowerOfTwo(double x) const{
		
			return pow(2.0, ceil(log2(x)) );
		}
		
		
		inline long values(double eps1, int d) const{
		
			const double e1 = 1.0/(nextPowerOfTwo(1.0/eps1));
			
			const long buckets = 1/e1;
			
			long ret = (DIMS+1)*(buckets+1);
			
			if (d == 1)
				return ret;
			
			for (int i = 1; (1L << i) <= buckets; i++)
				ret += (buckets >> i) * values(e1 * (1L << i), d-1 );
			
			return ret;
		}
		
		inline double operator()(double eps1) const override{

			assert (DIMS != -1);
			
			return values(eps1, DIMS);
			
			
			
			const double e1 = 1.0/(nextPowerOfTwo(1.0/eps1));
			
			const long splitbudget = log2(1.0/e1);
			
			const long n = DIMS;
			const long k = splitbudget;
			
			// n multichoose k = number of multisets of cardinality k from a set with cardinality n
			
			return multichoose(n, k)*1.0/e1;
		}
	};
	
	class DyadicPrecisionFunction : public RealFunction{
	public:
		const int DIMS ;
	
		inline DyadicPrecisionFunction(int dims): DIMS(dims){
		
		}
		
		inline double nextPowerOfTwo(double x) const{
		
			return pow(2.0, ceil(log2(x)) );
		}
		
		inline double f(double eps1, int d) const{
		
			const double e1 = 1.0/(nextPowerOfTwo(1.0/eps1));
			
			if (e1 >= 1)
				return 1;
			
			if (d == 1)
				return 2;
			
			const long L = (1.0/e1)/2.0;
			
			long ret = 2;
			
			for (int p = 1; (1L << p) < L; p++)
				ret += f(e1 * (1L << p), d-1);
			
			ret *= 2;
			
			return ret;	
		}
		
		inline double operator()(double eps1) const override{

			assert (DIMS != -1);

			// log2(1/e1)
			
			const double e1 = 1.0/(nextPowerOfTwo(1.0/eps1));
			
			
			if (e1 >= 1)
				return 1;
			
			const double ret = f(e1, DIMS)*e1;
			
			if (ret > 1)
				return 1;
			
			return ret;
			
			
			const long l = 2*ceil( log2(1.0/eps1) )-2;
			
			long ret = 0;
			
			for (int i = 0; i <= DIMS-2; i++)
				ret += 4*ceil(pow(l*1.0, i*1.0));
			
			ret += 2*pow(l, DIMS-1);
			
			if (ret*eps1 >= 1)
				return 1;
			
			return ret*eps1;
		}
	};

public:

	const int DIMS;

	inline DyadicHistCalc(int dims) : DIMS(dims){
	
	
	}

	inline double eps1ToKb(double eps1) const{
	
		DyadicStorageFunction f(DIMS);
		
		const double kbPerVal = 8.0/1024.0;
		
		return f(eps1)*kbPerVal;
	}
	
	
	inline double kbToEps1(double kb) const{
	
		DyadicStorageFunction f(DIMS);
		
		const double kbPerVal = 8.0/1024.0;
		
		return f.inv(0, 1, kb/kbPerVal, false);
	}
	
	inline double eps1ToEps(double eps1) const{
	
		DyadicPrecisionFunction f(DIMS);
		
		return f(eps1);
	}
	
	inline double epsToEps1(double eps) const{
	
		DyadicPrecisionFunction f(DIMS);
		
		return f.inv(0, eps, eps, true);
	}
	
	inline double epsToKb(double eps) const{
	
		return eps1ToKb(epsToEps1(eps) );
	}
	
	inline double kbToEps(double kb) const{
	
		return eps1ToEps(kbToEps1(kb) );
	}
};*/


typedef RankSummary GKSummary;




class DyadicNode{

typedef unique_ptr<DyadicNode> Node;	
	
private:

	unique_ptr<GKSummary> quantiles;
	
	int state = -1;
	
	long quantileInserts = 0;
	long quantileLookups = 0;
	
	const int dim;
	const int DIMS;
	
	const double eps1n;
	const double countLB;
	const double countUB;
	const double eps1;
	
	Box boundingBox;
	
	// children, first listing all children spanning whole space
	vector<Node> children;
	
	int maxLevel = -1;
	int minLevel = 0;
	
	
	inline long getChildIndex(int level, int c) const{
	
		// (1L << level)-1 = 1+2+4+8+16... where 1, 2 etc are the number of children in each level
		return (1L << level)-1+c;
	}
	
	inline Node& getChild(int level, int c){
	
		assert (children.size() > 0);
		
		
		return children.at(getChildIndex(level, c));
	}
	
	inline const Node& getChild(int level, int c) const{
	
		assert (children.size() > 0);
		
		return children.at(getChildIndex(level, c));
	}
	
	inline bool hasChild(int level, int c) const{
		
		if (quantiles){
		
		} else {
		
			return false;
		}
		
		if (level == maxLevel)
			return false;
		
		if ( UtilDyadic::getFirst(c, maxLevel-level) >= quantiles->getQuantileNum() )
			return false;
			
		if ( getChildIndex(level, c) >= children.size())
			return false;
		
		const Node& n = children.at( getChildIndex(level, c) );
		
		if (n)
			return true;
		
		return false;
	}
	
	inline bool isLeaf() const{
	
		return dim == (DIMS-1);	
	}
	
	long fastParam;

public:

	
	inline DyadicNode(int dims, int mainDim, double eps1n, double countLB, double countUB, int fastParam=0) : fastParam(fastParam), DIMS(dims), dim(mainDim), eps1n(eps1n), countLB(countLB), countUB(countUB), eps1(eps1n/countUB){
		
		assert (eps1n == eps1n);
		assert( dim < dims);
	}
	
	
	
	inline long getQuantileLookups() const{
	
		long ret = quantileLookups;
		
		for (auto it = children.begin(); it != children.end(); it++){
		
			const Node& n = (*it);
			
			if (n)
				ret += (*it)->getQuantileLookups();
		}
		
		return ret;
	}
	
	inline long getQuantileInserts() const{
	
		long ret = quantileInserts;
		
		for (auto it = children.begin(); it != children.end(); it++){
		
			const Node& n = (*it);
			
			if (n)
				ret += (*it)->getQuantileInserts();
		}
		
		
		
		return ret;
	}
	
	
	
	inline bool isEmpty() const{
	
		if (quantiles)
			return false;
			
		return true;
	}
	
	inline void print() const{
	
		if (isEmpty())
			return;
		
		cout << UtilString::repeat("\t", dim) << "eps1 " << eps1 << " count " << quantiles->getCount() << " quantiles " << quantiles->getQuantileNum() << endl;
		
		for (int level = 0; level < maxLevel; level++){
		
			for (int j = 0; j < (1L << level); j++){
			
				const Node& n = getChild(level, j);
				
				cout << UtilString::repeat("\t", dim+1) << "child(" << level << "," << j << ") ";
				
				if (n){
				
					n->print();
				
				} else {
				
					//cout << UtilString::repeat("\t", dim+1) << "empty" << endl;
				}
				
				cout << endl;
			}
		}
	}
	
	
	inline long getSize() const{
		
		if (isEmpty())
			return 0;
		
		long ret = quantiles->getQuantileNum()+1;
		
		for (auto it = children.begin(); it != children.end(); it++){
		
			const Node& n = (*it);
			
			if (n)
				ret += n->getSize();
			
			ret++;
		}
	
		return ret;
	}
	
	inline void createChildren(){
		
		if (isEmpty())
			return;
		
		const double count = quantiles->getCount();
	
		assert (count > 0);
		
		assert (dim < (DIMS-1));
	
		const long qnum = quantiles->getQuantileNum();
		
		assert (qnum > 0);
		assert (qnum <= (1L << maxLevel) );
		assert (qnum >= (1L << (maxLevel-1) ) );
		
		for (int level = minLevel; level < maxLevel; level++){
		
			//cout << "level " << level << endl;
		
			const int shift = maxLevel-level;
			const long childNum = (1L << level);
			
			// j is the child
			for (int j = 0; j < childNum; j++){
				
				// when level == maxLevel, there are the most nodes
				
				const long qmin = UtilDyadic::getFirst(j, shift); // first covered quantile
				const long qmax = UtilDyadic::getLast(j, shift); // last covered quantile
				
				assert ( (qmin >> shift) == (qmax >> shift));				
				assert ( ( (qmax+1) >> shift) == 1+(qmax >> shift));
				
				Node n;
				
				if (qmin < qnum){
				
					const double countLB = quantiles->getQuantileCount(qmin, UtilMath::minVal<long>(qnum-1, qmax), QueryMode::LB);
					const double countUB = quantiles->getQuantileCount(qmin, UtilMath::minVal<long>(qnum-1, qmax), QueryMode::UB);
					
					n.reset(new DyadicNode(DIMS, dim+1, eps1n, countLB, countUB, fastParam) );
				}
				
				assert (children.size() == (j-1+(1L << level) ) );
				
				children.push_back( std::move(n) );
				
				//const Node& m = getChild(level, j);
				
				//assert (!isReachable || (m));
			}
		}
	}
	
	inline void setState(int newState){
	
		//if (isEmpty() )
		//	return;
			
		assert (newState == (state+1));
		
		state = newState;
		
		if (state == 1){
		
			if (quantiles){
		
				quantiles->finalise();
				
				quantiles->startCounting();
			
				const long qnum = quantiles->getQuantileNum();
			
				maxLevel = UtilDyadic::getBinaryPowUB(qnum);
				
				//if (EQUIDEPTH)
				//	minLevel = maxLevel;
			
				
				
				if (!isLeaf())
					createChildren();
					
			}
		}
		
		
		if (state == 2){
		
			if (quantiles)
				quantiles->endCounting(true);
		}
		
		if (!isLeaf() && state > 0){
			
			
			for (auto it = children.begin(); it != children.end(); it++){
			
				if (*it)
					(*it)->setState(state-1);
			}
		}
		
		
	}
	
	inline void add(const Point& p){
	
		assert (state >= 0);
		
		if (state == 0){
		
			if (isEmpty() ){
		
				if (true){
				
					const Point domainMin(DIMS, std::nexttoward(0.0, -1.0) );
					const Point domainMax(DIMS, std::nexttoward(1.0, 1.1) );
					const PointComparator domainComp(dim, true);
					
					RankSummaryParameters params;
				
					params.eps = eps1*0.5;
					
					params.datasizeLB = countLB;
					params.datasizeUB = countUB;
					
					params.datasize = -1;
					
					params.merge = 5;
					params.naive = 10;
					params.gk = 5;
					params.equi = 100;
					params.sketch = 0;
					
					quantiles.reset( new GKSummary(dim, domainMin, domainMax, domainComp) );
				
					quantiles->init(params);
				}
				
			}
			
			quantileInserts++;
			
			if (quantiles)
				quantiles->add(p);
			boundingBox.enclose(p);
		} 
		
		if ( ( !isLeaf() ) && ( state > 0 ) ){
			
			//assert (!isEmpty());
			
			if (quantiles){
			
				quantileLookups++;
			
				long q = quantiles->getQuantile(p);
			
				// add point for all aggregation granularities (lower level -> less granules)
			
				for (int level = maxLevel; level >= minLevel; level--){
				
					if (hasChild(level, q))
						getChild(level, q)->add(p);
				
					q >>= 1;
				}
			}
		}
		
		if (state == 1){
			
			if (quantiles)
				quantiles->countPoint(p); // add point for counts
		}
	}
	
	
	
private:	
	inline double childBoxCount(int level, long q, const Point& qmin, const Point& qmax, QueryMode mode) const{
	
		assert (state >= 1);
	
		if (isEmpty() )
			return 0.0;
		
		if (q >= quantiles->getQuantileNum() )
			return 0;
		
		if (hasChild(level, q)){
		
			const Node& n = getChild(level, q);
			return n->boxCount(qmin, qmax, mode);
		}
		
		if (mode == QueryMode::LB)
			return 0;
		
		if (mode == QueryMode::UB)
			return quantiles->getQuantileCount(q, mode);
		
		
		Point a = qmin;
		Point b = qmax;
			
		a[dim] = UtilMath::maxVal<double>(a[dim], quantiles->getQuantileInfimum(q)[dim]);
		b[dim] = UtilMath::minVal<double>(b[dim], quantiles->getQuantileMax(q)[dim]);
		
		Box bq(a,b);
			
		bq.intersect(boundingBox);
			
		if (bq.getVolume() == 0)
			return 0;
			
		return bq.getVolume()/boundingBox.getVolume()*quantiles->getQuantileCount(q, mode);
	}

public:


	inline double getMaxError() const{
	
		return 0;
	}
	
	

	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode) const{
	
		assert (state >= 1);
		
		if (isEmpty() )
			return 0.0;
			
		
		if (dim == (DIMS-1))
			return quantiles->intervalCount(qmin, qmax, mode);
		
		assert (quantiles->getCount() > 0);
		
		const long q1 = quantiles->getQuantile(qmin);
		const long q2 = quantiles->getQuantile(qmax);
		
		double ret = 0;
		
		ret += childBoxCount(maxLevel, q1, qmin, qmax, mode);
		
		if (q2 > q1)
			ret += childBoxCount(maxLevel, q2, qmin, qmax, mode);
		
		const long q3 = q1+1;
		const long q4 = q2-1;
		
		if (q4 >= q3)
		for (DyadicIterator it(q3, q4); !it.endReached(); it.increment() )
			ret += childBoxCount(maxLevel-it.p, it.e, qmin, qmax, mode);
		
		return UtilMath::makeBetween<double>(0, quantiles->getCount(), ret);
		
	}


	inline void drawRandomPoint(RandomGenerator& gen, Point& p) const{
	
		if (quantiles){
		
			//
		
		} else { 
		
			return; 
		}
	
		const long shift = isLeaf() ? 0 : gen.randomInteger(minLevel, maxLevel-minLevel);
		const long qmax = quantiles->getQuantileNum()-1;
		const long level = maxLevel-shift;
		
		if (level == maxLevel){
			
			long q = gen.randomInteger(0, qmax);
			
			for (int i = 0; i < DIMS; i++){
			
				const double pmin = i == dim ? quantiles->getQuantileInfimum(q)[i] : boundingBox.min[i];
				const double pmax = i == dim ? quantiles->getQuantileMax(q)[i] : boundingBox.max[i];
				
				const double w = gen.randomDouble();
				
				p[i] = w*pmin+(1-w)*pmax;
			}
			
			return;
		}
		
		long q = gen.randomInteger(0, (1L << level)-1 );
		
		while (!hasChild(level, q))
			q = gen.randomInteger(0, q-1);
		
		getChild(level, q)->drawRandomPoint(gen, p);
	}
};


class DyadicHistSummary : public DataSummary{

public:

	const long datasize;
	
	
	unique_ptr<DyadicNode> root;
	Params params;
	
	double eps1;
	double eps;
	
	const int DIMS;
	
	TheoreticalDyadicHist calc;
	
	
	inline long getSize() const{
	
		if (root)
			return root->getSize();
		else
			return 0;
	}
	
	inline double getBytes() const override{
	
		return getSize()*8;
	}
	
	
	inline void drawRandomPoint(RandomGenerator& rg, Point& p) const override{
	
		if (root)
			root->drawRandomPoint(rg, p);
	}
	
	inline DyadicHistSummary( const DataSet& data, const string& paramStr) : datasize(data.getSize() ), DIMS(data.getDims() ), calc(DIMS), params(paramStr){
	
		eps = params.get("eps")/100.0;
		const double kb = params.get("kb");
		
		eps1 = 0;
		
		if (eps > 0){
			
			
			calc.setEps(eps);
			
			eps = calc.getEps();
			eps1 = calc.getEps1();
			
		}
		
		if (kb > 0){
			
			calc.setBytes(kb*1000);
			
			eps = calc.getEps();
			eps1 = calc.getEps1();
		}
		
		observer["eps"] = eps;
		observer["eps1"] = eps1;
		observer["modelbytes"] = calc.getBytes();
		observer["p"] = calc.getP();
		
		cout << toString() << endl;
		
		root.reset( new DyadicNode(data.getDims(), 0, eps1*datasize, datasize, datasize, params.get("f") ) );
		
		const int dataScans = DIMS;//+1; //params.get("extra") ? DIMS+1 : DIMS; 
		
		Timer t("dyadichist");
		
		t.start(10, 1000000, data.getSize()*dataScans);
		
		for (int i = 0; i < dataScans; i++){
			
			root->setState(i);
			
			const string key = "datascan"+to_string(i+1);
			
			observer.begin(key);
		
			for (auto it = data.begin(); it != data.end(); it++){
			
				t.tick();
			
				root->add(*it);
			}
				
			observer.end(key);
		}
		
		t.end();
		
		root->setState(dataScans);
		
		cout << "getSize() " << getSize() << endl;
		
		observer["bytes"] = getBytes();
		observer["quantileInserts"] = root->getQuantileInserts();
		observer["quantileLookups"] = root->getQuantileLookups();
	}
	
	inline DyadicHistSummary( const shared_ptr<DataSet>& data, const string& paramStr) : DyadicHistSummary( (*data), paramStr){
	
	}
	
	inline DyadicHistSummary( const unique_ptr<DataSet>& data, const string& paramStr) : DyadicHistSummary( (*data), paramStr){
	
	}
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override{
		
		if (root)
			return root->boxCount(qmin, qmax, mode);
		else
			return 0;
	}
	
	inline void print() const{
		
		if (root)
			return root->print();
		else
			return;
	}
	
	inline void test() const{
		
		//root->test(datasize);
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const DyadicHistSummary& m) { 
    	
    	return os << m.toString();
	}
	
	inline const string toString() const {
		
		return "dyadichist "+params.toString();
	}
};

#endif
