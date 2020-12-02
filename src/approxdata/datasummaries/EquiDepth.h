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

#ifndef SHEKELYAN_EQUIDEPTH_H
#define SHEKELYAN_EQUIDEPTH_H


class TheoreticalEquidepth : public TheoreticalSummaryModel{

public:
	inline const string toString() const override{
 	
 		stringstream s;
 	
 		s << "equidepth";
 		
 		s << "\t -size " << UtilString::removeAll(UtilString::bytesToString(bytes), ' ');
 		
 		s << "\t -eps " << UtilString::doubleToString(eps*100, 10) << "%";
 		s << "\t -eps1 " << UtilString::doubleToString(eps1*100, 10) << "%";;
 		
 		return s.str();
 	}

private:

	long double bytes = -1;
	long double eps = -1;
	long double eps1 = -1;

	long splits = -1;
	
	long datasize = 0;

	const int DIMS;
	
	const int BYTES_PER_VALUE = 8;
	const int VALUES_PER_BOUNDARY;
	
	inline virtual const string getParameterLatex() const override{
	
		return "$$\\begin{bmatrix} b = "+to_string(splits)+"\\end{bmatrix}$";
	};
	
	inline long double intersectedBuckets(int dims, long L) const{
		
		long double ret = 2;
		
		for (int i = 2; i <= dims; i++)
			ret = 2*pow(L, i-1)+(L-2)*ret;
		
		return ret;
	}
	
	inline long double values(int dims, long L) const{
		
		long double ret = (L-1)*VALUES_PER_BOUNDARY;
		
		for (int i = 2; i <= dims; i++)
			ret = (L-1)+L*ret;
		
		return ret;
	}
	
	inline long double getEps(int dims, long L) const{
	
		return intersectedBuckets(dims, L)/pow(L, DIMS);
	}
	
	inline long double getBytes(int dims, long L) const{
	
		return values(dims, L)*BYTES_PER_VALUE;
	}
	
public:
	
	
	inline TheoreticalEquidepth(int dims) : DIMS(dims), VALUES_PER_BOUNDARY(dims){
	
	}
	
	inline long double getEps() const override{
	
		return eps > 1 ? 1 : eps;
	}
	
	inline void reset() override{
	
		bytes = -1;
		eps = -1;
		eps1 = -1;

		splits = -1;
	}
	
	inline long double getAsymptotic(long double eps) const override{
	
		return pow(1.0/eps, DIMS);
	}
	
	inline long double getBytes() const override{
	
		return bytes;
	}
	
	inline long getSplits() const{
	
		return splits;
	}
	
	
	inline void setSplits(long L){
	
		splits = L;
		bytes = getBytes(DIMS, L);
		eps1 = 1.0/L;
		eps = getEps(DIMS, L);
	}
	
	inline void setEps1(double e1){
	
		setSplits( ceil(1.0/e1) );
	}
	
	inline long double getEps1() const{
	
		return eps1;
	}
	
	inline void setBucketSelectivity(long double e){
	
		setEps1( pow(e, 1.0/DIMS) );
	}
	
	inline long getQuantileOps() const override{
		
		long n = datasize;
		
		long ret = 0;
		
		for (int i = 1; i < DIMS; i++){
		
			for (int j = 1; j < i; j++)
				ret += n; // rank-transform
			
			if (i < DIMS)
				ret += n; // add to rank-transform of last level
		}
		
		return ret;
	}
	
	inline void setBytes(long double size) override{
	
		long maxP = 60/DIMS;
		
		long L = 0;
		
		while (getBytes(DIMS, 1L << maxP) <= size)
			maxP++;
		
		for (long x = 1L << maxP; x >= 1; x >>= 1)
			if (getBytes(DIMS, L|x) <= size)
				L = L|x;
		
		setSplits(L);
	}
	
	inline void setEps(long double EPS) override{
	
		long maxP = 60/DIMS;
	
		while (getEps(DIMS, 1L << maxP) > EPS)
			maxP++;
	
		long L = (1L << (maxP+1))-1L;
		
		for (long x = 1L << maxP; x >= 1; x >>= 1)
			if (getEps(DIMS, L^x ) <= EPS)
				L = L^x;
		
		setSplits(L);
	}
	
	inline void setDataSize(long n) override{
	
		datasize = n;
	}
	
	inline long getDataSize() const override{
	
		return datasize;
	}

	inline const string getName() const override{
	
		return "equidepth";
	}

	inline int getDims() const override{
	
		return DIMS;
	}	
};

class EquidepthNode{

typedef unique_ptr<EquidepthNode> Node;	
	
private:

	unique_ptr<RankSummary > quantiles;
	
	const int dim;
	const int DIMS;
	const double eps1;
	
	const bool invDims;
	
	// children, first listing all children spanning whole space
	vector<Node> children;
	
	int maxLevel = -1;
	
	Box boundingBox;
	
	const double countLB;
	const double countUB;
	
public:
	
	inline EquidepthNode(int dims, int mainDim, double eps1, double countLB, double countUB, bool invDims = false) : DIMS(dims), dim(mainDim), eps1(eps1), countLB(countLB), countUB(countUB),invDims(invDims){
		
		assert (eps1 == eps1);
		assert( dim < dims);
	}
	
	
	inline bool isLeaf() const{
	
		return invDims ? dim == 0 : dim == (DIMS-1);
	}
	
	
	inline void addBoxes(const Box& b, shared_ptr<RangeQueries>& boxes) const{
		
		for (long q = 0; q < quantiles->getQuantileNum(); q++){
		
			Box b2 = b;
			
			b2.min[dim] = UtilMath::maxVal( b2.min[dim], quantiles->getQuantileInfimum(q)[dim] );
			b2.max[dim] = UtilMath::minVal( b2.max[dim], quantiles->getQuantileMax(q)[dim] );
			
			if (isLeaf()){
			
				boxes->add(b2);	
				
			} else {
				
				children.at(q)->addBoxes(b2, boxes);
			}
		}
	}
	
	inline void startCounting(){
	
		quantiles->startCounting();
	
		for (auto it = children.begin(); it != children.end(); it++)
			(*it)->startCounting();
	}
	
	inline void countPoint(const Point& p){
	
		quantiles->countPoint(p);
		
		if (children.size() > 0){
			
			long q = quantiles->getQuantile(p);
			
			Node& n = children.at(q);
			
			if (n)
				n->countPoint(p);
		}
	}
	
	inline void endCounting(){
	
		quantiles->endCounting();
	
		for (auto it = children.begin(); it != children.end(); it++)
			(*it)->endCounting();
	}
	
	inline bool isEmpty() const{
		
		return quantiles ? false : true;
	}
	
	inline bool hasNoChildren() const{
	
		return children.size() == 0;
	}
	
	inline void print() const{
	
		if (isEmpty())
			return;
	
		cout << UtilString::repeat("\t", dim) << "eps1 " << eps1 << " count " << quantiles->getCount() << " quantiles " << quantiles->getQuantileNum() << endl;
		
		for (int j = 0; j < (1L << maxLevel); j++){
		
			const Node& n = children.at(j);
			
			cout << UtilString::repeat("\t", dim+1) << "child(" << maxLevel << "," << j << ") ";
			
			if (n){
			
				n->print();
			
			} else {
			
				//cout << UtilString::repeat("\t", dim+1) << "empty" << endl;
			}
			
			cout << endl;
		}
	}
	
	inline double intervalCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const{
	
		if (isEmpty())
			return 0.0;
	
		return quantiles->intervalCount(qmin, qmax, mode);
	}
	
	inline long getBytes() const{
		
		if (isEmpty())
			return 0;
		
		double ret = (quantiles->getQuantileNum()+1)*8*DIMS;
		
		for (auto it = children.begin(); it != children.end(); it++){
		
			const Node& n = (*it);
			
			if (n)
				ret += n->getBytes();
		}
		
		return ret;
	}
	
	inline void createChildren(){
		
		if (isEmpty())
			return;
			
		const double count = quantiles->getCount();
	
		assert (count > 0);
		assert (!isLeaf() );
		
		const long qnum = quantiles->getQuantileNum();
		
		assert (qnum > 0);
		assert (qnum <= (1L << maxLevel) );
		assert (qnum >= (1L << (maxLevel-1) ) );
		
		{
			//cout << "level " << level << endl;
			
			// j is the child
			for (int j = 0; j < qnum; j++){
				
				Node n;
				
				const double countLB = quantiles->getQuantileCount(j, j, QueryMode::LB);
				const double countUB = quantiles->getQuantileCount(j, j, QueryMode::UB);
					
				assert (countUB <= quantiles->getCount() );
				
				if (countUB > 0){
					
					//const double eps2 = eps1*count/UtilMath::makeBetween<double>(0, quantiles->getCount(), countUB);
						
					const int childDim = invDims ? dim-1 : dim+1;
					n.reset(new EquidepthNode(DIMS, childDim, eps1, countLB, countUB, invDims) );
				}
				
				children.push_back( std::move(n) );
				
				const Node& m = children.at(j);
			}
		}
	}
	
	inline void add(int d, const Point& p){
	
		if (d == dim){
		
			if (isEmpty() ){
		
				const Point domainMin(DIMS, std::nexttoward(0.0, -1.0) );
				const Point domainMax(DIMS, std::nexttoward(1.0, 1.1) );
				const PointComparator domainComp(dim, true);
					
				quantiles.reset( new RankSummary(dim, domainMin, domainMax, domainComp) );
				
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
				
				quantiles->init(params);
			}
			
			quantiles->add(p);
			boundingBox.enclose(p);
			
		} else {
			
			assert (!isEmpty());
			
			if (dim+1 == d)
				quantiles->countPoint(p);
			
			if (d < DIMS){
				long q = quantiles->getQuantile(p);
			
				Node& n = children.at(q);
			
				if (n)
					n->add(d, p);
			}
		}
	}
	
	inline void finalise(int d){
	
		if (isEmpty() )
			return;
		
		if (dim == d){
		
			quantiles->finalise();
			
			quantiles->startCounting();
			
			const long qnum = quantiles->getQuantileNum();
		
			maxLevel = UtilMath::getPowerOfTwoUB(qnum);
			
			if (!isLeaf() )
				createChildren();
			
			return;
		}
		
		if (dim+1 == d)
			quantiles->endCounting();
		
		if (dim == DIMS)
			return;
		
		for (auto it = children.begin(); it != children.end(); it++){
			
			Node& n = (*it);
				
			if (n)
				n->finalise(d);
		}
	}
	
private:	
	inline double childBoxCount(long q, const Point& qmin, const Point& qmax, QueryMode mode) const{
	
		if (isEmpty() )
			return 0.0;
		
		if (!isLeaf() ){
		
			const bool qminIsAligned = !quantiles->compare(quantiles->getQuantileInfimum(q), qmin);
			const bool qmaxIsAligned = !quantiles->compare(qmax, quantiles->getQuantileMax(q));
			
			if (qminIsAligned && qmaxIsAligned)
				return children.at(q)->boxCount(qmin, qmax, mode);
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
			
		return bq.getVolume() / boundingBox.getVolume() * quantiles->getQuantileCount(q, mode);
	}
	
public:	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode) const{
		
		if (isEmpty() )
			return 0.0;
		
		//cout << qmin << " " << qmax << endl;
		
		assert (quantiles->getCount() > 0);
		
		if (isLeaf() )
			return intervalCount(qmin, qmax, mode);
		
		const long q1 = quantiles->getQuantile(qmin);
		const long q2 = quantiles->getQuantile(qmax);
		
		double ret = 0;
		
		for (long j = q1; j <= q2; j++)
			ret += childBoxCount(j, qmin, qmax, mode);
		
		assert(quantiles);
		assert (ret <= quantiles->getCount());
		
		return ret;
	}
	
	inline long test(long datasize){
	
		if (isLeaf() ){

			assert (quantiles->getCount() <= datasize*pow(eps1, DIMS-1)*1.0);
			assert (quantiles->getMaxQuantileCount() <= datasize*pow(eps1, DIMS)*1.0);
			
			return quantiles->getCount();
		}
		
		long ret = 0;
		
		for (auto it = children.begin(); it != children.end(); it++){
		
			ret += (*it)->test(quantiles->getCount() );
		}
		
		return ret;
	}


	inline void drawRandomPoint(RandomGenerator& gen, Point& p) const{
	
		const long qmax = quantiles->getQuantileNum()-1;
		const long q = gen.randomInteger(0, qmax);
		
		if ( isLeaf() ){
		
			for (int i = 0; i < DIMS; i++){
			
				const double pmin = i == dim ? quantiles->getQuantileInfimum(q)[i] : boundingBox.min[i];
				const double pmax = i == dim ? quantiles->getQuantileMax(q)[i] : boundingBox.max[i];
				
				const double w = gen.randomDouble();
				
				p[i] = w*pmin+(1-w)*pmax;
			}
			
			return;
		}
		
		children.at(q)->drawRandomPoint(gen, p);
	}
};


class EquidepthSummary : public DataSummary{

public:
	
	unique_ptr<EquidepthNode> root;
	Params params;
	
	const int DIMS;	
	
	inline double getBytes() const override{
	
		if (root)
			return root->getBytes();
		else
			return 0;
	}
	
	const long datasize;
	
	inline void test(){
	
		const long ret = root->test(datasize);
	
		assert(ret == datasize);
		
		cout << "test successful" << endl;
	}
	
	TheoreticalEquidepth model;
	
	
	inline EquidepthSummary( const DataSet& data, const string paramStr) : model(data.getDims() ), datasize(data.getSize() ), DIMS(data.getDims() ), params(paramStr){
		
		bool invDims = params.get("rev");
		
		bool noCount = params.get("nc");
		
		if (params.get("sel") > 0)
			model.setBucketSelectivity(params.get("sel")/100.0);
		
		if (params.get("eps") > 0)
			model.setEps(params.get("eps")/100.0);
		
		if (params.get("kb") > 0)
			model.setBytes(params.get("kb")*1000);
		
		
		observer["eps"] = model.getEps();
		observer["eps1"] = model.getEps1();
		observer["modelbytes"] = model.getBytes();
		
		cout << "constructing " << model.toString() << endl;
		
		
		
		//cout << toString() << endl;
		
		const int initDim = invDims ? DIMS-1-0 : 0;
		
		root.reset( new EquidepthNode(data.getDims(), initDim, model.getEps1(), data.getSize(), data.getSize(), invDims) );
		
		int scans = noCount ? DIMS-1: DIMS;
		
		for (int i = 0; i <= scans; i++){
		
			const int ii = invDims ? (DIMS-1-i) : i;
			
			//cout << "start " << i << endl;
		
			for (auto it = data.begin(); it != data.end(); it++)
				root->add(ii, (*it) );
			
			root->finalise(ii);
			
			//cout << "end " << ii << endl;
		}
		
		/*
		if (!noCount){
			
		
		root->startCounting();
		
		for (auto it = data.begin(); it != data.end(); it++)
			root->countPoint( (*it) );
		
		root->endCounting();
		
		}*/
		
		observer["bytes"] = getBytes();
		
		//cout << "getBytes() " << UtilString::bytesToString(getBytes()) << endl;
	}
	
	
	inline EquidepthSummary( shared_ptr<DataSet> data, const string paramStr) : EquidepthSummary( (*data), paramStr){
	
		
	}
	
	
	inline void addBoxes(shared_ptr<RangeQueries> boxes) const{
	
		Point dmin(DIMS, 0);
		Point dmax(DIMS, 1);
	
		const Box b(dmin, dmax);
	
		root->addBoxes(b, boxes);
	}
	
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override{
	
		return root->boxCount(qmin, qmax, mode);
	}
	
	inline void print() const{
	
		return root->print();
	}
	
	inline void drawRandomPoint(RandomGenerator& rg, Point& p) const{
	
		root->drawRandomPoint(rg, p);
	}
	
	inline const string toString2() const{
	
		stringstream s;
		
		const double kb = getBytes()*8.0;
		
		s << "eps " << model.getEps()*100 << "% ";
		s << "eps1 " << model.getEps1()*100 << "% ";
		
		s << "psiz " << UtilString::bytesToString(model.getBytes() ) << " ";
		
		s << "siz " << UtilString::bytesToString(getBytes() ) << " ";
		//s << "E[siz] " << UtilString::bytesToString(calc.eps1ToKb(eps1)*1024) << " ";
		//s << "E[eps] " << calc.kbToEps(getBytes())*100 << "% ";
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const EquidepthSummary& m) { 
    	
    	return os << m.toString();
	}
	
	inline const string toString() const override{
		
		return "equidepth "+params.toString();
	}
};


#endif