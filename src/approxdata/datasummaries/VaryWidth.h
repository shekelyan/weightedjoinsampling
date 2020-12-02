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

#ifndef SHEKELYAN_VARYWIDTH_H
#define SHEKELYAN_VARYWIDTH_H

#include <approxdata/datasummaries/datasummaries.h>


enum GridHistType{

	AVI, EQUIAVI, EQUI, VARY, EQUIVARY, EQUIVARYAVI, DYADIC
};

class GridDimensions{

public:
	
	long grid;
	
	vector<long> sorted;
	vector<long> unsorted;
	
	inline GridDimensions(long grid, vector<long> coords) : grid(grid){
	
		
		for (int i = 0; i < coords.size(); i++){
			
			sorted.push_back( coords[i] );
			unsorted.push_back( coords[i] );
		}
		
		std::sort(sorted.begin(), sorted.end() );
		
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		/*for (int i = 0; i < sorted.size(); i++)
			s << sorted[i] << " ";
			
		s << " ::: ";*/
		for (int i = 0; i < unsorted.size(); i++)
			s << log2(unsorted[i]) << " ";
			
		return s.str();
	}

	inline bool operator<(const GridDimensions& other) const{
    
    
    	for (int i = 0; i < sorted.size(); i++){
    	
    		if (sorted[i] < other.sorted[i]) 
    			return false;
    			
    		if (sorted[i] > other.sorted[i])
    			return true;
    	}
    	
    	for (int i = 0; i < unsorted.size(); i++){
    	
    		if (unsorted[i] > other.unsorted[i])
    			return true;
    		
    		if (unsorted[i] < other.unsorted[i])
    			return false;
    	}
    	
    	return false;
	}
};


class VaryGrids{

public:

	LongCoords lc;
	vector<long> superWidths;
	
	vector<long> minWidths;
	
	
	
	vector< shared_ptr<DenseGridHistogram<double>> > grids;
	
	
	vector<long> samplingOrder;
	
	
	
	inline void setAllCounts(long n){
	
		for (auto it = grids.begin(); it != grids.end(); it++){
			(*(*it)).setAllCounts(n);
		}	
	}
	
	inline long allCountsEqualTo(){
	
		long ret = -1;
	
		for (auto it = grids.begin(); it != grids.end(); it++){
			
			const long m = (*(*it)).allCountsEqualTo();
			
			if (ret != -1 && ret != m)
				return -1;
				
			ret = m;
		}
		
		return ret;
	}
	
	inline bool allCountsEqualTo(long n){
	
		for (auto it = grids.begin(); it != grids.end(); it++){
			
			
			if (!(*(*it)).allCountsEqualTo(n))
				return false;
		}	
		
		return true;
	}
	
	long totalCount = 0;
	const int DIMS;
	
	long totalBins = 0;
	
	
	bool isFinalised = false;
	
	inline long getBins() const{
	
		return totalBins;
	}
	
	inline void addQuadNoise(){
	
		double privacyBudget = 1;
		
		double lambda0 = 2; //grids.size() * 2;
		
			
		long minWidth = -1;
	
		for (int i = 0; i < DIMS; i++){
		
			if (minWidth == -1 || minWidths[i] < minWidth)
				minWidth = minWidths[i];
		}
	
		
		const long levels = ceil(log2(minWidth));
		
		lambda0 = (levels+grids.size() );
		
		for (auto it = grids.begin(); it != grids.end(); it++){
			
			(*it)->addLaplacianNoise( lambda0 );
			
			privacyBudget -= 1.0/lambda0;
		}
		
		cout << "lambda0 " << lambda0 << endl;
		cout << "priv " << privacyBudget << endl;
	
		RandomGenerator rndGen("");
		double lambda = 1.0/privacyBudget;
		
		for (int l = levels-1; l >= 0; l--){
		
			const long s = 0;
			const long t = lc.setAllCoords( (1L << l)-1 );
			
			vector<long> widths(DIMS, 1L << l);
			
		 	lambda *= 2;
		 	
		 	//lambda = (levels+grids.size() );
		 	
		 	if (l == 0)
		 		lambda = 1.0/privacyBudget;
		 	
			privacyBudget -= 1.0 / lambda;
			
			cout << "lambda " << lambda << endl;
			cout << "priv " << privacyBudget << endl;
			
			for (long c = s; c <= t; c = lc.getNext(c, s, t) ){
				
  	 	 		const double rnd = rndGen.laplace(lambda);
				
				for (auto it = grids.begin(); it != grids.end(); it++){
					
					DenseGridHistogram<double>& h = (*(*it));
					
					const long ss = getMin( widths, h.ind.getLengths(), c );
					const long tt = getMax( widths, h.ind.getLengths(), c );
					
					if (ss == tt){
					
						if (h.contains(ss))
							h.setCount(ss, rnd);
						
					} else {
						
						double sum = 0;
						long num = 0;
						
						for (long cc = ss; cc <= tt; cc = lc.getNext(cc, ss, tt) ){
						
							if (h.contains(cc)){
								sum += h.sum(cc, cc);
								num++;
							}
						}
						
						assert (num > 0);
						
						if (true){
						
							const double add = (rnd-sum)/num;
						
							for (long cc = ss; cc <= tt; cc = lc.getNext(cc, ss, tt) ){
						
								if (h.contains(cc))
									h.addCount(cc, add);
							}
						} else {
						
						
							const double mul = rnd/sum;
						
							for (long cc = ss; cc <= tt; cc = lc.getNext(cc, ss, tt) ){
						
								if (h.contains(cc))
									h.mulCount(cc, mul);		
							}
						}
						
					}
				}
			}
		} 	
		
		cout << "priv " << privacyBudget << endl;
		assert (privacyBudget >= 0);
	}
	
	inline void print(){
	
		for (auto it = grids.begin(); it != grids.end(); it++){
					
			DenseGridHistogram<double>& h = (*(*it));
			
			h.print();	
			
			cout << "sum " << h.sum() << endl;	
		}
		
		for (int i = 0; i < DIMS; i++)
			cout << "minWidths[" << i << "] " << minWidths[i] << endl;
		
		for (int i = 0; i < DIMS; i++)
			cout << "superWidths[" << i << "] " << superWidths[i] << endl;
		
	}
	
	inline void finalise(bool sorted= false){
		
		if (isFinalised)
			return;
	
		//cout << "finalise" << endl;
		
		
		for (auto it = grids.begin(); it != grids.end(); it++){
		
			if (minWidths.size() == 0){
			
				for (int i = 0; i < DIMS; i++)
					minWidths.push_back( (*it)->ind.getLengths()[i] );
			} else {
			
				for (int i = 0; i < DIMS; i++)
					minWidths[i] = UtilMath::gcd(minWidths[i], (*it)->ind.getLengths()[i]);
			
			}
		}
		
		
		vector<GridDimensions> v;
		
		int j = 0;
		for (auto it = grids.begin(); it != grids.end(); it++){
		
			v.push_back(GridDimensions(j, (*it)->ind.getLengths()));
			j++;
		}
		
		if (!sorted)
			std::sort(v.begin(), v.end() );
		
		for (auto it = v.begin(); it != v.end(); it++){
		
			//cout << it->toString() << endl;
			
			samplingOrder.push_back(it->grid);
		}
		
		isFinalised = true;
	}
	
	inline VaryGrids(int dims) : DIMS(dims), lc(dims), superWidths(dims, 0){
	
		
	}
	
	inline void addGrid(long x, long y, long z, long u){
	
		vector<long> v = {1L << x,1L << y,1L << z, 1L << u};
		
		addGrid(v);
	}
	
	
	inline void addGrid(long x, long y, long z){
	
		vector<long> v = {1L << x,1L << y,1L << z};
		
		addGrid(v);
	}
	
	inline void addGrid(long x, long y){
	
		vector<long> v = {1L << x,1L << y};
		
		addGrid(v);
	}
	
	inline void addGrid(const vector<long>& v){
		
		shared_ptr<DenseGridHistogram<double>> d( new DenseGridHistogram<double>(v, true) );
		
		
		for (int i = 0; i < DIMS; i++)
			superWidths[i] = UtilMath::maxVal(superWidths[i], v[i] );
		
		grids.push_back(d);
		
		totalBins += d->getSize();
	
	}
	
	
	inline void addDyadicGrids(long p, int b=2){
	
		long min = lc.setAllCoords(0);
		long max = lc.setAllCoords(p);

		vector<long> v(DIMS);

		for (long l = min; l <= max; l = lc.getNext(l, min, max) ){

			long sum = 0;
	
			for (int i = 0; i < DIMS; i++){
	
				const int k = lc.getDimCoord(l, i);	
		
				sum += k;
				v[i] = pow(b, k);
			}
			
			if (sum == p)
				addGrid(v);
		}
		
		finalise();
	}
	
	
	
	inline void addGrids(GridHistType type, long bins){
	
		if (type == GridHistType::EQUI){
		
			const long l = pow(bins*1.0, 1.0/(DIMS) );
		
			vector<long> v(DIMS, l);
		
			addGrid(v);
			
			finalise();
			
			return;
		}
		
		if (type == GridHistType::AVI){
		
			const long l = bins/DIMS;	
		
			for (int i = 0; i < DIMS; i++){
		
				vector<long> v(DIMS, 1);
		
				v[i] = l;
				
				addGrid(v);
			}
			
			finalise();
			
			return;
		}
		
		if (type == GridHistType::EQUIAVI){
			
			const long bins2 = bins/(DIMS+1);
			
			const long l = pow(bins*1.0, 1.0/(DIMS) );
		
			vector<long> v(DIMS, l);
		
			addGrid(v);
			
			for (int i = 0; i < DIMS; i++){
		
				vector<long> v(DIMS, 1);
		
				v[i] = l;
			
				addGrid(v);
			}
			
			finalise();
			
			return;
		}
		
		
		if (type == GridHistType::VARY){
		
			long l = pow(bins*1.0/DIMS, 1.0/(DIMS+1) );
		
			for (int i = 0; i < DIMS; i++){
		
				vector<long> v(DIMS, l);
		
				v[i] *= l;
			
				addGrid(v);
			}
			
			finalise();
			
			return;
		}
		
		if (type == GridHistType::EQUIVARY){
		
			long l = pow(bins*1.0/(DIMS+1), 1.0/(DIMS+1) );
			long l2 = 1;
			long l3 = 1;
		
			while (true){
		
				const long ll = l+1;
			
				const long ll2 = pow(ll*1.0, 1.0/DIMS);
		
				const long ll3 = pow(ll2, DIMS);
			
				if ( DIMS*pow(ll, DIMS)*ll3 + pow(ll*ll2, DIMS) > bins)
					break;
			
				l = ll;
				l2 = ll2;
				l3 = ll3;
			}
		
			vector<long> v2(DIMS, l*l2);
		
			addGrid(v2);
		
			for (int i = 0; i < DIMS; i++){
		
				vector<long> v(DIMS, l);
			
				v[i] *= l3;
			
				addGrid(v);
			}
			
			finalise();
			
			return;
		}
			
		if (type == GridHistType::EQUIVARYAVI){
		
		
			long l = pow(bins*1.0/(2*DIMS+1), 1.0/(DIMS+1) );
			long l2 = 1;
			long l3 = 1;
		
			while (true){
		
				const long ll = l+1;
			
				const long ll2 = pow(ll*1.0, 1.0/DIMS);
		
				const long ll3 = pow(ll2, DIMS);
			
				if ( 2*DIMS*pow(ll, DIMS)*ll3 + pow(ll*ll2, DIMS) > bins)
					break;
			
				l = ll;
				l2 = ll2;
				l3 = ll3;
			}
		
			vector<long> v2(DIMS, l*l2);
		
			addGrid(v2);
		
			for (int i = 0; i < DIMS; i++){
		
				vector<long> v(DIMS, l);
			
				v[i] *= l3;
			
				addGrid(v);
			}
			
			for (int i = 0; i < DIMS; i++){
		
				vector<long> v(DIMS, 1);
			
				v[i] = pow(l,DIMS)*l3;
				
				addGrid(v);
			}
			
			finalise();
			
			return;
		}
		
		if (type == GridHistType::DYADIC){
		
			long p = 0;
		
			while ( UtilMath::choose( (p+1)+DIMS-1, (p+1) )*(1L << (p+1) ) < bins)
				p++;
		
			addDyadicGrids(p);
			
			return;
		}
	}
	
	
	
	
	inline long getBinNum(const vector<long>& v) const{
	
		long ret = 1;
		for (int i = 0; i < DIMS; i++)
			ret *= v[i];
			
		return ret;
	}

	inline long translateMin(const vector<long>& widths, long l, vector<long>& v) const{
	
		
		for (int i = 0; i < DIMS; i++){
		
			const int p = superWidths[i];
			const int k = widths[i];
			
			v[i] = lc.getDimCoord(l, i) * (p/k);
		}
	}
	
	inline long getMin(const vector<long>& widths1, const vector<long>& widths2, long l) const{
	
		long ret = 0;
		
		for (int i = 0; i < DIMS; i++){
		
			const int a = widths1[i];
			const int b = widths2[i];
			
			assert (b >= a);
			
			ret = lc.setDimCoord(ret, i, lc.getDimCoord(l, i) * (b/a));
		}
		
		return ret;
	}
	
	inline long getMax(const vector<long>& widths1, const vector<long>& widths2, long l) const{
	
		long ret = 0;
		
		for (int i = 0; i < DIMS; i++){
		
			const int a = widths1[i];
			const int b = widths2[i];
			
			assert (b >= a);
			
			ret = lc.setDimCoord(ret, i, (lc.getDimCoord(l, i)+1) * (b/a)-1);
		}
		
		return ret;
	}
	
	inline long translateMax(const vector<long>& widths, long l, vector<long>& v) const{
		
		
		for (int i = 0; i < DIMS; i++){
		
			const int p = superWidths[i];
			const int k = widths[i];
			
			v[i] = ((lc.getDimCoord(l, i)+1) * (p/k))-1;
		}
	}
	
	inline bool intersect(vector<long>& v1, vector<long>& v2, const vector<long>& v3, const vector<long>& v4) const{
	
		bool ret = false;
	
		for (int i = 0; i < DIMS; i++){
		
			if (v3[i] > v1[i]){
				v1[i] = v3[i];
				ret = true;
			}
				
			if (v4[i] < v2[i]){
				v2[i] = v4[i];
				ret = true;
			}
		}
		
		return ret;
	}
	
	inline long translate(vector<long>& v, const vector<long>& widths) const{
	
		long ret = 0;
	
		for (int i = 0; i < DIMS; i++){
			
			const int p = superWidths[i];
			const int k = widths[i];
			
			ret = lc.setDimCoord(ret, i, v[i] / (p/k));
		}
		
		return ret;
	}
	
	bool withoutReplacement = false;
	
	bool repairInconsistencies = false;
	
	
	
	
	
	inline bool sample(RandomGenerator& gen, Point& point){
	
		vector<long> min(DIMS);
		vector<long> max(DIMS);
		
		vector<long> min2(DIMS);
		vector<long> max2(DIMS);
		
		vector<long> min3(DIMS);
		vector<long> max3(DIMS);
		
		
		for (int i = 0; i < DIMS; i++){
		
			min[i] = 0;
			max[i] = superWidths[i]-1;
		}
		
		bool first = true;
		
		for (auto itt = samplingOrder.begin(); itt != samplingOrder.end(); itt++){
		
			shared_ptr<DenseGridHistogram<double>> it = grids[(*itt)];
		
			//cout << UtilString::toString( it->ind.getLengths() ) << endl;
		
			const long s = translate(min, it->ind.getLengths());
			const long t = translate(max, it->ind.getLengths());
			
			translateMin(it->ind.getLengths(), s, min2);
			translateMax(it->ind.getLengths(), t, max2);
			
			//cout << UtilString::toString( min2 ) << endl;
			//cout << UtilString::toString( max2 ) << endl;
			
			for (int i = 0; i < DIMS; i++){
			
				assert (min2[i] <= min[i]);
				assert (max2[i] >= max[i]);
			}
			
			
			long c = -1;
			
			if (it->sum(s,t) <= 0){
			
				if ( first || (!repairInconsistencies))
					return false;
				
				if (repairInconsistencies){
				
					c = it->getRandom(gen, s, t);
					
					if (withoutReplacement)
						it->addCount(c, 1);
				}	
			} else {
			
				c = it->getRandomNonEmpty(gen, s, t ); //it->getMax(s, t); //
			}
			
			assert (c != -1);
			
			translateMin(it->ind.getLengths(), c, min2);
			translateMax(it->ind.getLengths(), c, max2);
			
			assert( translate(min2, it->ind.getLengths()) == c);
			assert( translate(max2, it->ind.getLengths()) == c);
			
			if (withoutReplacement)
				it->subtractCount(c, 1);
			
			intersect(min, max, min2, max2);
			
			first = false;
		}
		
		for (int i = 0; i < DIMS; i++){
			
			const double aa = gen.randomDouble();
			const double bb = 1.0-aa;
			
			point[i] = aa*UtilHist::bucketMin( min[i], superWidths[i])+ bb*UtilHist::bucketMax( max[i], superWidths[i]);
		}
		
		return true;
	}
	
	inline void addNoise(double l){
	
		for (auto it = grids.begin(); it != grids.end(); it++)
			(*it)->addLaplacianNoise( l*grids.size() );
	}
	
	inline void removeNoise(double l){
	
		for (auto it = grids.begin(); it != grids.end(); it++){
			
			DenseGridHistogram<double>& h = (*(*it));
			
			for (long j = 0; j < h.getSize(); j++){
			
				h.counts.at(j) -= l*grids.size();
			}
		}
	}
	
	
	inline shared_ptr<DataSet> sample(long s=-1){
	
		Point pt(DIMS);
		
		RandomGenerator gen("");
		
		shared_ptr<vector<Point>> v(new vector<Point>() );
		
		/*
		cout << "SAMPLING ..." << endl;
		
		Timer t("t");
		
		t.start();*/
		
		for (bool b = sample(gen, pt); b ; b = sample(gen, pt)){
		
			(*v).push_back(pt);
			
			s--;
			
			if (s == 0)
				break;
		}
		
		/*
		while (s > 0){
			
			if (sample(gen, pt)){
				
				cout << "s " << s << endl;
				(*v).push_back(pt);
				s--;
			}
			
			
		}*/
			
		//t.end();
		
		Point pmin(DIMS, 0);
		Point pmax(DIMS, 1);
		
		Box b(pmin, pmax);	
		shared_ptr<DataSet> ret( new VectorDataSet(b, v, DIMS ));
		
		ret->selectAllDims();
		
	
		return ret;
	}

	inline void count(const Point& p){
	
		totalCount++;
		
		for (auto it = grids.begin(); it != grids.end(); it++)
			(*it)->count(p);
		
		
	}
	
	
	inline void countSub(const Point& p){
	
		totalCount--;
		
		for (auto it = grids.begin(); it != grids.end(); it++)
			(*it)->countSub(p);
		
		
	}
	
	inline void count(const DataSet& data){
	
		for (auto it = data.begin(); it != data.end(); it++)
			count(*it);
		
	}
	
	inline void countSub(const DataSet& data){
	
		for (auto it = data.begin(); it != data.end(); it++)
			countSub(*it);
		
	}
};





class VarywidthSummary : public DataSummary{

public:

	const int DIMS;
	Params params;
	
	GridHistType type = GridHistType::VARY;
	
	VaryGrids grids;
	
	inline double getBytes() const override{
	
		
		return reconstructed->getSize()*8;
	}
	
	shared_ptr<DataSet> reconstructed;
	
	inline VarywidthSummary(const DataSet& data, const string s) : DIMS(data.getDims() ), params(s), grids(data.getDims() )  {
		
		double buckets = 1;
		
		if (params.get("kb") > 0){
		
			double bytes = params.get("kb")*1024;
			double bytesPerBucket = 8;
		
		
			buckets = bytes/bytesPerBucket;
		}
		
		if (params.get("bins") > 0){
		
			buckets = params.get("bins");
		}
		
		const bool AVI = params.get("avi") == 1;
		
		const bool VARY = params.get("vary") == 1;
		
		const bool DYAD = params.get("dyad") == 1;
		
		const bool EQUI = params.get("equi") == 1;
		
		
		if (DYAD){
		
			type = GridHistType::DYADIC;
		
		} else if (AVI & VARY & EQUI){
		
			type = GridHistType::EQUIVARYAVI;
			
		} else if (EQUI & VARY){
		
			type = GridHistType::EQUIVARY;
			
		} else if (EQUI & AVI){
		
			type = GridHistType::EQUIAVI;
			
		} else if (EQUI){
		
			type = GridHistType::EQUI;
		
		} else if (AVI){
		
			type = GridHistType::AVI;
			
			
		}  else {
		
			type = GridHistType::VARY;
		}
		
		grids.addGrids(type, buckets);
		
		grids.count(data);
		
		if (params.get("dp") > 0)
			grids.addNoise(	params.get("dp") );	
		
		grids.withoutReplacement = true;
		grids.repairInconsistencies = true;
		
		reconstructed = grids.sample();
	}
	
	inline VarywidthSummary(shared_ptr<DataSet>& data, const string s) : VarywidthSummary( (*data), s) {
		
		
	}
	
	inline shared_ptr<DataSet> getReconstructed() const{
	
		return reconstructed;
	}
	
	inline virtual double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override{
	
		return reconstructed->boxCount(qmin, qmax);
	}
	
	inline const string toString() const override{
		
		return "varywidth "+params.toString();
	}
	
};


#endif