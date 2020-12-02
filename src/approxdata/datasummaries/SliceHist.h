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

#ifndef SHEKELYAN_SLICEHIST_H
#define SHEKELYAN_SLICEHIST_H

#include <approxdata/datasummaries/datasummaries.h>





class TheoreticalSliceHist : public TheoreticalSummaryModel{
public:
	inline const string toString() const override{
 	
 		stringstream s;
 	
 		s << "slicehist";
 		
 		s << "\t -size " << UtilString::removeAll(UtilString::bytesToString(bytes), ' ');
 		
 		s << "\t -eps " << UtilString::doubleToString(getEps()*100, 10) << "%";
 		//s << "\t -eps1 " << UtilString::doubleToString(eps1*100, 10) << "%";;
 		
 		s << "\t -c " << getC();
 		
 		s << "\t -c1 " << c1;
 		
 		
 		s << "\t -gridbytes " << gridCountBytes;
 		
 		
 		s << "\t -k " << getK();
 		
 		return s.str();
 	}
 	
 	const int ANCHORED_QUERIES = 1;
 	const int ALL_QUERIES = 2;
	
private:

	long double bytes = -1;
	long double eps = -1;
	long double unifeps = -1;
	long double eps1 = -1;
	long double eps0 = -1;
	int k = -1;	
	long double c = -1;
	long double c1 = -1;
	long double c2 = -1;
	long double c3 = -1;

	int queryDims;

	int supportedQueries = ALL_QUERIES;
	
 	const int DIMS;
 	
 	int gridCountBytes = 8;
 	const int rankValueBytes = 8;
 	
 	
 	
	inline long getBoundaryValues(int i) const {
	
		return DIMS; //;i == 1 ? DIMS : 1;
	}
 	
 	inline const string getParameterLatex() const override{
 	
 	
 		
 		const string S_EPS1 = UtilString::doubleToString(getEpsI(eps, k, 1),5);
 		const string S_EPS2 = UtilString::doubleToString(getEpsI(eps, k, 2),5);
 		const string S_EPS3 = UtilString::doubleToString(getEpsI(eps, k, 3),5);
 		const string S_EPS4 = UtilString::doubleToString(getEpsI(eps, k, 4),5);
 		
 		const string S_A1 = UtilString::doubleToString(getAlpha(1),5);
 		const string S_A2 = UtilString::doubleToString(getAlpha(2),5);
 		const string S_A3 = UtilString::doubleToString(getAlpha(3),5);
 		const string S_A4 = UtilString::doubleToString(getAlpha(4),5);
 		
 		if (k == 1){
 		
 			return "$\\begin{bmatrix} \\veps_1 = "+S_EPS1+" \\\\ \\alpha_1 = "+S_A1+"\\end{bmatrix}$";
 			
 		} else if (k == 2){
 		
 			return "$\\begin{bmatrix} \\veps_1 = "+S_EPS1+" \\\\\\veps_2 = "+S_EPS2+" \\\\\\alpha_1 = "+S_A1+" \\\\ \\alpha_2 = "+S_A2+"\\end{bmatrix}$";
 			
 		} else if (k == 3){
 		
 			return "$\\begin{bmatrix} \\veps_1 = "+S_EPS1+"  \\\\ \\veps_2 = "+S_EPS2+" \\\\ \\veps_3 = "+S_EPS3+" \\\\ \\alpha_1 = "+S_A1+" \\\\ \\alpha_2 = "+S_A2+" \\\\ \\alpha_3 = "+S_A3+"\\end{bmatrix}$";
 		} else {
 		
 			return "$\\begin{bmatrix} \\veps_1 = "+S_EPS1+"  \\\\ \\veps_2 = "+S_EPS2+" \\\\ \\veps_3 = "+S_EPS3+" \\\\ \\veps_4 = "+S_EPS4+" \\\\ \\alpha_1 = "+S_A1+" \\\\ \\alpha_2 = "+S_A2+" \\\\ \\alpha_3 = "+S_A3+" \\\\ \\alpha_4 = "+S_A4+"\\end{bmatrix}$";	
 		}
 	}
 	
	inline long getIntersectedSlices(int K) const{
	
		if (K == 0)
			return 1;
	
		if (supportedQueries == ANCHORED_QUERIES)
			return pow(queryDims, K);
		
		long ret = 2*queryDims;
		
		for (int i = 2; i < K; i++)
			ret *= (2*queryDims-1);
		
		return ret;
	}
 	
 	inline long double getEpsI(const long double EPS, const int K, const int I) const{
		
		
		
		
		//if (EPS >= 1)
		//	return 1;
			
		if (EPS <= 0)
			throw 0;
		
		assert (EPS >= 0);
		assert (EPS != 0);
		//assert (EPS < 1);
		
		const long double D = DIMS;
		
		long double d1 = 1;
		
		for (int j = 1; j <= (K-1); j++)
			d1 *= pow(D-1, 1 - pow(1-1/D, j) );
		
		const long double p1 = pow(1-1/D, I-1);
		const long double p2 = D * ( 1 - pow(1-1/D, K) );
		
		long double d2 = 1;
		
		for (int j = 0; j <= (I-2); j++)
			d2 *= pow(D-1, 1/D * pow(1-1/D, j) );
		
		const double ret = d2 * pow( EPS/(d1 * getIntersectedSlices(K) ), p1/p2);
	
		if (ret <= 0)
			throw 0;
			
		if (ret >= 1)
			throw 0;
	
		assert (ret > 0);
		assert (ret < 1);
		
		return ret;
	}
	
	
	inline long getSlicesI(const long double EPS, const int K, const long double C, const int I) const{
		
		//return ceil( (1.0+2.0/c) / getEpsI(EPS, K, I) );
		
		const long double A= 1.0L+1.0L/C;
		
		return ceil( A / getEpsI(EPS, K, I) );
		
		//return ceil( 1.0L / getEpsI(EPS, K, I) );
	}
	
	
	
	inline long getGridCellsI(const long double EPS, const int K, const long double C, const int I) const{
	
		return pow(getSlicesI(EPS, K, C, I), DIMS);
	}
	
	inline long double getRankEpsI(const long double EPS, const int K, const long double C, const int I) const{
	
		//return getEpsI(EPS, K, I)/C;
		
		const long double A= 1.0L+1.0L/C;
		
		return ((A-1.0L)/A) * getEpsI(EPS, K, I);
	}
	
	
	inline long getRankValuesI(const long double EPS, const int K, const long double C, const int I) const{
		
		return DIMS*getBoundaryValues(I)*ceil( 3.0L/getRankEpsI(EPS, K, C, I) );
	}
 	
 	inline long double getAlpha(int i) const{
 	
 		if (i == 1)
 			return 1.0L+1.0L/c1;
 		if (i == 2)
 			return 1.0L+1.0L/c2;
 		if (i == 3)
 			return 1.0L+1.0L/c3;
 		else
 			return 1.0L+1.0L/c;
 	}
 	
 	
 	/*
	inline vector<long double> getBestC(const long double EPSI) const{
	
		
		const long double minC = 1;
		const long double maxC = 1000000;
		
		long double min = -1;
		long double bestC = -1;
		
		for (long double c = minC; c <= maxC; c++){
		
			const long double A= 1.0L+1.0L/c;
			
			const double values = ceil(A/EPSI)+DIMS*DIMS*ceil( (A/(A-1.0L))*(3.0L/EPSI) );
			
			if ( (bestC == -1) || values < min){
				
				bestC = c;
				min = values;
			}
		}
		
		return bestC;
	
	}*/
 	
 	
 	
 	inline long double getBytesI(const long double EPS, const int K, const long C1, const long C2, const long C3, const long double C, const long double I) const{
	
		//long double rankValueBytes = ceil(log2(1.0/EPS)/8.0);
		//long double gridValueBytes = (ceil(log2(datasize*1.0L))/8.0); // ceil(log2(1.0/EPS)/8);
		
		/*
		long double rtErr = 0;
		
		long double mul = 1;
		
		for (int i = 1; i <= K; i++){
		
			const long c = i == 1 ? C1 : C;
		
			rtErr += getRankEpsI(EPS, K, c, i) * mul;
			
			mul *= getEpsI(EPS, K, i);
			mul *= (2*DIMS);
		}*/
		
		//cout << "EPS " << EPS << " rtErr " << rtErr << " C " << C << endl;
		
		const long double EPS2 = EPS;
		
		if (EPS2 <= 0)
			return numeric_limits<long double>::max();
		
		long double slices = 1;
		long double ret = 0;
		
		for (int i = I; i <= K; i++){
			
			if (slices < 0)
				return numeric_limits<long double>::max();
				
			const long c = i == 1 ? C1 : i == 2 ? C2 : i == 3 ? C3 : C;
			
			const long double a1 = slices * getRankValuesI(EPS2, K, c, i) * rankValueBytes;
			const long double a2 = slices * getGridCellsI(EPS2, K, c, i) * gridCountBytes;
			
			if (a1 < 0 || a2 < 0 || slices < 0)
				return numeric_limits<long double>::max();
				
			ret += a1;
			
			if (ret < 0)
				return numeric_limits<long double>::max();
			
			ret += a2;
			
			if (ret < 0)
				return numeric_limits<long double>::max();
			
			
			slices *= DIMS * getSlicesI(EPS2, K, c, i);
		}
		
		if (ret < 0)
			return numeric_limits<long double>::max();
		
		return ret;
	}
	
	inline long double getBytes(const long double EPS, const int K, const long C1, const long C2, const long C3, const long double C) const{
	
		return getBytesI(EPS, K, C1, C2, C3, C, 1);
	}
	
public:

	inline TheoreticalSliceHist(int dims) : DIMS(dims), queryDims(dims){
	
		
	}

	inline bool setQueryDimensionality(int dims) override{
	
		queryDims = dims;
		return true;
	}
	
	
	inline long double getAsymptotic(long double eps) const override{
	
		return pow( 1/eps, DIMS/(2.0-1.0/DIMS));
	}
	
	inline long double getAsymptotic2(long double eps) const override{
	
		return pow( 1/eps, DIMS/(2.0-1.0/DIMS));
	}
	
	inline long double getEps0(const long double EPS, const long double K, const long double C1, const long double C){
	
		return EPS;
	
		/* DEPRECATED
		long double rtErr = 0;
		
		long double mul = 2*DIMS;
		
		for (int i = 1; i <= K; i++){
		
			const long double c = i == 1 ? C1 : C;
		
			rtErr += getRankEpsI(EPS, K, c, i) * mul;
			
			mul *= getEpsI(EPS, K, i);
			mul *= (2*DIMS);
		}
		
		return EPS-rtErr; */
	}
	
	long datasize = 1L << 62;
	
	inline void setDataSize(long n) override{
	
		assert (n > 0);
	
		gridCountBytes = ceil(ceil(log2(n))/8.0);
		
		datasize = n;
	}
	
	inline void setGridCountBytes(int n){
	
		assert (n > 0);
	
		gridCountBytes = n;
	}
	
	inline long getDataSize() const override{
	
		return datasize;
	}
	
	inline long getDecompositionCardinality() const{
		
		long ret = 0;
		long a = 1;
		
		for (int i = 1; i <= (k-1); i++){
		
			ret += a;
			a *= supportedQueries == 1 ? queryDims : (i == 1) ? (2*queryDims) : (2*queryDims-1);
		}
		
		return ret+a;
	}
	
	
	inline void setSupportedQueries(int type){
	
		if (supportedQueries == type)
			return;
		
		
		eps /= getIntersectedSlices(k);
		
		supportedQueries = type;
		
		eps *= getIntersectedSlices(k);
	}
	
	inline void reset() override{
	
		bytes = -1;
		eps = -1;
		unifeps = -1;
		eps1 = -1;
		k = -1;	
		c = -1;
		c1 = -1;
		c2 = -1;
		c3 = -1;
		eps0 = -1;
		queryDims = DIMS;
		supportedQueries = ALL_QUERIES;
 		gridCountBytes = 8;
 		//rankValueBytes = 8;
	}
	

	
	inline TheoreticalSliceHist(const TheoreticalSliceHist& tsh) : DIMS(tsh.DIMS){
	
		k = tsh.k;
		eps = tsh.eps;
		eps1 = tsh.eps1;
		c = tsh.c;
		bytes = tsh.bytes;
		gridCountBytes = tsh.gridCountBytes;
		//rankValueBytes = tsh.rankValueBytes;
		unifeps = tsh.unifeps;
		datasize = tsh.datasize;
		queryDims = tsh.queryDims;
		supportedQueries = tsh.supportedQueries;
		eps0 = tsh.eps0;
		c1 = tsh.c1;
		c2 = tsh.c2;
		c3 = tsh.c3;
	}
	
	
	inline int getDims() const override{
	
		return DIMS;
	}
	
	
	inline void setK(int kk){
	
		k = kk;
	}
	
	inline int getK() const{
		
		return k;
	}
	
	inline long getSlices(int i) const{
	
		return getSlicesI(eps0, k, i == 1 ? c1 : i == 2 ? c2 : i == 3 ? c3 : c, i);
		//return getSlices(eps1, c, i);
	}
	
	
	inline double getRankEps(int i) const{
	
		return getRankEpsI(eps0, k, i == 1 ? c1 : i == 2 ? c2 : i == 3 ? c3 : c, i);
	
		//return getRankEps(eps1, c, i);
	}
	
	
	inline double getC() const{
		
		return c;
	}
	
	inline long getQuantileOps() const override{
	
		long ret = 0;
		
		const long n = datasize;
		
		for (int i = 1; i <= k; i++){
		
			for (int j = 1; j < i; j++)
				ret += n*DIMS; // rank-transform
			
			if (i <= k)
				ret += n*DIMS; // add to rank-transform of last level
		}
		
		return ret;
	}
	
	inline void setC(long double cc){
		
		c = cc;
	}
	
	
	
	inline double getEps(int i) const{
		
		return getEpsI(eps, k, i);//getEpsI(eps1, i);
	}
	
	
	 
	
	inline void setEps(long double e) override{
	
		//cout << "setEps(" << e << ")" << endl;
	
		unifeps = 1;
	
		if ( (e <= 0) ){
			
			if (k <= 0)
				k = 1;
			
			if (c <= 0)
				c = 1;
				
			c1 = c;
			c2 = c;
			c3 = c;
			
			
			eps = 1;
			eps1 = 1;
			bytes = numeric_limits<double>::max();
			
			return;
		}
		
		eps = e;
		
		
		Finder<int> bestK;
		
		const int minK = k > 0 ? k : 1;
		const int maxK = k > 0 ? k : 4;
		
		
		const long double minC = c > 0 ? c : 1;
		const long double maxC = c > 0 ? c : 1000000;
		
		Finder<long double> bestC4;
		Finder<long double> bestC3;
		Finder<long double> bestC2;
		Finder<long double> bestC1;
		
		vector<long double> cs;
			
		for (int j = 0; j <= 10; j++)
			cs.push_back(0);
			
		
		for (k = minK; k <= maxK; k++){
			
			for (int j = k; j >= 1; j--){
			
				Finder<long double> bestC;
				
				for (long double cc = minC; cc <= maxC; cc*=10){
				
					cs[j] = cc;
				
					bestC.feed(cc, getBytesI(eps, k,cs[1],cs[2],cs[3], cs[4], j) );
				}
				
				long double sofar = bestC.getMin();
				
				for (long double cc = sofar/10; cc <= sofar*10; cc++){
				
					if (cc <= 1)
						continue;
						
					cs[j] = cc;
					
					bestC.feed(cc, getBytesI(eps, k,cs[1],cs[2],cs[3], cs[4], j) );
				}
				
				long double sofar2 = bestC.getMin();
				
				for (long double cc = sofar2-1; cc <= sofar2+1; cc+=0.01){
				
					if (cc <= 1)
						continue;
						
					cs[j] = cc;
					
					bestC.feed(cc, getBytesI(eps, k,cs[1],cs[2],cs[3], cs[4], j) );
				}
				
				cs[j] = bestC.getMin();
			}
			
			const long double size = getBytes(eps, k,cs[1],cs[2],cs[3], cs[4]); //getBytes(eps1, c, k);
			
			bestC1.feed(cs[1], size);
			bestC2.feed(cs[2], size);
			bestC3.feed(cs[3], size);
			bestC4.feed(cs[4], size);
			bestK.feed(k, size);
		}
		
		/*
		
		const long double minC = c > 0 ? c : 3;
		const long double maxC = c > 0 ? c : 100000;
		
		
		Finder<long double> bestC;
		Finder<long double> bestC1;
		Finder<long double> bestEps0;
		Finder<int> bestK;
		
		for (k = minK; k <= maxK; k++){
		
			//eps1 = getEps1(k, eps);
			
			for (c1 = minC; c1 <= maxC; ){
			
				for (c = minC; c <= maxC; ){
				
					try {
						const long double EPS0 = getEps0(eps, k, c1, c);
						const long double size = getBytes(EPS0, k, c1, c); //getBytes(eps1, c, k);
						bestK.feed(k, size);
						bestC.feed(c, size);
						bestC1.feed(c1, size);
						bestEps0.feed(EPS0, size);
						
					} catch (...){
					
					}
					
					if (c >= 10000)
						c += 5000;
					else if (c >= 1000)
						c += 500;
					else if (c >= 100)
						c += 50;
					else if (c >= 10)
						c += 5;
					else
						c += 1;
				}
				
				if (c1 >= 10000)
					c1 += 5000;
				else if (c1 >= 1000)
					c1 += 500;
				else if (c1 >= 100)
					c1 += 50;
				else if (c1 >= 10)
					c1 += 5;
				else
					c1 += 1;
			}
		}*/
		
		if (!bestK.hasMin() )
			return setEps(0);
		
		k = bestK.getMin();	
		
		c1 = bestC1.getMin();
			
		c2 = bestC2.getMin();
		c3 = bestC3.getMin();
			
		c = bestC4.getMin();
			
		
		//c = bestC.getMin();
		//c1 = bestC1.getMin();
		eps0 = eps; //bestEps0.getMin();
		
		eps1 = getEps(1);
		bytes = getBytes(eps0, k, c1, c2, c3, c);
		
		//cout << "bytes " << bytes << endl;
		
		if (bytes == numeric_limits<long double>::max())
			return setEps(0);
		
		unifeps = getIntersectedSlices(k) * 4.0*sqrt( datasize*(eps/getIntersectedSlices(k))*0.25)/datasize;
	}
	
	
	
	
	inline long double getEps() const override{
		
		return eps; // eps > 1 ? 1 : eps;
	}
	
	inline long double getUnifEps() const{
		
		return unifeps > 1 ? 1 : unifeps;
	}
	
	inline const string getName() const override{
	
		return "slicehist";
	}
	
	inline void setBytes(long double size) override{
		
	
		const long maxP = 45;
		
		long L = 0;
		
		const double long X = (1.0L * (1L << 15) );
		
		for (long bit = 1L << maxP; bit > 1; bit >>= 1){
			
			const int oldK = getK();
			const long double oldC = getC();
			
			setEps( X/(L|bit) );
			
			if (bytes <= size)
				L |= bit;
			
			setK( oldK );
			setC( oldC );
		}
		
		if (L > 3){
		
			setEps(X/L);
			
		} else {
		
			setEps(0);
		}
		
		
		
		/*
		long double epsmin = 0L;
		long double epsmax = 1L;
		
		long double prev = 1L;
		
		for (long double x = 1.0L-1.0e-10L; x > 1e-100L; x *= 0.1L){
		
			const int oldK = getK();
			const long double oldC = getC();
		
			setEps(x);
			
			if (bytes == size)
				return;
			
			setK( oldK );
			setC( oldC );
			
			if (bytes > size){
			
				epsmin = x;
				epsmax = prev;
				break;
			}
			
			prev = x;
		}
		
		{
			BinarySearcher s;
		
			s.setOutput(size);
		
			s.inputGreaterThan(epsmin);
			s.inputLessThan(epsmax);
			
			for (int i = 0; i < 100; i++){
		
				const long double x = s.getInput();
			
				const int oldK = getK();
				const long double oldC = getC();
				
				setEps(x);
				
				setK( oldK );
				setC( oldC );
				
				const long double y = bytes;
				
				s.feed(x,y);
			}
			
			setEps(s.getInputLB() );
		}
		
		long double preveps = eps;
		
		for (long j = 0; j < 100; j++){
		
			const int oldK = getK();
			const long double oldC = getC();
		
			const double x = preveps*0.99;
		
			setEps(x);
			
			if (bytes == size)
				return;
			
			setK( oldK );
			setC( oldC );
			
			if (bytes > size){
			
				setEps(preveps);
				break;
			}
			
			preveps = x;
		}
		
		if (bytes == numeric_limits<long double>::max()){
		
			eps = 1;
			bytes = 0;
		}*/
	}
	
	inline long double getBytes() const override{
		
		return bytes;
		
	}
	
	
};

/*
class SliceHistParams : public Params{

public:
	int k = 0;
	int DIMS = 2;
	double eps = 0.01;
	int level = 1;
	
	double eps1 = -1;
	
	double c = 5;

	inline double getEps() const{
	
		return eps > 1 ? 1 : eps;
	}

	inline int getK() const{
	
		return k;
	}

	inline SliceHistParams(const SliceHistParams& pars) : Params(pars){
	
		k = pars.k;
		DIMS = pars.DIMS;
		setEps(pars.eps);
		level = pars.level;
	}
	
	inline double f(double e1, int i) const{
	
		double ret = e1;
		
		for (int ii = 0; ii < (i-1); ii++)
			ret = pow(ret, 1-1.0/DIMS)*pow(DIMS-1, 1.0/DIMS);
		
		return ret;
	}
	
	inline double fprod(double e1) const{
	
		double ret = 1.0;
		
		for (int i = 1; i <= k; i++)
			ret *= f(e1, i);
		
		return ret;
	}

	inline double getEps(int i) const{
		
		return f(eps1, i);
		
		
		// e[k+1] = ( e[k]^(1-1/d) ) (d-1)^(1/d)
		
		
		// e1 * e1^(1-1/d) * (d-1)^(1/d) * (e1^(1-1/d) * (d-1)^(1/d))(1-1/d) * (d-1)^(1/d)
		
		// prod[i] e[i] = (2d)^k e
		
		
		
		
		
		// e0 = e / (2d)^k
		
		
		
		const double a1 = eps; ///pow(2.0*DIMS, k);
		//const double a2 = pow(2.0*DIMS, k)*DIMS*pow( DIMS-1.0, (k-1.0)/DIMS );
		
		const double a2 = pow(2.0*DIMS, k)*pow( DIMS-1.0, (k-1.0)/DIMS );
		
		const double a3 = pow(1.0-1.0/DIMS, i-1.0);
		const double a4 = DIMS*(1.0 - pow(1.0-1.0/DIMS, k) );
		
		const double a5 = pow(DIMS-1.0, (i-1.0)/DIMS);
		
		return pow(a1/a2, a3/a4)*a5;
	}
	
	inline double setEps(double e){
	
		this->eps = e;
		
		{
			BinarySearcher s;
		
			s.setOutput( e / pow(2*DIMS, k) );
			
			s.inputGreaterThan(0);
			s.inputLessThan(1);
			
			for (int i = 0; i < 1000; i++){
			
				const double x = s.getInput();
				const double y = fprod(x);
				
				s.feed(x,y);
			}
			
			eps1 = s.getInputLB();
		}
		
		
		optimizeC();
		
	}
	
	inline void setK(int kk){
	
		setEps(eps);
		k = kk;
	}
	
	inline void optimizeC(){
	
		double cbytes = getBytes();
		
		for (double cc = 2; cc <= 100; cc += 0.5){
		
			double oldc = c;
			
			c = cc;		
			const double ccbytes = getBytes();
			
			if (ccbytes < cbytes)
				cbytes = ccbytes;
			else
				c = oldc;
		}
	}
	
	inline SliceHistParams(int dims, const string& s) : Params(s){
	
		DIMS = dims;
		
		if (get("k") != 0)
			setK( get("k") );
		
		if (get("kb") != 0){
			
			long bytes = get("kb")*1024;
			
			if (k == 0){
			
				setK(1);
				
				setBytes(bytes);
				
				int bestK = k;
				double bestEps = eps;
			
				for (int i = 2; i < 4; i++){
				
					setK(i);
				
					setBytes(get("kb")*1024);
					
					if (eps < bestEps){
						bestK = i;
						bestEps = eps;
					}	
				}
				
				setK(bestK);
			}
			
			setBytes(bytes);
			
		} else if (get("eps") != 0){
		
			setEps(get("eps")/100.0);
			
			if (getK() == 0){
			
				setK(1);
				
				int bestK = getK();	
				long bestBytes = getBytes();
			
				for (int i = 2; i < 4; i++){
				
					setK(i);
				
					const long bytes = getBytes();
					
					if (bytes < bestBytes){
						
						bestK = i;
						bestBytes = bytes;
					}
				}
				
				setK(bestK);
			}
			
		} else {
		
			assert (false);
		}
		
		
		
		
		double check = 2*DIMS;
		
		for (int i = 1; i <= k; i++)
			check *= getEps(i);
		
		cout << "check " << check << endl;
		
	}
	
	inline bool isLeaf() const{
	
		return level == k;
	}
	
	inline void print(){
	
		cout << DIMS << "d-SliceHist k " << k << " eps " << eps*100 << "% size " << UtilString::bytesToString(getBytes() ) << endl;
	}
	
	inline void setBytes(long bytes){
		
		{
			BinarySearcher s;
		
			s.setOutput(bytes);
		
			s.inputGreaterThan(0);
			s.inputLessThan(1);
			
			for (int i = 0; i < 100; i++){
		
				const double x = s.getInput();
			
				setEps(x);
			
				const double y = getBytes();
			
				s.feed(x,y);
			}
		
			setEps(s.getInputLB() );
		}
		
		if (eps >= 1){
		 
			BinarySearcher s;
		
			s.setOutput(bytes);
		
			s.inputGreaterThan(0.99999);
			s.inputLessThan(10000000);
			
			for (int i = 0; i < 100; i++){
		
				const double x = s.getInput();
			
				setEps(x);
				
				const double y = getBytes();
			
				s.feed(x,y);
			}
			
			setEps( s.getInputLB() );
		}
	}

	inline long getSlices(int i) const{
	
		const double e = getEps(i)/(1+2/c);///(2.0*DIMS);
		
		return ceil(1.0/e);
	}
	
	inline long getSlices() const{
	
		return getSlices(level);
	}
	
	inline long getGridCells(int i) const{
	
		return pow(getSlices(i), DIMS);
	}
	
	inline double getRankEps(int i) const{
	
		return getEps(i)/c;
	}
	
	inline long getRankValues(int i) const{
	
		return DIMS*(DIMS+1)*ceil( 3.0/getRankEps(i) );
	}
	
	inline double getRankEps() const{
	
		return getRankEps(level);
	}
	
	inline long getBytes() const{
		
		long num = 1;
		
		long values = 0;
		
		for (int i = 1; i <= k; i++){
			
			values += num * ( getGridCells(i) + getRankValues(i) );
			num *= DIMS * getSlices(i);
		}
		
		return values*8;
	}
	
	inline long getMegabytes() const {
	
		return getBytes()/1000000.0;
	}
};

*/

class SliceHistParams : public Params{

private:
	int level = 1;

public:
	
	
	TheoreticalSliceHist model;

	inline double getEps() const{
	
		return model.getEps();
	}


	inline double getC() const{
	
		return model.getC();
	}

	inline int getK() const{
	
		return model.getK();
	}
	
	inline int getLevel() const{
	
		return level;
	}
	
	inline int incLevel(){
	
		this->level++;
	}

	inline int getDims() const{
		
		return model.getDims();
	}
	
	inline long getFastParam() const{
	
		return get("f");
	}

	inline SliceHistParams(const SliceHistParams& pars) : model(pars.model), Params(pars){
		
		level = pars.level;
	}
	
	inline SliceHistParams(int dims, const DataSet& data, const string& s) : model(dims), Params(s){
	
		model.setDataSize(data.getSize() );
	
		if (get("k") != 0)
			model.setK( get("k") );
		
		if (get("c") != 0)
			model.setC( get("c") );
		
		if (get("kb") != 0){
			
			long bytes = get("kb")*1024;
			
			model.setBytes(bytes);
			
		} else if (get("eps") != 0){
		
			model.setEps(get("eps")/100.0);
			
		} else {
		
			assert (false);
		}
	}
	
	inline bool isLeaf() const{
	
		return getLevel() == getK();
	}
	
	inline void print(){
	
		cout << model.getDims() << "d-SliceHist k " << model.getK() << " eps " << model.getEps()*100 << "% size " << UtilString::bytesToString(model.getBytes() ) << endl;
	}
	

	
	inline long getSlices() const{
	
		return model.getSlices(level);
	}
	
	inline double getRankEps() const{
	
		return model.getRankEps(level);
	}
	
	inline long getBytes() const{
		
		return model.getBytes();
	}
	
	inline long getMegabytes() const {
	
		return getBytes()/1000000.0;
	}
};

class SliceHist1 : public DataSummary{

public:
	
	const int DIMS;
	
	RankTransform rt;
	DenseGrid g;
	
	unique_ptr<DenseGrid> ref; // temporary object related to verification of samples
	
	
	int state = -1; // 0 = create rank transform, 1 = count in grid, 2 = finished
	
	Point temp;
	
	long quantileInserts = 0;
	long quantileLookups = 0;
	
	vector<Point> mins;
	vector<Point> maxs;
	
	inline double getBytes() const{
	
		return g.getBytes()+rt.getBytes();
	}

	inline SliceHist1(int dims, long slices, double rankEps, bool isRoot=true, long size=-1, long fastBudget=0) : 
	
		DIMS(dims),
		rt(DIMS ),
		g(DIMS, slices, false),
		temp(DIMS)
	{
	
		RankSummaryParameters params;
		
		params.eps = rankEps;
		params.datasize = size;
		
		params.merge = 10;
		params.naive = 10;
		params.gk = 5;
		params.equi = 100;
		params.sketch = 0;
		
		params.keepEqui = false;
		
		params.oneDim = !isRoot;
		
		rt.init(params);
	}
	
	inline void drawRandomPoint(RandomGenerator& rg, Point& p) const{
	
		assert (DIMS == 2);
		
		const int dim = rg.randomInteger(0, DIMS-1);
		
		const long slice = rg.randomInteger(0, g.getCellsPerDim()-1);
		
		const int dim2 = 1-dim;
		
		long a = 0;
		long b = 0;
		
		for (int i = 0; i < DIMS; i++){
		
			a = g.lc.setDimCoord(a, i, 0);
			b = g.lc.setDimCoord(b, i, g.getCellsPerDim()-1);
		}
		
		a = g.lc.setDimCoord(a, dim, slice);
		b = g.lc.setDimCoord(b, dim, slice);
		
		long sum = g.sum(a,b);
		
		const double target = sum * rg.randomDouble();
		
		long s = 0;
		long t = g.getCellsPerDim()-1;
		
		while (s < t){
		
			long k = (s+t)/2;
			
			if ( UtilMath::isIn(k, s, t) ){
			
				s = k;
				t = k;
				break;
			}
			
			const long sum1 = g.sum(a, g.lc.setDimCoord(b, dim2, k));
			
			if (sum1 >= target)
				t = k;
			
			if (sum1 <= target)
				s = k;
		}
		
		long cell = g.lc.setDimCoord(a, dim2, s);
		
		for (long s2 = s; ( s2 >= 0 ) && (g.sum(cell,cell) == 0 ); s2--)
			cell = g.lc.setDimCoord(cell, dim2, s2);
		
		for (long s2 = g.lc.getDimCoord(cell, dim2); (s2 < g.getCellsPerDim() ) && (g.sum(cell,cell) == 0); s2++)
			cell = g.lc.setDimCoord(cell, dim2, s2);
		
		p[dim] = (g.lc.getDimCoord(cell, dim)+rg.randomDouble())/g.getCellsPerDim();
		p[dim2] = (g.lc.getDimCoord(cell, dim2)+rg.randomDouble())/g.getCellsPerDim();
			
		rt.normRankTransformInv(p, p);
	}
	
	inline long getGridCoord(const Point& p) const{
		
		long ret = 0;
		
		for (int i = 0; i < DIMS; i++)
			ret = g.lc.setDimCoord(ret, i, UtilHist::discretize(p[i], g.getCellsPerDim()) );
		
		return ret;
	}
	
	
	
	inline long getGridCoord(const Point& p, QueryMode mode) const{
		
		assert (false);
		
		long ret = 0;
		
		long skipDims = 0;
		
		for (int i = 0; i < DIMS; i++){
			ret = g.lc.setDimCoord(ret, i, UtilHist::discretize(p[i], g.getCellsPerDim() ) );
			
			if (mode == QueryMode::LB && (g.lc.getDimCoord(ret, i) == 0) )
				skipDims |= 1L << i;
			
			
			if (mode == QueryMode::UB && (g.lc.getDimCoord(ret, i) >= (g.ind.getLength(i)-1)) )
				skipDims |= 1L << i;
		}
		
		if (mode == QueryMode::LB)
			return g.lc.decCoords(ret, skipDims);
		
		if (mode == QueryMode::UB)
			return g.lc.incCoords(ret, skipDims);
		
		return ret;
	}


	inline void add(const Point& p){
	
		assert (state >= 0);
		assert (state <= 1);
	
		if (state == 0){
		
			rt.add(p);
			
			
			
			quantileInserts += DIMS;

		} else if (state == 1){
		
			if (mins.size() == 0){
				
				for (int i = 0; i < DIMS; i++)
					mins.push_back(p);
					
				for (int i = 0; i < DIMS; i++)
					maxs.push_back(p);
			}
			
			for (int i = 0; i < DIMS; i++){
				
				if ( rt.less(i, p, mins[i]) )
					mins[i] = p;
					
				if ( rt.less(i, maxs[i], p) )
					maxs[i] = p;
			}
		
			rt.normRankTransform(p, temp);
			
			quantileLookups += DIMS;
			
			g.addCount( getGridCoord(temp), 1);
		}
	}
	
	inline void distInit(){
	
		if (ref){
		
			//
			
		} else {
	
			ref = std::move( unique_ptr<DenseGrid>( new DenseGrid(DIMS, g.getCellsPerDim(), false) ) );
			
		}
		
		ref->allocate();
	}
	
	inline void distAdd(const Point& p){
		
		assert (ref);
		
		ref->addCount( getGridCoord(p), 1);
	}
	
	inline void distEnd(){
	
		ref->finalise();
	}
	
	inline double distMax(double d1, double d2){
		
		return g.distMax( (*ref), d1, d2);
	}
	
	inline double getCellMin( long c, int dim) const{
	
		return UtilHist::bucketMin( g.lc.getDimCoord(c, dim), g.getCellsPerDim() );
	}
	
	inline double getCellMax( long c, int dim) const {
	
		return UtilHist::bucketMax( g.lc.getDimCoord(c, dim), g.getCellsPerDim() );
	}
	
	inline double tboxCount(const QueryBox& box) const{
		
		assert (false);
		
		assert (state == 2);
		
		return g.getCount(box);
	}
		
	
	inline double tboxCount(const Point& bmin, const Point& bmax, QueryMode mode=QueryMode::EST, int amin = 0, int amax = 0, int emin = 0, int emax = 0) const{
		
		assert (false);
		
		assert (state == 2);
		
		return g.getCount(bmin, bmax, mode, -1, -1, amin, amax, emin, emax);
	}
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override{
	
		assert (false);
	
		assert (state == 2);
	
		Box	b(qmin, qmax);
		
		Box tempBox(qmin, qmax);
		
		rt.normRankTransform(b, tempBox);
		
		return tboxCount(tempBox.min, tempBox.max, mode);
	}
	
	
	inline void setState(int newState){
	
		assert (newState == state+1);
		
		state = newState;
		
		assert(state >= 0);
		assert(state <= 2);
		
		if (state == 1){
		
			rt.finalise();
			g.allocate();
			
		} else if (state == 2){
			
			g.finalise();
		}
	}

	inline const string toString() const override{
		
		return "rankhist";
	}
};


class SliceHistSummary : public DataSummary{

public:

	const int DIMS;
	SliceHistParams params;

	SliceHist1 root;
	Point temp;
	
	long quantileInserts = 0;
	long quantileLookups = 0;
	
	Box boundingBox;
	
	
	
	vector<vector<unique_ptr<SliceHistSummary>>> children;
	int state = -1; // index of datascan
	
	const long datasize;
	
	inline SliceHistSummary(long siz, const SliceHistParams& pars) : datasize(siz), DIMS(pars.getDims() ), boundingBox(pars.getDims()), params(pars), root(DIMS, params.getSlices(), params.getRankEps(), params.getLevel() == 1, siz, params.getFastParam() ), temp(DIMS) {
		
		
	}
	
	inline long getQuantileInserts() const{
	
		long ret = root.quantileInserts+quantileInserts;
		
		for (auto it1 = children.begin(); it1 != children.end(); it1++)
				for (auto it = it1->begin(); it != it1->end(); it++)
					ret += (*it)->getQuantileInserts();
		
		return ret;
	}
	
	inline long getQuantileLookups() const{
	
		long ret = root.quantileLookups+quantileLookups;
		
		for (auto it1 = children.begin(); it1 != children.end(); it1++)
				for (auto it = it1->begin(); it != it1->end(); it++)
					ret += (*it)->getQuantileLookups();
		
		return ret;
	}
	
	inline double getBytes() const override{
	
		double ret = root.getBytes();
		
		if (!params.isLeaf() ){
		
			for (auto it1 = children.begin(); it1 != children.end(); it1++)
				for (auto it = it1->begin(); it != it1->end(); it++)
					ret += (*it)->getBytes();		
		}
					
		return ret;
	}
	
	inline void drawRandomPoint(RandomGenerator& rg, Point& p) const override{
	
		if (params.isLeaf() ){
		
			root.drawRandomPoint(rg, p);
			
		} else {
		
			rg.randomElement(rg.randomElement(children))->drawRandomPoint(rg, p);
			root.rt.normRankTransformInv(p, p);
		}
		
		
	}
	
	
	
	inline SliceHistSummary(const DataSet& data, const string& s) : datasize(data.getSize() ), DIMS(data.getDims() ), params(DIMS, data, s), root(DIMS, params.getSlices(), params.getRankEps(), params.getLevel() == 1, data.getSize(), params.getFastParam() ), temp(DIMS) {
		
		params.print();
		
		observer["eps"] = params.getEps();
		
		
		if (params.model.getBytes() >= pow(1000.0, 5))
			throw 0;
		
		//observer["eps1"] = params.model.getEps1();
		observer["modelbytes"] = params.model.getBytes();
		observer["k"] = params.model.getK();
		observer["c"] = params.model.getC();
		
		cout << "constructing " << params.model.toString() << endl;
		
		Timer t("slicehist");
		
		t.start(10, 1000000, data.getSize()*(params.getK()+1) );
		
		// (k+1) data scans
		for (int j = 0; j <= params.getK(); j++){
		
			cout << "state " << j << endl;
		
			setState(j);
			
			const string key = "datascan"+to_string(j+1);
			
			observer.begin(key);
			
			for (auto it = data.begin(); it != data.end(); it++){
				
				const Point& p = (*it);
				
				/*
				for (int i = 0; i < DIMS; i++){
				
					const double d = p[i];
					
					if (d > 1)
						cout << "1 < d == " << d << endl;
					
					if (d < 0)
						cout << "0 > d == " << d << endl;
				}*/
				
				add(p);
				
				
				t.tick();
			}
			
			observer.end(key);
		}
		
		
		cout << "state " << (1+params.getK()) << endl;
		setState(1+params.getK());
		
		t.end();
		
		observer["bytes"] = getBytes();
		
		observer["quantileInserts"] = getQuantileInserts();
		observer["quantileLookups"] = getQuantileLookups();
		
	}
	
	inline SliceHistSummary(const shared_ptr<DataSet>& data, const string& s) : SliceHistSummary( (*data), s) {
		
		
	}
	
	inline SliceHistSummary(const unique_ptr<DataSet>& data, const string& s) : SliceHistSummary( (*data), s) {
		
		
	}
	
	inline long getChild(long c, int dim) const{
		
		return root.g.lc.getDimCoord(c, dim);
	}
	
	// input rank-transformed p
	inline long getChild(const Point& r, int dim) const{
		
		return getChild(root.getGridCoord(r), dim);
	}
	
	inline void setState(int s){
	
		assert (s == state+1);
		assert (s <= (params.getK()-params.getLevel()+2) );
		
		state = s;
		
		if (s <= 2)
			root.setState(s);
			
		if (params.isLeaf())
			return;
		
		if (s == 0)
			return;
			
		if (s == 1){
			
			SliceHistParams childParams(params);
			
			childParams.incLevel();
			
			for (int i = 0; i < DIMS; i++){
	
				vector<unique_ptr<SliceHistSummary>> v;
	
				for (long c = root.g.getCellsPerDim(); c >= 1 ; c--){
					
					unique_ptr<SliceHistSummary> sh(new SliceHistSummary(0.5*datasize*params.model.getEps( params.getLevel() ), childParams));					
					v.push_back( std::move(sh) );
				}
				
				children.push_back( std::move(v) );
			}
		}
		
		for (auto it1 = children.begin(); it1 != children.end(); it1++)
			for (auto it = it1->begin(); it != it1->end(); it++)
				(*it)->setState(s-1);
	}
	
	inline void add(const Point& p){
	
		assert (state >= 0);
		
		if (root.state == 0)
			boundingBox.enclose(p);
		
		if (root.state < 2)		
			root.add(p);
		
		if (params.isLeaf())
			return;
		
		if (state == 0)
			return;
		
		root.rt.normRankTransform(p, temp);
		
		quantileLookups += DIMS;
		
		for (int i = 0; i < DIMS; i++)
			children[i][getChild(temp, i)]->add(temp);
	}
	
	
	inline void distInit(){
	
		root.distInit();
	
		if (!params.isLeaf() ){
		
			for (auto it1 = children.begin(); it1 != children.end(); it1++)
				for (auto it = it1->begin(); it != it1->end(); it++)
					(*it)->distInit();		
		}
	}
	
	inline void distAdd(const Point& p){
	
		root.rt.normRankTransform(p, temp);
		
		root.distAdd(temp);
		
		quantileLookups += DIMS;
		
		if (!params.isLeaf() ){
			
			for (int i = 0; i < DIMS; i++)
				children[i][getChild(temp, i)]->distAdd(temp);
		}
	}
	
	inline void distEnd(){
	
		root.distEnd();
	
		if (!params.isLeaf() ){
		
			for (auto it1 = children.begin(); it1 != children.end(); it1++)
				for (auto it = it1->begin(); it != it1->end(); it++)
					(*it)->distEnd();		
		}
	}
	
	inline double distMax(double w1, double w2){
	
		double ret = root.distMax(w1,w2);
		
		if (!params.isLeaf() ){
		
			for (auto it1 = children.begin(); it1 != children.end(); it1++)
				for (auto it = it1->begin(); it != it1->end(); it++){
				
					const double d = (*it)->distMax(w1, w2);
					
					if (d > ret)
						ret = d;
				}
		}
		
		return ret;
	}
	
	inline double getEpsilonDistance(const DataSet& d){
		
		TheoreticalSliceHist tsh(params.model);
		
		tsh.setSupportedQueries(1);
		
		distInit();
			
		for (auto it = d.begin(); it != d.end(); it++)
			distAdd( (*it) );
		
		distEnd();
		
		const double n1 = datasize;
		const double n2 = d.getSize();
		
		const double nmax = UtilMath::maxVal<double>(n1,n2);
		
		const double w1 = nmax/n1;
		const double w2 = nmax/n2;
		
		const double dm = distMax(w1, w2);
		
		cout << "n1 " << n1 << endl;
		cout << "n2 " << n2 << endl;
		
		cout << "sheps " << tsh.getEps() << endl;
		cout << "nmax " << nmax << endl;
		cout << "distMax " << dm << endl;
		cout << "decomp " << tsh.getDecompositionCardinality() << endl;
		
		const double ret = pow(2.0, DIMS) * ( tsh.getEps()*nmax + tsh.getDecompositionCardinality() * dm );
		
		//data.distClear();
		
		return UtilMath::makeBetween<double>(0, nmax, ret);
	}
	
	inline bool contains(const QueryBox& q, const Point& p) const{
	
		for (int i = 0; i < DIMS; i++){
		
			if (!q.isMinUnbounded(i)){
			
				if ( root.rt.less(i, p, q.min ) )
					return false;
			}
				
			if (!q.isMaxUnbounded(i)){
			
				if ( root.rt.less(i, q.max,p ) )
					return false;
			}
		}
		
		return true;
	}
	
	inline double childCount(QueryBox& tempBox, int i, long c, bool alignLeft=false, bool alignRight=false) const{
	
		const long flags = tempBox.getFlags();
		
		tempBox.setUnbounded(i, ( alignLeft ||tempBox.isMinUnbounded(i) ) , (alignRight || tempBox.isMaxUnbounded(i) ) );
		
		const double ret = children.at(i).at(c)->boxCount(tempBox);
		
		tempBox.setFlags(flags);
		
		return ret;
	}
	
	inline double boxCount(const QueryBox& box) const {
	
		assert (state == (params.getK()-params.getLevel()+2) );
		
		QueryMode mode = box.mode;
		
		QueryBox tempBox(box);
		
		root.rt.normRankTransform(box, tempBox);
		
		/////////////////
		
		const std::tuple<long, long, long, long> gridbox = root.g.getGridBox(tempBox);
		
		const long ubmin = std::get<0>(gridbox);
		const long ubmax = std::get<1>(gridbox);
		
		const long UBAREA = root.g.lc.getArea(ubmin, ubmax);
		
		if (UBAREA <= 0)
			return 0;
		
		const long lbmin = std::get<2>(gridbox);
		const long lbmax = std::get<3>(gridbox);
		
		const long LBAREA = root.g.lc.getArea(lbmin, lbmax);
				
		if (params.isLeaf() ){
			
			const double LB = LBAREA > 0 ? root.g.sum(lbmin, lbmax) : 0;
			const double UB = UBAREA > 0 ? root.g.sum(ubmin, ubmax) : 0;
			
			if (LB == UB){
				
				return LB;
			
			} else if (mode == QueryMode::LB){
			
				return LB;
				
			} else if (mode == QueryMode::UB){
			
				return UB;
				
			} else {
				
				double vol = 1.0;
			
				for (int i = 0; i < DIMS; i++){
				
					const long LEN = root.g.ind.getLength(i);
				
					vol *= (tempBox.max[i]-tempBox.min[i])*(LEN-1);
				}
			
				vol -= root.g.lc.getArea(lbmin, lbmax);
			
				vol /= root.g.lc.getArea(ubmin, ubmax)-root.g.lc.getArea(lbmin, lbmax);
				
				vol = UtilMath::makeBetween<double>(0.0, 1.0, vol);
				
				return UtilMath::makeBetween<double>(LB, UB, LB+(UB-LB)*vol );
			}
		}
		
		// CONTAINED CELLS
		
		// [  |  ]
		// [  |  ]
		// [  |  ]
		// [u1|u2]
		// [  |  ]
		// [  |  ]
		// [  |  ]
		
		// SLICE-ALIGNED QUERIES ARE PASSED ON TO CHILDREN COMPLETELY
		for (int i = 0; i < DIMS; i++){
		
			const long u1 = root.g.lc.getDimCoord(ubmin, i);
			const long u2 = root.g.lc.getDimCoord(ubmax, i);
		
			if (  (u2-u1) <= 1){
				
				assert (LBAREA <= 0);
				
				double count = 0;
				
				if (u1 == u2){
				
					// SINGLE SLICE
					count += childCount(tempBox, i, u1, false, false  );
					
				} else	{
				
					// DOUBLE SLICE
					count += childCount(tempBox, i, u1, false, true );
					count += childCount(tempBox, i, u2, true, false  );
				}
				
				return count;
			}
			
		}
		
		assert (LBAREA > 0);
		
		double count = root.g.sum(lbmin, lbmax);
		
		
		// [  |  u1  |  ]
		// [  |------|  ]
		// [  |      |  ]
		// [u1|l1..l2|u2]
		// [  |      |  ]
		// [  |------|  ]
		// [  |  u2  |  ]
		
		for (int i = 0; i < DIMS; i++){
			
			
			
			const long l1 = root.g.lc.getDimCoord(lbmin, i);
			const long l2 = root.g.lc.getDimCoord(lbmax, i);
			
			const long u1 = root.g.lc.getDimCoord(ubmin, i);
			const long u2 = root.g.lc.getDimCoord(ubmax, i);
			
			assert ( (u2-u1) > 1);
			
			// [u1|l1..
			if ( (u1+1) == l1 ){
				
				count += childCount(tempBox, i, u1, false, true);
				tempBox.setMin(i, UtilHist::bucketMin(l1, root.g.getCellsPerDim() ));
				
				
			} else {
			
				// [l1...]
			}
			
			// ..l2|u2]
			if ( (l2+1) == u2 ){
			
				count += childCount(tempBox, i, u2, true, false);
				tempBox.setMax(i, UtilHist::bucketMax(l2, root.g.getCellsPerDim() ));
				
			} else {
				
				// ..l2]
			}
			
		}
		
		return count;
	}
	
	inline void debug() override{
	
		/*
		root.rt.debug();
		
		if (!params.isLeaf() ){
		
			for (auto it1 = children.begin(); it1 != children.end(); it1++){
				for (auto it = it1->begin(); it != it1->end(); it++){
				
					(*it)->debug();
				}
			}
		}*/
	}
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override{
	
	
		QueryBox box(qmin, qmax, mode);
		
		return boxCount(box);
	}
	
	inline const string toString() const override{
		
		return params.getK() == 1 ? "rankhist "+params.toString() : "slicehist "+params.toString();
	}
};


#endif
