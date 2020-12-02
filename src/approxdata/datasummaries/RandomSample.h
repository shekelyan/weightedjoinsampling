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

#ifndef SHEKELYAN_RANDOM_SAMPLE_H
#define SHEKELYAN_RANDOM_SAMPLE_H

#include <approxdata/datasummaries/datasummaries.h>


class TheoreticalRandomSample : public TheoreticalSummaryModel{

public:

	long double f = -18;
	
	const int DIMS;
	
	const int BYTES_PER_SAMPLE;
	
	long datasize = 1L << 63;
	
	long bytes = -1;
	
	long double eps = -1;
	
	inline TheoreticalRandomSample(int dims) : DIMS(dims), BYTES_PER_SAMPLE(dims*8){
	
		
	}

	inline void setFailureProbability(long ff){
	
		f = ff;
	}

	inline void setDataSize(long k) override{
	
		this->datasize = k;
		
		
	}
	
	inline void setEps(long double eps) override{
	
	}

	inline int getDims() const override{
	
		return DIMS;
	}
	
	inline long getDataSize() const override{
	
		return datasize;
	}
	
	inline long getSampleSizeForEps(long double eps) const{
	
		// if we expect to test all samples, we expect to 
	
		const long double n = datasize;
		
		
		
		
		
		const long maxP = 10;
		
		long x = (1L << (maxP+1))-1;
		
		for (long bit = 1L << maxP; bit >= 1; bit >>= 1){
		
			const long f = -(x^bit);
			
			const long samplesize = ceil( (1.0L/(eps*eps))*( DIMS*log(n)-0.5L*f*log(10.0L) ) );
			
			// n = ceil( (1.0L/(eps*eps))*( DIMS*log(n)-0.5L*f*log(10.0L) ) );
			
			
			
			
			if (samplesize < 0)
				return -1;
			
			if (samplesize > n){
				x ^= bit;
				continue;
			}
			
			const long double probf = UtilMath::log10choose<long double>(n, samplesize) ;
			
			//cout << "datasize " << datasize <<" eps " << eps << " samplesize " << samplesize << " probf " << probf << " f " << f << endl;
			
			if ( (-f) > probf)
				x ^= bit;
			
		}
		
		const long ret = ceil( (1.0L/(eps*eps))*( DIMS*log(n)-0.5L*x*log(10.0L) ) );
		
		if (ret < 0 || ret > n)
			return n;
		
		return ret;
		
		
		//assert(false);
		//return ceil( (0.5L*DIMS/(eps*eps))*(2.0L*DIMS*log((long double)datasize)-(f*log(10) ) ) );
	}
	
	
	
	inline void setSampleSize(long k){
	
		
		const long double n = datasize;
		bytes = k*BYTES_PER_SAMPLE;
		
		
		//const long double f = k > 1000000 ? -k*(log10(n)-log10(k)) : -UtilMath::log10choose<long double>(n, k) ;
		
		//f = -307.6527;
		
		//k = 1+( (0.5L*DIMS/(eps*eps))*(2.0L*DIMS*log((long double)datasize)-(f*log(10) ) );
		//return ceil( (0.5L*DIMS/(eps*eps))*(2.0L*DIMS*log((long double)datasize)-(f*log(10) ) ) );
		
		// 1/eps^2 *  0.5*( 2*DIMS*log(n)-f*log(10.0L) = k
		
		// exp( 2/eps^2 * DIMS*log(n) - 2k) = pf
		
		eps = sqrt( 0.5*( 2*DIMS*log(n)-f*log(10.0L) )/(k-1) );
		
		if (eps > 1)
			eps = 1;
	}
	
	
	inline void setBytes(long double k) override{
	
		setSampleSize( ceil( k/BYTES_PER_SAMPLE ) );
	}

	inline long double getBytes() const override{
	
		return bytes;
	}
	
	inline const string getName() const override{
	
		return "sampling (cf ="+UtilString::doubleToString(100*(1.0-pow(10.0, -f)))+")";
	}
	
	inline long double getEps() const override{
	
		return eps;
	}
	
	inline long double getMegabytes() const{
	
		return bytes/1000000.0L;
	}
};

class RandomSampleSummary : public DataSummary{
	
	
public: 
	
	const int DIMS;
	
private:

	const Params params;
	shared_ptr<DataSet> sample;
	const long datasize;
	
	double samplePointWeight;
	
	double maxAbsErr = -1;
	
public:
	
	inline double getSamplePointWeight() const{
	
		return samplePointWeight;
	}
	
	inline double getSampleSize() const{
	
		return sample->getSize();
	}
	
	inline double getBytes() const override{
	
		return getSampleSize()*8*DIMS;
	
	}
	
	inline double getDataSize() const{
		
		return datasize;
	}
	
	inline void verify(SliceHistSummary& sh){
	
		long quantileInserts = 0;
		long quantileLookups = 0;
		
		const double epsn = sh.getEpsilonDistance( *sample);
			
		cout << "maxAbsErr " << maxAbsErr << " epsn " << epsn << endl;	
			
		if ( (maxAbsErr < 0) || (epsn < maxAbsErr) ){
			
			maxAbsErr = epsn;
			
			observer["eps"] = maxAbsErr/getDataSize();
			observer["k"]  = sh.observer["k"];
			observer["verifybytes"] = sh.getBytes();
		}
		
		quantileInserts += sh.getQuantileInserts();
		quantileLookups += sh.getQuantileLookups();
		
		if (observer.count("quantileInserts") > 0)
			observer["quantileInserts"] = observer["quantileInserts"]+quantileInserts;	
		else
			observer["quantileInserts"] = quantileInserts;
		
		
		if (observer.count("quantileLookups") > 0)
			observer["quantileLookups"] = observer["quantileLookups"]+quantileLookups;	
		else
			observer["quantileLookups"] = quantileLookups;
		
		
		if (observer.count("constructionhours") > 0)
			observer["constructionhours"] = observer["constructionhours"]+sh.observer["constructionhours"];
		else
			observer["constructionhours"] = sh.observer["constructionhours"];
		
	}
	
	inline RandomSampleSummary(const DataSet& data, string paramStr) : DIMS(data.getDims()), params(paramStr), datasize(data.getSize() ){

		const double kbPerSample = DIMS*8.0/(1000.0);
		const double sampleSize = UtilMath::minVal<long>(datasize, params.get("kb")/kbPerSample );
		
		assert (sampleSize > 0);
		
		sample = data.getRandomSubset(sampleSize, to_string(params.get("seed")) );
		
		sample->setVerbose(false);
		samplePointWeight = 1.0/getSampleSize()*getDataSize();
		
		
		long quantileInserts = 0;
		long quantileLookups = 0;
		
		if (params.get("verifykb") > 0){
			
			for (int k = 1; k <= 2; k++){
			
				if (DIMS == 2 && k == 2)
					continue;
				
				SliceHistSummary sh(sample, "slicehist -k "+to_string(k)+" -kb "+to_string(params.get("verifykb") ));
		
				const double epsn = sh.getEpsilonDistance(data);
				
				if (k == 1 || epsn < maxAbsErr){
					
					maxAbsErr = epsn;
					
					observer["eps"] = maxAbsErr/getDataSize();
					observer["k"]  = k;
					observer["verifybytes"] = sh.getBytes();
					
				}
				
				quantileInserts += sh.getQuantileInserts();
				quantileLookups += sh.getQuantileLookups();
			}
		}
		
		observer["quantileInserts"] = quantileInserts;
		observer["quantileLookups"] = quantileLookups;
	}
	
	inline RandomSampleSummary(shared_ptr<DataSet> data, string paramStr) : RandomSampleSummary( (*data), paramStr){
		
	}
	
	inline RandomSampleSummary(unique_ptr<DataSet> data, string paramStr) : RandomSampleSummary( (*data), paramStr){
		
	}
	
	inline void drawRandomPoint(RandomGenerator& rg, Point& p) const override{
	
		p = (*sample)[rg.randomInteger(0, sample->getSize() )];
	}
	
	~RandomSampleSummary() override {};
	
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode = QueryMode::EST) const override{
		
		if (getDataSize() == 0)
			return 0;
			
		double ret = 0;
		
		if (mode == QueryMode::LB){
		
			if (maxAbsErr >= 0 && maxAbsErr <= getDataSize() )
				ret = sample->boxCount(qmin, qmax)*samplePointWeight-maxAbsErr;
			else 
				ret = sample->boxCount(qmin, qmax);
				
		} else if (mode == QueryMode::UB){
		
			if (maxAbsErr >= 0 && maxAbsErr <= getDataSize() )
				ret = sample->boxCount(qmin, qmax)*samplePointWeight+maxAbsErr;
			else
				ret = datasize-(getSampleSize()-sample->boxCount(qmin, qmax));
				
		} else {
			
			ret = sample->boxCount(qmin, qmax)*samplePointWeight;
		}
		
		return UtilMath::makeBetween<double>(0, getDataSize(), ret);
	}
	
	inline const string toString() const override{
		
		if (params.get("verifykb") > 0)
			return "verifiedsampling "+params.toString();
		
		if (params.get("vkb") > 0)
			return "vsampling "+params.toString();
		
		
		return "randomsampling "+params.toString();
	}
};

#endif