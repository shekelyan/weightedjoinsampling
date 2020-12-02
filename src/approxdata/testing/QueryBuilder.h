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

#ifndef SHEKELYAN_QUERYBUILDER_H
#define SHEKELYAN_QUERYBUILDER_H

#include <approxdata/testing/testing.h>

enum QueryDistribution{

	BIASED, UNBIASED, EQUIDEPTH

};

enum QueryAggregation {

	NONE, MINIMUM, AVERAGE, MEDIAN, PERCENTILE_90, PERCENTILE_95, MAXIMUM
};


class QueryBuilder {


public:
	Sorter sorter;
	
	bool verbose = true;
	
	double selmin = 0.5;
	double selmax = 1.5;
	
	QueryDistribution distribution = BIASED;
	
	long smallSize = 100*pow(10,3);
	long mediumSize = pow(10,6);
	long bigSize = 100*pow(10,6);
	
	shared_ptr<const string> seed = std::make_shared<const string>("abc");
	
	inline void setSelectivityRange(double selmin, double selmax){
	
		this->selmin = selmin;
		this->selmax = selmax;
	}
	
	inline void setQueryDistribution(QueryDistribution d){
	
		distribution = d;
	}
	
	inline void calcCounts(RangeQueries& queries, const DataSet& data){
		
		const double mode0ops = data.getSize()*queries.size();
		const double mode1ops = (queries.size()*queries.size())+data.getSize()*(0.01*queries.size());
		
		queries.discretize();
		
		RangeQueries q;
		
		q.addFromString(queries.toString() );
		
		BoxIndex bi(data.getDims(), q.ranges, ((mode0ops < mode1ops) || (q.size() < 10)) ? 0 : 1);	
		
		for (long j = 0; j < q.ranges.size(); j++)
			assert (q.ranges[j].count == 0);
		
		for (auto it = data.begin(); it != data.end(); it++ )
			bi.count( (*it) );
			
		for (long j = 0; j < queries.size(); j++)
			queries.setCount(j, q[j].count);
	}
	
	inline void scaleBoxToMatchCount(Box& b, long targetedCount, const DataSet& data){
	
		if (data.getSize() <= targetedCount){
		
			b.makeEmpty();
			
			for (auto it = data.begin(); it != data.end(); it++)
				b.enclose( (*it) );
				
			return;
		}
		
		sorter.clear();
		
		{
			long pos = 0;
			for (auto it = data.begin(); it != data.end(); it++, pos++)
				sorter.add( pos, b.volumeDistance( (*it) ));
		}
		
		b.makeEmpty();
		
		for (long j = 0; j < targetedCount; j++){
		
			auto it = data.begin();
			
			it += sorter.getSortedIndex(j);
			
			const Point p = (*it);
			
			b.enclose(p);
		}
	} 
	
	
	inline shared_ptr<RangeQueries> createEquidepthQueries(const DataSet& data, long target){
	
		shared_ptr<RangeQueries> ret(new RangeQueries() );
		
		const double sel = selmin*0.5+selmax*0.5;
		
		shared_ptr< RangeQueries > cand(new RangeQueries() );
			
		EquidepthSummary h(data, "-nc 1 -rev 1 -sel "+to_string(sel));
		
		if (verbose)
			cout << "ed created" << endl;

		h.addBoxes( cand );
		
		if (data.getSize() > 5000000){

			shared_ptr<DataSet> sample = data.getRandomSubset(1000000);
			calcCounts( (*cand), (*sample) );
			
			cand = cand->getWithCounts( sample->getCount(selmin), sample->getCount(selmax));
		}
		
		assert (cand->size() >= target);
		
		cand = cand->getRandomSubset(target*2);
		
		calcCounts( (*cand), data);

		cand = cand->getWithCounts( data.getCount(selmin), data.getCount(selmax));

		assert (cand->size() >= target);

		cand = cand->getRandomSubset(target);
		
		if (verbose)
			cand->print();
	
		ret->add( cand );
		
		return ret;
	}
	
	inline shared_ptr<RangeQueries> createQueries(const DataSet& data, long target){
		
		const string seed0 = (*seed);
		
		if (distribution == QueryDistribution::EQUIDEPTH)
			return createEquidepthQueries(data, target);
		
		const string distributionName = distribution == QueryDistribution::EQUIDEPTH ? "equidepth" : distribution == QueryDistribution::BIASED ? "biased" : distribution == QueryDistribution::UNBIASED ? "unbiased" : "unknown";
		
		
		
		//data.print();
		
		const double sel = selmin*0.5+selmax*0.5;
		
		assert (UtilMath::isIn(distribution, QueryDistribution::BIASED, QueryDistribution::UNBIASED));
		
		const bool isBiased = distribution == QueryDistribution::BIASED;
		
		// sel/100.0 * mediumSize >= 2 follows mediumSize >= 200.0 / sel
		
		const long smallSize = this->smallSize;
		const long mediumSize = UtilMath::maxVal<long>(this->mediumSize, ceil(200.0/sel) );
		const long bigSize = UtilMath::maxVal<long>(mediumSize, this->bigSize);
		
		const long tt = sel/100.0 * smallSize;
		const long tt2 = sel/100.0 * mediumSize;
		
		cout << "tt " << tt << " tt2 " << tt2 << endl;
		
		assert (tt2 >= 2);
		
		const int DIMS = data.getDims();
		
		if (verbose)
			cout << "create big sample of size " << bigSize << endl;
		
		shared_ptr<DataSet> big = data.getRandomSubset(UtilMath::minVal<long>(data.getSize(), bigSize), seed0);
	
		const double sel1 = smallSize*100.0 / big->getSize();
		
		
				
		if (verbose)
			cout << "create medium sample of size " << mediumSize << endl;
			
		shared_ptr<DataSet> medium = big->getRandomSubset(UtilMath::minVal<long>(big->getSize(), mediumSize), seed0);
		
		if (verbose)
			cout << "create small sample of size " << smallSize << endl;
			
		
		shared_ptr<DataSet> small = medium->getRandomSubset(UtilMath::minVal<long>(medium->getSize(), smallSize), seed0);
		
		
		
		Sorter sorter;
		
		shared_ptr< RangeQueries > fullCands(new RangeQueries() );
		
		for (int k = 0; k < 10; k++){
			
			const string seed1 = seed0+" "+to_string(k);
			
			shared_ptr<DataSet> centers = isBiased ? big->getRandomSubset(2*target, seed1+" 1") : big->getRandomPoints(2*target, seed0+" 1");
		
			shared_ptr<DataSet> unbiased1 = big->getRandomPoints(2*target, seed1+" 3");
			shared_ptr<DataSet> unbiased2 = big->getRandomPoints(2*target, seed1+" 4");
			
			shared_ptr< RangeQueries > bigCands(new RangeQueries() );
			
			auto it1 = centers->begin();
			auto it2 = unbiased1->begin();
			auto it3 = unbiased2->begin();
			
			for (long j = 0; j < (2*target); j++, it1++, it2++, it3++){
			
				// Create box spanned by two points drawn where each coordinate is drawn from U(0,1)
				const Point& p1 = (*it2);
				const Point& p2 = (*it3);
			
				Box b;
			
				b.enclose(p1);
				b.enclose(p2);
				
				cout << "b1 " << b << endl;
				
				// move box to a center that is either drawn from same distribution, or dataset
			
				const Point& center = (*it1);
				
				//cout << "center " << center << endl;
				
				for (int i = 0; i < DIMS; i++){
				
					double delta = center[i]-(b.max[i]*0.5+b.min[i]*0.5);
					
					if (b.min[i]+delta < 0)
						delta = -b.min[i];
					else if (b.max[i]+delta > 1)
						delta = 1.0-b.max[i];
					
					b.min[i] = UtilMath::truncate(0, 1, b.min[i]+delta);
					b.max[i] = UtilMath::truncate(0, 1, b.max[i]+delta);
				}
				
				//cout << "b2 " << b << endl;
				
				// scale box to desired count
				
				if (tt >= 10){
					
					scaleBoxToMatchCount(b, tt, (*small) );
					
					//cout << "b3 " << b << endl;
					
					bigCands->add(b);
			
				} else {
					
					scaleBoxToMatchCount(b, 10, (*small) );
					
					//cout << "b3 " << b << endl;
					
					shared_ptr<vector<Point>> v(new vector<Point>() );
					
					v->reserve(1.5*ceil(10.0/smallSize*mediumSize));
					
					for (auto it = medium->begin(); it != medium->end(); it++)
						if (b.contains( (*it) ))
							v->push_back(*it);
					
					shared_ptr<DataSet> selection = medium->getCopy(v);
					
					//scaleBoxToMatchCount(b, tt2, (*small) );
					
					scaleBoxToMatchCount(b, tt2, (*selection) );
					
					//cout << "b4 " << b << endl;
					
					b.min.id = std::numeric_limits<long>::lowest();
					b.max.id = std::numeric_limits<long>::max();
					
					bigCands->add( b);	
				}
			}
			
			if (verbose)
				cout << "created " << bigCands->size() << " candidates" << endl;
			
			
			calcCounts( (*bigCands), (*big) );
			
			if (verbose)
				cout << "calculated counts of " << bigCands->size() << " candidates" << endl;
			
			
			
			fullCands->add( bigCands->getWithCounts(big->getCount(selmin), big->getCount(selmax) ) );
			
			if (verbose)
				cout << "now at " << fullCands->size() << " strong candidates" << endl;
				
			if (fullCands->size() == 0){
			
				bigCands->print();
			
			}
			
			if (fullCands->size() >= target*2)
				break;
		}
		
		cout<< "start creating subset " << endl;
		fullCands = fullCands->getRandomSubset(target*2);
		
		
		cout<< "stop creating subset" << endl;
		
		
		cout<< "calc counts " << endl;
		
		calcCounts( (*fullCands), data);
		
		cout<< "stop calc counts " << endl;
		
		fullCands = fullCands->getWithCounts(data.getCount(selmin), data.getCount(selmax) );
		fullCands = fullCands->getRandomSubset(target);
		
		
		shared_ptr<RangeQueries> ret(new RangeQueries() );
		
		ret->add( fullCands );	
		
		ret->setName(data.getFileName()+"_"+distributionName+"_num"+UtilString::integerToShortString(target)+"_sel"+UtilString::doubleToString(sel)+"pct");
		
		return ret;
	}
	
	

};



inline shared_ptr<RangeQueries> newWorkload(const DataSet& data, const string s){

	const string prefix = "file ";

	if (UtilString::startsWith(s, prefix)){
	
		cout << "load file workload " << s.substr(prefix.length() ) << endl;
	
		shared_ptr<RangeQueries> ret(new RangeQueries);
		
		ret->add(s.substr(prefix.length() ));
		
		return ret;
	}
	
	cout << "new synthetic workload " << s << endl;

	Params p( s);

	QueryBuilder b;
	
	assert (p.get("selectivity") > 0 && p.get("selectivity") <= 1);
		
	b.setSelectivityRange( p.get("selectivity")*50.0, p.get("selectivity")*150.0);
	
	
	if (UtilString::equals(p.name, "unbiased"))
		b.setQueryDistribution(QueryDistribution::UNBIASED);
	else if (UtilString::equals(p.name, "biased"))
		b.setQueryDistribution(QueryDistribution::BIASED);
	
		

	const shared_ptr<RangeQueries> ret = b.createQueries(data, p.get("num"));
	
	//ret->setName(s);
	
	return ret;
}


#endif