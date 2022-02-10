
#ifndef HEADER_JOINSAMPLING_LIB
#define HEADER_JOINSAMPLING_LIB

//#include <approxdata/utils/utils.h>
//#include <approxdata/external/pq.h>
#include <approxdata/external/csv.h>
#include <iostream>
//#include <libpq-fe.h>
//#include <sql.h>
//#include <sqlext.h>
#include <queue>
#include <map>

#include <set>






#include <approxdata/utils/utils.h>


#include <joinsampling_lib_rows.h>



class SortedTable{
public:

	shared_ptr<Rows> rows;

	vector<Row> v;
	
	const Tuple cols;
	
	
	inline SortedTable(const Tuple& t) : cols(t){
	
		rows = std::make_shared<Rows>(t.size() );
	}
	
	
	template<typename T>
	inline void add(const T& t, vector<int> sel){
	
		long byt = 0;
		
		for (auto it = t.begin(); it != t.end(); ++it)
			byt += (*it).getBytes();
			
		rows->reserve("SortedTable", byt);
		
		if (sel.size() == 1 && sel.at(0) == -1){
		
			for (auto it = t.begin(); it != t.end(); ++it)
				rows->push_back(*it);
			
		} else {
		
			for (auto it = t.begin(); it != t.end(); ++it)
				rows->push_back(*it, sel);
		}
		
		rows->getRows(v);
		
		//for (auto it = rows->begin(); it != rows->end(); ++it)
		//	v.push_back(*it);
			
		//cout << "vsize " << v.size() << endl;
	}
	
	template<typename T>
	inline void add(const T& t, const Tuple& tt){
	
		if (UtilTuple::equals(tt, cols)){
			add(t, vector<int>{-1});
			return;
		}
		
		cout << tt << " != " << cols << endl;
		
		vector<int> v = UtilTuple::find(tt, cols);
		
		for (int j = 0; j < v.size(); ++j)
			assert(v.at(j) != -1);
		
		/*
		for (int j = 0; j < cols.size(); ++j){
		
			const int i = UtilTuple::find(tt, cols[j]);
			
			assert (i != -1);
			//if (i != -1)
			v.push_back(i);
		}*/
		
		assert (UtilTuple::equals(cols, UtilTuple::select(tt, v) ) );
		
		add(t, v);
		
	}
	
	inline void add(const Rows& t){
	
		add(t, cols);
	}
	
	inline void add(const WeightedRowSample& t){
	
		add(t, t.cols);
	}
	
	inline void add(shared_ptr<WeightedRowSample> t){
	
		add(*t);
	}
	
	inline void add(const DatabaseTable& t){
	
		add(t, t.getColNames() );
	}
	
	inline bool reserve(long n){
	
		bool ret = false;
		if (rows){
		
			v.reserve(rows->size() );
		
		} else {
		
			rows = std::make_shared<Rows>(cols.size() );
			rows->reserve("SortedTable", n * 8 * cols.size() );
			ret = true;
			v.reserve(n );
		}
		
		return ret;
	}
	
	
	DatabaseTable* tabptr = NULL;
	
	inline SortedTable(DatabaseTable& tab) : cols(tab.getColNames() ){
	
		tabptr = &tab;
		
		if (tab.getRows())
			rows = tab.getRows();
		
		if (reserve(tab.getSize() )){
		
			for (auto it = tab.begin(); it != tab.end(); ++it)
				rows->push_back(*it);
		}
		
		for (auto it = rows->begin(); it != rows->end(); ++it)
			v.push_back(*it);
			
		//cout << "vsize1 " << v.size() << endl;
	}
	
	inline void sort(const RowComparator& cmp ){
	
		std::sort(v.begin(), v.end(), cmp);
	}
	
	
	inline MainMemoryIndex* getIndex(int k){
	
		if (tabptr)
			return tabptr->getIndex(k);
		
		return NULL;
	}
};



class LazyJoin{

public:
	vector< shared_ptr<SortedTable> > tables;
	
	
	
	vector< vector<Row>::const_iterator > it1; // begin
	vector< vector<Row>::const_iterator > it2; // end of block
	vector< vector<Row>::const_iterator > it3; // total end
	bool reachedEnd = false;

	int top = -1;
	
	
	inline bool joinEqual(const vector< vector<Row>::const_iterator >& iit1, int j1, const vector< vector<Row>::const_iterator >& iit2, int j2) const{
	
		return iit1.at(j1)->equals(projsJoin.at(j1), *(iit2.at(j2)), projsJoin.at(j2) );
	}
	
	inline bool joinLess(const vector< vector<Row>::const_iterator >& iit1, int j1, const vector< vector<Row>::const_iterator >& iit2, int j2) const{
	
		return iit1.at(j1)->less(projsJoin.at(j1), *(iit2.at(j2)), projsJoin.at(j2) );
	}
	
	// check if aligned
	inline bool allEqual(const vector< vector<Row>::const_iterator >& iit){
	
	
	
		for (int i = iit.size()-1; i >= 0; --i)
			if (!joinEqual(iit, i, iit, top))
				return false;
		
		return true;
	}
	
	inline void begin(){
	
		assert (initialised);
	
		reachedEnd = false;
		it1.clear();
		it2.clear();
		it3.clear();
		
		for (int i = 0; i < this->tables.size(); ++i){
			
			it1.push_back( tables.at(i)->v.cbegin() );
			it2.push_back( tables.at(i)->v.cbegin() );
			it3.push_back( tables.at(i)->v.cend() );
			
			if ( it1.at(i) == it3.at(i) )
				reachedEnd = true;
		}
	}
	
	inline int smallest() const{
	
		int ret = top;
		
		for (int i = 0; i < it1.size(); ++i){
		
			if ( joinLess(it1, i, it1, ret))
				ret = i;
		}
		
		return ret;
	}
	
	
	inline bool next(){
	
		if (reachedEnd)
			return false;
		
		long m = 0;
		
		// progress at least one
		it1.at(top) = it2.at(top);
		
		m |= 1L << top;
		
		for (int i = 0; i < it1.size(); ++i){
		
			if (it1.at(i) == it3.at(i))
				reachedEnd = true;
		}
			
		if (reachedEnd)
			return false;
		
		// alignment
		while (!allEqual(it1) ){
		
			const int sm = smallest();
			
			//for (int i = 0; i < it2.size(); ++i)
			//	cout << "1234 " << *(it1.at(i)) << endl;
			
			if (m & (1L << sm) ){
			
				++it1.at(sm);
				++it2.at(sm);
				
			} else {
			
				it1.at(sm) = it2.at(sm); // jump to end of group
				
				m |= (1L << sm);
			}
			
			if (it1.at(sm) == it3.at(sm))
				reachedEnd = true;
			
			if (reachedEnd)
				return false;
		}
		
		//it2 = it1;
		
		// move to end of blocks
		for (int i = 0; i < it2.size(); ++i){
		
			/*if (m & (1L << i))*/{	
				
				while ( (it2.at(i) != it3.at(i)) && joinEqual(it2, i, it1, i) ){
					++it2.at(i);
				}
			}
			
			
		}
		
		return true;
		
	}
	
	
	
	inline double weight(int j1, const vector<Row>::const_iterator& iit1){
	
		return (*weighters.at(j1))(*iit1);
	}
	
	
	shared_ptr<Rows> joinResult;
	
	
	
	shared_ptr<Rows> inverted;
	const vector<double>* toInvert = NULL;
	long toInvertPos = -1;
	
	vector<long double> dstat1;
	vector<long double> dstat2;
	vector<long double> dstat3;
	
	Tuple joincand;
	
	inline void setJoinColumns(const Tuple& t){
	
		joincand = t;
	}
	
	WeightingFunction* fptr = NULL;
	
	inline shared_ptr<Rows> join(){
	
		const long jb = joinbytes();
		
		init();
		
		WeightingFunction f;
		
		weights(f);
		
		/*
		if (tables.size() == 2){
		
		
			MainMemoryIndex* ind1 = a.getIndex(projsJoin.at(0) );
		
		}*/
		
		joinResult = std::make_shared<Rows>(joinCols.size() );
		
		joinResult->reserve(jb);
		
		internalJoin(true, true);
		
		shared_ptr<Rows> ret = joinResult;
		
		joinResult.reset();
		
		return ret;
	}
	
	inline shared_ptr<Rows> invert(const vector<double>& v, WeightingFunction& f){
			
		init();
		
		weights(f);
		
		toInvert = &v;
		toInvertPos = 0;
		
		inverted = std::make_shared<Rows>(joinCols.size() );
		
		inverted->reserve(v.size()*16 * joinCols.size()  );
		
		
		const long double joinweight = internalJoin(false, false);
		
		shared_ptr<Rows> ret = inverted;
		
		assert (ret->size() == v.size() );
		
		inverted.reset();
		toInvert = NULL;
		toInvertPos = -1;
		
		cout << endl;
		cout << "join weight " << ((long) joinweight) << endl;
		cout << endl;
		
		return ret;
	}
	
	inline shared_ptr<Rows> invert(const vector<double>& v){
	
		WeightingFunction f;
		
		return invert(v, f);
	}
	
	
	// m == n => full row
	inline bool proc(long cp1, long cp2, long double wp1, long double wp2,  vector<const Row*>& cursor, int m, int n){
	
	
		if (joinResult){
		
			if (m == n){
			
				joinResult->push_back(cursor, projsSel);
				return true;
			}
			
			return false;
		}
		
		if (toInvertPos != -1){
		
			if (toInvertPos >= toInvert->size() )
				return true;
				
			if (wp2 < toInvert->at(toInvertPos) )
				return true;
		
			if (m == n){
			
				assert (wp1 <= toInvert->at(toInvertPos) );
				assert (toInvert->at(toInvertPos) <= wp2);
				
				
				
				while (toInvertPos < toInvert->size() && toInvert->at(toInvertPos) <= wp2){
			
					inverted->push_back(cursor, projsSel);
					
					if (toInvertPos % 1000 == 0){
					
						cout << toInvertPos << endl;
					}
					
					++toInvertPos;
					
				}
				
				return true;
			}
		
			return false;
		}
		
		
		if (toRankPos != -1){
		
			const long nn = toRank.at(0)->v.size();
		
			if (toRankPos >= nn)
				return true;
		
			if (m == n){
			
				bool br = false;
				
				long count = 0;
			
				while (toRankPos < nn ){
		
					for (int i = 0; i < m; ++i){
				
						const Row& a = *cursor.at(i);
						const Row& b = toRank.at(i)->v.at(toRankPos);
					
						if ( cmp(i)(a,b) || cmp(i)(b,a) ){
						
							br = true;
							break;
						}
					}
					
					if (br)
						break;
					
					if ( (toRankPos == 0) || (toRankPos % 10000 == 0) || (toRankPos == (nn-1)) ){
					
						cout << "toRankPos " << toRankPos << " rank " << cp1 << " weightrank " << wp1 << " to " << wp2 << " nn " << nn << endl;
					}
					ranks.push_back(wp1);
					++toRankPos;
					++count;
				}
				
				if (count > 0){
				
					dstat1.push_back(wp1);
					dstat2.push_back(wp2-wp1);
					dstat3.push_back(count);
					
				}
				
				return true;
			}
		
			// skip if previous part of row is not in sample
			for (int i = 0; i < m; ++i){
				
				const Row& a = *cursor.at(i);
				const Row& b = toRank.at(i)->v.at(toRankPos);
				
				if ( cmp(i)(a,b) || cmp(i)(b,a) )
					return true;
			}
			
			return false;
		}
		
		return true;
	}
	
	vector< shared_ptr<SortedTable> > toRank;
	long toRankPos = -1;
	vector<long double> ranks;
	
	template<typename T>
	inline vector<long double> rank(const T& t, WeightingFunction& f, string seed ="abc123"){
	
		assert (toRank.size() == 0);
		assert (toRankPos == -1);
		assert (ranks.size() == 0);
		
		assert (dstat1.size() == 0);
		assert (dstat2.size() == 0);
		assert (dstat3.size() == 0);
		
		
		cout << "init" << endl;
		init();
		
		cout << "weights" << endl;
		
		weights(f);
		
		cout << "sortedtable" << endl;
	
		SortedTable st(joinCols);
		
		long siz = st.v.size();
		
		cout << "sortedtable add" << endl;
		
		st.add(t);
		
		dstat1.reserve(siz);
		dstat2.reserve(siz);
		dstat3.reserve(siz);
		
		st.sort( cmp() );
		
		Row r;
		
		int pos = 0;
		
		// check for sorted
		for (auto it = st.v.begin(); it != st.v.end(); ++it){
		
			if (pos > 0)
				assert (!cmp()(*it, r));
			
			r = (*it);
			++pos;
		}
		
		toRankPos = 0;
		
		for (int i = 0; i < tables.size(); ++i){
		
			cout << "tab " << i << " " << tables.at(i)->cols << endl;
		
			auto tab = std::make_shared<SortedTable>( tables.at(i)->cols );
			
			assert (UtilTuple::equals(UtilTuple::select(joinCols, projsTab.at(i)), tables.at(i)->cols));
			
			cout << "add " << i << " " << UtilTuple::select(joinCols, projsTab.at(i)) <<endl;
			
			tab->add( st.v, projsTab.at(i));
			
			for (auto it = tab->v.cbegin(); it != tab->v.cend(); ++it){
			
				cout << "compweight(" << (*it) << ") = " << weight(i, it) << endl;
			
				break;
			}
			
			toRank.push_back(tab);
		}
		
		ranks.reserve(st.v.size() );
		
		cout << "internaljoin" << endl;
		
		const long double joinweight = internalJoin(false, false);
		
		cout << "joinweight floor(" << joinweight << ") == " << ((long)joinweight) << endl;
		
		vector<long double> ret = ranks;
		
		
		cout << "dstatcomp" << endl;
		
		ret.push_back(dstat(joinweight, seed));
		
		ranks.clear();
		toRankPos = -1;
		toRank.clear();
		
		dstat1.clear();
		dstat2.clear();
		dstat3.clear();
		
		return ret;
	}
	
	inline long double dstat(const long double joinweight, const string seed = "test123") const{
	
		
		RandomNumbers rnd(seed);
	
		double ret = 0;
		
		long c = 0;
		
		vector<double> smp;
		
		smp.reserve (dstat1.size() );
		
		assert (dstat1.size() > 0);
	
		long smpsize = 0;
	
		for (long j = 0; j < dstat3.size(); ++j)
			smpsize += dstat3.at(j);
	
		for (long j = 0; j < dstat1.size(); ++j){
		
			const long cnt = dstat3.at(j);
			
			smp.clear();
			
			smp.push_back(0);
			
			for (int k = 0; k < cnt; ++k)
				smp.push_back( rnd.randomDouble() );
			
			std::sort(smp.begin(), smp.end() );
			
			//cout << "cnt " << cnt << endl;
			
			for (int k = 0; k <= cnt; ++k){
				
				const long double est = (c+k)*1.0L/smpsize;
				const long double tru = ( dstat1.at(j)+dstat2.at(j)*smp.at(k) )/joinweight;
				
				const long double d = est-tru;
			
				if (d < 0.0){
				
					if ( (-d) > ret)
						ret = -d;
				} else {
				
					if (d > ret)
						ret = d;
				}
				//cout << "k " << k << " d " << d << endl;
				
				
			}
			
			c += cnt;
		}
		
		return ret;
	}
	
	
	template<typename T>
	inline vector<double> rank(const T& t, string seed ="abc123"){
	
		WeightingFunction f;
		
		return rank(t, f, seed);
	}
	
	
	inline double joinweight(WeightingFunction& f){
	
		init();
		
		weights(f);
		
		
		bool weightsDisabled = true;
		
		for (int i = 0; i < tables.size(); ++i)
			weightsDisabled &= weighters.at(i)->isDisabled();
		
		if (weightsDisabled)
			return joinsize();
			
		
		return internalJoin(false, false);
	}
	
	/*
	inline long joinsize(){
	
		WeightingFunction f;
		
		return joinweight(f);
	}*/
	
	
	inline long double internalJoin(const bool skipcounts, const bool skipweights){
	
		bool weightsDisabled = true;
		
		for (int i = 0; i < tables.size(); ++i)
			weightsDisabled &= weighters.at(i)->isDisabled();
		
		if (weightsDisabled && (!skipweights))
			return internalJoin(skipcounts, true);
		
		begin();
		
		const int n = tables.size();
		
		vector<const Row*> cursor;
		
		cursor.clear();
			
		for (int i = 0; i < n; ++i)
			cursor.push_back(NULL);
		
		long countpos = 0;
		long double weightpos = 0;
		
		vector<long double> weightsums;
		
		vector<long> countsums;
		
		weightsums.clear();
		
		for (int i = n; i >= 0; --i){
		
			weightsums.push_back(1);
			countsums.push_back(1);
		}
		
		vector<shared_ptr<vector<long double>>> weights;
		
		for (int i = n-1; i >= 0; --i)
			weights.push_back( std::make_shared<vector<long double>>() );
		
		Timer t("internaljoin");
		
		t.start();
		
		long h = 0;
		
		while ( next() ){
		
			++h;
			
			if (h % 1000000 == 0){
			
				cout << h << endl;
			}
			
			for (int i = n; i >= 0; --i){
			
				weightsums.at(i) = 1;
				countsums.at(i) = 1;
			}
			
			if (!skipcounts){
			
				for (int i = n-1; i >= 0; --i){
			
					countsums.at(i) = 0;
			
					for (auto it = it1.at(i); it != it2.at(i); ++it){
				
						countsums.at(i) += countsums.at(i+1);
					}
				}
			}
			
			for (int i = n-1; i >= 0; --i){
		
				weightsums.at(i) = 0;
				
				vector<long double>& ws = *weights.at(i);
				
				ws.clear();
		
				for (auto it = it1.at(i); it != it2.at(i); ++it){
			
					const long double w = skipweights ? 1 : weight(i, it);
					
					if (ws.capacity() == ws.size() )
						ws.reserve( ws.size()*2 );
					
					ws.push_back(w);
					
					weightsums.at(i) += w * weightsums.at(i+1);
				}
			}
			
			for (int iii = 0; iii < 1; ++iii){
			
				const long cp = skipcounts ? 0 : countpos+countsums.at(0);
				const long double wp = skipweights ? cp : weightpos+weightsums.at(0);
			
				if (proc(countpos, cp, weightpos, wp, cursor, 0, n) ){
					weightpos = wp;
					countpos = cp;
					continue;
				}
				
				auto it_0w = weights.at(0)->begin();
					
				for (auto it_0 = it1.at(0); it_0 != it2.at(0); ++it_0, ++it_0w){
			
					cursor.at(0) = &(*(it_0));
					
					const long double w0 = (*it_0w);
			
					const long cp0 = skipcounts ? 0 : countpos+countsums.at(1);
					const long double wp0 = skipweights ? cp0 : weightpos+w0 * weightsums.at(1);
					
					if (proc(countpos, cp0, weightpos, wp0, cursor, 1, n) ){
						weightpos = wp0;
						countpos = cp0;
						continue;
					}
					
					auto it_1w = weights.at(1)->begin();
					
					for (auto it_1 = it1.at(1); it_1 != it2.at(1); ++it_1, ++it_1w){
			
						cursor.at(1) = &(*(it_1));
						
						const long double w1 = w0 * (*it_1w);
			
						const long cp1 = skipcounts ? 0 : countpos+countsums.at(2);
						const long double wp1 = skipweights ? cp1 : weightpos+w1 * weightsums.at(2);
						
				
						if (proc(countpos, cp1, weightpos, wp1, cursor, 2, n) ){
							weightpos = wp1;
							countpos = cp1;
							continue;
						}
						
						auto it_2w = weights.at(2)->begin();
						
						for (auto it_2 = it1.at(2); it_2 != it2.at(2); ++it_2, ++it_2w){
						
							cursor.at(2) = &(*(it_2));
		
							const long double w2 = w1 * (*it_2w);
		
							const long cp2 = skipcounts ? 0 : countpos+countsums.at(3);
							const long double wp2 = skipweights ? cp2 : weightpos+w2 * weightsums.at(3);
							
							if (proc(countpos, cp2, weightpos, wp2, cursor, 3, n) ){
								weightpos = wp2;
								countpos = cp2;
								continue;
							}
						
							auto it_3w = weights.at(3)->begin();
						
							for (auto it_3 = it1.at(3); it_3 != it2.at(3); ++it_3, ++it_3w){
				
								cursor.at(3) = &(*(it_3));
								
								const long double w3 = w2 * (*it_3w);
	
								const long cp3 = skipcounts ? 0 : countpos+countsums.at(4);
								const long double wp3 = skipweights ? cp3 : weightpos+w3 * weightsums.at(4);		
						
								if (proc(countpos, cp3, weightpos, wp3, cursor, 4, n) ){
									weightpos = wp3;
									countpos = cp3;
									continue;
								}
								
								assert (false);
								weightpos = wp3;
								countpos = cp3;
							} // it_3
						
							weightpos = wp2;
							countpos = cp2;
						} // it_2
					
						weightpos = wp1;
						countpos = cp1;
					} // it_1
					
					weightpos = wp0;
					countpos = cp0;
				} // it_0
				weightpos = wp;
				countpos = cp;
			} // iii
		} // while
		
		t.end();
		
		return weightpos;
	}
	
	
	inline long joinbytes(){
	
		init();
	
		begin();
		
		long ret = 0;
		
		const int n = tables.size();
		
		while ( next() ){
		
			long m = 1;
			
			for (int i = 0; i < n; ++i){
			
				long bytes = 0;
				
				for (auto it = it1.at(i); it != it2.at(i); ++it)
					bytes += it->getBytes();
			
				m *= bytes;
				
			}
				
			ret += m;
		}
	
		return ret;
		
	}
	
	
	
	
	
	inline double cyclicSelectivity(WeightingFunction& f, long smps, string seed="abc123"){
	
	
		assert (false);
	
		//if (true)
		//	return joinsizeEstimate()*1.0L/(tables.at(0)->tabptr->getSize()*tables.at(1)->tabptr->getSize());
		
		cout << "A4B1" << endl;
	
		init(false);
		
		cout << "A4B2" << endl;
		
		assert (tables.size() == 2);
		
		const DatabaseTable& r1 = *tables.at(0)->tabptr;
		const Tuple& t1 = r1.getColNames();
		
		
		const DatabaseTable& r2 = *tables.at(1)->tabptr;
		const Tuple& t2 = r2.getColNames();
			
		auto smp1_ = withrepl(r1, t1, smps, seed);
		auto smp2_ = withrepl(r2, t2, smps, seed+"2");
		
		auto smp1 = smp1_->getRowVector();
		auto smp2 = smp2_->getRowVector();
		
		RandomNumbers rnd1(seed);
		
		rnd1.randomPermutation(*smp1);
		rnd1.randomPermutation(*smp2);
		
		assert (smp1->size() == smps);
		assert (smp2->size() == smps);
		
		const vector<int> proj1 = UtilTuple::without(projsJoin.at(0), -1);
		const vector<long> hash1(proj1.size(), UtilHash::HASHMAX);
		
		const vector<int> proj2 = UtilTuple::without(projsJoin.at(1), -1);
		const vector<long> hash2(proj2.size(), UtilHash::HASHMAX);
		
		assert (proj1.size() > 0);
		
		assert (UtilTuple::equals( UtilTuple::select( t1, proj1), UtilTuple::select( t2, proj2)));
		
		auto it1 = smp1->begin();
		
		long cnt = 0; 
		
		long k = 0;
		
		for (auto it2 = smp2->begin(); it1 != smp1->end() && it2 != smp2->end(); ++it1, ++it2){
		
			const Row& row1 = (*it1);
			const Row& row2 = (*it2);
			
			
				
				
			++k;
			
			if (row1.equals(proj1, row2, proj2)){
			
				if (k%10000 == 0)
					cout << "tab1 " << t1 << " tab2 " << t2 << " join " << t1[proj1[0]] << " row1 " << row1[proj1[0]] << " == row2 " << row2[proj2[0]] << endl;
			
				++cnt;
			} else {
			
				if (k%10000 == 0)
					cout << "tab1 " << t1 << " tab2 " << t2 << " join " << t1[proj1[0]] << " row1 " << row1[proj1[0]] << " != row2 " << row2[proj2[0]] << endl;
			}
			
			
		}
		
		cout << "CNT " << cnt << endl;
		
		if (cnt > 0)
			return cnt*1.0/smps;
		
		return 0.1/smps;
		
		//shared_ptr<WeightedRowSample> s1 = std::make_shared<WeightedRowSample>(t1, seed, smps );	
		/*LongRowHashMap map1;
		
		
		for (auto it = r1.begin(); it != r1.end(); ++it)
			map1.addCount(*it, proj1, hash1);
		
		*/
		
		//double ret = 0;
		/*
		shared_ptr<WeightedRowSample> s2 = std::make_shared<WeightedRowSample>(t2, seed, smps );	
	
		Weighter weight(f, t2);
		
		double wsum = 0;

		for (auto it = r2.begin(); it != r2.end(); ++it)
			wsum += weight(*it);			
		
		
		
		s2->initOffline(wsum, true);
	
		for (auto it = r2.begin(); it != r2.end(); ++it)
			s2->add(*it, weight(*it));
	
		s2->withReplacement();
		
		
		
		
		double div = 0;
		
		for (auto it = s2->begin(); it != s2->end(); ++it){
			ret += map1.getCount(*it, proj2, hash2)*1.0/r1.size();
			
			++div;
		}
		
		ret /= div;*/
		
		/*
		cout << "A4B2" << endl;
		
		for (auto it = r2.begin(); it != r2.end(); ++it)
			ret += map1.getCount(*it, proj2, hash2)*1.0/r1.size();
		
		
		ret /= r2.size();
		
		cout << "A4B3" << endl;
		
		return ret;*/
		
		/*
		
		LongRowHashMap map1;
		{
			const Rows& r1 = *tables.at(0)->rows;
			const vector<int> proj1 = UtilTuple::without(projsJoin.at(0), -1);
			const vector<long> hash1(proj1.size(), UtilHash::HASHMAX);
		
			for (auto it = r1.begin(); it != r1.end(); ++it)
				map1.addCount(*it, proj1, hash1);
		}
		
		const vector<int> proj2 = UtilTuple::without(projsJoin.at(1), -1);
		const vector<long> hash2(proj2.size(), UtilHash::HASHMAX);
		
		double wsum = 0;
		
		const Rows& r2 = *tables.at(1)->rows;
		const Tuple& t2 = tables.at(1)->cols;
		
		Weighter weight(f, t2);

		for (auto it = r2.begin(); it != r2.end(); ++it)
			wsum += weight(*it);			
		
		shared_ptr<WeightedRowSample> s2 = std::make_shared<WeightedRowSample>(t2, seed, smps );	
			
		s2->initOffline(wsum, true);
		
		for (auto it = r2.begin(); it != r2.end(); ++it)
			s2->add(*it, weight(*it));
		
		s2->withReplacement();
		
		long k = 0;
		long c = 0;
		
		for (auto it = s2->begin(); it != s2->end(); ++it){
		
			if (map1.getCount(*it, proj2, hash2) > 0)
				++c;
				
			++k;
		}
		
		assert (k == smps);
		
		return c*1.0/smps;*/
		
	}
	
	inline long joinsize(){
	
	
		if (Global::THREADS_ENABLED == false)
		if (!sorted && tables.size() == 2){
		
			init(false);
			
			LongRowHashMap map1;
			{
				const Rows& r1 = *tables.at(0)->rows;
				const vector<int> proj1 = UtilTuple::without(projsJoin.at(0), -1);
				const vector<long> hash1(proj1.size(), UtilHash::HASHMAX);
			
				for (auto it = r1.begin(); it != r1.end(); ++it)
				map1.addCount(*it, proj1, hash1);
			}
			
			long ret = 0;
			{
				const Rows& r2 = *tables.at(1)->rows;
				const vector<int> proj2 = UtilTuple::without(projsJoin.at(1), -1);
				const vector<long> hash2(proj2.size(), UtilHash::HASHMAX);
			
				for (auto it = r2.begin(); it != r2.end(); ++it)
					ret += map1.getCount(*it, proj2, hash2);
			}
			
			return ret;
		}
		
		init();
	
		begin();
		
		long ret = 0;
		
		const int n = tables.size();
		
		while ( next() ){
		
			long m = 1;
			
			for (int i = 0; i < n; ++i){
			
				//cout << tables.at(i)->cols << " " << *(it1.at(i)) << " ... " << *(it2.at(i)) << endl;
				
				const long g = it2.at(i)-it1.at(i);
				
				m *= g;
				
			}
				
			ret += m;
		}
	
		return ret;
		
	}
	
	
	inline long joinsizeSort(){
	
		init();
	
		begin();
		
		long ret = 0;
		
		const int n = tables.size();
		
		while ( next() ){
		
			long m = 1;
			
			for (int i = 0; i < n; ++i){
			
				//cout << "123 " << *(it1.at(i)) << " " << *(it2.at(i)) << endl;
				
				const long g = it2.at(i)-it1.at(i);
				
				m *= g;
				
			}
				
			ret += m;
		}
	
		return ret;
		
	}
	
	
	inline void add(shared_ptr<SortedTable> tab){
	
		assert (!initialised);
		
		//cout << tab->cols << endl;
		tables.push_back(tab);
		
		
	}
	
	
	inline void add(DatabaseTable& tab){
	
		add( std::make_shared<SortedTable>(tab));
		
	}
	
	
	bool initialised = false;
	
	bool sorted = false;
	
	inline void init(bool sort = true){
	
		if (!initialised){
		
			update();
			initialised = true;
		
		}
		
		if (!sorted){
		
			if (sort){
			
				for (int i = 0;i < this->tables.size(); ++i)
					this->tables.at(i)->sort( cmp(i) );
					
				sorted = true;
			}
		}
	}
	
	
	// baseTableCols[{projsJoin}] = baseTableCols that are join columns
	vector<vector<int>> projsJoin;
	
	// baseTableCols[{projsNonJoin}] = baseTableCols that are non-join columns
	vector<vector<int>> projsNonJoin;
	
	// baseTableCols[{projsSel}] = baseTableCols that appear in join rows
	vector<vector<int>> projsSel;
	
	// joinRows[{projsTab[baseTableIndex]}] = baseTableCols
	vector<vector<int>> projsTab;
	
	// joinRows[{projsTab[baseTableIndex]}] = baseTableCols
	vector<vector<int>> projsLex;
	
	vector<int> projLex;
	
	
	
	
	vector<shared_ptr<Weighter>> weighters;
	
	inline void weights(WeightingFunction& f){
	
		assert (initialised);
	
		weighters.clear();
		
		assert (projsTab.size() == tables.size() );
	
		for (int i = 0; i < tables.size(); ++i){
		
			//cout << joinCols << endl;
			//cout << UtilTuple::select(joinCols, projsTab.at(i)) << endl;
			cout << "weighter " << UtilTuple::select( UtilTuple::select(joinCols, projsTab.at(i)), projsSel.at(i)) << endl;
		
			weighters.push_back( std::make_shared<Weighter>(f, UtilTuple::select(joinCols, projsTab.at(i)), UtilTuple::select( UtilTuple::select(joinCols, projsTab.at(i)), projsSel.at(i))) );
		}
	}
	
	
	
	
	shared_ptr<RowComparator> cmp0;
	
	vector<shared_ptr<RowComparator>> cmps;
	
	inline const RowComparator& cmp(int n=-1){
	
		if (n == -1)
			return *cmp0;
			
		return *cmps.at(n);
	}
	
	Tuple joinCols;
	
	inline void update(){
	
		top = -1;
	
		joinCols = Tuple();
	
		projsSel.clear();
		projsJoin.clear();
		projsNonJoin.clear();
		
		projsTab.clear();
		
		projsLex.clear();
		
		projLex.clear();
	
		vector<Tuple> tables;
		
		cmp0.reset();
		cmps.clear();
		
		for (int i = 0;i < this->tables.size(); ++i){
			tables.push_back(this->tables.at(i)->cols);
		}
		
		
		Tuple jcs = UtilTuple::nonUnique(tables);
		
		if (joincand.size() > 0)
			jcs = UtilTuple::intersection(jcs, joincand);
		
		/*
		//cout << tables.size() << endl;
	
		std::set<string> setJoinCols;
		
		for (int i = 0; i < tables.size(); ++i){
		
			const Tuple& ti = tables.at(i);
			
			for (int j = i+1; j < tables.size(); ++j){
			
				const Tuple& tj = tables.at(j);
			
				for (int k = 0; k < tj.size(); ++k){
				
					if ( UtilTuple::find(ti, tj[k]) != -1){
					
						if (setJoinCols.count(tj[k]) == 0){
						
							//cout << "add tj[j] " << tj[k] << endl;
						
							assert (i == 0); // join columns shared by all, otherwise ambigious!
							setJoinCols.insert(tj[k]);
							
							
						}
					}
				}
			}
		}
		
		Tuple jcs;
		
		for (auto it = setJoinCols.begin(); it != setJoinCols.end(); ++it){
		
			jcs.push_back(*it);
		}*/
		
		cout << "join columns " << jcs << endl;
		
		
		
		
		for (int i = 0; i < tables.size(); ++i){
		
			const Tuple& ti = tables.at(i);
			
			vector<int> projSel;
			vector<int> projJoin;
			vector<int> projNonJoin;
			
			bool topCand = true;
			
			for (int k = 0; k < jcs.size(); ++k){
			
				int f = UtilTuple::find(ti, jcs[k]);
				
				if (f == -1)
					topCand = false;
			
				projJoin.push_back( f );
			}
			
			if (top == -1 && topCand)
				top = i;
			
			for (int k = 0; k < ti.size(); ++k){
			
				const bool isJoin = UtilTuple::find(jcs,ti[k]) != -1;
			
				// first appearance of join column || non-join column
				if ( UtilTuple::find(joinCols, ti[k]) == -1 || !isJoin){ 
				
					projSel.push_back(k);	
					
					joinCols.push_back(ti[k]);
					
				}
				
				if (isJoin)
					;//projJoin.push_back(k);
				else
					projNonJoin.push_back(k);
			}
			
			assert (projJoin.size() == jcs.size() );
			
			
			cout << "join" << i << " " << UtilTuple::select(ti, projJoin) << endl;
			
			projsSel.push_back(projSel);
			projsJoin.push_back(projJoin);
			
			
			projsNonJoin.push_back(projNonJoin);
		}
		
		cout << "joinedCols " << joinCols << endl;
		
		assert (top != -1);
		
		assert (projsSel.size() == tables.size() );
		
		for (int i = 0; i < tables.size(); ++i){
		
			const Tuple& ti = tables.at(i);
			
			
			const Tuple& tnonjoin = UtilTuple::select(tables.at(i), projsNonJoin.at(i));
			
			const Tuple& tjoin = UtilTuple::select(tables.at(i), projsJoin.at(i));
			
			if (i == 0){
			
				for (int k = 0; k < tjoin.size(); ++k)
					projLex.push_back( UtilTuple::find(joinCols, tjoin[k]));
			
			}
			
			for (int k = 0; k < tnonjoin.size(); ++k)
				projLex.push_back( UtilTuple::find(joinCols, tnonjoin[k]));
			
		
			vector<int> projTab;
			
			vector<int> projLex1;
			
			for (int k = 0; k < ti.size(); ++k)
				projTab.push_back( UtilTuple::find(joinCols, ti[k]));
			
			for (int k = 0; k < projsJoin.at(i).size(); ++k)
				projLex1.push_back(projsJoin.at(i).at(k));
			
			for (int k = 0; k < projsNonJoin.at(i).size(); ++k)
				projLex1.push_back(projsNonJoin.at(i).at(k));
			
			projsLex.push_back(projLex1);
			
			cout << "lex" << i << " " << UtilTuple::select(ti, projLex1) << endl;
			
			projsTab.push_back(projTab);
		}
		
		cmp0 = std::make_shared<RowComparator>(projLex);
		
		cout << "lex " << UtilTuple::select(joinCols, projLex) << endl;
		
		for (int i = 0; i < tables.size(); ++i){
		
			const Tuple& ti = tables.at(i);
			
			cmps.push_back(std::make_shared<RowComparator>(projsLex.at(i)));
		
			//cout << ti << " == " << UtilTuple::select(joinCols, projsTab[i]) << endl;
			
			//cout << ti << " == " << UtilTuple::select(joinCols, projsTab[i]) << endl;
			
			assert( UtilTuple::equals( ti, UtilTuple::select(joinCols, projsTab.at(i)) ) );
			
		}
		
		
		
		//cout << "lex " << UtilTuple::select(joinCols, projLex) << endl;
		
		
	}
	
	
	

};


class DatabaseTables{

public:
	vector<DatabaseTable*> tabs;
	
	Tuple joincand;
	
	inline void add(DatabaseTable& t){
	
		tabs.push_back(&t);
	}

	inline void setJoinColumns(const Tuple& t){
	
		assert (t.size() > 0);
	
		joincand = t;
	}
	
	inline int joincolumns(){
	
		assert (joincand.size() > 0);
	
		Tuple t = joincand.size() > 0 ? joincand : tabs.at(0)->getColNames();
		
		for (int i = 0; i < tabs.size(); ++i)
			t = UtilTuple::intersection(t, tabs.at(i)->getColNames());
		
		return t.size() > 0 ? 1 : -1;
		
		/*
		int ret = 0;
		
		std::unordered_map<string,int> cols;
		
		for (int i = 0; i < tabs.size()-1; ++i){
		
			const Tuple& t1 = tabs.at(i)->getColNames();
		
		
			for (int ii = 0; ii < t1.size(); ++ii){
			
				const string s = t1[ii];
			
				for (int j = i+1; j < tabs.size(); ++j){
		
					const Tuple& t2 = tabs.at(j)->getColNames();
					
					
					if (UtilTuple::find(t2, s) != -1)
						cols[s]++;
				}
			
			}
		}
		
		for (int i = 0; i < tabs.size(); ++i){
		
			const Tuple& t1 = tabs.at(i)->getColNames();
			
			for (auto it = cols.begin(); it != cols.end(); ++it){
			
				if (UtilTuple::find(t1, it->first) == -1)
					return -1;
			}
		}
		
		return 1;*/
	}

	inline DatabaseTables(){
	
	
	}
	
	inline DatabaseTables(shared_ptr<DatabaseTable> t1){
	
		add(t1);
	}
	
	inline DatabaseTables(vector<shared_ptr<DatabaseTable>> ts){
	
		for (int i = 0; i < ts.size(); ++i)
			add(ts.at(i));
	}
	
	inline DatabaseTables(shared_ptr<DatabaseTable> t1, shared_ptr<DatabaseTable> t2){
	
		add(t1);
		add(t2);
	}
	
	inline DatabaseTables(shared_ptr<DatabaseTable> t1, shared_ptr<DatabaseTable> t2, shared_ptr<DatabaseTable> t3){
	
		add(t1);
		add(t2);
		add(t3);
	}
	
	inline DatabaseTables(DatabaseTable& t1){
	
		add(t1);
	}
	
	inline DatabaseTables(DatabaseTable& t1, DatabaseTable& t2){
	
		add(t1);
		add(t2);
	}
	
	inline DatabaseTables(DatabaseTable& t1, DatabaseTable& t2, DatabaseTable& t3){
	
		add(t1);
		add(t2);
		add(t3);
	}
	
	LazyJoin sj_;
	
	bool initialised = false;
	
	inline LazyJoin& sj(){
		
		
		if (initialised)
			return sj_;
			
		if (joincand.size() > 0)
			sj_.setJoinColumns(joincand);	
		
	
		for (auto it = tabs.begin(); it != tabs.end(); ++it)
			sj_.add( *(*it) );
		
		
		initialised = true;
		
		return sj_;
	}

	
	inline long double cyclicSelectivity(WeightingFunction& f, long smps, string seed="abc123"){
	
		return joinsizeEstimate()/cartesiansize();
		
		//return sj().cyclicSelectivity(f,smps, seed);
	}

	inline void add(shared_ptr<DatabaseTable> t){
	
		add(*t);
	}

	inline void add(DatabaseTable* t){
	
		add(*t);
	}
	
	inline void samplingActivate(double x, double d = -1){
	
		if (x >= 0.9)
			return;
	
		int k = 0;
		for (auto it = tabs.begin(); it != tabs.end(); ++it){
		
			if (d < 0)
				(*it)->disableSampling();
			else
				(*it)->enableSamplingFrequency(d, 10000, to_string(k) );
		}
	}
	
	inline double joinweight(WeightingFunction& f){
	
		return sj().joinweight(f);
	}
	
	inline long joinsize(){
	
		return sj().joinsize();
	}
	
	template<typename T>
	inline double dstat(T& t, WeightingFunction& f, string seed="abc123"){
	
		auto v = sj().rank(t, f, seed);
		
		return v.back();
	}
	
	template<typename T>
	inline double dstat(T& t, string seed="abc123"){
	
		auto v = sj().rank(t, seed);
		
		return v.back();
	}
	
	inline long joinbytes(){
	
		return sj().joinbytes();
	}
	
	
	inline shared_ptr<DatabaseTable> join(string name){
	
		shared_ptr<Rows> jr = sj().join();
	
		return std::make_shared<DatabaseTable>(name, sj().joinCols,  jr);
	}
	
	inline shared_ptr<DatabaseTable> join(WeightingFunction& f){
	
		shared_ptr<Rows> jr = sj().join();
	
		return std::make_shared<DatabaseTable>(UtilTuple::toString(sj().joinCols), sj().joinCols,  jr);
	}
	
	inline shared_ptr<DatabaseTable> join(){
	
		shared_ptr<Rows> jr = sj().join();
	
		return std::make_shared<DatabaseTable>(UtilTuple::toString(sj().joinCols), sj().joinCols,  jr);
	}
	
	
	
	
	/*
	inline shared_ptr<DatabaseTable> sample(long sampleSize, WeightingFunction& f, string seed = "abc123"){
	
		const double wsum = sj().joinweight(f);
		
		cout << "wsum " << wsum << endl;
		
		RandomNumbers rnd(seed);
		
		vector<double> v;
		
		v.reserve(sampleSize);
		
		for (long j = 0; j < sampleSize; ++j)
			v.push_back( rnd.randomDouble() * wsum );
		
		std::sort(v.begin(), v.end() );
		
		return std::make_shared<DatabaseTable>("noname", sj().joinCols, sj().invert(v, f) );
	}*/
	
	inline shared_ptr<WeightedRowSample> sample(long sampleSize, WeightingFunction& f, string seed = "abc123"){
	
		const double wsum = sj().joinweight(f);
		
		cout << "wsum " << wsum << endl;
		
		RandomNumbers rnd(seed);
		
		vector<double> v;
		
		v.reserve(sampleSize);
		
		for (long j = 0; j < sampleSize; ++j)
			v.push_back( rnd.randomDouble() * wsum );
		
		std::sort(v.begin(), v.end() );
		
		auto smp = sj().invert(v, f);
		
		Weighter weight(f, sj().joinCols );
		
		shared_ptr<WeightedRowSample> ret = std::make_shared<WeightedRowSample>(sj().joinCols, seed, smp->size() );
		
		ret->initOnlineWOR();
			
		for (auto it = smp->begin(); it != smp->end(); ++it)
			ret->add(*it, weight(*it));
			
		ret->withoutReplacement();
		
		ret->setPopulationWeight(wsum);
		
		return ret;
	}
	
	inline long double cartesiansize(){
	
		long ret = 1; 
	
		for (auto it = tabs.begin(); it != tabs.end(); ++it)
			ret *= (*it)->getSize();
			
		return ret;
	}
	
	inline double joinsizeEstimate(){
	
		double d = 0.1;
	
		for (int i = 0; i < 10; ++i){
		
			cout << "joinsizeEstimate " << d << endl;
		
			if (d >= 0.9)
				return joinsize();
		
			samplingActivate(d, d*d);
			
			LazyJoin sj1;
			
			for (auto it = tabs.begin(); it != tabs.end(); ++it)
				sj1.add( *(*it) );
			
			const long siz1 = sj1.joinsize();
			
			samplingActivate(d);
			
			if (siz1 == 0){
			
				if (d >= 0.9)
					return 0;
				
				d *= 100;
				continue;
			}
			
			LazyJoin sj2;
			
			samplingActivate(d, d);
			
			for (auto it = tabs.begin(); it != tabs.end(); ++it)
				sj2.add( *(*it) );
				
			const long siz2 = sj2.joinsize();
			
			samplingActivate(d);
			
			return (long) (siz2*1.0/siz1*siz2);
		}
	
		assert (false);
	
	}

};

class Extensions{

typedef long_double Weight;

private: 

	
	long cost = 0;
	
	bool simpleExtension = true;
	

public:
	vector<Row> samples;

	vector<Weight> randomValues;
	//RowHashMap<shared_ptr<vector<long>>> index;
	
	
	
	
	vector<long> empty;
	
	string tempString;
	
	
	RandomNumbers rnd;
	
	vector<int> smpProj;
	vector<long> smpHash;
	vector<int> extProj;
	vector<long> extHash;
	
	Tuple smpCols;
	Tuple extCols;
	
	Tuple retCols;
	
	
	vector<int> extToSmp;
	
	shared_ptr<WeightingFunction> f;
	
	vector<int> extIds;
	
	JoinColumns rel2;
	
	JoinColumns rel3;
	
	Weight populationWeight;
	long populationSize;
	long maxSize;
	long invalid;
	
	bool fulljoin;
	
	
	/*
	RowHashMap<long> index1;
	vector<long> index;*/
	
	RowHashMultiset index;
	
	/*
	inline RowHashMap<shared_ptr<vector<long>>>& getIndex(){
	
		return index;
	}*/
	
	inline const RowHashMap<long>& getIndex() const{
	
		return index.getHeadCatalog();
	}
	
	unique_ptr<Rows> smpRows;
	unique_ptr<Rows> extRows;
	
	vector<long> extensions;
	
	vector<Weight> weights;
	//vector<Weight> checks;
	
	long idfactory = 0;
	
	
	vector<vector<int> > sources;
	vector<vector<int> > positions;
	
	
	vector<Weight> weightCursors;
	vector<int> skips;
	
	vector<Row> tempRow;
	unique_ptr<Rows> tempRows;
	
	inline void compress(){
	
	
		tempRows->clear();
		
		{
			
			while (tempRow.size() < extRows->size() )
				tempRow.push_back(Row() );
			
			long jj = 0;
			
			for (auto it = extRows->begin(); it != extRows->end(); ++it){
			
				assert (tempRow.at(jj).empty() );
			
				tempRow.at(jj).set( (*it) );
				++jj;
			}
			
			//assert (jj == extensions.size() );
			
			tempRows->reserve(extRows->getBytes() );
			
			for (long j = 0; j < extensions.size(); ++j){
		
				if (extensions.at(j) >= 0){
				
					tempRows->push_back( tempRow.at(extensions.at(j)) );
					extensions.at(j) = tempRows->size()-1;
				}
			}
			
			for (auto it = tempRow.begin(); it != tempRow.end(); ++it)
				it->set();
		}
		
		unique_ptr<Rows> temp = std::move(extRows);
		extRows = std::move(tempRows);
		tempRows = std::move(temp);
	}
	
	
	
	shared_ptr<WeightedRowSample> preparedSample;
	
	inline Extensions(const WeightedRowSample& smp, const JoinColumns& rel, const RowHistogram<Weight>& wt, const string seed = "", MainMemoryIndex* memIndex_ = NULL) : rnd(seed){
		
		
		assert (smp.getSize() > 0);
		
		MainMemoryIndex* memIndex = memIndex_;
		
		tempString.reserve(100);
		
		const long sampleSize = smp.isPopulation() ? smp.getPopulationSize() : smp.getSize();
		
		smpCols = smp.cols;
		extCols = rel.getMainTuple();
		
		tempRows = make_unique<Rows>(extCols.size() );
		
		
		
		fulljoin = smp.isPopulation();
		
		populationWeight = smp.getPopulationWeight();
		populationSize = smp.getPopulationSize();
		maxSize = smp.maxSize;
		invalid = smp.getInvalidSize();
		
		
		
		f = rel.getWeightingFunction();
		
		rel2 = rel;
		
		rel2.setParent(smpCols);
		
		rel3.setMain(extCols);
		rel3.addChild(smpCols);
		
		retCols = rel3.getNovelTuple(); // rel2.getMainTuple(); //
		retCols.push_back(smpCols);
		
		cout << smp.json() << endl;
		cout << "extension " << smpCols << " by " << extCols << " into " << retCols << " sampleSize " <<  sampleSize << " isPop " << smp.isPopulation() << endl;
		
		
		cout << "REL 1 REL 1 REL 1" << endl;
		cout << rel.toString() << endl;
		cout << "REL 1 REL 1 REL 1" << endl;
		
		cout << smpCols << endl;
		
		cout << "REL 2 REL 2 REL 2" << endl;
		cout << rel2.toString() << endl;
		cout << "/REL 2 REL 2 REL 2" << endl;
		
		extIds = f->getVarIds(extCols);
		
		extToSmp = rel2.getTranslationTable(JoinColumn::MAIN, JoinColumn::PARENT);
		
		smpProj = rel2.getMap(JoinColumn::TARGET_MAIN, JoinColumn::PARENT);
		smpHash = rel2.getHashings(JoinColumn::TARGET_MAIN, JoinColumn::PARENT);
		
		extProj = rel2.getMap(JoinColumn::TARGET_MAIN, JoinColumn::MAIN);
		extHash = rel2.getHashings(JoinColumn::TARGET_MAIN, JoinColumn::MAIN);
		
		
		assert (UtilTuple::equals( UtilTuple::select(smpCols, smpProj),  UtilTuple::select(extCols, extProj)));
		
		assert (smpHash.size() == smpProj.size() );
		assert (extHash.size() == extProj.size() );
		assert (smpProj.size() == extProj.size() );
		
		for (int j = 0; j < smpProj.size(); j++){
		
			assert( UtilString::equals( smpCols[smpProj[j]], extCols[extProj[j]]) );
			
			const long k = smpProj[j];
			const string key = smpCols[k];
			const long h = rel2.getHashings()->count(key) > 0 ? (*rel2.getHashings())[key] : std::numeric_limits<long>::max();
		
			assert (smpHash[j] == h);
		}
		
		for (int j = 0; j < extProj.size(); j++){
		
			const long k = extProj[j];
			const string key = extCols[k];
			const long h = rel2.getHashings()->count(key) > 0 ? (*rel2.getHashings())[key] : std::numeric_limits<long>::max();
		
			assert (extHash[j] == h);
		}
		
		
		//const vector<int>& smpToTarget = rel2.getMap(JoinColumn::PARENT, JoinColumn::TARGET);
		
		if (memIndex){
			
			cout << "EXTENSION USING INDEX " << endl;
			
			assert (memIndex->hasJoiner() );
			assert (smpProj.size() == 1);
		
			const long column = smpProj.at(0);
			
			long elemInd = 0;
			
			const vector<int>& novel = rel3.getMap(JoinColumn::NOVEL, JoinColumn::MAIN);
			
			auto ret = std::make_shared<WeightedRowSample>(retCols, "", maxSize );
			
			ret->init(smp.maxSize-smp.getInvalidSize() );
			ret->setInvalidSize(invalid);
			
			for (auto it = smp.begin(); it != smp.end(); ++it){
	
				const long freq = it.getFrequency(); //smp.isPopulation() ? it.getWeight() : it.getFrequency();

				if (freq == 0)
					continue;

				const Row& smpRow = (*it);
	
				assert (freq > 0);
				
				//const Weight weight = wt.getCount(smpRow, smpProj); 
			
				const Row& extRow = *memIndex->get( smpRow, column, rnd.randomDouble() );
			
				const bool valid = validExtension(smpRow, extRow);
				
				for (long i = freq; i >= 1; --i){
					
					if (valid)
						ret->push_back(extRow, novel, smpRow, 1 );
					else
						ret->push_back_invalid(1);
				}
			}
			
			ret->setPopulationWeight(populationWeight);
			ret->setPopulationSize(populationSize);
			
			ret->finalise();
			
			preparedSample = ret;
			
			return;
		
		} else {
		
			//cout << "X20" << endl;
		
			{
				smpRows = make_unique<Rows>(smp.cols.size() );
				extRows = make_unique<Rows>(rel.getMainTuple().size() );
		
				smpRows->reserve(smp.getBytes() ); // smpCols.size()*8*smp.getBytes()
				extRows->reserve(smp.getBytes() ); // extCols.size()*8*smp.getBytes()*2
		
				assert (weightCursors.size() == 0);
		
				weightCursors.reserve(sampleSize);
				extensions.reserve(sampleSize);
				index.reserve(sampleSize);
				weights.reserve(sampleSize);
			
				long entryInd = 0;
				long elemInd = 0;
				
				//cout << "X200" << endl;
				
				assert(!smp.isPopulation());
		
				for (auto it = smp.begin(); it != smp.end(); ++it){
	
					const int freq = smp.isPopulation() ? it.getWeight() : it.getFrequency();
	
					if (freq == 0)
						continue;
	
					const Row& smpRow = (*it);
		
					assert (freq > 0);
				
					//cout << "X201" << endl;
				
					// add entryInd at beginning of sample group
				
					if ( index.insert(smpRow, smpProj, smpHash, entryInd) ){
					
						weightCursors.push_back(0);
						skips.push_back(0);
					
						
						entryInd++;
						
						assert (weightCursors.size() == entryInd);
						assert (skips.size() == entryInd);
					}
					
					
					//cout << "X202" << endl;
					
					//cout << "smprow " << smpRow << " " << smpRow[smpProj[0]] << endl;
				
					const Weight weight = wt.getCountConst(tempString, smpRow, smpProj); 
					
					if (weight == 0){
					
						wt.print(100);
					}
						
					assert (weight > 0);
					
					//cout << "X203" << endl;
				
				
					// add elemInd for elements of sample group
				
					for (int i = freq; i >= 1; --i){
					
						index.push_back(smpRow, smpProj, smpHash, elemInd);
						smpRows->push_back(smpRow);
						weights.push_back(weight);
						extensions.push_back(-1);
						randomValues.push_back( smp.isPopulation() ? 0 : weight * rnd.randomDouble()); //weight * rnd.randomDouble()
						
						
						elemInd++;
						
						assert (randomValues.size() == elemInd);
						assert (extensions.size() == elemInd);
						assert (weights.size() == elemInd);
						
					}
					
					//cout << "X204" << endl;
				}
				
				//cout << "X205" << endl;
				
				
				smpRows->shrink_to_fit();

				if (smpRows->getSize() != sampleSize)
					cout << smpRows->getSize() << " smpRows->getSize() != sampleSize " << sampleSize << endl;

				assert (smpRows->getSize() == sampleSize);
				
				index.finalise(randomValues);
				
				cout << "X206" << endl;
			}
			
			cout << "X21" << endl;
			
			smpRows->getRows(samples);
			/*
			samples.reserve(sampleSize);
		
			for (auto it = smpRows->begin(); it != smpRows->end(); ++it)
				samples.push_back( (*it) );
			*/
		
			assert (samples.size() == sampleSize);
			assert (randomValues.size() == extensions.size());
		}
		
		cout << "X22" << endl;
		
		for (int i = 0; i < rel2.getChildNum(); ++i){
		
		
			const Tuple& t = rel2.getTuple(JoinColumn::CHILD+i);
			
			
			vector<int> sources1;
			vector<int> positions1;
			
			for (int j = 0; j < t.size(); j++){
			
				int k = UtilTuple::find(extCols, t[j]);
				
				if (k >= 0){
					
					assert (UtilString::equals(extCols[k], t[j]));
					sources1.push_back(JoinColumn::MAIN);
					positions1.push_back(k);
					continue;
				}
				
				k = UtilTuple::find(smpCols, t[j]);
				
				if (k >= 0){
					
					assert (UtilString::equals(smpCols[k], t[j]));
					sources1.push_back(JoinColumn::PARENT);
					
					simpleExtension = false;
					
					positions1.push_back(k);
					continue;
				}
				
				assert (false);
			}
			
			
			sources.push_back(sources1);
			positions.push_back(positions1);
			
			
		}
		
		cout << "X23" << endl;
		
		assert (samples.size() == extensions.size());
		
		
		
		for (int i = 0; i < rel2.getChildNum(); ++i){
			
			const vector<int>& sources1 = sources.at(i);
			const vector<int>& positions1 = positions.at(i);
			
			const int last = sources1.size()-1;
			
			for (int j = 0; j <= last; j++){
			
				const Tuple& t = sources1.at(j) == JoinColumn::MAIN ? extCols : smpCols; 
				const int k = positions1.at(j);
				
				assert ( UtilString::equals(t[k], rel2.get(JoinColumn::CHILD+i)[j]));
			}
		}
	
		cout << "extensionbegin " << smpCols << " by " << extCols << endl;
	}
	
	inline bool validExtension(const Row& smpRow, const Row& extRow) const{
		
		for (int i = 0; i < extCols.size(); i++){
			
			const int j = extToSmp[i];
		
			if (j >= 0){
		
				if (!extRow.equals(i, smpRow, j)){
				
					//cout << "[" << smpCols << "]{" << smpRow << "} != [" << extCols << "]{" << extRow << "}" << endl;
				
					return false;
				}
			}
		}
		
		return true;
	}
	
	inline Weight getJoinWeight(string& temp, const Row& smpRow, const Row& extRow, const vector<const RowHistogram<Weight>*>& children) const{
		
		Weight w = 1;
		
		for (int i = 0; i < children.size(); ++i){
			
			if (w == 0)
				break;
				
			if (children.at(i) == NULL)
				continue;
			
			const RowHistogramNode<Weight>* m = children.at(i)->getConst();
			
			assert (m != NULL);
			
			const vector<int>& sources1 = sources.at(i);
			const vector<int>& positions1 = positions.at(i);
			
			const int last = sources1.size()-1;
			
			for (int j = 0; j <= last; j++){
			
				if (w == 0)
					break;
			
				const Row& row = sources1.at(j) == JoinColumn::MAIN ? extRow : smpRow; 
				const int k = positions1.at(j);
				
				if (!m->containsConst(temp, row, k)){
				
					w = 0;
					
				} else if (j == last){
				
					w *= m->getCountConst(temp, row, k);
					
				} else {
				
					m = m->getConst(temp, row, k);
					assert (m != NULL);
				}
			}
		}
		
		//cout << "getJoinWeight " << smpRow << " " << extRow << " weight " << w << endl;
		
		return w;
	}
	
	inline bool has(const Row& extRow){
	
		return index.hasKey(extRow, extProj, extHash);
	}
	
	
	inline void add(const long& entryInd, const long& elemInd, const Row& smpRow, const Row& extRow){
	
		if (!validExtension(smpRow, extRow)){
				
			extensions.at(elemInd) = -2;

		} else {

			extRows->push_back(extRow);
			extensions.at(elemInd) = extRows->size()-1;
		}
		
		if(entryInd != -1)
			skips.at(entryInd) += 1;
	}
	
	vector<Row> temp;
	
	
	template<typename T>
	inline void addExtensions(const T& buffer, const vector<const RowHistogram<Weight>*>& children){
	
		cout << "add extensions" << endl;
		
		long res = 0;
		for (auto it = buffer.begin(); it != buffer.end(); ++it)
			++res;
		
		
		temp.clear();
		temp.reserve(res);
		
		for (auto it = buffer.begin(); it != buffer.end(); ++it)
			temp.push_back(*it);
		
		for (int i = 0; i < extProj.size(); ++i){
		
			for (long j = 0; j < temp.size(); ++j)
				temp.at(j).setSortKey( extProj[i] );
				
			std::stable_sort(temp.begin(), temp.end() );
		}
		
		for (long j = 0; j < temp.size(); ++j)
			temp.at(j).setSortKey(-1);
		
		Row last;
		long head = -2;
		
		const long* begin = NULL;
		const long* end = NULL;
		
		const long* ptr = NULL;
		
		long entryInd = -1;
		
		for (auto it = temp.begin(); it != temp.end(); ++it){
		
			const Row& extRow = (*it);
		
			if ( (head == -2) ||!extRow.equals(last, extProj) ){
			
				head = index.hasKey(extRow, extProj, extHash) ? index.get(extRow, extProj, extHash) : -1;
				last = extRow;
				
				if (head != -1){
					
					entryInd = index.getEntry(head);
					
					begin = index.begin(head)+skips.at(entryInd);
					end = index.end(head);
					
					ptr = begin;
				}
			}
			
			if (head == -1)	
				continue;
			
			
			if (children.size() > 0 && !simpleExtension){
			
				addExtension(head, *it, children);
				continue;
			}	
			
			if (ptr == end)
				continue;
			
			const Weight novelWeight = rel2.getNovelWeightInMain(extRow);
			
			const Weight w = children.size() == 0 ? novelWeight : novelWeight * getJoinWeight(tempString, extRow, extRow, children);
			
			if (w == 0)
				continue;
				
				
			Weight& weightCursor = weightCursors.at(entryInd);
		
			weightCursor += w;
			
			while (ptr != end && randomValues.at(*ptr) <= weightCursor){
			
				add(entryInd, (*ptr), samples.at(*ptr), extRow);
				++ptr;
			}
				
			
			if (extRows->size() > (smpRows->size() * 2) )
				compress();
			
		}
		
		cout << "/add extensions" << endl;
	}
	
	inline void addExtension(const Row& extRow, const vector<const RowHistogram<Weight>*>& children){
	
		if (!has(extRow))
			return;
		
		const long head = index.get(extRow, extProj, extHash);
		
		addExtension(head, extRow, children);
	}
	
	
	inline void addExtension(const long& head, const Row& extRow, const vector<const RowHistogram<Weight>*>& children){
		
		const long entryInd = index.getEntry(head);
		Weight& weightCursor = weightCursors.at(entryInd);
		
		const long* begin = index.begin(head)+skips.at(entryInd);
		const long* end = index.end(head);
		
		if (begin >= end)
			return;
		
		const Weight novelWeight = rel2.getNovelWeightInMain(extRow);
			
		if (children.size() > 0 && !simpleExtension){
		
			for (const long* ptr = begin; ptr != end; ++ptr){
			
				const long& elemInd = (*ptr);
				
				const Row& smpRow = samples.at(elemInd);
				const Weight w = novelWeight * getJoinWeight(tempString, smpRow, extRow, children);
				
				randomValues.at(elemInd) -= w;
				
				if (randomValues.at(elemInd) < 0){
				
					add(entryInd, elemInd, smpRow, extRow);
					index.swap(begin, ptr);
				}
			}
			
			return;
		}
		
		const Weight w = children.size() == 0 ? novelWeight : novelWeight * getJoinWeight(tempString, extRow, extRow, children);
		
		if (w == 0)
			return;
		
		weightCursor += w;
		
		for (const long* ptr = begin; ptr != end; ++ptr){
		
			const long& elemInd = (*ptr);
			
			if (randomValues.at(elemInd) > weightCursor)
				return;
			
			add(entryInd, elemInd, samples.at(elemInd), extRow);
		}
		
		if (extRows->size() > (smpRows->size() * 2) )
			compress();
	}
	
	
	
	inline shared_ptr<WeightedRowSample> getSample(){
	
	
		if (preparedSample)
			return preparedSample;
		//cout << "extensionend " << smpCols << " by " << extCols << endl;
	
		long fail = 0;
	
		for (long j = 0; j < extensions.size(); ++j){ 
		
			assert (fail < 5);
			
			/*
			if (checks.at(j) != 0){
			
				cout << "checks.at(j) " << checks.at(j) << endl;
				fail++;
			}*/
		
			if (extensions.at(j) == -1){
			
				cout << smpCols << " " << samples.at(j) << endl;
				
				//cout << "key " << index.getKey(samples.at(j), smpProj, smpHash) << endl;
				
				cout << smpCols << " " << UtilString::toString(smpProj) << " " << UtilString::toString(smpHash) << endl;
				
				cout << extCols << " " << UtilString::toString(extProj) << " " << UtilString::toString(extHash) << endl;
				
				//cout << "weight left " << priorities[j] << " out of " << weights[j] << endl;
				
				fail++;
			}
		}
		
		assert (fail == 0);
		
		//cout << "cost " << cost << endl;
		
		compress();
		
		auto ret = std::make_shared<WeightedRowSample>(retCols, "", maxSize );
		
		ret->initPop();
		
		ret->setInvalidSize(invalid);
		
		
		auto it2 = extRows->begin();
		
		
		const vector<int>& novel = rel3.getMap(JoinColumn::NOVEL, JoinColumn::MAIN);
		
		for (long j = 0; j < samples.size(); ++j){
			
			if (extensions.at(j) == -2){
			
				ret->push_back_invalid();
				continue;
			}
			
			const Row& smpRow = samples.at(j);
			
			const Row& extRow = (*it2);
			
			assert (extRow.size() > 0);
			
			//cout << extRow << endl;
			//cout << smpCols << " " << smpRow << " " << extCols << " " << extRow << endl;
			//cout << "#" << endl;
			
			const long_double w = weights.size() > 0 ? weights.at(j) : 1;
			
			if (validExtension(smpRow, extRow) ){
			
				ret->push_back(extRow, novel, smpRow, w );
				
				//ret->push_back(extRow, smpRow, weights.at(j) );
				
			} else {
			
				ret->push_back_invalid(w);
			}
			
			++it2;
		}
		
		ret->setPopulationWeight(populationWeight);
		ret->setPopulationSize(populationSize);
		
		ret->finalise();
		
		//cout << "extensionsample " << smpCols << " by " << extCols << endl;
		
		return ret;
	}
	
};




class JoinParameters{

public:

	const map<string, string> stringParams;
	
	double accrate = 1.0;
	
	long maxvals = 0;

	JoinParameters(const std::map<string, string>& m) : stringParams(m){
	
	
		param_gb = real("maxgb", 999.9);	
		
		weights = enabled("weights") || enabled("weight");
		
		maxvals = (134217728.0L/k)*real("maxgb", 999.9);
	}
	
	inline bool defined(const string& s) const{
	
		return stringParams.count(s) > 0;
	}
	
	inline bool undefined(const string& s) const{
	
		return stringParams.count(s) == 0;
	}
	
	inline const string get(const string& s, string defaultString ="NULL") const {
	
		return undefined(s) ? defaultString : stringParams.at(s);
	}
	
	
	inline long_double _real(const string& s) const{
	
		return UtilString::stringToLongDouble(stringParams.at(s));
	}
	
	inline bool isReal(const string& s) const{
	
		if (undefined(s))
			return false;
	
		const long_double d = _real(s);
	
		return d >= 0 && d != std::numeric_limits<long_double>::max();
	}
	
	inline long_double real(const string& s, long_double def = -1) const{
	
		return !isReal(s) ? def : _real(s);
	}
	
	inline bool enabled(const string& s) const{
	
		return defined(s) && real(s) != 0;
	}
	
	inline bool disabled(const string& s) const{
	
		return defined(s) && real(s) == 0;
	}
	
	long_double param_gb = -1;
	
	
	bool weights = false;
	
	bool breakCycles = true;
	bool preload = false;
	bool buildIndex = true;
	
	int sampler = WeightedSamplers::ONLINE_WOR;
	
	bool cyclic = false;
	
	bool hash = false;
	bool debugcolnames = false;
	
	
	bool validate = false;
	
	string query;
	
	long sampleSize = 0;
	
	long batchSize = 0;
	
	vector<long> tableSizes;
	vector<long> mainCols;
	vector<long> targetCols;
	
	
	long maxSamples = 0;
	long_double n = 0;
	long_double m = 0;
	long_double k = 0;
	
	long_double w = 0;
	
	
	inline void print(){
	
		cout << endl << endl;
		
		cout << "JOINSAMPLING" << endl;
	
		cout << "join                 " << (cyclic ? "cyclic" : "acyclic") << " "<< k << "-way join" << endl;
	
		cout << "number of samples    " << sampleSize << (weights ? " weighted" : "") << endl;
		
		cout << "maximal table size   " << ((long) n) << endl;
		
		cout << "batch size   " << batchSize << endl;
		
		cout << "maximal sample size   " << maxSamples << endl;
	
		cout << "memory budget        " << (param_gb == 999.9 ? "unlimited" : UtilString::bytesToString(param_gb*1024*1024*1024)) << endl;
	
		cout << "strategies           ";
		
		if (buildIndex)
			 cout << "[build index] ";
		
		if (breakCycles)
			 cout << "[break cycles] ";
			
		if (hash)
			 cout << "[hashing trick m=" << m << "] ";
		
		if (maxSamples > sampleSize)
			 cout << "[superset sample s=" << maxSamples << "] ";
		
		if (validate)
			cout << "[ks validation] ";
		
		cout << endl;	
			
		cout << "query           " << query;
		
		cout << endl << endl;
		
	}
	
	
	inline long getMul() const{
	
		if (enabled("index"))
			return 1;
	
		long mul = breakCycles ? (real("mul1", 10.0)/(accrate*accrate) )+0.5 :1;
		
		
		
		if (mul > real("mul2", 50.0))
			mul = real("mul2", 50.0);
		
			
		return mul; //mul;
	}
	
	
	inline void update(){
		
		sampler = real("sampler", sampler);
		
		
		
		breakCycles =  cyclic; // enabled("breakcycles") &&
		
		buildIndex = enabled("index"); // || (!disabled("index") && breakCycles);
		
		hash = enabled("hash");
		
		const long mul = getMul();
		
		m = n;
		
		if (hash){
		
			BasicAggregator aggk;
			BasicAggregator aggm;
			
			aggm.add(m,m);
			aggk.add(1,m);
			
			vector<long> sizes = tableSizes;
			std::sort(sizes.begin(), sizes.end() );
			
			for (int kk = 2; kk <= k; ++kk){
			
				BasicAggregator agg;
				
				// agg.max() = limit such that we only hash kk tables
				for (long j = 0; j <= sizes.size()-kk; ++j)
					agg.add( sizes[j] );
				
				
				long mm = pow( 2.0L*powl(n*1.0L, kk-1)*mul*sampleSize, 1.0/kk );
				
				cout << "kk " << kk << " max " << agg.max() << " mm " << mm << endl;
				
				if (mm <= agg.max() ){
					aggm.add(mm, mm);
					aggk.add(kk, mm);
				}
			}
			
			
			
			m = aggm.argMin();
			maxSamples = getSampleSize(m, sampleSize, aggk.argMin());
			
			cout << "m " << m << " maxSamples " << maxSamples << endl;
		}
		
		if (m >= maxvals)
			maxvals = m;
		
		if (maxSamples >= maxvals)
			maxSamples = maxvals;
		
		
		if (m >= n){
		
			m = n;
			hash = false;
			maxSamples = getSampleSize(m, sampleSize);
		}
		
		if (m != n){
		
			if (maxSamples <= 1000000)
				maxSamples = UtilMath::minVal<long>(n,1000000);
				
			if (m <= 1000000)
				m = UtilMath::minVal<long>(n,1000000);
		
		}
		
		batchSize = ((long)maxSamples) >> 3;
		
		if (batchSize <= 1000000)
			batchSize = UtilMath::minVal<long>(n,1000000);
		
		/*
		batchSize = m/10;
		
		if (batchSize <= 10000)
			batchSize = 10000;
			
		if (batchSize >= 1000000)
			batchSize = 1000000;*/
	}
	
	inline void addTable(long size, long main, long targ = 1, long nov = -1){
	
		tableSizes.push_back(size);
		mainCols.push_back(main);
		targetCols.push_back(targ);
		
		n = UtilMath::maxVal<long>(n, size);
		
		++k;
		
		if (nov == -1)
			w += main;
		else
			w += nov;
			
		update();
	}
	
	inline long getSampleSize(long m, long sampSize, long kk = -1) const{
	
		if (kk = -1)
			kk = k;
		// m = max table size
	
		
		const long mul = getMul();
		
		if (m == n)
			return mul*sampSize;
		
		const long_double p = UtilMath::maxVal<long_double>(n/m, 1.0L);
	
		return UtilMath::minVal<long>(m/2, mul*sampSize * 2.0L*powl(p, kk-1) );
	}
	
	inline long_double sampleGB() const {
	
		return w*stringGB();
	}
	
	inline long_double stringGB() const {
	
		return 16.0L/powl(1024.0L, 3);
	}
	
	inline long_double longGB() const {
	
		return 4.0L/powl(1024.0L, 3);
	}
	
	/*
	inline long_double expectedGB(long m, long sampSize) const {
	
		if (n == 0)
			return 0;
	
		long_double gb = 0;
		
		long_double gbdiv = 1.0L/powl(1024.0L, 3);
	
		for (int j = 0; j < tableSizes.size(); j++){
		
			const long size = tableSizes.at(j);
			const long targ = targetCols.at(j);
			const long main = mainCols.at(j);
			
			gb += pow(m, targ) * (stringGB()*targ+longGB() );
		}
		
		gb += getSampleSize(m, sampSize) * sampleGB();
		
		return gb;
	}*/
	
	
	
	inline void print(long_double gb, long sampSize){
	
		//long m = getM(gb, sampSize);
		long s = getSampleSize(m, sampSize);
		
		if (s < sampSize)
			cout << "gb " << gb << " \tsmp " << sampSize << " not enough space!" << endl;
		else
			cout << "gb " << gb << " \tsmp " << sampSize << " \tm " << m << " \ts " << s << endl;
	}

};

class JoinTree{

private:


	
	long randomId = 1;
	
	JoinTree* parent = NULL;
	
	//const bool LEGACY = false;
	
	shared_ptr<RowHistogram<Weight>> batchesJoinWeights; // JOIN WEIGHTS OF ALL BATCHES, AS MUCH AS FITS INTO MEMORY	
	shared_ptr<RowHistogram<Weight>> scanJoinWeights; // JOIN WEIGHTS OF ALL ROWS, USING HASHINGS TO FIT INTO MEMORY
	
	//unique_ptr<Joiner<long>> scanJoinWeightsJoiner;
	//unique_ptr<RowHeapSampler> sampler;
	
	//shared_ptr<WeightedRowSample> weightedSampler;
	
	//shared_ptr< RowHashMap<TupleSamplerPtr> > samplers;
	//shared_ptr< RowHashMap<shared_ptr<WeightedRowSample>> > weightedSamplers;
	
	shared_ptr<DatabaseTable> table;
	vector<shared_ptr<JoinTree>> children;
	
	bool semijoin = false;
	
	bool fk = false;
	
	shared_ptr< unordered_map<string, long> > hashings;
	
	
	// ABC parent
	// BC joint attributes of parent and table
	// BCD table attributes
	// BCA partial attributes
	
	/*
	vector<int> jointInParent; // coords of BC in ABC
	vector<int> jointInTable; // coords of BC in BCD
	
	
	vector<int> jointInSample; // coords of BC in BCD
	
	
	vector<int> jointInParentHashings; // coords of BC in ABC
	vector<int> jointInTableHashings; // coords of BC in BCD*/
	
	RowHashMap<bool> obligatory;
	
	bool complete = false;
	
	bool completedFirstPass = false;
	
	
	bool exactJoinSize = true;
	
	static const int PHASE_ZERO = 0;
	static const int PHASE_JOINWEIGHTS = 1;
	static const int PHASE_ROOT = 2;
	
	
	
	static const int PHASE_EXTEND = 3;
	static const int PHASE_EXTENDCHILDREN = 4;
	
	static const int STATE_STREAMPASS = 5;
	static const int STATE_PREPCHILDREN = 6;
	static const int STATE_PROCESSBATCH = 7;
	static const int STATE_PREPFORPARENT = 8;
	
	static const int STATE_INDEXPASS = 9;
	
	static const int PHASE_GETSAMPLE = 10;
	
	vector<clock_t> beginTimes;
	
	bool verbose = true;
	
	long maxEntries = 1000000000;
	long maxSamples = 1000000000;
	long maxBatchSize = 1000000000;
	
	uint64_t state = 0;
public:

	
	inline void setFK(bool b = true){
	
		fk = b;
	}
	
	
	inline bool getFK() const{
	
		return fk;
	}
	
	inline void setSemiJoin(bool b = true){
	
		semijoin = b;
	}
	
	
	inline bool getSemiJoin() const{
	
		return semijoin;
	}
	
	shared_ptr<JoinColumns> rel;

	Weight joinSize = 0;
	Weight joinSizeEst = 0;
	
	int scan = 1;

	vector<string> verboseprefix;
	
	Rows buffer;
	
	int batchNum = 0;
	
	inline const vector<shared_ptr<JoinTree>>& getChildren() const {
	
		return children;
	}
	
	
	
	inline static const string explain(const int s){
	
		stringstream ss;
	
		if (s == PHASE_ZERO)
			ss << "PHASE_ZERO";
		
		if (s == PHASE_JOINWEIGHTS)
			ss << "PHASE_JOINWEIGHTS";
			
		if (s == PHASE_ROOT)
			ss << "PHASE_ROOT";
			
		if (s == PHASE_EXTEND)
			ss << "PHASE_EXTEND";
			
		if (s == PHASE_GETSAMPLE)
			ss << "PHASE_GETSAMPLE";
			
		if (s == STATE_STREAMPASS)
			ss << "PASS";
			
		if (s == STATE_INDEXPASS)
			ss << "INDEX";
		
		if (s == STATE_PREPCHILDREN)
			ss << "PREP";
			
		if (s == STATE_PROCESSBATCH)
			ss << "BATCH";
			
		if (s == STATE_PREPFORPARENT)
			ss << "PARENTPREP";
		
		return ss.str();
	}
	
	
	
	inline JoinTree& setLimits(long m, long s, long b){
	
		if (m == -1){
		
			reserve(100000);
		
			return (*this);
		}
		
		maxEntries = m;
		maxSamples = s;
		maxBatchSize = b;
		
		if (!exactJoinSize)
			aprioriShrink(maxEntries);
		
		for (int j = 0; j < children.size(); j++)
			children.at(j)->setLimits(m,s, b);
		
		reserve(100000);
		
		
		return (*this);
	}
	
	inline long getLargestTableSize() const{
	
		long max = table->getSize();
	
		for (int j = 0; j < children.size(); j++)
			max = UtilMath::maxVal<long>(max, children[j]->getLargestTableSize() );
			
		return max;
	}
	
	inline long getLargestChildTableSize() const{
	
		long max = 0;
	
		for (int j = 0; j < children.size(); j++)
			max = UtilMath::maxVal<long>(max, children[j]->getLargestTableSize() );
		
		return max;
	}
	
	const bool wasCyclic;
	const bool isCyclic;
	
	inline JoinTree(shared_ptr<DatabaseTable> table, const Tuple& joinAttributes, shared_ptr<WeightingFunction> f, bool wasCyclic, bool isCyclic) : buffer(table->colNames.size() ), wasCyclic(wasCyclic), isCyclic(isCyclic){
	
		this->hashings = std::make_shared< unordered_map<string, long> >();
		
		rel = std::make_shared<JoinColumns>();
		
		rel->set(f);
		
		
		rel->setMain(table->colNames);
		
		rel->setTarget(joinAttributes);
		rel->setHashings(hashings);
		
		this->table = table;
		
		
		verboseprefix.clear();
		verboseprefix.push_back("");
		
		
		joinSize = 0;
		//sampler.reset();
		
		//weightedSampler.reset();
		
		//allJoinWeightsJoiner.reset();
		//samplers.reset();
		batchesJoinWeights.reset();
		scanJoinWeights.reset();
		
		obligatory.clear();
		
		batchesJoinWeights = std::make_shared<RowHistogram<Weight>>( rel->getTargetTuple() );
		
		if (batchesJoinWeights->update(*hashings))
			rel->updateHashings();
		
		const Tuple& tableAttributes = rel->getMainTuple();
		
		if (batchesJoinWeights->update(*hashings))
			rel->updateHashings();
	
		beginState(PHASE_ZERO);	
	}
	
	inline void reserve(long buf){
	
		const Tuple& tableAttributes = rel->getMainTuple();
	
		if (buf < table->getSize() ){
		
			buffer.reserve(buf*tableAttributes.size()*8);
			
		} else{
		
			buffer.reserve(table->getSize()*tableAttributes.size()*8);
		}
		
		if (maxEntries < table->getSize() ){
		
			batchesJoinWeights->reserve(maxEntries);			
			
		} else{
		
			batchesJoinWeights->reserve(table->getSize() );
		}
	
	}
	
	
	double secsPhase1 = -1;
	double secsPhase1and2 = -1;
	
	inline void beginState(const int s){
	
		if (inState(s)){
			cout << explain(s) << endl;
			assert (false);
		}
	
		state = UtilLongSet::add(state, s);
		stateChange(s, -1);
		
		
		const long ind = s; //UtilLongSet::getIndex(s);
		
		assert (ind != -1);
		
		if (ind != -1){
		
			while (beginTimes.size() <= ind)
				beginTimes.push_back(0);	
		}
		
		if (ind != -1 && beginTimes.size() > ind){
		
			if (verbose){
			
				const string time = UtilTime::dateTimeString(UtilTime::dateTimeNow() );
				
				cout << verboseprefix[0] << table->getTableName() << "#" << scan << " BEGIN " << explain(s) << " AT " << time << endl;
			}
		
			assert (beginTimes[ind] == 0);
			
			beginTimes[ind] = clock();	
		}
	}
	
	inline void endState(const int s){
		
		if (!inState(s)){
			cout << explain(s) << endl;
			assert (false);
		}
	
		state = UtilLongSet::remove(state, s);
		stateChange(-1, s);
		
		const long ind = s; //UtilLongSet::getIndex(s);
		
		if (ind != -1 && beginTimes.size() > ind){
		
			assert (beginTimes[ind] != 0);
			
			if (verbose){
			
				clock_t earlier = beginTimes[ind];
				clock_t later = clock();
				
				//const string duration = UtilString::msToString( UtilMath::makeMultipleOf(((long)(later-earlier))/(CLOCKS_PER_SEC/1000.0), 1000) );
				
				const double secs = ((long)(later-earlier))/(1.0*CLOCKS_PER_SEC);
				
				if (s == PHASE_ROOT){
				
					secsPhase1 = secs;
				}
				
				if (s == PHASE_GETSAMPLE){
				
					secsPhase1and2 = secs;
				}
				
				const string duration = to_string( secs );
				
				const string time = UtilTime::dateTimeString(UtilTime::dateTimeNow() );
			
				cout << verboseprefix[0] << table->getTableName() << "#" << scan << " END " << explain(s) << " AT " << time << " DURATION " << duration << " secs" << endl;
			}
			
			beginTimes[ind] = 0;
		}
		
		
	}
	
	inline bool inState(const int s) const{
	
		return UtilLongSet::contains(state, s);
	}
	
	inline void setExact(bool b){
	
		exactJoinSize = b;
		
		for (auto it = children.cbegin(); it != children.cend(); ++it)
			(*it)->setExact(b);
	}
	

	inline bool isRoot() const{
	
		return parent == NULL;
	}
	
	inline const string toString(int tabs = 0) const{
	
	
		stringstream s;
		
		for (int i = 0; i < tabs; ++i)
			s << "   ";
		
		s << "============================" << table->getTableName() << endl;
		
		if (semijoin){
		
			for (int i = 0; i < tabs; ++i)
				s << "   ";
				
			s << "semijoin" << endl;
		}
		
		if (fk){
		
			for (int i = 0; i < tabs; ++i)
				s << "   ";
				
			s << "fk" << endl;
		}
		
		
		s << rel->toString(tabs);
		
		for (int i = 0; i < tabs; ++i)
			s << "   ";
		
		s << "============================" << endl;
		
		for (int i = 0; i < children.size(); ++i){
		
			s << children[i]->toString(tabs+1);
		}
		
		return s.str();
	
	}
	
	
	inline void setHashings(shared_ptr< unordered_map<string, long> > h){
	
		for (int j = 0; j < children.size(); j++)
			children.at(j)->setHashings(h);
		
		this->hashings = h;
		rel->setHashings(h);
	}
	
	inline void setParent(JoinTree* parent){
	
		this->parent = parent;
		
		setHashings(parent->hashings);
		
		if (parent){
			
			rel->setParent(parent->rel->getMainTuple() );
			parent->rel->addChild(rel->getTargetTuple() );
		}
		
		stringstream s;
		
		for (JoinTree* p = parent; p != NULL; p = p->parent){
			
			s << "\t";
		}
		
		verboseprefix.clear();
		verboseprefix.push_back( s.str() );
		
	}
	
	inline void addChild(shared_ptr<JoinTree> j){
	
		children.push_back(j);
		j->setParent(this);	
	}
	
	
	inline long getJoinSize(){
	
		
		vector<JoinTree*> queue = {this};
		vector<JoinTree*> newQueue = {};
		
		shared_ptr<DatabaseTable> ret = table;
		
		while (queue.size() > 0){
		
			for (auto it1 = queue.cbegin(); it1 != queue.cend(); ++it1){
		
				JoinTree* t = (*it1);
				
				const vector<shared_ptr<JoinTree>>& children = t->getChildren();
				
				for (auto it2 = children.cbegin(); it2 != children.cend(); ++it2){
			
					JoinTree& c = *(*it2);
					
					
					if (newQueue.size() == 0 && (it2+1) == children.cend() )
						return UtilTables::getJoinSize( *(ret), *(c.table) );
					
					ret = UtilTables::getJoin( *(ret), *(c.table), -1, -1, -1, -1, true);
					
					
					
					newQueue.push_back(&c);
				}
			}
			
			queue = newQueue;
			newQueue.clear();
		}
		
		return ret->getSize();
	}
	
	inline shared_ptr<DatabaseTable> getJoin(){
	
		
		vector<JoinTree*> queue = {this};
		vector<JoinTree*> newQueue = {};
		
		shared_ptr<DatabaseTable> ret = table;
		
		while (queue.size() > 0){
		
			for (auto it1 = queue.cbegin(); it1 != queue.cend(); ++it1){
		
				JoinTree* t = (*it1);
				
				const vector<shared_ptr<JoinTree>>& children = t->getChildren();
				
				for (auto it2 = children.cbegin(); it2 != children.cend(); ++it2){
			
					JoinTree& c = *(*it2);
					
					ret = UtilTables::getJoin( *(ret), *(c.table), ret->getSize(), ret->getSize()*64, c.table->getSize(), c.table->getSize()*64, true);
					
					ret->init();
					ret->buildIndices();
					
					newQueue.push_back(&c);
				}
			}
			
			queue = newQueue;
			newQueue.clear();
		}
		
		return ret;
	}
	
	inline void load(JoinParameters& p){
	
		p.addTable(table->getSize(), rel->getMainTuple().size(), rel->getTargetTuple().size(), rel->getNovelTuple().size() );
		
		for (int j = 0; j < children.size(); j++)
			children.at(j)->load(p);
		
	}
	
	inline string getSecondsString(const string s, double secs) const{
	
		stringstream ss;
		
		ss << to_string(secs) << " seconds";
		
		ss << " for " << s;
		
		if (secs >= 3600*24){
		
			ss << " ( " << to_string(secs/(24*3600.0)) << " days )";
			
		} else if (secs >= 3600){
		
			ss << " ( " << to_string(secs/3600.00) << " hours )";
			
		} else if (secs >= 60){
		
			ss << " ( " << to_string(secs/60.0) << " minutes )";
		}
		
		return ss.str();
	}
	
	
	inline shared_ptr<WeightedRowSample> getSample(JoinParameters& p, const string seed = "abc"){
		
		beginState(PHASE_GETSAMPLE);
		
		p.update();
		
		load(p);
		
		p.print();
		
		cout << toString() << endl;
		
		shared_ptr<WeightedRowSample> sample;
		
		//if (p.buildIndex)
		//	p.maxSamples = 1000*1000;
		
		
		
		for (int i = 0; (sample == nullptr) || (sample->getSize() < p.sampleSize); i++){
		
			shared_ptr<WeightedRowSample> samp1 = getSample1(p, !p.hash, seed+"_"+to_string(++randomId));
			
			cout << "join size " << to_string( ((long) joinSize) ) << endl;
			
			if (samp1->getSize() > 0){
			
				if (sample){
				
					sample->push_back(*samp1);
		
				} else {
			
					sample = samp1;
				}
			
			}
			
			if (sample){
				cout << sample->cols << endl;
				cout << sample->toString(5) << endl;
			} else {
			
				cout << samp1->cols << endl;
				cout << samp1->toString(5) << endl;
			}
			
		}
		
		endState(PHASE_GETSAMPLE);
		
		cout << endl << endl;
		cout << endl << endl;
		
		if (joinSizeEst != joinSize){
			
			cout << " superset join size = " << to_string( ((long) joinSize) ) << endl;
			cout << "estimated join size ~ " << to_string( ((long) joinSizeEst) ) << endl;
		} else {
		
			cout << "exact join size = " << to_string( ((long) joinSize) ) << endl;
		}
		
		cout << endl;
		
		cout << getSecondsString("join weights", secsPhase1) << endl;
		
		cout << getSecondsString("sample extensions", secsPhase1and2-secsPhase1) << endl;
		
		cout << getSecondsString("both together", secsPhase1and2) << endl;
		
		
		cout << endl << endl;
		cout << endl << endl;
		//if (sample)
		//	cout << sample->toString(5) << endl;
		
		return sample;
	}
	
	
	inline void aprioriShrink(long maxEntries){
	
		const long tableSize = table->getSize();
		const Tuple& t = rel->getTargetTuple();
		
		for (int i = 0; i < t.size(); ++i){
			
			const long estSize = t.size() == 1 ? tableSize : exp(log(maxEntries)*t.size() );
			const long v = t.size() == 1 ? maxEntries : exp(log(maxEntries)/t.size() );
			
			const long oldval = hashings->count(t[i]) == 0 ? std::numeric_limits<long>::max() : (*hashings)[ t[i] ];
			const long newval = estSize <= maxEntries ? std::numeric_limits<long>::max() : v;
			
			
			if (newval < oldval){
			
				(*hashings)[ t[i] ] = newval ;
				
				cout << t[i] << " => " << (*hashings)[ t[i] ] << endl;
			}
			
		}
		
		for (int j = 0; j < children.size(); j++)
			children[j]->aprioriShrink(maxEntries);
	}
	
	
	
	
	std::shared_ptr<WeightedRowSample> pop;
	
	
	inline bool hasIndex() const {
		
		if (rel->getJoinTuple().size() != rel->getTargetTuple().size() )
			return false;
			
		//const vector<int>& jointInParent = rel->getMap(JoinColumn::JOIN, JoinColumn::PARENT);
		
		const vector<int>& jointInMain = rel->getMap(JoinColumn::JOIN, JoinColumn::MAIN);
	
		if (jointInMain.size() > 0){
		
			if (jointInMain.size() > 1)
				return false;
				
			if (!table->hasIndex(jointInMain.at(0)))
				return false;
		}
		
		/*
		for (int j = 0; j < children.size(); j++)
			if (!children.at(j)->f)
				return false;*/
				
		return true;
	}
	
	inline bool hasChildrenWithIndex() const {
	
		for (int j = 0; j < children.size(); j++)
			if (!children.at(j)->hasIndex())
				return false;
				
		return true;
	}
	
	
	inline shared_ptr<WeightedRowSample> getSample1(const JoinParameters& p,  bool exact=true, const string seed="abc123"){
		
		cout << "*** getSample1(sampleSize=" << p.maxSamples << ", maxEntries=" << p.m << ", exact=" << exact << ") ***" << endl;
		
		setExact(exact);
		setLimits(fk ? 10000 : exact ? 1000000000 : p.m, p.maxSamples, p.batchSize);
		
		
		
		const bool buildpop = hasChildrenWithIndex() || p.sampleSize == 0;
		
		const bool indexpop = hasChildrenWithIndex() && table->hasIndices();
		
		if (inState(PHASE_ZERO))
			endState(PHASE_ZERO);
		
		joinSize = 0;
		
		shared_ptr<WeightedRowSample> sample;
		
		const bool twoscan = p.sampler == WeightedSamplers::OFFLINE_WR;
		
		if ( (!buildpop) && (!twoscan) ){
		
			cout << "X1" << endl;
			sample = std::make_shared<WeightedRowSample>(rel->getMainTuple(), seed+"_"+to_string(++randomId), p.sampleSize > 0 ? p.maxSamples: table->getSize() );
			
			sample->initOnlineAuto(p.sampler);
			
		}
		
		
		if (pop){
		
			// nothing to do here
			
			cout << "X2" << endl;
		
		} else {
		
			cout << "X3" << endl;
		
			WeightedRowSample* smp = NULL;
		
			if (buildpop){
			
				if (!indexpop){
			
					pop = std::make_shared<WeightedRowSample>(rel->getMainTuple(), seed+"_"+to_string(++randomId), table->getSize(), true);
					pop->initPop();
					smp = &(*pop);
				}
				
			} else {
				
				if (sample)
					smp = &(*sample);
					
				cout << "X4" << endl;
			}
			
			beginState(PHASE_ROOT);
			
			root1(smp);
			
			if (indexpop){
		
				shared_ptr<RowHistogram<Weight>> empty = std::make_shared<RowHistogram<Weight>>(Tuple());
				auto joiner = std::make_shared<Joiner>(empty, rel, fk);
	
				for (int j = 0; j < children.size(); j++)
					joiner->addChild( children.at(j)->batchesJoinWeights, children.at(j)->getSemiJoin(), children.at(j)->getFK() );
			
				table->getIndex()->setJoiner(joiner);
				
			} else if (pop){
			
				cout << "X5" << endl;
				pop->finalise();
			
			} else {
			
				cout << "X6" << endl;
				if (sample)
					sample->withReplacement(p.maxSamples);
				
				//cout << "JOIN SIZE " << joinSize << " (" << ((long) joinSize) << ")" << endl;
				//cout << "SAMPLE " << (*sample->begin()) << endl;
			}
			
			endState(PHASE_ROOT);
		}
		
		
		vector<JoinTree*> extensionList;
		
		vector<JoinTree*> queue = {this};
		vector<JoinTree*> newQueue = {};
		
		while (queue.size() > 0){
		
			for (auto it1 = queue.cbegin(); it1 != queue.cend(); ++it1){
		
				JoinTree* t = (*it1);
				
				const vector<shared_ptr<JoinTree>>& children = t->getChildren();
				
				if (t->getSemiJoin() )
					continue;
				
				for (auto it2 = children.cbegin(); it2 != children.cend(); ++it2){
			
					JoinTree& c = *(*it2);
					
					if (c.getSemiJoin() )
						continue;
					
					extensionList.push_back(&c);
					
					newQueue.push_back(&c);
				
					
				}
			}
			
			queue = newQueue;
			newQueue.clear();
		}
		
		
		
		if (indexpop || buildpop){
		
			if (!indexpop){
				
				cout << "X7" << endl;
		
				cout << "pop " << pop->cols << endl;
				cout << pop->toString(5) << endl;
			
				if (p.sampleSize == 0)
					sample = pop;
				
			}
			/* else
				sample = pop->getWithReplacementSampleFromPopulation(p.maxSamples, seed+"_"+to_string(++randomId) ); */
				
			// STREAMLINED INDEX-SUPPORTED SAMPLING CODE
			if (p.sampleSize > 0){
			
			
				RandomNumbers rnd(seed+"_"+to_string(++randomId));
				
				vector<string> j_st; // j_st[1] is the column with identity 1
				
				map<string, int> st_j;  // inverse of j_st using a hash map
				map<string, vector<int> > st_ks; // st[ks][0] the first k-id of a table that has column st
				
				vector<JoinTree*> k_t; // k_t[0] = pointer to table with k-id = 0
				
				k_t.push_back(this);
				for (auto it = extensionList.begin(); it != extensionList.end(); ++it){
				
					//if ( (*it)->getSemiJoin() )
						//continue;
				
					k_t.push_back(*it);
				}
				
				for (int k = 0; k < k_t.size(); ++k){
					
					JoinTree& t = *k_t.at(k);
					
					for (int col = 0; col < t.table->getColNames().size(); ++col){
					
						const string& st = t.table->getColNames()[col];
					
						
						if (st_j.count(st) == 0){ // unseen column
						
							const int j = j_st.size();
						
							st_j.insert( {st, j } );
							j_st.push_back(st);
							
							cout << "j_st[" << j << "] = " << st << endl;
							
							vector<int> v;
							
							st_ks.insert( {st, v } );
						}
						
						st_ks.at(st).push_back(k); // add to tables for column
					}
				}
				
				
				// COLUMN IDS OF EACH TABLE
				vector<vector<int>> k_js;
				
				for (int k = 0; k < k_t.size(); ++k){
					
					JoinTree& t = *k_t.at(k);
					
					vector<int> js;
					
					for (int col = 0; col < t.table->getColNames().size(); ++col){
					
						const string& st = t.table->getColNames()[col];
						
						assert (st_j.count(st) > 0);
						
						const int j = st_j[st];
						
						js.push_back(j);
					}
					
					k_js.push_back(js);
				}
				
				// TABLE IDS FOR EACH COLUMN ID
				vector< vector<int> > j_ks;
				vector< vector<int> > j_kcols;
				
				for (auto it = j_st.begin(); it != j_st.end(); ++it){
				
					const string& st = (*it);
					
					assert (st_ks.count(st) > 0);
					
					vector<int>& ks = st_ks[st];
					
					std::sort( ks.begin(), st_ks[st].end() );
					
					const int j = j_ks.size();
					
					j_ks.push_back( ks );
					
					vector<int> kcols;
					
					for (int l = 0; l < ks.size(); ++l){
					
						const int k = ks.at(l);
						
						const vector<int>& js = k_js.at(k);
						
						int kcol = -1;
						
						for (int col = 0; col < js.size(); ++col){
						
							if (js.at(col) == j){
								kcol = col;
								break;
							}		
						}
						
						assert (kcol != -1);
						kcols.push_back(kcol);
					}
					
					j_kcols.push_back( kcols);
				}
				
				
				vector<vector<int>> k_projdown_ks;
				vector<vector<int>> k_projdown_kcols;
				
				vector<vector<int>> k_projup;
				
				vector<vector<int>> k_across_k;
				vector<vector<int>> k_across_col;
				vector<vector<int>> k_across_kcol;
				
				
				// RELATIONSHIPS BETWEEN TABLES
				
				for (int k = 0; k < k_t.size(); ++k){
					
					JoinTree& t = *k_t.at(k);
					
					cout << t.table->getColNames() << endl;
					vector<int> projup;
					
					vector<int>& js = k_js.at(k);
					
					for (int col2 = 0; col2 < t.rel->getTargetTuple().size(); ++col2){
					
						const string& st = t.rel->getTargetTuple()[col2];
					
						assert (st_j.count(st) > 0);
						
						const int j = st_j[st];
						const int col = std::find(js.begin(), js.end(), j)-js.begin();
						
						projup.push_back(col);
						
						const vector<int> ks = j_ks.at(j);
						const vector<int> kcols = j_kcols.at(j);
						
						if (k > 0){
						
							for (int l = 0; l < ks.size(); ++l){
					
								const int kk = ks.at(l);
								const int kcol = kcols.at(l);
							
								if (kk < k){
									
									bool alreadyContained = false;
									
									for (int i = 0; i < k_projdown_ks.at(k-1).size(); ++i){
									
										if ( k_js.at( k_projdown_ks.at(k-1).at(i) ).at( k_projdown_kcols.at(k-1).at(i) ) == j){
											
											alreadyContained = true;
										} 
									}
									
									if (!alreadyContained){
									k_projdown_ks.at(k-1).push_back(kk);
									k_projdown_kcols.at(k-1).push_back(kcol);
									}
								}
							}
						}	
					}
					
					vector<int> across_col;
					vector<int> across_k;
					vector<int> across_kcol;
					
					for (int col = 0; col < js.size(); ++col){
					
						const int j = js.at(col);
						const vector<int> ks = j_ks.at(j);
						
						const vector<int> kcols = j_kcols.at(j);
						
						for (int kki = 0; kki < ks.size(); ++kki){
						
							int kk = ks.at(kki);
							int kcol = kcols.at(kki);
						
							if (kk == k)
								continue;
								
							if (kk < k){
							
								across_col.push_back(col);
								across_k.push_back(kk);
								across_kcol.push_back(kcol);
								
								assert (k_js.at(k).at(col) == k_js.at(kk).at(kcol) );
							}
						}
					}
					
					k_across_col.push_back(across_col);
					k_across_k.push_back(across_k);
					k_across_kcol.push_back(across_kcol);
					
					k_projup.push_back(projup);
					
					k_projdown_ks.push_back( vector<int>() );
					k_projdown_kcols.push_back( vector<int>() );
				}
				
				vector<MainMemoryIndex*> k_ind;
				
				for (int k = 0; k < k_t.size(); ++k){
					
					JoinTree& t = *k_t.at(k);
					
					if ( k == 0){
					
						k_ind.push_back(NULL);
					
					} else {
					
						assert (k_projup.at(k).size() == 1);
						
						assert (t.table->hasIndex(k_projup.at(k).at(0)));
						
						MainMemoryIndex* ind = t.table->getIndex( k_projup.at(k).at(0) );
						
						if (!ind->hasJoiner() ){
			
							shared_ptr<RowHistogram<Weight>> empty = std::make_shared<RowHistogram<Weight>>(Tuple());
							auto joiner = std::make_shared<Joiner>(empty, t.rel, fk);
						
							for (int j = 0; j < t.children.size(); j++)
								joiner->addChild( t.children.at(j)->batchesJoinWeights );
							
							ind->setJoiner(joiner);
						}
					
						k_ind.push_back(ind);
					}
				}
				
				vector<vector<int>> k_sels;
				
				std::bitset<1000> bs;
				
				Tuple smpCols;
				
				for (int k = 0; k < k_t.size(); ++k){
				
					const vector<int>& js = k_js.at(k);
					
					vector<int> sels;
				
					for (int col = 0; col < js.size(); ++col){
					
						const int j = js.at(col);
						
						if (!bs.test(j)){
						
							smpCols.push_back( j_st.at(j) );
							sels.push_back(col);
							bs.set(j);
						}
					}
					
					k_sels.push_back(sels);
				}
				
				cout << "sample " << smpCols << endl;
				
				for (int k = 0; k < k_t.size(); ++k){
				
					cout << "===========" << endl;
					cout << "table " << UtilString::toString(k_js.at(k), j_st) << endl;
					cout << "projup " << UtilString::toString(k_projup.at(k), j_st) << endl;
					
					const vector<int>& projdown_ks = k_projdown_ks.at(k);
					const vector<int>& projdown_kcols = k_projdown_kcols.at(k);
					
					cout << "projdown";
					
					for (int i = 0; i < projdown_kcols.size(); ++i){
						
						const int projdown_kcol = projdown_kcols.at(i);
						const int projdown_k = projdown_ks.at(i) ;
					
						cout << " " << j_st.at( k_js.at(projdown_k).at(projdown_kcol) );
					}
					
					cout << endl;
				
				}
				
				
				
				beginState(PHASE_EXTEND);
				
				long validsamples = 0;
				
				shared_ptr<WeightedRowSample> sample = std::make_shared<WeightedRowSample>(smpCols, "", p.sampleSize);
				
				sample->initPop();
				
				vector<const Row* > k_row;
				
				k_row.reserve(k_t.size() );
				
				long invalid = 0;
				
				MainMemoryIndex* ind0 = table->getIndex();
				
				while (validsamples < p.sampleSize){
				
					k_row.clear();
				
					for (int k = 0; k < k_t.size(); ++k){
					
						if (k == 0){
						
							if (ind0){
							
								k_row.push_back(ind0->get( ind0->begin(), ind0->end(), rnd.randomDouble() ));
							
							} else {
								
								k_row.push_back(&pop->getRow(rnd.randomDouble() ));
							}
							
							continue;
						}
					
						//JoinTree& t = *(k_t.at(k));
						
						//assert (k_projup.at(k).size() == 1);
						
					
						MainMemoryIndex* ind = k_ind.at(k);
						
						const int kk = k_projdown_ks.at(k-1).at(0);
						const int kcol = k_projdown_kcols.at(k-1).at(0);
						
						if (ind->begin( *k_row.at(kk), kcol ) == NULL){
						
							k_row.clear();
							break;
						}
						
						const Row* ext = ind->get(*k_row.at(kk), kcol, rnd.randomDouble() );
						
						const vector<int>& across_col = k_across_col.at(k);
						const vector<int>& across_k = k_across_k.at(k);
						const vector<int>& across_kcol = k_across_kcol.at(k);
						
						for (int l = 0; l < across_col.size(); ++l){
						
							const int col = across_col.at(l);
							const int kk = across_k.at(l);
							const int colk = across_kcol.at(l);
							
							assert (kk < k);
							
							if (!ext->equals(col, *k_row.at(kk), colk) ){
							
								k_row.clear();
								break;
							}
						}
						
						if (k_row.size() == 0)
							break;
						
						k_row.push_back(ext);
					}
					
					if (k_row.size() == 0){
					
						++invalid;
						if (invalid % 1000000 == 0){
						
							const long double acceptanceRate = validsamples == 0 ? 0 : ((long double) validsamples) / (validsamples+invalid);
							
							const long double joinSize = this->joinSize;
							
							if (acceptanceRate == 0){
							
								const long double acceptanceRateUpperBound = 1.0 / p.sampleSize;
							
								cout << "collected " << validsamples << " dropped " << invalid << " join size estimate between " << to_string( ((long) (joinSize*acceptanceRateUpperBound+0.5)) ) << " and " << to_string( ((long) joinSize) ) << endl;	
								
							} else {
							
								cout << "collected " << validsamples << " dropped " << invalid << " join size estimate " << to_string((long) (joinSize*acceptanceRate+0.5) ) << endl;	
							}
						}
					
						continue;
					}
				
					assert (k_row.size() == k_t.size() );
					
					//cout << validsamples << endl;
				
					++validsamples;
					
					sample->push_back(k_row, k_sels, 1);
				}
				
				
				{
					const long double acceptanceRate = validsamples == 0 ? 0 : ((long double) validsamples) / (validsamples+invalid);
							
					const long double joinSize = this->joinSize;
					
					joinSizeEst = joinSize*acceptanceRate+0.5;
					
					if (acceptanceRate == 0){
							
						const long double acceptanceRateUpperBound = 1.0 / p.sampleSize;
							
							cout << "collected " << validsamples << " dropped " << invalid << " join size estimate between " << to_string( ((long) (joinSize*acceptanceRateUpperBound+0.5)) ) << " and " << to_string( ((long) joinSize) ) << endl;	
								
						} else {
							
							cout << "collected " << validsamples << " dropped " << invalid << " join size estimate " << to_string((long) (joinSize*acceptanceRate+0.5) ) << endl;	
					}
				
				
				}
			
				sample->finalise();
			
				endState(PHASE_EXTEND);
				
				
				
				return sample;
			
			}	
			
			
			//cout << sample->toString(5) << endl;
		} else {
		
			cout << "X8" << endl;
		
			if (twoscan){
			
				cout << "X9" << endl;
				
				shared_ptr<WeightedRowSample> sample2 = std::make_shared<WeightedRowSample>(rel->getMainTuple(), seed+"_"+to_string(++randomId), p.maxSamples);
				
				sample2->initOfflineWR(joinSize);
				
				beginState(PHASE_ROOT);
				
				root1( &(*sample2) );
				
				endState(PHASE_ROOT);
				
				sample2->withReplacement(joinSize);
				
				sample = sample2;
			}
			
		}
		
		cout << "X10" << endl;
		
		cout << "MAIN TABLE SAMPLE" << endl;
		cout << sample->toString(5) << endl;
		cout << "/MAIN TABLE SAMPLE" << endl;
		
		/*
		if (exact)
			cout << "EXACT: ";
		
		cout << "SUPERSET SAMPLE SIZE " << sample->getSize() << endl;*/
		
		beginState(PHASE_EXTEND);
		
		for (auto it2 = extensionList.cbegin(); it2 != extensionList.cend(); ++it2){
		
			if (sample->getSize() == 0){
			
				endState(PHASE_EXTEND);
				return sample;
			}
			
			JoinTree& c = *(*it2);
			
			sample = c.extend(sample, seed+"_"+to_string(++randomId));	
		}
		
		endState(PHASE_EXTEND);
		
		const long validsamples = sample->getSize();
		const long invalid = sample->getInvalidSize();
		
		const long double acceptanceRate = validsamples == 0 ? 0 : ((long double) validsamples) / (validsamples+invalid);
							
		const long double joinSize = sample->getPopulationWeight();
		
		joinSizeEst = joinSize*acceptanceRate+0.5;
		
		
		//cout << sample->toString(5) << endl;
		
		return sample;
	}
	
	typedef shared_ptr<WeightedRowSample> WeightedRowSamplePtr;
	
	
	inline shared_ptr<WeightedRowSample> extend(shared_ptr<WeightedRowSample> sample, const string seed="abc123"){
		
		assert (getSemiJoin() == false);
		
		if (sample->getSize() == 0 && !sample->isPopulation()){
		
			cout << "X11" << endl;
		
			return sample;
		}
		
		beginState(PHASE_EXTEND);
		
		shared_ptr<WeightedRowSample> ret = sample;
		
		cout << "X12" << endl;
		
		if (batchesJoinWeights->update(*hashings))
			rel->updateHashings();
		
		
		
		
		if (hasIndex() ){
		
			cout << "X13" << endl;
		
			const vector<int>& jointInTable = rel->getMap(JoinColumn::JOIN, JoinColumn::MAIN);
		
			MainMemoryIndex* ind = table->getIndex(jointInTable[0]);
		
			if (!ind->hasJoiner() ){
			
				shared_ptr<RowHistogram<Weight>> empty = std::make_shared<RowHistogram<Weight>>(Tuple());
				auto joiner = std::make_shared<Joiner>(empty, rel, fk);
		
				for (int j = 0; j < children.size(); j++)
					joiner->addChild( children.at(j)->batchesJoinWeights, children.at(j)->getSemiJoin(), children.at(j)->getFK() );
				
				ind->setJoiner(joiner);
			}
			
			ret = Extensions(*sample, *rel, *batchesJoinWeights, seed+"_"+to_string(++randomId), ind ).getSample();
			
		} else {
		
			cout << "X14" << endl;
			
			Extensions ext(*sample, *rel, *batchesJoinWeights, seed+"_"+to_string(++randomId) );
		
			cout << "X15" << endl;
		
			extend1(ext);
			
			cout << "X16" << endl;
			
			ret = ext.getSample();
			
			cout << "X17" << endl;
		}
		
		
	
		
		
		//cout << ret->toString(5) << endl;
		
		endState(PHASE_EXTEND);
		
		//beginState(PHASE_EXTENDCHILDREN);
		
		//for (int j = 0; j < children.size(); j++)
		//	ret = children.at(j)->extend(ret, seed);
		
		//endState(PHASE_EXTENDCHILDREN);
	
		return ret;		
	}
	
	
	
	inline void stateChange(int state_begin, int state_end){
	
		if (state_begin == PHASE_ZERO){			
		}
	
		if (state_begin == PHASE_EXTEND){
		
			//requestedJoinWeights.reset();
			
			//if (!complete)
			//	scanJoinWeights.reset();
		}
	}
	
	inline void root2(Rows* buf, WeightedRowSample* sample = NULL){
	
		
		for (auto it = children.begin(); it != children.end(); ++it)
			(*it)->joinweights1(buf);
		
		shared_ptr<RowHistogram<Weight>> empty = std::make_shared<RowHistogram<Weight>>(Tuple());
		Joiner joiner(empty, rel, fk);
	
		for (int j = 0; j < children.size(); j++){
		
			//children.at(j)->batchesJoinWeights->print(5);
			
			if (children.at(j)->batchesJoinWeights->update(*hashings))
				rel->updateHashings();
			
			//cout << "root2 add child " << children.at(j)->batchesJoinWeights->getSize() << endl;
			
			joiner.addChild( children.at(j)->batchesJoinWeights, children.at(j)->getSemiJoin(), children.at(j)->getFK() );
		}
		
		for (auto it = buf->begin(); it != buf->end(); ++it){
		
			const Row& r = (*it);
		
			const Weight f = joiner.getJoinSize(r);
			
			//cout << r << " " << f << endl;
			
			if (sample)
				sample->add(r, f);
				
			joinSize += f;
		}
	}
	
	inline void root1(WeightedRowSample* sample = NULL){
	
		joinSize = 0;
	
		beginState(STATE_STREAMPASS);
	
		Rows* buf = &buffer;
		
		buf->clear();
		
		const long bufsize = maxBatchSize;
	
		for (auto it = table->begin(); it != table->end(); ++it){
		
			const Row& row = (*it);
			
			buf->push_back(row);
			
			//cout << "ROW " << row << endl;
			
			if (buf->size() == bufsize){
				
				root2(buf, sample);
				buf->clear();
			}
		}
		
		if (buf->size() > 0){
				
			root2(buf, sample);
			buf->clear();
		}
		
		
		
		
		cout << "root1, joinweight = " << joinSize << " joinsize "<< ( (long) joinSize) << endl;
		
		endState(STATE_STREAMPASS);
	}
	
	
	// updates hashings
	inline void joinweights2(Rows* buf, RowHashMap<bool>* obligatory = NULL){

	

		assert (buf != NULL);
	
		// GET WEIGHTS FROM CHILDREN
		for (auto it = children.begin(); it != children.end(); ++it)
			(*it)->joinweights1(buf);

		shared_ptr<RowHistogram<Weight>> batchJoinWeights = std::make_shared<RowHistogram<Weight>>(rel->getTargetTuple());
		
		
		batchJoinWeights->reserve(buf->size() );
	
		// GET REQUESTED JOIN WEIGHTS	
		Joiner joiner1(batchJoinWeights, rel, fk);
		
		if (batchJoinWeights->update(*hashings))
			rel->updateHashings();
	
		for (int j = 0; j < children.size(); j++){
	
			if (children.at(j)->batchesJoinWeights->update(*hashings))
				rel->updateHashings();
		
			joiner1.addChild( children.at(j)->batchesJoinWeights, children.at(j)->getSemiJoin(), children.at(j)->getFK() );
		}
		
		{
			const vector<int>& jointInTable = rel->getMap(JoinColumn::JOIN, JoinColumn::MAIN);
			const vector<long>& jointInTableHashings = rel->getHashings(JoinColumn::JOIN, JoinColumn::MAIN);
			
			//const vector<int> jointInTableHashings(jointInTable.size(), 0);
			
			for (auto it = buf->begin(); it != buf->end(); ++it){
		
				const Row& row = (*it);
		
				if ( (obligatory == NULL) || obligatory->hasKey(row, jointInTable, jointInTableHashings))
					joiner1.join(row);
			}
		}
		
		/*
		cout << "batchesJoinWeights" << endl;
		batchesJoinWeights->print(5);
		
		cout << "batchJoinWeights" << endl;
		batchJoinWeights->print(5);*/
		
		// STORE WITH OLDER REQUESTED JOIN WEIGHTS WHEN POSSIBLE
		
		const bool index = hasIndex();
		
		//const bool nonEmpty = batchesJoinWeights->getSize() > 0;
		//const bool noSpace = (maxEntries) > 0 && (batchesJoinWeights->getSize()+batchJoinWeights->getSize() < maxEntries);
		
		
		//cout << "preBatchesJoinWeights" << endl;
		//batchesJoinWeights->print(5);
				
		//cout << "batchJoinWeights" << endl;
		//batchJoinWeights->print(5);
		
		
		if (batchesJoinWeights->getSize() == 0){
		
			batchesJoinWeights.reset();
			batchesJoinWeights = std::move( batchJoinWeights );
			batchJoinWeights.reset();
			
		} else {
		
			batchesJoinWeights->push_back(*batchJoinWeights);
			batchJoinWeights.reset();
		}
		
		/*
		if ( nonEmpty && (index || noSpace) ){
			
			if (batchJoinWeights->getSize() > 0){
			
				//cout << batchesJoinWeights->getSize() << " " << batchJoinWeights->getSize() << endl;
			
				//batchesJoinWeights->set(*hashings);
				
				if (batchesJoinWeights->update(*hashings))
					rel->updateHashings();
				
				batchesJoinWeights->push_back_overwriting(*batchJoinWeights);
				
				
				//cout << "//" << batchesJoinWeights->getSize() << " " << batchJoinWeights->getSize() << endl;
			}
			
			batchJoinWeights.reset();
		
		} else {
			
			batchesJoinWeights.reset();
			batchesJoinWeights = std::move( batchJoinWeights );
			batchJoinWeights.reset();
		}*/
		
		//cout << "postBatchesJoinWeights" << endl;
		//batchesJoinWeights->print(5);
				
		
		if (!index){
		
			if (scanJoinWeights){
		
				Joiner joiner2(scanJoinWeights, rel, fk);
				
				if (scanJoinWeights->update(*hashings))
					rel->updateHashings();
		
				for (int j = 0; j < children.size(); j++){
				
					if (children.at(j)->batchesJoinWeights->update(*hashings))
						rel->updateHashings();
				
					joiner2.addChild( children.at(j)->batchesJoinWeights, children.at(j)->getSemiJoin(), children.at(j)->getFK() );		
				}
			
				for (auto it = buf->begin(); it != buf->end(); ++it){
		
					const Row& row = (*it);
					
					joiner2.join(row);
					
					/* NO SHRINK 
					if ( joiner2.shrinkToSize(maxEntries) ){
				
						for (int j = 0; j < children.size(); j++){
							
							if (children.at(j)->batchesJoinWeights->update(*hashings))
								rel->updateHashings();
						}
						
						rel->updateHashings();
					}*/
				}
			}
		}
	
	}
	
	inline void joinweights1(const Rows* batch){
	
		
	
		//const Tuple& parentAttributes = rel->getParentTuple();
	
		// ALREADY PREPARED ALL JOIN WEIGHTS
		if (complete){
		
			
			//endState(PHASE_JOINWEIGHTS);
			
			/*
			beginState(STATE_INDEXPASS);
			
			cout << "SKIP" << endl;
			
			endState(STATE_INDEXPASS);
			*/
			
		
			return;
		}
		//inState(PHASE_JOINWEIGHTS)
		
		//cout << "NOSKIP" << endl;
		
		if (!inState(PHASE_JOINWEIGHTS))
			beginState(PHASE_JOINWEIGHTS);
		
		
		if (fk){
		
			cout << table->getColNames() << " FK JOINWEIGHTS" << endl;
			complete = true;
	
			return;
		}
		
		
		
		scanJoinWeights = std::make_shared<RowHistogram<Weight>>(rel->getTargetTuple() );
		scanJoinWeights->update(*hashings);
		
		batchesJoinWeights->update(*hashings);

					
		if (maxEntries < table->getSize() )
			scanJoinWeights->reserve(maxEntries);
		else
			scanJoinWeights->reserve(table->getSize());
		
		Rows* buf = &buffer;
		
		const bool index = hasIndex();
		
		const long bufsize = maxBatchSize; //index ? 1000000 : maxEntries;
		
		// IF THERE IS AN INDEX ONE CAN USE
		
		
		const vector<int>& jointInTable = rel->getMap(JoinColumn::JOIN, JoinColumn::MAIN);
		const vector<int>& jointInParent = rel->getMap(JoinColumn::JOIN, JoinColumn::PARENT);
		
		const vector<long>& jointInParentHashings = rel->getHashings(JoinColumn::JOIN, JoinColumn::PARENT);
		//vector<int> jointInParentHashings(jointInParent.size(), 0);
		
		//cout << rel->toString() << endl;
		
		//cout << "parent(" << rel->getParentTuple() << ") join(" << rel->getJoinTuple() << ") [" << UtilString::toString(jointInParent) << "] [" << UtilString::toString(jointInParentHashings)<< "]" << endl;
		
		if (index){
		
			beginState(STATE_INDEXPASS);
			
			obligatory.clear();
		
			// JUMP TO ROWS JOINING WITH BATCH
			
			MainMemoryIndex* ind = table->getIndex(jointInTable[0]); 
			
			assert (ind != NULL);
			
			for (auto pit = batch->begin(); pit != batch->end(); ++pit){
			
				const Row& parentrow = (*pit);
				
				
				
				if (obligatory.hasKey(parentrow, jointInParent, jointInParentHashings))
					continue;
				
				// ALREADY HAVE JOINWEIGHT FOR THIS ONE
				if (batchesJoinWeights->contains(parentrow, jointInParent) ){
				
					//cout << "contains parentrow " << parentrow << endl;
				
					//batchJoinWeights->set(parentrow, jointInParent, batchesJoinWeights->get(parentrow) );
					continue;
				}
				
				obligatory.insert(parentrow, jointInParent, jointInParentHashings, true);
				
				//cout << "parentrow " << parentrow << endl;
				
				auto begin = ind->begin( parentrow, jointInParent[0]);
				auto end = ind->end(parentrow, jointInParent[0]);
			
				assert (begin <= end);		
				
				for (auto it = begin; it != end; ++it){
				
					const Row& row = (*it);
					
					//cout << "parentrow " << parentrow << " childrow " << row << endl;
					
					buf->push_back(row);
					
					if (buf->size() == bufsize){
						joinweights2(buf);
						buf->clear();
					}
				}
			}
			
			if (buf->size() > 0){
				joinweights2(buf);
				buf->clear();
			}
			
			endState(STATE_INDEXPASS);
		
		} else {
		
			// PASS OVER TABLE
			
			beginState(STATE_STREAMPASS);
			
			obligatory.clear();
			
			{
				//const vector<int>& jointInParentHashings = rel->getHashings(JoinColumn::MAIN, JoinColumn::PARENT);
		
				for (auto pit = batch->begin(); pit != batch->end(); ++pit){
			
					const Row& parentrow = (*pit);
			
					if (!obligatory.hasKey(parentrow, jointInParent, jointInParentHashings))
						obligatory.insert(parentrow, jointInParent, jointInParentHashings, true);
				}
			}
			
			Timer t("scan "+table->getTableName()+" to get join weights" );
			
			t.start(10, UtilMath::maxVal<long>(1000000, table->getSize()/10), table->getSize() );
			
			for (auto it = table->begin(); it != table->end(); ++it){
			
				const Row& row = (*it);
				buf->push_back(row);
			
				t.tick();
			
				if (buf->size() == bufsize){
				
					joinweights2(buf, &obligatory);
					buf->clear();
				}
			}
			
			if (buf->size() > 0){
				
				joinweights2(buf, &obligatory);
				buf->clear();
			}
			
			t.end();
			
			obligatory.clear();
				
			completedFirstPass = true;
			
			//allJoinWeightsJoiner.reset();
			
			//if ( (!exactJoinSize) || (!scanJoinWeights->usesHashing()) ){
			
			if (true){
	
				assert (!complete);
	
				batchesJoinWeights = scanJoinWeights;
				scanJoinWeights.reset();
				complete = true;
				
				//if (batchesJoinWeights->update(*hashings))
				//	batchesJoinWeights->updateHashings();
				
				if (verbose)
					cout << verboseprefix[0] << "COMPLETED " << table->getTableName() << " SIZE " << batchesJoinWeights->getSize() << endl;
					
				batchesJoinWeights->lock();
			}
			
			endState(STATE_STREAMPASS);
		}
	}
	
	inline void extend2(const Rows* buffer, Extensions& ext){
	
		for (auto it = children.begin(); it != children.end(); ++it)
			(*it)->joinweights1(buffer);
		
		shared_ptr<RowHistogram<Weight>> empty = std::make_shared<RowHistogram<Weight>>(Tuple());
		Joiner joiner(empty, rel, fk);
		
		for (int j = 0; j < children.size(); j++){
		
			if (children.at(j)->batchesJoinWeights->update(*hashings))
				rel->updateHashings();
			
			joiner.addChild( children.at(j)->batchesJoinWeights, children.at(j)->getSemiJoin(), children.at(j)->getFK() );
		}
		
		Timer t("extend2");
		
		t.start(30, ceil(buffer->getSize()/10.0), buffer->getSize() );
		
		ext.addExtensions(*buffer, joiner.getChildren2() );
		
		/*
		for (auto it = buffer->begin(); it != buffer->end(); ++it){
			
			const Row& extRow = (*it);

			t.tick();
			
			ext.addExtension(extRow, joiner.getChildren() );
			
			
			//assert (map.hasKey(extRow, extProj, extHash) );
			
			//const W f = joiner.getJoinSize(extRow);
			
			//if (f == 0)
			//	continue;
			
			//shared_ptr<WeightedRowSample> smp = map.get(extRow, extProj, extHash);
			
			//assert (smp->maxSize > 0);
			
			//smp->add(extRow,f);	
		}*/
		
		t.end();
	}
	
	inline void extend1(Extensions& ext){
		
		
		for (int j = 0; j < children.size(); j++){
			
			if (children.at(j)->batchesJoinWeights->update(*hashings))
				rel->updateHashings();
		}
		
		const long bufsize = maxBatchSize;
		
		// IF THERE IS AN INDEX ONE CAN USE
		
		Rows* buf = &buffer;
		
		const vector<int>& jointInTable = rel->getMap(JoinColumn::JOIN, JoinColumn::MAIN);
		
		if (hasIndex() ){
		
			assert (false);
			/*
			MainMemoryIndex* ind = table->getIndex(jointInTable[0]);
		
			if (!ind->hasJoiner() ){
			
				auto joiner = std::make_shared<Joiner>(empty, rel);
		
				for (int j = 0; j < children.size(); j++)
					joiner->addChild( children.at(j)->batchesJoinWeights );
				
				ind->setJoiner(joiner);
			}
			
			ext.addExtensions(ind);
			*/
		
			/*
			beginState(STATE_INDEXPASS);
		
			// JUMP TO ROWS JOINING WITH BATCH
			
			auto& map = ext.getIndex();
			
			MainMemoryIndex* ind = table->getIndex(jointInTable[0]);
			
			for (auto pit = map.longMap.begin(); pit != map.longMap.end(); ++pit){
			
				auto begin = ind->begin(pit->first);
				auto end = ind->end(pit->first);
				
				assert(begin <= end);
				
				for (auto it = begin; it != end; ++it){
				
					const Row& extRow = (*it);
					
					buf->push_back(extRow);
					
					if (buf->size() == bufsize){
					
						extend2(buf, ext);
						buf->clear();
					}
				}
			}
			
			if (buf->size() == bufsize){
					
				extend2(buf, ext);
				buf->clear();
			}
			
			for (auto pit = map.stringMap.begin(); pit != map.stringMap.end(); ++pit){
			
				auto begin = ind->begin(pit->first);
				auto end = ind->end(pit->first);
				
				for (auto it = begin; it != end; ++it){
				
					const Row& extRow = (*it);
					
					buf->push_back(extRow);
					
					if (buf->size() == bufsize){
					
						extend2(buf, ext);
						buf->clear();
					}
				}
			}
			
			if (buf->size() > 0){
					
				
					
				extend2(buf, ext);
				buf->clear();
			}
			
			endState(STATE_INDEXPASS);*/
			
		
		} else {
		
			beginState(STATE_STREAMPASS);
		
			// PASS OVER TABLE
			
			Timer t("extend1");
			
			t.start(30, ceil(table->getSize()/100.0), table->getSize() );
			
			for (auto it = table->begin(); it != table->end(); ++it){
			
				const Row& extRow = (*it);
				
				t.tick();
				
				if (ext.has(extRow) ){
					
					buf->push_back(extRow);
					
					if (buf->size() == bufsize){
				
						extend2(buf, ext);
						buf->clear();
					}
				}
			}
			
			if (buf->size() > 0){
				
				extend2(buf, ext);
				buf->clear();
			}
			
			t.end();
			
			endState(STATE_STREAMPASS);
		}
	}
	
};

struct subgraph1{

	class iterator{
	
		private: 
		
			long cursor = 0;
			const subgraph1* data = NULL;
			
			std::set<int>::const_iterator it;
			std::set<int>::const_iterator end;
			
		public:

			iterator(){
			
			
			}
			
			inline int pos(int id) const{
	
				return id < 0 ? -id-1 : id;
			}
			
			iterator(const std::set<int>::const_iterator begin, const std::set<int>::const_iterator end) : it(begin), end(end){
			
				
			}
			
			iterator(const subgraph1* s) : data(s){
			
				while (cursor < (data->capacity() ) && data->count(getIdFromPos(cursor) ) == 0)
					++cursor;
			}
			
			inline int nodeIdFromPos(int x) const {
	
				return x;
			}
	
			inline int edgeIdFromPos(int x) const {
	
				return -x-1;
			}

			inline char operator()() const{
			
				return 'x';
			}
			
			inline const int operator*() const{
			
				if (data == NULL)
					return *it;
		
				return getIdFromPos(cursor);
			}
			
			inline bool operator ==(char b) const{
	
				assert (b == 'x');
	
				if (data == NULL)
					return it == end;
	
				return cursor == data->capacity();
			}
			
			inline bool operator !=(char b) const{
	
				assert (b == 'x');
				
				if (data == NULL)
					return it != end;
				
				return cursor != data->capacity();
			}
	
			inline int getIdFromPos(int x) const{
			
				return x < data->edges.size() ? edgeIdFromPos(data->edges.size()-1-x) : nodeIdFromPos(x-data->edges.size());
			}
			
			inline iterator& operator++(){
	
				if (data == NULL){
					
					if ( it == end)
						return (*this);
					
					++it;
					return (*this);
				}
				
				++cursor;
				
				while ( cursor < data->capacity() && data->count(getIdFromPos(cursor) ) == 0)
					++cursor;
				
				return (*this);
			}
	
			
		};

	const iterator begin() const{
	
		if (big)
			return iterator( big->cbegin(), big->cend() );
	
		return iterator(this);
	}
	
	const iterator end;

	
	const iterator cbegin() const{
	
		return iterator(this);
	}
	
	inline subgraph1(){
	
		vertices.reset();
		edges.reset();
	
	}
	
	inline subgraph1(const subgraph1& copy){
	
		vertices = copy.vertices;
		edges = copy.edges;
		
		if (copy.big)
			big.reset(new std::set<int>( *copy.big ));
	}
	
	const iterator cend;

	std::bitset<64> vertices;
	std::bitset<64> edges;
	
	unique_ptr<std::set<int> > big;
	
	inline int pos(int id) const{
	
		return id < 0 ? -id-1 : id;
	}
	
	inline int capacity() const {
	
		return vertices.size()+edges.size();
	}
	
	inline int size() const {
		
		return big ? big->size() : vertices.count()+edges.count();
	}
	
	inline void clear(){
	
		if (big)
			big->clear();
	
		vertices.reset();
		edges.reset();
	}
	
	inline int count(int n) const{
	
		return big ? big->count(n) : pos(n) >= vertices.size() ? 0 :  n < 0 ? edges.test( pos(n) ) : vertices.test(pos(n) );
	}
	
	inline void switchToBig(){
	
		//cout << "switch to big" << endl;
	
		std::unique_ptr<std::set<int> > uniq;
		
		uniq.reset(new std::set<int>() );
		
		for (auto it = begin(); it != end(); ++it)
			uniq->insert(*it);
		
		big = std::move(uniq);
	}
	
	inline void insert(int n){
	
		if (big){
		
			big->insert(n);
			
		} else {
		
			if (pos(n) >= vertices.size() ){
			
				switchToBig();
				big->insert(n);
				return;
			}
			
			if (n < 0)
				edges.set( pos(n) );
			else
				vertices.set( pos(n) );
		}
	}

};

//typedef std::set<int> subgraph;
typedef subgraph1 subgraph;

class UndirectedHypergraph{

	typedef vector<int> X;
	
	typedef vector<int> E;
	typedef vector<int> V;
	
	typedef vector<int>::const_iterator Iter;
	
	vector<E> edges;
	vector<V> vertices;

public:		

	inline void clear(){
	
		edges.clear();
		vertices.clear();
	}	
	
	
	inline long edgeCount() const {
	
		return edges.size();
	}

	inline long vertexCount() const {
	
		return vertices.size();
	}

	inline int id(const X& e) const{
	
		return e.at(0);
	}
	
	inline X& edge(int id){
	
		return edges.at( pos(id) );
	}
	
	inline X& vertex(int id){
	
		return vertices.at( pos(id) );
	}
	
	inline bool isVertex(int id) const{
	
		return id >= 0;
	}
	
	inline bool isEdge(int id) const{
	
		return id < 0;
	}
	
	inline bool isVertex(const X& x) const{
	
		return isVertex( id(x) );
	}
	
	inline bool isEdge(const X& x) const{
	
		return isEdge( id(x) );
	}
	
	inline int pos(int id) const{
	
		return id < 0 ? -id-1 : id;
	}
	
	inline int vertexIdFromPos(int pos) const {
	
		return pos;
	}
	
	inline int edgeIdFromPos(int pos) const {
	
		return -pos-1;
	}
	
	inline X& get(int id){
	
		return id < 0 ? edge(id) : vertex(id);
	}
	
	inline const X& getConst(int id) const{
	
		return id < 0 ? edges.at( pos(id) ) : vertices.at( pos(id) );
	}
	
	inline int addVertex(){
		
		const int id1 = vertices.size();
		
		vertices.push_back( {id1} );
		
		return id1;
	}
	
	inline Iter adjBegin(int id1) const{
	
		if ( getConst(id1).size() == 1)
			return getConst(id1).cend();
	
		return getConst(id1).cbegin()+1;
	}
	
	inline Iter adjEnd(int id1) const {
	
		return getConst(id1).cend();
	}
	
	
	inline Iter adjEdgesBegin(int id1) const{
		
		return adjBegin(id1);
	}
	
	inline Iter adjEdgesEnd(int id1) const{
	
		return std::lower_bound(adjBegin(id1), adjEnd(id1), 0);
	}
	
	
	inline int degree(int id1) const{
	
		return adjEdgesEnd(id1)-adjEdgesBegin(id1);
	}
	
	inline Iter adjVerticesBegin(int id1) const{
		
		return std::lower_bound(adjBegin(id1), adjEnd(id1), 0);
	}
	
	inline Iter adjVerticesEnd(int id1) const{
	
		return adjEnd(id1);
	}
	
	inline int addEdge(int id1, int id2){
	
		assert (id1 != id2);
		
		if (id1 > id2){
		
			int old1 = id1;
			int old2 = id2;
			id2 = old1;
			id1 = old2;
		}
		
		const int id3 = edgeIdFromPos(edges.size() );
	
		edges.push_back( {id3, id1, id2} );
		
		connect(id1, id2);
		connect(id1, id3);
		connect(id2, id3);
		
		return id3;
	}
	
	
	inline bool adjacent(int id1, int id2) const{
		
		if (id1 == id2)
			return true;
		
		return std::binary_search( adjBegin(id1), adjEnd(id1), id2 );
	}
	
	inline int connectingEdge(int id1, int id2) const {
	
		auto it1 = adjEdgesBegin(id1);
		auto end1 = adjEdgesEnd(id1);
	
		auto it2 = adjEdgesBegin(id2);
		auto end2 = adjEdgesEnd(id2);
	
		while (it1 != end1 && it2 != end2){
	
			if ( (*it1) == (*it2) )
				return (*it1);
		
			if ( (*it1) < (*it2) ){
		
				++it1;
		
			} else {
		
				++it2;
			}
		}
		
		return 0;
	}
	
	/*
	inline void mergeNodes(int id1, int id2){
	
		int e = connectingEdge(id1, id2);
		
		assert (isEdge(e) );
		
		remove(id1, e);
		remove(id2, e);
		
		remove(e);
		
		remove(id2, id1);
		
		for (auto it = adjNodesBegin(id2); it != adjNodesEnd(id2); ++it)
			connect(id1, *it);
		
		for (auto it = adjEdgesBegin(id2); it != adjEdgesEnd(id2); ++it)
			connect(id1, *it);
		
		remove(id2);	
		
		
	}*/
	
	inline bool adjacent2(int id1, int id2) const {
	
		if ( isVertex(id1) && isVertex(id2) ){
		
			if ( isEdge(connectingEdge(id1, id2) ))
				return true;
				
		} else {
		
			return adjacent(id1, id2);
		}
		
		return false;
	}
	
	inline void remove(int id1, int id2){
	
		X& x = get(id1);
		
		assert (x.at(0) != id2);
	
		auto it = std::lower_bound( x.begin()+1, x.end(), id2);
		
		if (it == x.end() )
			return;
		
		if ( (*it) == id2)
			get(id1).erase(it);
	}
	
	inline void remove(int id1){
		
		auto st = adjBegin(id1);
		auto en = adjEnd(id1);
	
		for (auto it = st; it != en; ++it)
			remove ( *it, id1 );
			
		if (isEdge(id1) ){
		
			auto vst = adjVerticesBegin(id1);
			auto ven = adjVerticesEnd(id1);
			
			for (auto it = vst; it != ven; ++it){
			
				const int v1 = (*it);
				
				if (v1 == id1)
					continue;
			
				for (auto it2 = vst; it2 != ven; ++it2){
			
					const int v2 = (*it2);
					
					if (v2 == id1 || v2 == v1)
						continue;
					
					if (!adjacent2(v1, v2) ){
					
						remove(v1, v2);
						remove(v2, v1);
					}
				}
			}
		
		}
		
		 get(id1).resize(1);
		//x.shrink_to_fit();
	}
	
	inline void connect(int id1, int id2){
	
		add(id1, id2);
		add(id2, id1);
	}
	
	inline int otherEnd(int edge, int vertex) const{
	
		for (auto it = adjVerticesBegin(edge); it != adjVerticesEnd(edge); ++it){
		
			if (*it != vertex)
				return (*it);
		}
		
		return vertex;
	}
	
	inline void add(int id1, int k){
	
		X& x = get(id1);
	
		if (k == x[0])
			return;
	
		if (k == x.back() )
			return;
	
		if (x.size() == 1 || k > x.back() ){
			x.push_back(k);
			return;
		}
	
		auto it = std::lower_bound(x.begin()+1, x.end(), k );
		
		assert (it != x.end() );
		
		if ( (*it) == k)
			return;
		
		x.insert( it, k);
	}
	
	template<typename T>
	inline const T label(const vector<T>& vertexLabels, const vector<T>& edgeLabels, int id){
	
		if (isEdge(id))
			return edgeLabels[pos(id)];
		
		if (isVertex(id))
			return vertexLabels[pos(id)];
		
		assert (false);
	}
	
	
	
	subgraph component( int id1, bool addVertices = true, bool addEdges = false, subgraph removed = subgraph() ) const{
	
		subgraph found;
		
		subgraph ret;
		
		if (addVertices)
			ret.insert(id1);
		
		std::queue<int> q;
		
		q.push(id1);
		
		while (q.size() > 0){
			
			int u = q.front();
			q.pop();
			
			if (removed.count(u) > 0)
				continue;
			else
				removed.insert(u);
				
			auto e_beg = adjEdgesBegin(u);
			auto e_end = adjEdgesEnd(u);
			
			for (auto it = e_beg; it != e_end; ++it){
			
				const int e = (*it);
				
				if (removed.count(e) > 0)
					continue;
					
				auto v_beg = adjVerticesBegin(e);
				auto v_end = adjVerticesEnd(e);
				
				for (auto vit = v_beg; vit != v_end; ++vit){
				
					const int v = (*vit);
				
					if (v != u && removed.count(v) == 0){
					
						if ( found.count(v) == 0 ){
				
							if (addEdges && found.count(e) == 0){
								found.insert(e); 
								ret.insert(e);
							}
							
							if (addVertices)
								ret.insert(v);
							
							found.insert(v);
							q.push(v);
						}
					}
				}
			}
		}
		
		if ( (!addVertices) && (!addEdges) ){
		
			assert (ret.size() == 0);
		
			ret.insert( found.size() );
		}
		
		return ret;
	}
	
	int componentSize( int id1, subgraph removed = subgraph() ) const{
	
		subgraph ret = component(id1, false, false, removed);
		
		assert (ret.size() == 1);
		
		return *(ret.begin());
	}
	
	vector<int> shortestCycle(int id1, bool addVertices = true, bool addEdges = false, subgraph removed = subgraph() ) const {
	
		vector<int> ret;
		
		auto e_beg = adjEdgesBegin(id1);
		auto e_end = adjEdgesEnd(id1);
		
		for (auto it = e_beg; it != e_end; ++it){
		
			const int e = (*it);
			
			if (removed.count(e) > 0)
				continue;
			
			const int v = otherEnd(e, id1);
		
			subgraph r2 = removed;
			r2.insert(e);
			
			vector<int> sp = shortestPath(id1, v, addVertices, addEdges, r2);
			
			if (sp.size() == 0)
				continue;
			
			if (addEdges)
				sp.push_back(e);
			
			if (addVertices)
				sp.push_back(v);
			
			if (sp.size() < ret.size() || ret.size() == 0)
				ret = sp;
			
		}
		
		return ret;
	}
	
	vector<int> shortestPath( int id1, int id2 , bool addVertices = true, bool addEdges = false, subgraph removed = subgraph() ) const{
	
		map<int, int> dist;
		map<int, int> predEdge;
		map<int, int> pred;
		
		std::priority_queue< std::pair<int, int> > q;
		
		q.push( std::make_pair(0, id2) );
		
		dist[id2] = 0;
		
		int mindist = vertices.size();
		
		while (q.size() > 0){
			
			int u = q.top().second;
			q.pop();
			
			assert (isVertex(u) );
			
			if ( (dist[u]+1) >= mindist)
				break;
			
			if (removed.count(u) > 0)
				continue;
			else
				removed.insert(u);
				
			auto e_beg = adjEdgesBegin(u);
			auto e_end = adjEdgesEnd(u);
		
			for (auto it = e_beg; it != e_end; ++it){
			
				const int e = (*it);
				
				if (removed.count(e) > 0)
					continue;
				
				assert (isEdge(e) );
				
				assert (adjacent(e, u));
				
				auto v_beg = adjVerticesBegin(e);
				auto v_end = adjVerticesEnd(e);
				
				for (auto vit = v_beg; vit != v_end; ++vit){
				
					const int v = (*vit);
					
					assert (isVertex(v) );
					
					assert (adjacent(e, v));
					
					assert (adjacent(u, v));
					assert (adjacent2(u, v));
					
				
					if (v != u && removed.count(v) == 0){
					
						const int length_uv = 1;
						
						if ( (dist.count(v) == 0) ||((dist[u]+length_uv) < dist.at(v) ) ){
				
							dist[v] = dist[u]+length_uv;
							predEdge[v] = e;
							pred[v] = u;
							
							if (v == id1)
								mindist = dist[v] < mindist ? dist[v] : mindist;
							else
								q.push( std::make_pair(-dist[v], v) );
						}
					}
				}
			}
		}
		
		if (!addVertices && !addEdges){
		
			if (pred.count(id1) == 0)
				return {-1};
			else
				return {dist[id1]};
		}
		
		if (pred.count(id1) == 0)
			return vector<int>();
		
		vector<int> ret;
		
		ret.reserve(dist.at(id1));
		
		for (int x = id1; x != id2; x = pred.at(x)){
		
			if (addVertices)
				ret.push_back(x);
			
			if (addEdges)
				ret.push_back(predEdge[x]);
		}
		
		if (addVertices)
			ret.push_back(id2);
		
		return ret;
	}
	
	inline long distance(int id1, int id2, const subgraph& removed) const{
	
		return shortestPath(id1, id2, false, false, removed).at(0);
	}
	
	inline long distance(int id1, int id2) const{
	
		return shortestPath(id1, id2, false, false).at(0);
	}
	
	inline bool connected(int id1, int id2) const{
	
		return distance(id1, id2) >= 0;
	}
	
	inline bool connected(int id1, int id2, const subgraph& removed) const{
	
		return distance(id1, id2, removed) >= 0;
	}
	
	inline int countCycles() const{
		
		if (vertices.size() == 0)
			return 0;
		
		subgraph removed;
		removed.insert( vertexIdFromPos(0) );
		
		return countCycles( vertexIdFromPos(0), removed);
	}
	
	
	inline int countCyclesWithout(int edge) const{
		
		if (vertices.size() == 0)
			return 0;
		
		subgraph removed;
		removed.insert( vertexIdFromPos(0) );
		removed.insert(edge);
		
		return countCycles( vertexIdFromPos(0), removed);
	}

private:
	inline int countCycles(int root, subgraph& removed, int level=0) const{
		
		vector<int> children;
		
		int cycles = 0;
		
		const UndirectedHypergraph& g = *this;
		
		while (true){
			
			BasicAggregator ag1;
			BasicAggregator ag2;
		
			for (auto it = g.adjEdgesBegin(root); it != g.adjEdgesEnd(root); ++it){
		
				if (removed.count(*it) > 0)
					continue;
		
				int e = (*it);
				int v = g.otherEnd(*it, root);
				
				assert (g.isVertex(v) && (v != root) && g.adjacent(v, root) );
			
				if (removed.count(v) > 0)
					continue;
				
				ag1.add( v, 0.0);
				ag2.add( e, 0.0);

			}
			
			if (ag1.count() == 0)
				break;
			
			const int child = ag1.argMinLast();
			const int childEdge = ag2.argMinLast();
			
			removed.insert(childEdge);
			
			const int oldcycles = cycles;
			
			for (auto it = children.begin(); it != children.end(); ++it){
			
				const int otherChild = (*it);
				
				assert (child != otherChild);
			
				if (g.connected(child, otherChild, removed)){
					
					/*
					auto sp = g.shortestPath(child, otherChild, true, false, removed);
					
					cout << (*removed.begin() ) << endl;
					
					for (auto it = sp.begin(); it != sp.end(); ++it){
					
						cout << "--" << (*it) << endl;
					}
					*/
					
					++cycles;
					break;
				}
			}
			
			if (oldcycles == cycles){
			
				removed.insert(childEdge);
				children.push_back(child);
			}
		}
		
		for (int k = 0; k < children.size(); ++k)
			cycles += countCycles(children.at(k), removed, level+1);
		
		return cycles;
	}
};

class JoinGraph1{

public:
	vector<shared_ptr<DatabaseTable>> tables;
	
	vector< vector<string> > edgeLabels;
	
	vector<Tuple> joinAttrs;
	
	shared_ptr<WeightingFunction> f;
	shared_ptr<SQLExpression> parsed;
	
	UndirectedHypergraph g;
	
	bool spoiled = false;
	
	bool registerWeights = false;
	
	inline JoinGraph1(){
	
		
	}
	
	inline JoinGraph1( const JoinGraph1& gr){
	
		f = gr.f;
		parsed = gr.parsed;
		
		for (long j = 0; j < gr.tables.size(); ++j){
		
			if (gr.joinAttrs.size() > 0)
				add( gr.tables.at(j), gr.joinAttrs.at(j) );
			else
				add( gr.tables.at(j) );
		}
		
		cout << "copy joingraph" << endl;
	}
	
	
	inline subgraph getCyclebreaks(double& accrate, int rank=0){
	
		subgraph removed;
	
		for (long j = 0; j < tables.size(); ++j){
			
			const int v = g.vertexIdFromPos(j);
		
			if (removed.count(v) > 0)
				continue;
			
			for (int k = 0; k < 123; ++k){
			
				assert (k < 100);
			
				vector<int> cy = g.shortestCycle(v, false, true, removed);
			
				if (cy.size() == 0)
					break;
					
				//cout << cy.size() << endl;
				
				BasicAggregator agg1;
				BasicAggregator agg1_a;
				BasicAggregator agg1_b;
				
				vector<std::pair<long double, int>> cands;
				
				for (int i = 0; i < cy.size(); ++i){
				
					const int ee = cy.at(i);
					
					auto e = g.get(ee);
					
					auto a = g.pos(e[1]);
					auto b = g.pos(e[2]);
					
					const Tuple& t = UtilTuple::intersection(tables.at(a)->colNames, tables.at(b)->colNames);
					
					//cout << i << " " << t << endl;
					
					vector<int> v = {a,b};
					
					cout << "A4A" << endl;
					
					DatabaseTables tabs1(*tables.at(a),*tables.at(b));
					
					cout << "A4B" << endl;
					
					//DatabaseTables tabs2(*tables.at(b),*tables.at(a));
					
					const long double cyclicSel = tabs1.cyclicSelectivity(*f, 1000000, "abc");
					
					cout << "A4C" << endl;
					
					/*
					
					std::set<string> s;
					
					
					
					
					
					
					for (int j = 0; j < v.size(); ++j){
					
						DatabaseTable& tab = *tables.at(v.at(j) );
						
						const int ind = UtilTuple::find(tab.colNames, t[0]);
						
						long cnt = 0;
						
						for (auto it = tab.begin(); it != tab.end(); ++it){
						
							if (cnt > 1000)
								break;
								
							s.insert( (*it)[ind] );
						
							cnt++;
						}
					}
					
					const int distinct = s.size();
					*/
					
					
					/*
					if (UtilString::equals(t[0], "n2"))
						agg1.add(cy.at(i), 20);
					else
						agg1.add(cy.at(i), 1000);*/
						
					
					cout << "[" << label(a) << "]--" << label(ee) << "--["<< label(b)<< "] ACCEPTRATE " << cyclicSel << endl;
					
					cands.push_back( std::make_pair(-cyclicSel, ee) );
					
					//agg1.add(ee, cyclicSel);
					//agg1_a.add(a, cyclicSel);
					//agg1_b.add(b, cyclicSel);
				}
				
				std::sort(cands.begin(), cands.end() );
				
				auto pair = cands.at(rank % cands.size() );
				
				const int ee = pair.second;
					
				auto e = g.get(ee);
					
				auto a = g.pos(e[1]);
				auto b = g.pos(e[2]);
				
				const long double cyclicSel = -pair.first;
				
				/*
				const int a = agg1_a.argMax();
				const int b = agg1_b.argMax();
				const int ee = agg1.argMax();
				
				const double cyclicSel = agg1.max();*/
				
				cout << "REMOVED EDGE ";
				
				cout << "[" << label(a) << "]--" << label(ee) << "--["<< label(b)<< "] ACCEPTRATE " << cyclicSel << endl;
				
				accrate = cyclicSel;
				
				removed.insert( ee );
			}
			
		}
		
		return removed;
	}
		
	
	inline JoinGraph1( const JoinGraph1& gr, long sample, const string& seed){
	
		f = gr.f;
		parsed = gr.parsed;
		
		for (long j = 0; j < gr.tables.size(); ++j){
		
			auto tab = gr.tables.at(j)->getSample(sample, seed);
		
			if (gr.joinAttrs.size() > 0)
				add( tab, gr.joinAttrs.at(j) );
			else
				add( tab );
		}
		
		cout << "copy joingraph" << endl;
	}
	
	
	
	inline int addLink(int i, int j, const string& s){
	
		const int vi = g.vertexIdFromPos(i);
		const int vj = g.vertexIdFromPos(j);
	
		subgraph sg; // subgraph with edges not containing s removed
		
		for (long j = 0; j < g.edgeCount(); ++j){
		
			int e = g.edgeIdFromPos(j);
			
			int pos = g.pos(e);
		
			bool remove = true;
			
			for (int k = 0; k < edgeLabels.at(pos).size(); ++k)
				if (UtilString::equals(edgeLabels.at(pos).at(k), s))
					remove = false;
			
			if (remove)
				sg.insert(e);
		}
		
		// avert transitivity redundancies
		if (g.connected(vi, vj, sg))
			return 0;
		
		
		
		int e = g.connectingEdge(vi, vj);
		
		if (e != 0){
		
			vector<string>& vs = edgeLabels.at( g.pos(e) );
			vs = UtilTuple::stringvector( UtilTuple::merge(UtilTuple::tuple(vs), Tuple(s) ) );
			
			return e;
			
		} else {
		
			const int ret = g.addEdge(i, j);
		
			edgeLabels.push_back(vector<string>() );	
			edgeLabels.back().push_back(s);
			
			return ret;
		}
		
		//cout << "addLink(" << tables.at(i)->tableName << ", " << tables.at(j)->tableName << ", " << s << ")" << endl;
	}
	
	inline void setWeights(shared_ptr<SQLExpression> p){
	
		parsed = p;
		
		auto m = p->getAliases();
		f = std::make_shared<WeightingFunction>(p->parsed, m);
	}
	
	map<string, vector<int>> attrmap;
	
	inline long joinSize(int a, int b, bool approx=false){
	
		if (true){
			
			/*	
			auto ta = tables.at(a);
			auto tb = tables.at(b);
		
		
			WeightingFunction f;
			
			
			JoinProcessor jproc(Tuple(), f, joinAttrs.at(a), joinAttrs.at(b) );
		
			
			ta->join(*tb, jproc );
		
			return jproc.getJoinSize();*/
			
			DatabaseTables tabs(tables.at(a), tables.at(b));
			
			const vector<int> vs = { g.vertexIdFromPos(a), g.vertexIdFromPos(b)};
		
			const int e1 = g.connectingEdge(vs.at(0), vs.at(1) );
			
			Tuple t = UtilTuple::tuple(edgeLabels.at( g.pos( e1 ) ) );
			
			tabs.setJoinColumns( t );
			
			return approx ? tabs.joinsizeEstimate() : tabs.joinsize();
		}
		
		return UtilTables::getJoinSize(*tables.at(a), joinAttrs.at(a), *tables.at(b), joinAttrs.at(b), tables.at(a)->getSize(), tables.at(a)->estimatedBytes(), tables.at(b)->getSize(), tables.at(b)->estimatedBytes() );
	}
	
	
	
	inline int join(int a, int b){
		
		/*
		cout << "INPUT 1: " << endl;
		cout << tables.at(a)->toString(3) << endl;
		cout << "INPUT 2: " << endl;
		cout << tables.at(b)->toString(3) << endl;*/
	
		
		
		
		
	
		DatabaseTables tabs(tables.at(a), tables.at(b));
		
		const vector<int> vs = { g.vertexIdFromPos(a), g.vertexIdFromPos(b)};
		
		const int e1 = g.connectingEdge(vs.at(0), vs.at(1) );
	
		Tuple t = UtilTuple::tuple(edgeLabels.at( g.pos( e1 ) ) );
		
		tabs.setJoinColumns( t );
		
		cout << tables.at(a)->getColNames() << " join " << tables.at(b)->getColNames() << " along " << t << endl;
		
		Tuple t2 = t;
		
		for (int ii = 0; ii < 2; ++ii){
		
		
			for (auto it = g.adjEdgesBegin(vs.at(ii)); it != g.adjEdgesEnd(vs.at(ii)); ++it){
		
				int k = (*it);
			
				if (k != e1){
			
					Tuple tk = UtilTuple::tuple( edgeLabels.at( g.pos(k) ) );
			
					t2 = UtilTuple::merge(t2, tk );
				}
			}
		}
	
		auto jtab = registerWeights ? tabs.join(*f) :  tabs.join();
		
		//UtilTables::getJoin(*tables.at(a), joinAttrs.at(a), *tables.at(b), joinAttrs.at(b), tables.at(a)->getSize(), tables.at(a)->estimatedBytes(), tables.at(b)->getSize(), tables.at(b)->estimatedBytes() );
		
		//cout << "OUTPUT 1: " << endl;
		//cout << jtab->toString(3) << endl;
	
		jtab->init(10000);
		
		/*
		// union of join attributes
		for (long j = 0; j < joinAttrs.at(b).size(); ++j){
		
			if (UtilTuple::find(t, joinAttrs.at(b)[j] ) == -1)
				t.push_back(joinAttrs.at(b)[j]);
		}
		*/
		
		
		tables.at(a) = jtab;
		joinAttrs.at(a) = t2;
		
		tables.at(b) = tables.at(tables.size()-1);
		joinAttrs.at(b) = joinAttrs.at(tables.size()-1) ;
		
		tables.pop_back();
		joinAttrs.pop_back();
		
		cout << "after joining: " << UtilTuple::getString(jtab->getColNames(), ' ') << " size " << jtab->getSize() << " join columns " << t2 << endl;
		
		updateGraph();
		
		return a;
	}
	
	
	
	inline void addToGraph(int jj){
	
	
		long newInd = g.addVertex();
		
		Tuple& joinColNames = joinAttrs.at(jj);
		
		//cout << "add " << tables.at(jj)->tableName << " [" << tables.at(jj)->getColNames() << "] joined along " << joinColNames << endl;
	
		
		for (int j = 0; j < joinColNames.size(); ++j){
		
			const string col = joinColNames[j];
			
			if (attrmap.count(col) == 0)
				attrmap.insert( std::make_pair(col, vector<int>() ) );
			
			vector<int>& tabs = attrmap.at(col);
			
			for (auto it = tabs.begin(); it != tabs.end(); ++it){
			
				if (UtilTuple::find(joinAttrs.at(*it), col) != -1)
					addLink((*it), newInd, col);
			}
			
			tabs.push_back(newInd);
		}
		
		cout << "add to graph" << endl;
	}
	
	inline void updateGraph(){
	
		g.clear();
		attrmap.clear();
		edgeLabels.clear();
		
		for (int j = 0; j < tables.size(); ++j)
			addToGraph(j);
	
		cout << "updated graph" << endl;
	}
	
	inline int add(shared_ptr<DatabaseTable> table, Tuple joinColNames=Tuple()){
	
		
		
		tables.push_back(table);
		joinAttrs.push_back(joinColNames);
		
		addToGraph(tables.size()-1);
		
		return tables.size()-1;
	}
	
	
	inline int addLink(int i, shared_ptr<DatabaseTable> table, Tuple joinColNames=Tuple()){
	
		tables.push_back(table);
		
		const int j = g.addVertex();
		
		Tuple t = joinColNames;
		t = UtilTuple::intersection(t, tables.at(i)->getColNames());
		t = UtilTuple::intersection(t, tables.at(j)->getColNames());
		
		joinAttrs.push_back(t);
		
		
		const int vi = g.vertexIdFromPos(i);
		const int vj = g.vertexIdFromPos(j);
		
		const int ret = g.addEdge(vi, vj);
		
		edgeLabels.push_back( UtilTuple::stringvector(t) );	
		
		return ret;
	}
	
	
	// DEPRECATED
	inline long getJoinSize(subgraph removed = subgraph()) const{
		
		removed.insert(0);
		
		vector<double> edgeScores(edgeLabels.size(), 0);
		vector<double> vertexScores(tables.size(), 0);
		
		return JoinGraph1(*this).generateTree(vertexScores, edgeScores, f, Tuple(),0, removed, false, subgraph() )->getJoinSize();
	}
	
	
	inline string label(int x) const {
	
		if (g.isEdge(x)){
			return UtilString::toString(edgeLabels.at(g.pos(x) ));
		} else {
		
		
			return UtilTuple::getString(tables.at(g.pos(x))->getColNames(), ' '); //tables.at(g.pos(x))->tableName;
		}
	}
	
	inline string toString(int v, subgraph& removed, subgraph& semi, bool multi=true, bool skipEmpty=true) const {
	
	
		stringstream ss;
		
		bool empty = true;
		
		for (auto it1 = g.adjVerticesBegin(v); it1 != g.adjVerticesEnd(v); ++it1){
			
			int b = (*it1);
		
			int e = g.connectingEdge(v,b);
		
			if (removed.count(e) > 0)
				continue;
				
			assert (b != v);
		
			assert (g.isEdge(e) );
			
			removed.insert(e);
		
			if (semi.count(e) > 0){
			
				ss  << "{  " << label(v)<< " }~~[" << label(e)<< "]~~" << toString(b, removed, semi, false, false);
			
			} else {
			
				ss  << "{  " << label(v)<< " }--[" << label(e)<< "]--" << toString(b, removed, semi, false, false);
			}
		
			
			
			empty = false;
			
			
			if (!multi)	
				return ss.str();
		}
		
		if (empty && !skipEmpty)
			ss << "{  " << label(v)<< " }" << endl;
		
		
		return ss.str();
	}
	
	inline string toString(subgraph removed= subgraph(), subgraph semi= subgraph()) const {
	
		//subgraph removed;
		
		stringstream ss;
		
		vector<std::pair<long, long> > v;
		
		for (int j = 0; j < tables.size(); ++j)
			v.push_back( std::make_pair( g.degree( g.vertexIdFromPos(j) ), g.vertexIdFromPos(j)) );
		
		std::sort(v.begin(), v.end() );
		
		for (int j = 0; j < v.size(); ++j)
			ss  << toString(v.at(j).second, removed, semi, true, true);
		
		return ss.str();
	}
	
	inline string toString(int x) const {
	
		if (g.isEdge(x)){
			
			const int e = x;
			auto vs = g.getConst(e);
			
			const int v1 = vs.at(1);
			const int v2 = vs.at(2);
		
			return toString(v1)+"--["+label(e)+"]--"+toString(v2);
		} else {
		
			return "{  "+label(x)+"  }";
		}
	}
	
	
	inline void simplify(){
	
		
		for (int a = 0; a < (tables.size()-1); ++a){
			for (int b = a+1; b < tables.size(); ++b){
			
				const int e = g.connectingEdge(a,b);
				
				const vector<string>& v = edgeLabels.at( g.pos(e) );
				
				if (v.size() > 1){
				
					Tuple t = UtilTuple::tuple(v);
					
					Tuple is = t;
					
					is = UtilTuple::intersection(is, tables.at(a)->getColNames() );
					is = UtilTuple::intersection(is, tables.at(b)->getColNames() );
					
					if (t.size() == is.size() ){
					
						tables.at(a) = tables.at(a)->addColumnUnion(t);
						tables.at(b) = tables.at(b)->addColumnUnion(t);
						
						vector<string>& vv = edgeLabels.at( g.pos(e) );
						
						vv.clear();
						vv.push_back(UtilTuple::getString(t, '_'));
					}
				}
			}
		}
	
	}
	
	
	
	inline void shrinkCycles2(long budget, int mode, bool connected = false, subgraph black= subgraph(), subgraph white=subgraph(), bool alsoAcyclic = false){
	
		map<int, long> joinSizes;	
		
		int compressedTable = -1;
		
		const bool approx = true; //tableNum = 1; //tableNum == 1; //budget >= 10000000000L;
		
		for (int j = 0; j < 1000; ++j){
		
			if (tables.size() == 1)
				return;
		
			cout << "budget " << budget << endl;
			cout << toString() << endl;
			joinSizes.clear();
			
			BasicAggregator agg1;
			BasicAggregator agg2;
			
			for (int a = 0; a < tables.size(); ++a){
			
				//cout << "a " << a << endl;
			
				if (compressedTable != -1 && (a != compressedTable))
					continue;
					
				if (white.size() > 0 && white.count(a) == 0)
					continue;
				
				for (auto it1 = g.adjVerticesBegin(a); it1 != g.adjVerticesEnd(a); ++it1){
			
					int b = (*it1);
					
					//cout << "b " << b << endl;
					
					if (white.size() > 0 && white.count(b) == 0)
						continue;
					
					int e = g.connectingEdge(a,b);
					
					if (black.size() > 0 && black.count(e) > 0)
						continue;
					
					assert (g.isEdge(e) );
				
					subgraph removed;
					removed.insert(e);
					
					if ( (g.distance(a,b, removed) > 0) || (alsoAcyclic) ){
					
						if (joinSizes.count(e) == 0){
						
							//cout << "computing join size " << endl;
						
							joinSizes.insert( std::make_pair(e, joinSize(a,b, approx)) );
							
							long cost = joinSizes.at(e); //-tables.at(a)->getSize()-tables.at(b)->getSize();
						
							agg1.add( a, cost);
							agg2.add( b, cost);
						
							//cout << endl;
							cout << tables.at(a)->colNames << " joined with " << tables.at(b)->colNames << " = " << joinSizes.at(e) << " cost " << cost << endl;
							//cout << endl;
						}
						
						
					}
				}
			}
			
			if (agg1.count() == 0){
			
				cout << "DONE" << endl;
			
				return;
			}
			
				
			const long cost = agg1.min();
			
			budget -= cost;
			
			if (budget < 0){
			
				cout << toString() << endl;
			
				return;
			}
			
			const int a = agg1.argMin();
			const int b = agg2.argMin();
			
			cout << tables.at(a)->colNames << " joining with " << tables.at(b)->colNames << endl;
			
			const int joined = join(a, b);
			
			if (connected)
				compressedTable = joined;
			
			cout << endl;	
			cout << "TABLE NUMBER " << tables.size() << endl;
			
			simplify();
			
			DatabaseTables tabs(tables);
			
			Tuple tt;
			
				
			
			
			
			for (int i = 0; i < joinAttrs.size(); ++i){
			
				tt = UtilTuple::merge(joinAttrs.at(i), tt);
			}
			
			tabs.setJoinColumns(tt);
			
			if ( (mode == 999) && !isCyclic() ){
			
				cout << "CYCLIC " << isCyclic() << endl;
				cout << toString() << endl;
				cout << "DONE" << endl;
			
				return;
			}
			
			// Inversion mode, TODO fix legacy code
			if (mode == 1 && tabs.joincolumns() == 1){
			
				cout << "CYCLIC " << isCyclic() << endl;
				cout << toString() << endl;
				cout << "DONE" << endl;
				
				return;
			}
			
			
			
			for (long j = 0; j < tables.size(); ++j){
			
				cout << "size(" << tables.at(j)->colNames << ") = " << tables.at(j)->getSize() << endl;
			}
			
			cout << "/TABLE NUMBER " << tables.size() << endl;
			cout << endl;
			
			joinSizes.clear();
		}
	}
	
	inline void shrinkCycles(long budget, bool connected = false, subgraph black= subgraph(), subgraph white=subgraph(), bool alsoAcyclic = false){
	
		shrinkCycles2(budget, 999, connected, black, white, alsoAcyclic);
	}
	
	
	inline bool containedIn(const string& s, const string& s2){
	
		const string s3 = UtilString::toLowerCase(s2);
		
		return s3.find( UtilString::toLowerCase(s)) !=std::string::npos;
	}
	
	inline vector<int> findVertices(const string& s){
	
		vector<string> vec = UtilString::getSplits(s, ' ', true);
	
		vector<int> ret;
		
		for (long j = 0; j < tables.size(); ++j){
		
			int v1 = g.vertexIdFromPos(j);
			
			if (containedIn(s, tables.at(v1)->tableName))
				ret.push_back(v1);
		}
		
		return ret;
	}
	
	inline vector<int> findEdges(const string& s){
	
		vector<string> vec = UtilString::getSplits(s, ' ', true);
		
		vector<int> ret;
		
		if (vec.size() == 0)
			return ret;
			
		for (long j = 0; j < edgeLabels.size(); ++j){
		
			int e = g.edgeIdFromPos(j);
			auto vs = g.getConst(e);
			
			int v1 = vs.at(1);
			int v2 = vs.at(2);
			
			const vector<string>& cols = edgeLabels.at(j);
			
			vector<string> v;
			
			v.push_back( tables.at(v1)->tableName);
			v.push_back( tables.at(v2)->tableName);
			
			for (auto it = cols.begin(); it != cols.end(); ++it)
				v.push_back(*it);
			
			int found = 0;
			
			for (auto it1 = vec.begin(); it1 != vec.end(); ++it1){
			
				bool found1 = false;
			
				for (auto it2 = v.begin(); it2 != v.end(); ++it2)
					found1 |= containedIn( (*it1), (*it2) );
				
				if (found1)
					++found;
			}
			
			if (found == vec.size())
				ret.push_back(e);
		}
	
		return ret;
	}
	
	inline shared_ptr<JoinTree> getJoinTree(int root, int edge, bool breakcycles, subgraph& removed, subgraph& semi, subgraph& fk, int level=0 ) const{
	
	
		Tuple t;
		
		if (edge != 0){
			t = UtilTuple::tuple( edgeLabels.at( g.pos(edge) ) );
			removed.insert(edge);
		}
	
		removed.insert(root);
			
		const bool cyclic = breakcycles;
		
		shared_ptr<JoinTree> ret = std::make_shared<JoinTree>( tables.at(g.pos(root)), t, f, cyclic, cyclic && !breakcycles );
		
		if (semi.count(edge) > 0)
			ret->setSemiJoin();
		
		if (fk.count(edge) > 0)
			ret->setFK();
		
		for (auto it = g.adjEdgesBegin(root); it != g.adjEdgesEnd(root); ++it){
		
			const int e = (*it);
			
			if (removed.count(e) > 0)
				continue;
			
			if (edgeLabels.at(g.pos(e)).size() == 0)
				continue;
				
			const int v = g.otherEnd(e, root);
			
			if (removed.count(v) > 0)
				continue;
				
			Tuple te = UtilTuple::tuple( edgeLabels.at( g.pos(e) ) );
			
			for (int i = 0; i < level; ++i)
				cout << "   ";
			cout << "root " << tables.at(g.pos(root))->getColNames() << " edge " << te << " child " << tables.at(g.pos(v))->getColNames() << endl;
				
			
			ret->addChild(getJoinTree(v, e, breakcycles, removed, semi, fk, level+1) );
		}
		
		return ret;
	}
	
	
	inline shared_ptr<JoinTree> getJoinTree(bool breakcycles, subgraph removed=subgraph(),subgraph semi=subgraph(),subgraph fk=subgraph()) const{
	
		BasicAggregator ag;
		
		
		
		
		for (long j = 0; j < tables.size(); ++j){
		
			if (removed.count( g.vertexIdFromPos(j)) > 0)
				continue;
				
			bool sem = true;
			for (auto it = g.adjEdgesBegin(j); it != g.adjEdgesEnd(j); ++it)
				sem &= semi.count(*it) > 0;
			
				
			if (!sem){
				ag.add(j, tables.at(j)->hasIndices() ? tables.at(j)->getSize() >> 5 : tables.at(j)->getSize()  );
			}
		}
		
		const int root = g.vertexIdFromPos(ag.argMaxLast() );
		
		return getJoinTree(root, 0, breakcycles, removed, semi, fk);
		
	}
	
	inline shared_ptr<JoinTree> generateTree(bool breakcycles = false, subgraph cyclebreaks=subgraph()) const{
		
		assert (false);
		
		assert (!spoiled);
		
		vector<double> edgeScores; // prefer to break cycles at edges with higher scores
		vector<double> vertexScores; // prefer to pick roots with higher scores
		
		
		
		//cout << "optimise " << optimise << " cyclic " << isCyclic() << " breakcycles " << breakcycles << endl;
		
		/*
		if (optimise && (isCyclic() && breakcycles) ){
		
			JoinGraph1 proto(*this, 1000, "abc");
			
			cout << "full join " << JoinGraph1(proto).getJoinSize() << endl;
			
			for (long j = 0; j < edgeLabels.size(); ++j){
		
				int e = g.edgeIdFromPos(j);
				
				auto vs = g.getConst(e);
				
				subgraph removed;
				
				removed.insert(e);
				
				long part = proto.getJoinSize(removed);
				
				int node1 = vs.at(1);
				int node2 = vs.at(2);
			
				auto t1 = tables.at( g.pos(node1) );
				auto t2 = tables.at( g.pos(node2) );
				
				cout << t1->getColNames() << "[" << t1->getSize() << "] join " << t2->getColNames() << "[" << t2->getSize() << "] removed from join = " << part << endl;
				
				edgeScores.push_back(part);
			}
			
		} else {*/
		
			JoinGraph1 proto(*this, 10000, "abc");
		
			for (long j = 0; j < g.edgeCount(); ++j){
		
				auto e = g.getConst(g.edgeIdFromPos(j));
			
				assert (e.size() == 3);
			
				int node1 = e.at(1);
				int node2 = e.at(2);
			
				auto t1 = proto.tables.at( g.pos(node1) );
				auto t2 = proto.tables.at( g.pos(node2) );
				
				//long joinSize = UtilTables::getJoinSize(*t1, *t2);
			
				//cout << t1->getColNames() << "[" << t1->getSize() << "] join " << t2->getColNames() << "[" << t2->getSize() << "] = " << joinSize << endl;
				
				edgeScores.push_back(0);
			}
		//}
		
		for (long j = 0; j < tables.size(); ++j)
			vertexScores.push_back(tables.at(j)->getSize() );
		
		BasicAggregator ag;
		
		for (long j = 0; j < tables.size(); ++j)
			ag.add(j, tables.at(j)->getSize()  );
		
		cout << ag.max() << endl;
		cout << tables.at(ag.argMaxLast() )->getSize() << endl;
		
		assert (ag.count() > 0);
		assert (ag.max() == tables.at(ag.argMaxLast() )->getSize() );
		
		const int root = g.vertexIdFromPos(ag.argMaxLast() );
		
		subgraph removed;
		removed.insert(root);
	
		return JoinGraph1(*this).generateTree(vertexScores, edgeScores, f, Tuple(), root, removed, breakcycles, cyclebreaks, 0);
	}
	
	inline bool isCyclic() const{
	
		return countCycles() > 0;
	}
	
	inline int countCycles() const{
		
		return g.countCycles();
		/*assert (!spoiled);
		
		std::set<int> removed;
		removed.insert(0);
		
		return countCycles(0, removed);*/
	}
	
private:
	
	/*
	inline int countCycles(int root, std::set<int>& removed, int level=0) const{
		
		vector<int> children;
		
		int cycles = 0;
		
		while (true){
			
			BasicAggregator ag1;
			BasicAggregator ag2;
		
			for (auto it = g.adjEdgesBegin(root); it != g.adjEdgesEnd(root); ++it){
		
				if (removed.count(*it) > 0)
					continue;
		
				int v = g.otherEnd(*it, root);
				
				assert (g.isVertex(v) && (v != root) && g.adjacent(v, root) );
			
				if (removed.count(v) > 0)
					continue;
				
				ag1.add( v, 0 );
				ag2.add( *it, 0);

			}
			
			if (ag1.count() == 0)
				break;
			
			int next = ag1.argMinLast();
			int nextedge = ag2.argMinLast();
			
			removed.insert(nextedge);
			
			for (auto it = children.begin(); it != children.end(); ++it){
				if (g.connected(*it, next, removed)){
					
					++cycles;
				}
			}
			
			removed.insert(next);
			children.push_back(next);
		}
		
		for (int k = 0; k < children.size(); ++k){
		
			const int child = children.at(k);
			
			cycles += countCycles(child, removed, level+1))
			
		}
		
		return cycles;
	}*/

	inline shared_ptr<JoinTree> generateTree(const vector<double>& vertexScores, const vector<double>& edgeScores, shared_ptr<WeightingFunction> f, const Tuple t, const int root, subgraph removed1, bool breakcycles = false, subgraph cyclebreaks = subgraph(), int level=0){
		
		assert (false);
		
		cout << level << " begin " << tables[root]->getColNames() << " " << t << endl;
		
		vector<std::tuple<double, int, int> > children;
		
		bool cyclic = false;

		subgraph removed = removed1;
		
		for (auto it = cyclebreaks.begin(); it != cyclebreaks.end(); ++it)
			removed.insert( *it);
		
		while (true){

			
			BasicAggregator ag1;
			BasicAggregator ag2;
		
			for (auto it = g.adjEdgesBegin(root); it != g.adjEdgesEnd(root); ++it){
		
				const int e = (*it);
		
				if (removed.count(e) > 0)
					continue;
					
				//if (cyclebreaks.count(e) > 0)
				//	continue;
				
				int v = g.otherEnd(e, root);
			
				if (removed.count(v) > 0)
					continue;
			
				const long double vertexScore = vertexScores.at( g.pos(v) );
				const long double edgeScore = edgeScores.at( g.pos(e) );
			
				ag1.add( v, vertexScore );
				ag2.add( e, vertexScore);
			}
			
			if (ag1.count() == 0)
				break;
			
			int child = ag1.argMaxLast();
			int childEdge = ag2.argMaxLast();
			
			removed.insert(childEdge);
			
			bool circlesBack = false;
			
			for (auto it = children.begin(); it != children.end(); ++it){
				
				const int otherChild = std::get<2>(*it);
				
				assert (otherChild != child);
				
				if (!g.connected(otherChild, child, removed))
					continue;
				
				auto p = g.shortestPath(otherChild, child, false, true, removed);
				
				if (p.size() > 0){
					
					circlesBack = true;
					
					vector<string>& reloc = edgeLabels.at( g.pos(childEdge) );
					
					for (auto e_it = p.cbegin(); e_it != p.cend(); ++e_it){
					
						vector<string>& dst = edgeLabels.at( g.pos(*e_it) );
					
						for (auto it2 = reloc.cbegin(); it2 != reloc.cend() ; ++it2)
							dst.push_back(*it2);
					}
					
					reloc.clear();
					break;
				}
			
			}
			
			if (circlesBack){
			
				cyclic = true;
			
			} else {
			
				removed.insert(child);
				children.push_back( std::make_tuple(-edgeScores.at(g.pos(childEdge)),  child, childEdge) );
				
				//std::stable_sort(children.begin(), children.end() );
			}
			
		
		}
		
		
		
		//cout << level << " extend " << tables[g.pos(root)]->getColNames() << " " << t << endl;
		
		shared_ptr<JoinTree> ret = std::make_shared<JoinTree>( tables[g.pos(root)], t, f, cyclic, cyclic && !breakcycles );
		
		
		for (int k = 0; k < children.size(); ++k){
		
			if (k > 0){	
			
				assert ( edgeScores.at( g.pos( std::get<2>(children.at(k))) ) <=  edgeScores.at( g.pos(std::get<2>(children.at(k-1)))  ) ); 
			}
			
			const int child = std::get<1>(children.at(k));
		
			const vector<string>& set = edgeLabels.at( g.pos(std::get<2>(children.at(k))) );
			
			std::set<string> included;
			
			Tuple attrs;
			
			for (auto it2 = set.cbegin(); it2 != set.cend(); it2++){
			
				const string& s = (*it2);
				
				if (included.count(s) > 0)
					continue;
				
				attrs.push_back(s);
				
				included.insert(s);
				
				if (breakcycles)
					break;
			}
			
			//cout << level << " addChild " << tables[g.pos(child)]->getColNames() << " " << attrs << endl;
			
			ret->addChild(generateTree(vertexScores, edgeScores, f, attrs, child, removed, breakcycles, cyclebreaks, level+1));
			
		}
		
		return ret;
	}
};


typedef JoinGraph1 JoinGraph;



#endif
