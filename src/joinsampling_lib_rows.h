#ifndef SHEKELYAN_JOINSAMPLING_LIB_ROWS
#define SHEKELYAN_JOINSAMPLING_LIB_ROWS




#include <approxdata/utils/utils.h>
//#include <approxdata/external/pq.h>
#include <approxdata/external/csv.h>
#include <iostream>
//#include <libpq-fe.h>
//#include <sql.h>
//#include <sqlext.h>

#include <external/robin_hood.h>

#include <queue>
#include <map>



const bool ONLYINTEGERS = false;

template <typename E, typename N>
using long_map = robin_hood::unordered_flat_map< E, N>;

//using long_map = std::unordered_map< E, N>;

#define CSV_IO_NO_THREAD



const vector<string> ParsedMathExpression_functions = {"", "abs","log", "exp", "sqrt", "ind"};
const vector<char> ParsedMathExpression_ops = {'0', 'a', 'l', 'e', 's','i'};

const vector<char> ParsedMathExpression_binops = {'+','-','*','/','^',';','?','<','>','=','!'};

const vector<string> ParsedMathExpression_binops_strings = {"+","-","*","/","**","<=",">=","<",">","==","!="};
		
const std::hash<std::string> HASHER_STRING;
	




const long TERM = std::numeric_limits<long>::max();

const long REMOVE = std::numeric_limits<long>::max()-1;




namespace UtilBytePacking{



	template<typename T>
	inline int packLength(T l){
	
		
	
		return l == 0 ? 1 : sizeof(T)+2;
	}
	
	template<typename T>
	inline T unpack(const char* chptr_){
	
		if ( (*chptr_) != 1 && (*chptr_) != 2)
			return -1;
	
		const uint8_t* chptr = reinterpret_cast<const uint8_t*>(chptr_);
		
		++chptr;
		
		T ret = 0;
		
		uint8_t* retbytes = reinterpret_cast<uint8_t*>(&ret);
		
		++chptr;
		
		//int i = 0;
		
		
		unsigned int sel = *(chptr-1);
		
		/**/
		if (sel == 0)
			return ret;
		
		if (sel & 1)
			(*retbytes) = (*chptr);
	
		sel >>= 1;
		++chptr;
		++retbytes;
		/**/
		/**/
		if (sel == 0)
			return ret;
		
		if (sel & 1)
			(*retbytes) = (*chptr);
	
		sel >>= 1;
		++chptr;
		++retbytes;
		/**//**/
		if (sel == 0)
			return ret;
		
		if (sel & 1)
			(*retbytes) = (*chptr);
	
		sel >>= 1;
		++chptr;
		++retbytes;
		/**//**/
		if (sel == 0)
			return ret;
		
		if (sel & 1)
			(*retbytes) = (*chptr);
	
		sel >>= 1;
		++chptr;
		++retbytes;
		/**//**/
		if (sel == 0)
			return ret;
		
		if (sel & 1)
			(*retbytes) = (*chptr);
	
		sel >>= 1;
		++chptr;
		++retbytes;
		/**//**/
		if (sel == 0)
			return ret;
		
		if (sel & 1)
			(*retbytes) = (*chptr);
	
		sel >>= 1;
		++chptr;
		++retbytes;
		/**//**/
		if (sel == 0)
			return ret;
		
		if (sel & 1)
			(*retbytes) = (*chptr);
	
		sel >>= 1;
		++chptr;
		++retbytes;
		/**//**/
		if (sel == 0)
			return ret;
		
		if (sel & 1)
			(*retbytes) = (*chptr);
	
		sel >>= 1;
		++chptr;
		++retbytes;
		/**//**/
		if (sel == 0)
			return ret;
		
		if (sel & 1)
			(*retbytes) = (*chptr);
		
		/*
		
		
		for (unsigned int sel = *(chptr-1); sel != 0 ; sel >>= 1, ++chptr, ++retbytes){
		
			//cout << sel << " " << i << endl;
			
			if (sel & 1)
				(*retbytes) = (*chptr);
				
			//++i;
		}*/
		
		return ret;
		
	}
	
	template<typename T>
	inline void pack( char* chptr, T l){
	
		char* sanity = chptr;
		
		assert (l > 0);
		
		{
	
			const int BYTES = sizeof(T);
		
			(*chptr) = BYTES >> 2;
		
			++chptr;
		
			uint8_t* selptr = reinterpret_cast<uint8_t*>(chptr);
		
			++chptr;
		
			T* lptr = reinterpret_cast<T*>(chptr);
		
			(*lptr) = l;
		
			uint8_t sel = 0;
		
			for (int i = 0; i < BYTES; ++i){
		
				if ( (*chptr) == 0)
					(*chptr) = 1;
				else
					sel |= 1 << i;
			
				++chptr;
			}
		
			(*selptr) = sel;
		
			uint8_t sel2 = *(selptr);
		
			assert (sel2 == sel);
		}
		
		for (int i = 0; i < packLength<T>(l); ++i){
		
			assert ( *(sanity+i) != 0);
		}
		
		assert (unpack<T>(sanity) == l);
	}


};


/*

Structure to maintain multiple sequences (of unknown length) in one vector (to avoid vector<vector<long>> as the vector<long> grows)

insert() starts a new sequence and gives the id of the sequence

next(id) moves towards the next sequence id (needs to be next(id+1) to guarantee moving forward)

get(id) gives sequence element at id

to add elements to the end of a sequence one can do:

long seqbeg = gv.insert();

long seqend = seqbeg;

seqend = gv.push_back(seqend, element1);
seqend = gv.push_back(seqend, element2);
seqend = gv.push_back(seqend, element3);

to read sequence one can do a loop such as:

long seqbeg = gv.insert();

for (long x = gv.next(seqbeg); x != -1; x = gv.next(x+1))
	cout << gv.get(x) << endl;

*/

class GrowingVector{

public:	
	vector<long> v;
	
	long size = 0;
	
	inline void reserve(long x){
	
		v.reserve(x*2);
	}
	
	inline void shrink_to_fit(){
	
		v.shrink_to_fit();
	}
	
	inline void print(){
	
		for (long j = 0; j < v.size(); ++j){
		
			cout << "[" << v[j] << "]";
		}
		cout << endl;
	
	}
	
	inline long getBack(long ind){
	
	
		while (v[ind] != TERM){
		
			if (v[ind] < 0)
				ind = -v[ind];
			else
				++ind;
		}
		return ind;
	}
	
	
	inline long getLast(long ind){
	
		long last = -1;
		
		for (long x = next(ind); x != -1 ; x = next(x+1))
			last = x;
		
		return last;
	}
	
	inline long get(long index) const {
	
		return v[index];
	}
	
	inline long next(long ind) const{
		
		assert (ind >= 0);
		
		while (v[ind] < 0)
			ind = -v[ind];
		
		if (v[ind] == TERM)
			return -1;
			
		return ind;
	}
	
	inline long insert(){
	
		v.push_back(TERM);
		
		return v.size()-1;
	}
	
	inline long getSize() const{
	
		return v.size();
	}
	
	inline long push_back(long ind, long val){
	
		assert( val >= 0);
	
		ind = getBack(ind);
		
		if (ind < (v.size()-1) ){
		
			v[ind] = -v.size();
			v.push_back(val);
			v.push_back(TERM);
			
			return v.size()-2;
			
		} else {
		
			v[ind] = val;
			v.push_back(TERM);
			return ind;
		}
	}
	
	inline void remove(long ind){
	
		assert (v[ind] >= 0);
		
		long last = getLast(ind);
		
		assert (last != -1);
		
		v[ind] = v[last];
		v[last] = TERM;
	}
	
	
	inline void markForRemoval(long ind){
	
		assert (v[ind] >= 0);
		
		v[ind] = REMOVE;
		--size;
	}
	
	inline void removeMarked(long ind){
			
		long r = next(ind);
		long w = next(ind);
		
		while (r != -1){
		
			while ( (r != -1) && (v[r] == REMOVE) )
				r = next(r+1);
		
			if (r == -1)
				break;
			
			v[w] = v[r];
			
			
			r = next(r+1);
			w = next(w+1);
		}
		
		if (w != -1)
			v[w] = TERM;

		
	}
	

};


class SampleEntry{
public:
	long_double u = -1;
	long id = -1;
	long_double w = -1;

	inline SampleEntry(double u, long id) : u(u), id(id){
	
		
	}
	
	inline bool operator<(const SampleEntry& e) const{
		
		assert (id != -1);
		assert (e.id != -1);
		
		return u < e.u;
	}
};

class RandomNumbers{

typedef long_double Weight;

public:

	std::uniform_real_distribution<double> doubleDistribution;

	std::seed_seq seed1;
	std::mt19937 generator;
	
	
	
	
	long_map<long,long> sparse;
	
	vector<double> pool;
	
	double* poolit = NULL;
	double* poolend = NULL;
	
	long counter = 0;
	
	inline RandomNumbers(const string seed) : doubleDistribution(0.0, 1.0), seed1(seed.begin(), seed.end()), generator(seed1){
	
		
	}
	
	

	inline long randomInteger(const long a, const long b){
		
		if (a == b)
			return a;
		
		std::uniform_int_distribution<long> distribution(a,b);
		
		return distribution(generator);
	}
	
	inline void loadPool(){
	
		if (pool.size() != 1000){
		
			pool.clear();
		
			pool.reserve(1000);
		
			for (long j = 0; j < 1000; ++j)
				pool.push_back(doubleDistribution(generator));
				
		} else {
		
			for (long j = 0; j < 1000; ++j)
				pool[j] = doubleDistribution(generator);
		}
			
		poolit = pool.data();
		poolend = poolit+pool.size();
		
		counter += 1000;
	}
	
	inline double randomDouble(){
		
		return doubleDistribution(generator);
		/*
		if (poolit == poolend)
			loadPool();
		
		const double ret = (*poolit);
		
		++poolit;
		return ret;*/
	}
	
	
	inline void randomIntegers(long u, vector<long>& v){
	
		for (long i = 0; i < v.size(); i++)
			v[i] = randomInteger(0, u-1);
	}
	
	
	inline void randomSortedIntegers(long u, vector<long>& v){
	
		randomIntegers(u, v);
		
		std::sort(v.begin(), v.end() );
	}
	
	inline void randomUniqueIntegers(long u, vector<long>& v){
	
		assert (u >= v.size() );
		
		const long n = v.size();
		
		// V = (1...u) represented by v (first n elements) and sparse (rest of elements);
		
		for (long i = 0; i < n; i++)
			v[i] = i;
		
		sparse.reserve(UtilMath::maxVal<long>(0, UtilMath::minVal<long>(v.size(), u-v.size() )));
		
		// if j is not contained in sparse, that means sparse[j] = j
		
		// obtain first n elements of random permutation of V using Knuth's shuffle
		for (long i = 0; i < v.size(); i++){
		
			const long j = randomInteger(i, u-1);
		
			if (j < v.size() ){
			
				// swap inside v
				const long vj = v[j];
				v[j] = v[i];
				v[i] = vj;
				
			} else {
				
				// swap between v and sparse
				const long sj = sparse.count(j) > 0 ? sparse[j] : j;
				sparse[j] = v[i];
				v[i] = sj;
			}
		}
		
		sparse.clear();
	}
	
	template<typename T>
	inline void randomPermutation(vector<T>& v){
	
		std::shuffle(v.begin(), v.end(), generator);
		
		/*
	
		for (long i = 0; i < v.size(); i++){
		
			const long j = randomInteger(i, v.size()-1);
				
			const T t = v[j];
			
			v[j] = v[i];
			v[i] = t;
		}*/
	}
	
	inline double maxOfRandomDoubles(long num){
	
		return exp( log( randomDouble() )/num);
	}
	
	inline double minOfRandomDoubles(long num){
	
		return 1.0-maxOfRandomDoubles(num);
	}
	
	
	
	// get smaller WOR sample from WR sample
	inline void getWORfromWR(vector<int>& freqs){
	
		for (long j = 0; j < freqs.size(); j++){
		
			assert (freqs[j] >= 0);
			
			if (freqs[j] > 0)
				freqs[j] = 1;
		}
	}
	
	
	// get WOR sample from WOR treating repetitions as distinct items
	inline void getWORfromWR(vector<long>& freqs, long m){
	
		if (freqs.size() == 0){
		
			assert (m == 0);
			return;
		}
		
		const long inputSize = UtilCollection::sum(freqs);
	
		assert (inputSize >= m);
		
		if (inputSize == m)
			return;
		
		vector<long> ordered;
		
		{
			vector<Weight> ws;
			vector<long> js;
			
			ws.reserve(inputSize);
			js.reserve(inputSize);
		
			for (long j = 0; j < freqs.size(); ++j){
			
				for (int freq = freqs.at(j); freq >= 1; --freq){
				
					ws.push_back(-1);
					js.push_back(j);
				}
			}
			
			randomSort(ws, js, ordered, true);
		}
		
		for (long j = 0; j < freqs.size(); j++)
			freqs.at(j) = 0;
		
		long count = 0;
		
		for (long j = 0; (j < ordered.size() ) && (count < m); j++){
		
			if (freqs.at(ordered.at(j)) == 0){
			
				++freqs.at(ordered.at(j));
				++count;
			}
		}
		
		assert (UtilCollection::sum(freqs) == m);
	}
	
	
	// get WOR sample from WOR treating repetitions as distinct items
	inline void getWORfromWOR(const Weight totalWeight, const vector<Weight>& weights, vector<int>& freqs, long m){
	
		if (freqs.size() == 0){
		
			assert (m == 0);
			return;
		}
		
		
		const long inputSize = UtilCollection::sum(freqs);
	
		assert (inputSize >= m);
		
		if (inputSize == m)
			return;
			
		assert (weights.size() == 0);
		
		
		vector<long> ordered;
		
		{
			
			vector<long> js;
			js.reserve(inputSize);
			
			for (long j = 0; j < freqs.size(); ++j)
					for (int freq = freqs.at(j); freq >= 1; --freq)
						js.push_back(j);
			
			if (weights.size() > 0){
			
				vector<Weight> ws;
				ws.reserve(inputSize);
			
				for (long j = 0; j < freqs.size(); ++j)
					for (int freq = freqs.at(j); freq >= 1; --freq)
						ws.push_back(weights.at(j) );
				
				//Weight smpWeight = UtilCollection::sum(ws);
				
				//for (long j = 0; j < ws.size(); ++j)
				//	ws.at(j) = exp( (log( ws.at(j) )/totalWeight)*smpWeight );
				
				randomSort(ws, js, ordered, true);
				
			} else {
			
				randomSort(js, ordered);
			}
			
		}
		
		for (long j = 0; j < freqs.size(); j++)
			freqs.at(j) = 0;
		
		for (long j = 0; j < m; j++)
			++freqs.at(ordered.at(j));
		
		assert (UtilCollection::sum(freqs) == m);
	}
	
	
	
	// WOR from multiset treating each repetition as distinct item (also get smaller WR from WR)
	inline void getWORfromWOR(vector<int>& freqs, long m){
	
	
		if (false){
		
			vector<Weight> v;
			getWORfromWOR(1, v, freqs, m);
			return;
		}
	
		long left = UtilCollection::sum(freqs);
		
		while (left > m){
		
			long k = left-m;
		
			long removed = 0;
			
			const long space = left;
			
			double g = exp( log( randomDouble() )/k );
			double ind = (1-g)*space;
			--k;
			
			double cursor = 0;
			
			for (long j = 0; j < freqs.size(); ++j){
		
				if (k < 0)
					break;
		
				const int f = freqs.at(j);
		
				for (int i = f; i >= 1; --i){
				
					cursor += 1;
					
					for (bool first = true; k >= 0 && ind <= cursor; first = false){
						
						if (first){
						
							--freqs.at(j);
							--left;
						}
						
						if (k > 0){
							g *= exp( log( randomDouble() )/k );
							ind = (1-g)*space;
						}
						
						--k;
					}
					
					
					if (k < 0)
						break;
					
				}
			}
		
		}
		
	}
	
	// WOR from multiset treating each repetition as distinct item (also get smaller WR from WR)
	inline void getWRfromWR(vector<int>& freqs, long m){
		
		getWORfromWOR(freqs, m);
	}
	
	
	
	
	inline void setCoords(const vector<Weight>& ws, vector<int>& freqs, const vector<Weight>& sortedCoords){
		
		assert (ws.size() == freqs.size() );
		
		Weight cursor = 0;
		
		long k = 0;
		
		for (long j = 0; j < freqs.size(); ++j){
		
			cursor += ws.at(j) * freqs.at(j);
		
			freqs.at(j) = 0;
			
			while (sortedCoords.at(k) <= cursor){
					
				++freqs.at(j);
				++k;
					
				if (k == sortedCoords.size()){
				
					for (long j2 = j+1; j2 < freqs.size(); ++j2)
						freqs.at(j2) = 0;
				
					return;
				}
			}
		}
		
		assert (false);
		
	}
	
	
	
	inline void randomSort(const vector<Weight>& ws, const vector<long>& js, vector<long>& ordered, const bool descending = false, bool special = false){
		
		vector<SampleEntry> ord;
	
		ord.reserve(ws.size() );
		
		for (long j = 0; j < ws.size(); ++j){
		
			if (descending)
				ord.push_back( SampleEntry( ws[j] == -1 ? -randomDouble() : -logl(randomDouble() )/ws.at(j), js.at(j)) );
			else
				ord.push_back( SampleEntry( ws[j] == -1 ? randomDouble() : logl( randomDouble() )/ws.at(j), js.at(j)) );
				
		}
		
		std::sort(ord.begin(), ord.end() );
		
		ordered.clear();
		ordered.reserve(ord.size() );
		
		for (long j = 0; j < ord.size(); ++j)
			ordered.push_back( ord.at(j).id);
	}
	
	
	inline void randomSort(const vector<long>& js, vector<long>& ordered){
	
		vector<SampleEntry> ord;
	
		ord.reserve(js.size() );
		
		for (long j = 0; j < js.size(); ++j)
			ord.push_back( SampleEntry( randomDouble(), js.at(j) ));
		
		std::sort(ord.begin(), ord.end() );
		
		ordered.clear();
		ordered.reserve(ord.size() );
		
		for (long j = 0; j < ord.size(); ++j)
			ordered.push_back( ord.at(j).id);
	}
	
	
	// get WR from population
	inline void getWRfromPOP(const vector<Weight>& weights, vector<int>& freqs, long m){
		
		vector<Weight> cmf;
		
		cmf.reserve(freqs.size() );
	
		cmf.push_back(0);
	
		for (long j = 0; j < freqs.size(); j++)
			cmf.push_back( weights.at(j) * freqs.at(j) );
		
		for (long j = 1; j < cmf.size(); j++)
			cmf.at(j) += cmf.at(j-1);
			
		for (long j = 0; j < freqs.size(); j++)
			freqs.at(j) = 0;
		
		for (long j = 0; j < m; j++){
		
			const long ind = (std::upper_bound(cmf.begin(), cmf.end(), randomDouble()*cmf.back() )-cmf.begin())-1;
		
			if (ind < freqs.size() )
				++freqs.at( ind );
			else
				++freqs.back();
		}
	}
	
	// get WOR from population
	inline void getWORfromPOP(const vector<Weight>& weights, vector<int>& freqs, long m){
		
		vector<Weight> cmf;
		
		cmf.reserve(freqs.size() );
	
		cmf.push_back(0);
	
		for (long j = 0; j < freqs.size(); j++)
			cmf.push_back( weights.at(j) * freqs.at(j) );
		
		for (long j = 1; j < cmf.size(); j++)
			cmf.at(j) += cmf.at(j-1);
			
		for (long j = 0; j < freqs.size(); j++)
			freqs.at(j) = 0;
		
		if (m >= freqs.size() ){
		
			for (long j = 0; j < freqs.size(); j++)
				freqs.at(j) = 1;
			
			return;
		}
		
		long s = 0;
		
		while (s != m){
		
			const long ind =  (std::upper_bound(cmf.begin(), cmf.end(), randomDouble()*cmf.back() )-cmf.begin())-1 ;
			
			if (ind < freqs.size() ){
				
				if (freqs.at(ind) == 0)
					s++;
		
				++freqs.at(ind);
				
			} else {
				
				if (freqs.back() == 0)
					s++;
		
				++freqs.back();
			}
			
			
		}
	}
	
	
	inline void getWRfromWOR(const Weight totalWeight, const vector<Weight>& weights, vector<int>& freqs, const vector<long>& ordered, const long m){
		
		if (weights.size() == 0){
			
			assert (m == 0);
			
			return;
		}
		
		/*
		if (UtilCollection::sum(weights) >= totalWeight){
		
			getWRfromPOP(weights, freqs, m);
			return;
		}*/
		
		
		//getWORfromWOR(totalWeight, weights, freqs, m);
		
		const long freqsum = UtilCollection::sum(freqs);
		
		//assert (m == freqsum );
		
		
		/*
		vector<long> coords;
		vector<Weight> ws;
		
		
		Weight smpWeight = 0;
		
		for (long j = 0; j < freqs.size(); ++j){
		
			assert (freqs.at(j) <= 1);
		
			const Weight w = weights.at(j);
		
			for (long freq = freqs.at(j); freq >= 1; --freq){
				
				coords.push_back(j);
				
				ws.push_back(w);
				
				smpWeight += w;
			}
		}
		
		for (long j = 0; j < ws.size(); ++j)
			ws.at(j) = exp( (log( ws.at(j) )/totalWeight)*smpWeight );
		
		cout << "ws_sum " << UtilCollection::sum(ws) << " pop " << totalWeight << " smpWeight " << smpWeight << endl;
		
		vector<long> ordered;
		
		ordered.reserve(ws.size() );
		
		for (long j = 0; j < ws.size(); ++j)
			ordered.push_back(j);
		
		//randomSort(coords, ordered);
		//randomSort(ws, coords, ordered, true);
		*/
		
		vector<Weight> cmf;
		
		cmf.reserve(ordered.size()+1);
		
		cmf.push_back(0);
		
		for (long j = 0; j < ordered.size(); j++)			
			cmf.push_back( weights.at(ordered.at(j)) );
		
		for (long j = 1; j < cmf.size(); j++)
			cmf.at(j) += cmf.at(j-1);
		
		for (long j = 0; j < freqs.size(); ++j)
			freqs.at(j) = 0;
		
		long comp = 0;
		
		for (long k = 0; k < m; ++k){
		
			long sel = k;
			
			if (comp > 0){
				
				const Weight u = randomDouble()*totalWeight;
			
				if (u < cmf.at(comp)){
					
					sel = (std::upper_bound(cmf.begin(), cmf.end(), u)-cmf.begin())-1;
				
					assert (sel < comp);
				}
			}
			
			if (sel < ordered.size() )
				++freqs.at(ordered.at(sel));
			else
				++freqs.at(ordered.back() );
			
			++comp;	
		}
	}

};	



/*
class SortedSample{



private:
	
	RandomNumbers rnd;
	
	// draw random integer between a and b, based on g(b+1-a)
	inline long realToInt(long_double logG, long a, long b) const{
	
		const long_double logW = logl(b+1-a);
		
		return UtilMath::minVal<long>(b, a+floorl(expl(logG+logW)) );
	}
	
	long big = -1;
	long small = -1;

	long prev = -1;
	
	long n = -1;
	long u = -1;
	
	long_double logG = -1;
	
	long k = -1;
	
	long u2 = -1;
	
	long next = -1;
	
	long nextPrev = -1;
	
	long pos = 0;

public:
	inline SortedSample(long size, long univ, const string seed="abc") : rnd(seed){
	
		n = size;
		u = univ;
		
		small = 0;
		big = 0;
		
		if (u > n){
		
			// [small=1..n][big=n+1...u]
			
			
		
			for (long i = 1; i <= n; ++i){
			
				if (rnd.randomInteger(i,u) <= n) // switching small with small => guaranteed to stay in small
					++small;	
				else
					++big;	// switching small with big => might go back to small
			}	
		}
		
		logG = 0;
	
		prev = u+1;
		
		u2 = n;
		
		next = u+1;
		
		++(*this);
	}
	
	inline SortedSample& operator++(){
	
		assert (next != -1);
		
		nextPrev = next;
		next = getNext();
	
		return (*this);
	}
	
	inline long operator*() const {
	
		return next == -1 ? -1 : (u-1)-(next-1);
	}
	
	inline long getSkip() const{
	
		assert (next != -1);
	
		return (nextPrev-next)-1;
	}
	
	inline operator bool() const {	 
	
    	return next != -1; 
	}

private:
	inline long getNext(){
	
		++pos;
	
		while (big > 0){
		
			logG += logl( rnd.randomDouble() )/big;
		
			--big;
		
			const long s = realToInt(logG, n+1, u);
			
			if (s < prev){
			
				prev = s;
				return s;
			}
			
			++small;
		}
		
		// draw 'small' balls from urn with 'u2' without replacement
		
		assert (small <= u2);
		
		while (small > 0){
		
			// pr of drawing any ball is number_of_draws / balls
			
			if (u2 == small || (rnd.randomDouble() <= (((double) small)/u2)) ){

				const long s = u2;
				
				--u2; // remove one ball
				--small; // move it into sample
				
				return s;
			
			} else {
			
				--u2; // remove one ball
			}
		}
		
		//assert (pos == n);
		
		assert (small == 0);
		
		return -1;
	}


};*/



/*
class SortedSample{



private:
	
	RandomNumbers rnd;
	
	// draw random integer between a and b, based on g(b+1-a)
	inline long realToInt(long_double logG, long a, long b) const{
	
		const long_double logW = logl(b+1-a);
		
		return UtilMath::minVal<long>(b, a+floorl(expl(logG+logW)) );
	}
	
	long h = -1;
	long l = -1;

	long prev = -1;
	
	long n = -1;
	long u = -1;
	
	long_double logG = -1;
	
	long next = -1;
	
	long nextPrev = -1;
	
	long pos = 0;

public:
	inline SortedSample(long size, long univ, const string seed="abc") : rnd(seed){
	
		n = size;
		u = univ;
		
		small = 0;
		big = 0;
		
		if (u > n){
		
			// [small=1..n][big=n+1...u]
			
			
		
			for (long i = 1; i <= n; ++i){
			
				if (rnd.randomInteger(i,u) <= n) // switching small with small => guaranteed to stay in small
					++small;	
				else
					++big;	// switching small with big => might go back to small
			}	
		}
		
		logG = 0;
	
		prev = u+1;
		
		u2 = n;
		
		next = u+1;
		
		++(*this);
	}
	
	inline SortedSample& operator++(){
	
		assert (next != -1);
		
		nextPrev = next;
		next = getNext();
	
		return (*this);
	}
	
	inline long operator*() const {
	
		return next == -1 ? -1 : (u-1)-(next-1);
	}
	
	inline long getSkip() const{
	
		assert (next != -1);
	
		return (nextPrev-next)-1;
	}
	
	inline operator bool() const {	 
	
    	return next != -1; 
	}

private:
	inline long getNext(){
	
		++pos;
	
		while (big > 0){
		
			logG += logl( rnd.randomDouble() )/big;
		
			--big;
		
			const long s = realToInt(logG, n+1, u);
			
			if (s < prev){
			
				prev = s;
				return s;
			}
			
			++small;
		}
		
		// draw 'small' balls from urn with 'u2' without replacement
		
		assert (small <= u2);
		
		while (small > 0){
		
			// pr of drawing any ball is number_of_draws / balls
			
			if (u2 == small || (rnd.randomDouble() <= (((double) small)/u2)) ){

				const long s = u2;
				
				--u2; // remove one ball
				--small; // move it into sample
				
				return s;
			
			} else {
			
				--u2; // remove one ball
			}
		}
		
		//assert (pos == n);
		
		assert (small == 0);
		
		return -1;
	}


};
*/

class SortedSample{

private:
	
	RandomNumbers rnd;
	
	long h = -1;
	long l = -1;

	long prev = -1;
	
	long n = -1;
	long u = -1;
	
	double logG = -1;
	
	long next = -1;
	
	long nextPrev = -1;
	
	long pos = 0;
	
	long sampleSize;

public:
	inline SortedSample(long size, long univ, const string seed="abc") : rnd(seed){
	
		sampleSize = size;
		n = size;
		u = univ;
		
		h = 0;
		l = n;
		
		if (u > n){
		
			// [small=1..n][big=n+1...u]
			
			double u_n = (u-n);
		
			for (long i = 0; i < n; ++i){
				if (rnd.randomDouble()*(u-i) <= u_n ){
					++h;
					--l;
				}
			}
		}
		
		logG = 0;
		
		prev = u+1;
		next = u+1;
		
		++(*this);
	}
	
	inline SortedSample& operator++(){
	
		assert (next != -1);
		
		nextPrev = next;
		next = getNext();
	
		return (*this);
	}
	
	inline long operator*() const {
	
		return next == -1 ? -1 : (u-1)-(next-1);
	}
	
	inline long getSkip() const{
	
		assert (next != -1);
	
		return (nextPrev-next)-1;
	}
	
	inline operator bool() const {	 
	
    	return next != -1; 
	}

private:
	inline long getNext(){
	
		++pos;
	
		while (h > 0){
		
			logG += logl( rnd.randomDouble() )/h;
		
			--h;
		
			const long s = n+UtilMath::minVal<long>(u-n-1, floor( exp(logG)*(u-n) ));
			
			if (s < prev){
			
				prev = s;
				return s;
			}
			
			++l;
		}
		
		
		while (l > 0){
		
			// pr of drawing any ball is number_of_draws / balls
			
			if (n == l || (rnd.randomDouble() <= (((double) l)/n)) ){

				const long s = n-1;
				
				--l; // remove one ball
				--n; // move it into sample
				
				return s;
			
			} else {
			
				--n; // remove one ball
			}
		}
	

		assert (pos == sampleSize+1);
		
		
		return -1;
	}


};




//template<typename T1, typename T2> using hash_map = robin_hood::unordered_map<T1,T2>;

//template<typename T1> using hash_set = robin_hood::unordered_set<T1>;

template<typename T1, typename T2> using hash_map = std::unordered_map<T1,T2>;
template<typename T1> using hash_set = std::unordered_set<T1>;

typedef constvector<string> Tuple;

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}





class Row{

private:
	const char* s = NULL;

	offarray off;
	
public:
	
	inline const char* begin() const{
	
		return s;	
	}
	
	inline Row(){
	
		off.fill(0);
		setSortKey(-1);
	}
	
	inline bool empty() const {
	
		return s == NULL;
	}
	
	
	inline void set(){
	
		s = NULL;
		off.fill(0);
		setSortKey(-1);
	}
	
	inline Row(const Row& r){// : s(r.s), off(r.off){
	
		s = r.s;
		off = r.off;
	}
	
	inline void setSortKey(int i = -1){
	
		off[off.size()-1] = (i == -1) ? std::numeric_limits<uint8_t>::max() : i;
	}
	
	inline const char* getSortKey() const{
	
		return get(off[off.size()-1]);
	}
	
	
	inline bool stringless(const char* p1, int l1, const char* p2, int l2) const{
	
		return l1 < l2 ? true : l1 > l2 ? false : (*p1) < (*p2) ? true : (*p1) > (*p2) ? false : strcmp(p1, p2) < 0;
	}
	
	inline bool stringeqPossible(const char* p1, const int l1, const char* p2, const int l2) const{
		
		return l1 == l2 ? (*p1) == (*p2) : false;
	}
	
	inline int stringcmp(const char* p1, const int l1, const char* p2, const int l2) const{
	
		return l1 < l2 ? -1 : l1 > l2 ? +1 : (*p1) < (*p2) ? -1 : (*p1) > (*p2) ? +1 : strcmp(p1, p2);
	}
	
	inline bool stringeq(const char* p1, const int l1, const char* p2, const int l2) const{
		
		return stringeqPossible(p1,l1,p2,l2) && strcmp(p1, p2) == 0;
	}
	
	
	inline bool less(const char* ptr1, const char* ptr2) const{
	
	
		if (true)
			return strcmp(ptr1, ptr2) < 0;
		
		const char cc1 = *ptr1;
		const char cc2 = *ptr2;
		
		if (cc2 == 0){
		
			return false;
			
		} else if (cc1 == 0){
		
			return true;
		
		} else if ( cc1 < cc2){
			
			return true;
			
		} else if ( cc2 < cc1){
		
			return false;
			
		}/* else if ( (cc1 == 0) && (cc2 == 0) ){
		
			return false;
			
		}*/ else {
		
			for (int q = 0; q < 1001; ++q){
			
				assert (q < 1000);
			
				++ptr1;
				++ptr2;	
			
				const char c1 = (*ptr1);
				const char c2 = (*ptr2);
		
				if (c2 == 0)
					return false;
		
				if (c1 == 0)
					return true;
				
				if (c1 < c2)
					return true;
				
				if (c2 < c1)
					return false;
				
				assert ( c1 != 0);
				assert ( c2 != 0);
			}
			
			return false;
		}
		
	}
	
	bool less(const vector<int>& sortKeys1, const Row& row, const vector<int>& sortKeys2) const{
	
		assert (sortKeys1.size() == sortKeys2.size() );
	
		bool eq = false;
	
		for (int i = 0; i < sortKeys1.size(); i++){
		
			const int i1 = sortKeys1.at(i);
			
		
			if (i1 == -1) // SKIPSKIP
				continue;
			
			const int i2 = sortKeys2.at(i);
			
			if (i2 == -1){
				eq = true;
				continue;
			}
			
			const char* ptr1 = get(i1);
			const char* ptr2 = row.get(i2);
			
			const int len1 = getLength(i1);
			const int len2 = row.getLength(i2);
			
			const int c = stringcmp(ptr1, len1, ptr2, len2);
			
			if (c < 0)
				return true;
				
			if (c > 0)
				return false;
			
			/*
			const char* ptr1 = get( sortKeys1.at(i) );
			const char* ptr2 = row.get( sortKeys2.at(i) );
			
			if ( less(ptr2, ptr1) )
				return false;
			
			if ( less(ptr1, ptr2) )
				return true;*/
		}
		
		return eq;
		
		/*
		setSortKey(-1);
		row.setSortKey(-1);
		
		return (*this) < row;*/
	}
	
	
	
	inline int cmp(vector<int>& sortKeys1, Row& row, vector<int>& sortKeys2) const{
	
		return equals(sortKeys1, row, sortKeys2) ? 0 : less(sortKeys1, row, sortKeys2) ? -1 : +1;
	}
	
	bool operator<(const Row& row) const {
	
		//assert (false);
	
		if (off[off.size()-1] == std::numeric_limits<uint8_t>::max() ){
		
			//for (int i = size()-1; i >= 0 ; --i){
			
			/*
			if ( (*s) < (*(row.s)))
				return true;
			
			if ( (*(row.s)) < (*s))
				return false;*/
			
			for (int i = 0; i < size(); ++i){
			
				/*
				const char* ptr1 = get(i);
				const char* ptr2 = row.get(i);
				
				const int len1 = get(i);
				const int len2 = row.get(i);
				
				const int c = stringcmp(ptr1, len1, ptr2, len2);
				
				if (c < 0)
					return true;
					
				if (c > 0)
					return false;
				*/
				
				const char* ptr1 = get(i);
				const char* ptr2 = row.get(i);
				
				if (less(ptr1, ptr2)){
				
					return true;
					
				} else if (less(ptr2, ptr1)){
					
					return false;
				}
			}
			
			return false;
		}
		
		const char* ptr1 = getSortKey();
		const char* ptr2 = row.getSortKey();
		
		return less(ptr1, ptr2);
    }
	
	inline long size() const{
	
		return off[0];
	}
	
	
	inline long getBytes() const{
	
		return end()-begin();
	}
	
	inline long getBytes(int i) const{
	
		return get(i+1)-get(i);
	}
	
	
	inline void set(const char* begin, offarray offset){
		
		s = begin;
		off = offset;
		setSortKey(-1);
	}
	
	
	inline const Tuple getTuple() const{
	
		Tuple ret;
		
		write(ret);
		
		return ret;
	}
	
	inline void set(const Row& row){
	
		s = row.s;
		off = row.off;
	}
	
	// n = number of columns
	inline const char* set(int n, const char* snew, const char* end){

		off.fill(0);
		off[0] = n;
		
		assert (n < off.size());
		
		s = snew;
		
		const char* ss = snew;
		
		for (int i = 1; i <= n; i++){
		
			assert (ss != end);
			
			while ((*ss) != 0){
			
				ss++;
				assert (ss != end);
			}
			
			assert ( (*ss) == 0);	
			ss++;
			
			assert ( (ss-s) < 255);
			
			off[i] = (ss-s);
		}
		
		return ss;
	}
	
	inline const char* get(int n) const{
	
		//cout << "off[" << n << "] = " << ((int) off[n]) << "#" << ((long)(s)) << "#";
	
		assert (n >= 0);
		
		assert (n <= size());
		
		return n == 0 ? s : s+off[n];
	}
	
	inline const char* end() const{
	
		return s+off[size()];
	}
	
	
	
	inline const string operator[](int n) const{
	
		if ( (*get(n) == 1) || ( (*get(n) == 2) ) )
			return to_string( getLong(get(n) ));
		
		return string(get(n));
	}

	inline bool equals(const Row& t2) const{
	
		if (size() != t2.size() )
			return false;
	
		for (int j = size()-1; j >= 0; --j){
		
			if (!equals(t2, j))
				return false;
		}
		
		return true;
	}
	
	
	inline bool equalsPossible(int i, const Row& t2, int j) const{	
		
		return stringeqPossible(get(i), getLength(i), t2.get(j), t2.getLength(j));
	}
	
	inline bool equals(int i, const Row& t2, int j) const{	
		
		assert (i >= 0);
		assert (i < size() );
		
		assert (j >= 0);
		assert (j <= t2.size() );
		
		return stringeq(get(i), getLength(i), t2.get(j), t2.getLength(j));
		
		/*
		if (true)
			return strcmp(ptr1, ptr2) == 0;
		
		while (true){
		
			if ( (*ptr1) != (*ptr2))
				return false;
		
			if ( (*ptr1) == 0 && (*ptr2) == 0 )
				return true;
			
			assert ( (*ptr1) != 0);
			assert ( (*ptr2) != 0);
			
			++ptr1;
			++ptr2;
		}
		
		return true;*/
	}
	
	// r[v1[j]] = r2[v2[j]]
	inline bool equals(const vector<int>& v1, const Row& t2, const vector<int>& v2) const{
	
		if (v1.size() != v2.size() )
			return false;
			
		for (int j = v1.size()-1; j >= 0; --j){
		
			if (v1.at(j) == -1 || v2.at(j) == -1)  // SKIPSKIP
				continue;
			
			if (!equalsPossible(v1.at(j), t2, v2.at(j) ))
				return false;
		}
		
		for (int j = v2.size()-1; j >= 0; --j){
		
			if (v1.at(j) == -1 || v2.at(j) == -1)  // SKIPSKIP
				continue;
			
			if (!equals(v1.at(j), t2, v2.at(j) ))
				return false;
		}
		
		return true; // no difference found
	}
	
	/*
		for (auto it = v.begin(); it != v.end(); ++it){
		
			if (!equals((*it), t2, (*it)))
				return false;
		}
		
		return true;
	}*/
	
	// r[j] = r2[v[j]]
	inline bool equals(const Row& t2, const vector<int>& v) const{
	
		for (auto it = v.begin(); it != v.end(); ++it){
		
			if (!equals((*it), t2, (*it)))
				return false;
		}
		
		return true;
	}

	// r[j] = r2[j]
	inline bool equals(const Row& t2, int j) const{
		
		return equals(j, t2, j);
		/*
		const char* ptr1 = get(j);
		const char* ptr2 = t2.get(j);
		
		while (true){
		
			if ( (*ptr1) != (*ptr2))
				return false;
		
			if ( (*ptr1) == 0 && (*ptr2) == 0 )
				return true;
				
			assert ( (*ptr1) != 0);
			assert ( (*ptr2) != 0);
			
			++ptr1;
			++ptr2;
		}
		
		return true;*/
	}
	
	inline int getOff(int i) const{
	
		return i == 0 ? 0 : off[i];
	}
	
	inline int getLength(int i) const{
	
		return getOff(i+1)-getOff(i);
	}
	
	// r[i] == r[j]
	inline bool equals(int i, int j) const{
		
		if (getLength(i) != getLength(j))
			return false;
		
		return strcmp(get(i), get(j)) == 0;
		
		/*
		
		const char* ptr1 = get(i);
		const char* end1 = get(i+1);
		
		const char* ptr2 = get(j);
		const char* end2 = get(j+1);
		
		for ( ; (ptr1 != end1) && (ptr2 != end2) ; ++ptr1, ++ptr2){
		
			if ( (*ptr1) != (*ptr2) )
				return false;
		}
		
		return true;*/
	}
	
	// map integer strings to integers and one-char strings to 2^62-char
	inline static long getLong(const char* s){
	
		const char c0 = (*s);
	
		if (c0 == 1)
			return UtilBytePacking::unpack<int32_t>(s);
	
		if (c0 == 2)
			return UtilBytePacking::unpack<int64_t>(s);
	
		
		if (c0 == 0)
  			return 1L << 62;
  		
  		long ret = 0;
  		
  		if (c0 > '9')
			ret -= (int) c0;
		else if (c0 < '0')
			ret -= (int) c0;
		else
			ret += (int)(c0-'0');
  		
  		++s;
  		
  		if ((*s) == 0){
  		
  			if (ret < 0)
  				return (1L << 62)+ret;
  			
  			return ret;
  		}
  		
  		if (ret <= 0)
  			return -1;
  		
  		for (; (*s) != 0; ++s){
    		
    		ret *= 10;
    		
    		if ((*s) > '9')
    			return -1;
    		else if ((*s) < '0')
    			return -1;
    		
    		ret += (int)((*s)-'0');
  		}
		
		return ret;
	}
	
	inline long getLong(int n) const{
		
		assert (n >= 0);
		assert (n < size());
		
  		return getLong(get(n));
	}
	
	inline bool isNull(const vector<int>& v) const{
	
		if (v.size() == 0)
			return true;
			
		for (auto it = v.begin(); it != v.end(); ++it){
		
			if (isNull(*it) )
				return true;
		}
		
		return false;
	}	
	
	inline bool isNull(int j) const {
	
		const char* ch = get(j);
		
		if ( (*ch) != 'N' &&(*ch) != 'n' )
			return false;
		
		++ch;
		
		if ( (*ch) != 'U' &&(*ch) != 'u' )
			return false;
			
		++ch;
		
		if ( (*ch) != 'L' &&(*ch) != 'l' )
			return false;
		
		++ch;
		
		if ( (*ch) != 'L' &&(*ch) != 'l' )
			return false;
			
		return true;	
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const Row& v) { 
    	
    	for (int i = 0; i < v.size(); i++){
			
			if (i > 0)
				os << ",";
			
			os << "'" << v[i] << "'"; // v.get(i)
			
			/*
			
			for (const char* ptr = v.get(i); (*ptr) != 0; ++ptr){
			
				os << (*ptr);
			}*/
					
		}
		
		return os;
	}
	
	inline void write(Tuple& t) const{
	
		for (int i = 0; i < size(); ++i){
			
			t.push_back(  string(get(i))  );
		}
	}
	
};



template<typename E>
class RowHashMap{

private:


	long toReserve = -1;
	
	long loadFactor = -1;
	
	string stringKey;
	string hashTarget;
	
	unique_ptr<std::bitset<100000>> filter;

public:
	
	unordered_map<string, E> stringMap;
	long_map<long, E> longMap;

	
	
	inline RowHashMap(long size = -1){
	
		
		if (size > 0)
			reserve(size);
			
	}
	
	
	inline long getSize() const {
	
		return stringMap.size()+longMap.size();
	}
	
	inline void print(){
	
		for (auto it = stringMap.begin(); it != stringMap.end(); ++it){
		
			cout << it->first << " " << it->second << endl;
		}
		
		for (auto it = longMap.begin(); it != longMap.end(); ++it){
		
			cout << it->first << " " << it->second << endl;
		}
	
	}
	inline void reserve(long m, long lf=-1){
	
		toReserve = m;
		loadFactor = lf;
		
		if (getSize() == 0){
			
			if (m > 500000){
		
				filter.reset( new std::bitset<100000>() );
				filter->reset();
			}
		}
		
	}
	
	inline void clear(){
	
		stringMap.clear();
		longMap.clear();
		onlyIntegers = true;
		
		if (filter)
			filter->reset();
	}
	
	bool onlyIntegers = true;
	
	
	inline const string& combine(const Row& r, const vector<int>& projection, const vector<long>& hashings){
		
		assert (!ONLYINTEGERS);
		assert (!onlyIntegers);
		
		stringKey.reserve(1024);
		hashTarget.reserve(1024);
		
		stringKey.clear();
		
		for (int i = 0; i < projection.size(); i++){
		
			if (i > 0)
				stringKey.push_back('`');
			
			const long key = getLongKey(r, projection[i], hashings[i]);
			
			if (key >= 0){
				
				const long mask = ((1L << 3)-1);
				
				for (long hv = key; hv != 0; hv >>= 3){
					
					stringKey.push_back( (char) (((int)'0')+(hv & mask)) );
				}
				
				//for (long hv = key; hv != 0; hv /= 10)
				//	stringKey.push_back( '0'+(hv % 10) );
					
			} else {
			
				for (const char* ptr = r.get(projection[i]); (*ptr) != 0; ++ptr)
					stringKey.push_back( (*ptr) );
			}
		}
		
		stringKey.push_back(0);
		
		//cout << "stringKey " << stringKey << endl;
		
		return stringKey;
	}
	
	inline long getLongKey(const Row& row, int proj, long h){
	
		const long key = row.getLong(proj);
			
		if (key >= 0)
			return h < UtilHash::HASHMAX ? UtilHash::hash(key, h): key;
					
		return h < UtilHash::HASHMAX ? UtilHash::hash( row.get(proj), h) : -1;
	}
	
	inline long getLongKey(const Row& row, const vector<int>& projection, const vector<long>& hashings){
	
		if (projection.size() == 1)
			return getLongKey(row, projection[0], hashings[0]);

		return -1;
	}
	
	inline void insert(const Row& row, const vector<int>& projection, const vector<long>& hashings, const E e){
	
		if (row.isNull(projection))
			return;
	
		//assert (!hasKey(row, projection, hashings));
		
		const long key = getLongKey(row, projection, hashings);
		
		
		if (key >= 0){
		
			if (toReserve != -1){
			
				//if (loadFactor != -1)
				//	longMap.max_load_factor(loadFactor);
					
				longMap.reserve(toReserve);
				toReserve = -1;
				loadFactor = -1;
			}
			
			//std::make_pair(key, e) 
			longMap.insert( {key,e});
			
			//longMap.insert(key, e);
			
			if (filter)
				filter->set(key % filter->size());
			
			return;
		}
		
		assert (!ONLYINTEGERS);
		onlyIntegers = false;
		
		if (toReserve != -1){
			
			//if (loadFactor != -1)
			//	stringMap.max_load_factor(loadFactor);
			
			stringMap.reserve(toReserve);
			toReserve = -1;
			loadFactor = -1;
		}
		
		const string& s = combine(row, projection, hashings);
		
		if (filter)
			filter->set(UtilHash::hash(s, filter->size()) );
		
		stringMap.insert( {s, e} );
	}
	
	
	inline bool hasKey(const Row& row, const vector<int>& projection, const vector<long>& hashings, long longkey=-1) {
	
		if (row.isNull(projection))
			return false;
	
		const long key = longkey != -1 ? longkey: getLongKey(row, projection, hashings);
		
		if (key >= 0){
			
			
			if (filter){
				if (!filter->test(key % filter->size()))
					return false;
			}
		
			return longMap.count(key) > 0;
		}
		
		if (ONLYINTEGERS || onlyIntegers)
			return false;
			
		
		const string& s = combine(row, projection, hashings);
		
		if (filter){
			if (!filter->test(UtilHash::hash(s, filter->size()) ))
				return false;
		}
		
		return stringMap.count(s) > 0;
	}
	
	inline E& get(const string& key){
		
		return stringMap[key];
	}
	
	inline E& get(long key){
		
		return longMap[key];
	}
	
	inline const string getKey(const Row& row, const vector<int>& projection, const vector<long>& hashings){
	
		if(!hasKey(row, projection, hashings))	
			return "NULL";
	
		const long key = getLongKey(row, projection, hashings);
		
		if (key >= 0)
			return "long("+ to_string(key)+")";
			
		assert (!ONLYINTEGERS);
		onlyIntegers = false;
		
		return "string("+combine(row, projection, hashings)+")";
	}
	
	
	inline E& get(const Row& row, const vector<int>& projection, const vector<long>& hashings, long longkey=-1){
	
		//assert (hasKey(row, projection, hashings));
	
		const long key = longkey != -1 ? longkey: getLongKey(row, projection, hashings);
	
		if (key >= 0)
			return longMap[key];
		
		assert (!ONLYINTEGERS);
		onlyIntegers = false;
		
		return stringMap[combine(row, projection, hashings)];
	}
	
	inline void set(const Row& row, const vector<int>& projection, const vector<long>& hashings, const E& val){
	
		if (row.isNull(projection))
			return;
	
		//assert (hasKey(row, projection, hashings));
	
		const long key = getLongKey(row, projection, hashings);
		
		if (key >= 0){
			longMap[key] = val;
			
			if (filter)
				filter->set(key % filter->size() );
			
			return;
		}
		
		assert (!ONLYINTEGERS);
		onlyIntegers = false;
		
		string s = combine(row, projection, hashings);
		stringMap[s] = val;
		
		if (filter)
			filter->set(UtilHash::hash(s, filter->size()) ); 
	}
};

class LongRowHashMap : public RowHashMap<long>{

public:
	inline LongRowHashMap(long size = -1) : RowHashMap<long>(size){
	
		
	}
	
	inline long getCount(const Row& row, const vector<int>& projection, const vector<long>& hashings){
	
		const long lk = this->getLongKey(row, projection, hashings);
	
		if (!this->hasKey(row, projection, hashings, lk))
			return 0;
		
		return this->get(row, projection, hashings, lk);
	}
	
	inline void addCount(const Row& row, const vector<int>& projection, const vector<long>& hashings, long v=1){
		
		if (v == 0)
			return;
			
		const long lk = this->getLongKey(row, projection, hashings);
					
		if (!this->hasKey(row, projection, hashings, lk) )
			this->insert(row, projection, hashings, v);
		else
			this->get(row, projection, hashings, lk) += v;
		
		/*
		if (!this->hasKey(row, projection, hashings))
			this->insert(row, projection, hashings, v);
		else
			this->get(row, projection, hashings) += v;*/
	}
};



class Rows{

	
	vector<char> v;
	long rows = 0;
	bool open = false;
	
	//Row r;
	
	bool locked = false;
	
public:	
	const int cols;

	inline Rows (int cols) : cols(cols){
	
	
	}
	
	inline void reserve(long m){
	
		
	
		if (m >= v.capacity() ){
		
			cout << "reserve " << m << " (" << UtilString::bytesToString(m) << ")" <<  endl;
		
			v.reserve(m);
			
			cout << "/reserve " << m << " (" << UtilString::bytesToString(m) << ")" <<  endl;
		}
	}
	
	inline void reserve(const string& s, long m){		
	
		if (m >= v.capacity() ){
		
			cout << "Rows " << s << " reserve " << m << " (" << UtilString::bytesToString(m) << ")" <<  endl;
		
			v.reserve(m);
			
			cout << "/Rows " << s << " reserve " << m << " (" << UtilString::bytesToString(m) << ")" <<  endl;
		}
	}
	
	inline long capacity() const{
	
		return v.capacity();
	}
	
	
	inline Rows (const Rows& r) : cols(r.cols), rows(r.rows){
	
		v = r.v;
	}
	
	inline void shrink_to_fit(){
	
		assert (!locked);
	
		if (v.capacity() > 1.5*v.size() ){
			
			cout << "shrink " << v.capacity() << " to " << 1.5*v.size() << endl;
			v.shrink_to_fit();
		}
	
		
	}
	
	inline void lock(){
	
		locked = true;
	}
	
	inline void unlock(){
	
		locked = false;
	}
	
	inline void permute(const vector<long>& posVec){
	
		assert (!locked);
	
		long check = v.size();
	
		Rows rows3(cols);
		
		rows3.reserve(v.size() );
		
		{
			vector<Row> rows2(rows);
		
			long k = 0;
			for (auto it = begin(); it != end(); ++it){
		
				rows2.at(posVec.at(k)) = (*it);
				++k;
			}
		
			for (auto it = rows2.begin(); it != rows2.end(); ++it)
				rows3.push_back(*it);
		}
		
		clear();
		
		reserve(rows3.size() );
		
		for (auto it = rows3.begin(); it != rows3.end(); ++it)
			push_back(*it);
		
		assert (v.size() == check);
	}
	
	
	inline void getRows(vector<Row>& v) {
	
		lock();
	
		v.reserve(rows);
		
		for (auto it = begin(); it != end(); ++it)
			v.push_back( *it); // Row(*it) 
	}
	
	inline shared_ptr<vector<Row>> getRows(){
	
		shared_ptr<vector<Row>> ret = std::make_shared<vector<Row>>();
		
		getRows(*ret);
	
		return ret;
	}
	
	inline int getCols() const{
	
		return cols;
	}
	
	inline void startRow(){
	
		assert (!open);
		
		assert (!locked);
		
		open = true;
	}
	
	inline void push_back_char(const char& c){
	
		assert (!locked);
		assert (open);
		//assert (c != 1);
		
		if (v.capacity() == v.size()){
		
			cout << "rows with " << cols << " columns reaches capacity " << (v.size()) << " (" << UtilString::bytesToString(v.size() ) << ")" << endl;
		
			cout << "reserve " << (v.size()*1.5) << " (" << UtilString::bytesToString(v.size()*1.5) << ")" << endl;
		
			v.reserve(v.size()*1.5);
			
			cout << "reserved " << (v.capacity() ) << " (" << UtilString::bytesToString(v.capacity() ) << ")" << endl;
			
			//assert (false);
		}
		
		v.push_back(c);
	}
	
	
	inline void push_back_string(const string& s){
	
		assert (!locked);
	
		assert (open);
	
		for (int i = 0; i < s.length(); i++)
			push_back_char(s[i]);
			
		push_back_char(0);
	}
	
	inline void push_back_string(const char* s){
	
		assert (!locked);
		
		assert (open);
	
		for ( const char* ptr = s ; (*ptr) != 0 ; ++ptr)
			push_back_char( (*ptr) );
		
		push_back_char(0);
	}
	
	
	inline void push_back_chars(const char* s){
	
		assert (!locked);
		
		assert (open);
	
		for ( const char* ptr = s ; (*ptr) != 0 ; ++ptr)
			push_back_char( (*ptr) );
	}
	
	
	
	inline void finishRow(){
	
		assert (open);
		++rows;
		open = false;
	}
	
	inline const char* beginPtr() const{
	
		if (v.size() == 0)
			return NULL;
	
		return v.data();
	}
	
	inline const char* endPtr() const{
	
		if (v.size() == 0)
			return NULL;
	
		return v.data()+v.size();
	}
	
	inline long size() const{
	
		return rows;
	}
	
	inline long getSize() const{
	
		return rows;
	}
	
	inline long getBytes() const{
	
		return v.size();
	}
	
	
	inline const char* front() const{
	
		if (v.size() == 0)
			return NULL;
	
		return v.data();
	}
	
	inline void front(Row& r) const{
	
		r.set(cols, front(), v.data()+v.size() );
	}
	
	inline void back(Row& r) const{
	
		r.set(cols, back(), v.data()+v.size() );
	}
	
	
	inline void pop_back(){
	
		const char* b = back();
	
		v.resize( b-v.data() );
		rows--;
	}
	
	inline const char* back() const{
	
		if (rows <= 1)
			return front();
		
		const char* ptr = v.data()+v.size()-1;
		
		assert ( (*ptr) == 0);
		
		for (int col = 0; col < cols; col++){
		
			--ptr;
		
			while ((*ptr) != 0)
				--ptr;
		}
		
		return ptr+1;
	}
	
	inline void remove(vector<long>& idsToRemove){
		
		remove(idsToRemove, idsToRemove, false);
	}
	
	inline void remove(vector<long>& ids, vector<long>& idsToRemove, const bool distinct = true){
	
		std::sort(idsToRemove.begin(), idsToRemove.end() );
	
		const char* readPtr = v.data();
		char* writePtr = v.data();
		
		int writej = 0;
		long newSize = 0;
		
		int m = 0;
		
		long lastId = -1;
		
		for (long readj = 0; readj < ids.size(); readj++){
			
			const long id = distinct ? ids.at(readj) : readj;
			
			assert (id >= lastId);
			
			lastId = id;
				
			if (m < idsToRemove.size() && lastId == idsToRemove.at(m)){
				
				for (int i = 0; i < cols; ){
				
					if ( (*readPtr) == 0)
						++i;
					
					++readPtr;
				}
				
				++m;
				
			} else {
			
				assert (m == idsToRemove.size() || lastId < idsToRemove.at(m));
				
				for (int i = 0; i < cols; ){
				
					if ( (*readPtr) == 0)
						++i;
					
					(*writePtr) = (*readPtr);
					
					++newSize;
					
					++readPtr;
					++writePtr;
				}
				
				if (distinct)
					ids.at(writej) = ids.at(readj);
				writej++;
			}
		}
		
		idsToRemove.clear();
		
		if (distinct)
			ids.resize(writej);
			
		v.resize(newSize);
		rows = writej;
	}
	
	
	
	
	inline void push_back_parse(const Row& r){
	
	
		assert (r.size() == cols);
	
		startRow();
	
		for (int i = 0; i < r.size(); ++i){
		
			const long l = r.getLong(i);
		
			/*
			if (l > 99999999L){
			
				const int len = UtilBytePacking::packLength<int64_t>(l);
			
				for (int ii = 0; ii < len; ++ii)
					push_back_char(0);				
				
				char* chptr = v.data()+(v.size()-len );
				
				UtilBytePacking::pack<int64_t>(chptr, l);
				
				//(*chptr) = 2;
				
				//assert ( (l > 0 && (*chptr) == 1) || (l == 0 && (*chptr) == 2) );
				
				for (int ii = 0; ii < len; ++ii){
				
					assert ( (*chptr) != 0);
					
					++chptr;
				}
				
				push_back_char(0);
			
			} else if (l > 9999L){
		
				const int len = UtilBytePacking::packLength<int32_t>(l);
		
				for (int ii = 0; ii < len; ++ii)
					push_back_char(0);
				
				char* chptr = v.data()+(v.size()-len );
				
				UtilBytePacking::pack<int32_t>(chptr, l);
				
				//(*chptr) = 1;
				//assert ( (l > 0 && (*chptr) == 1) || (l == 0 && (*chptr) == 2) );
				
				for (int ii = 0; ii < len; ++ii){
				
					assert ( (*chptr) != 0);
					
					++chptr;
				}
				
				push_back_char(0);
			
			} else*/ {
			
				const char* st = r.get(i);
				const char* en = r.get(i+1);
			
				for (const char* ptr = st; ptr != en; ++ptr){
				
					push_back_char(*ptr);
				}
			}
		}
		
		finishRow();
	}
	
	inline void push_back(const Row& r){
	
		assert (r.size() == cols);
		
		startRow();
	
		for (const char* ptr = r.begin(); ptr != r.end(); ++ptr){
		
			push_back_char( (*ptr) );
		}
		
		finishRow();
	}
	
	inline void push_back(vector<const Row*>& k_rows, vector<vector<int>>& k_sels){
	
		startRow();
	
		for (int k = 0; k < k_rows.size(); ++k){
		
			const vector<int>& sels = k_sels.at(k);
			
			const Row& row = *k_rows.at(k);
		
			for (int col = 0; col < sels.size(); ++col){
			
				push_back_string( row.get( sels.at(col) ));
			}
		}
	
		finishRow();
	}
	
	
	inline void push_back(const Row& r1, bool activate1, const vector<int>& coords1,  const Row& r2, bool activate2, const vector<int>& coords2){
		
		assert ( (activate1 ? coords1.size() : r1.size())+( activate2 ? coords2.size() : r2.size() ) == cols);
		
		startRow();
		
		if (!activate1){
		
			for (const char* ptr = r1.begin(); ptr != r1.end(); ++ptr)
				push_back_char( (*ptr) );
			
		} else {
		
			for (int j = 0; j < coords1.size(); ++j)
				push_back_string( r1.get( coords1.at(j) ));
		}
		
		if (!activate2){
		
			for (const char* ptr = r2.begin(); ptr != r2.end(); ++ptr)
				push_back_char( (*ptr) );
			
		} else {
		
			for (int j = 0; j < coords2.size(); ++j)
				push_back_string( r2.get( coords2.at(j) ));
		}
		
		finishRow();
	}
	
	inline void push_back(const Row& r1, const vector<int>& coords1, const Row& r2, const vector<int>& coords2){
	
		push_back(r1, true, coords1, r2, true, coords2);
	}
	
	inline void push_back(const Row& r1,  const Row& r2, const vector<int>& coords2){
	
		push_back(r1, false, coords2, r2, true, coords2);
	}
	
	inline void push_back(const Row& r1, const vector<int>& coords1,  const Row& r2){
	
		push_back(r1, true, coords1, r2, false, coords1);
	}
	
	
	inline void push_back(const Row& r1, const Row& r2){
	
		assert (r1.size()+r2.size() == cols);
	
		startRow();
	
		for (const char* ptr = r1.begin(); ptr != r1.end(); ++ptr)
			push_back_char( (*ptr) );
		
		for (const char* ptr = r2.begin(); ptr != r2.end(); ++ptr)
			push_back_char( (*ptr) );
		
		finishRow();
	}
	
	
	inline void push_back(const Tuple& t){
	
		assert (t.size() == cols);
		
		startRow();
	
		for (int i = 0; i < cols; i++){
		
			push_back(t[i]);
		}
		
		finishRow();
	}
	
	
	inline void push_back(const Row& r, const vector<int>& v){
	
		assert (v.size() == cols);
		
		startRow();
	
		for (int i = 0; i < cols; i++){
		
			push_back_string(r.get( v.at(i) ) );
		}
		
		finishRow();
	}
	
	inline void push_back(const string& s1, const string& s2, const string& s3){
		
		assert (cols == 3);
		
		startRow();
		
		push_back_string(s1);
		push_back_string(s2);
		push_back_string(s3);
		
		finishRow();
	}
	
	inline void push_back(const string& s1, const string& s2){
		
		assert (cols == 2);
		
		startRow();
		push_back_string(s1);
		push_back_string(s2);
		
		finishRow();
	}
	
	inline void clear(){
	
		v.clear();
		rows = 0;
	}
	
	class iterator{
	public:
		const int cols;
		const char* begin = NULL;
		const char* ending = NULL;
		
		const char* ptr = NULL;
		
		Row row;
	
		inline iterator() : cols(0){
		
			
		}
		
		inline iterator(const iterator& it) : cols(it.cols), begin(it.begin), ending(it.ending), ptr(it.ptr){
		
			
		}
		
		inline iterator(int cols, const char* begin, const char* ending) : cols(cols), begin(begin), ending(ending){
				
			
			if (begin != NULL && begin != ending)
				ptr = row.set(cols, begin, ending);
			
			//cout << "P" << (ptr-begin) << endl;	
			//cout << "E" << (ending-begin) << endl;	
		}
		
		inline bool operator()() const{
		
			return true;
		}
		
		inline bool operator!=(const bool& b) const{
		
			return row.size() > 0;
		}
		
		inline const Row& operator*() const{
		
			return row;
		}
		
		inline iterator& operator++(){
		
			//cout << "A" << (ptr-begin) << endl;
			//cout << "E" << (ending-begin) << endl;	
			
			if (ptr != ending)
				ptr = row.set(cols, ptr, ending);
			else
				row = Row();
				
			//cout << "B" << (ptr-begin) << endl;
			
			return (*this);
		}
	};
	
	
	
	inline void sort(int sortKey =-1){
	
		vector<Row> rows;
		
		rows.reserve(size() );
		
		for (auto it = begin(); it != end(); ++it){
			rows.push_back( (*it) );
			rows.back().setSortKey(sortKey);
		}
			
		std::sort(rows.begin(), rows.end() );
		
		Rows r2(cols);
		
		r2.reserve(v.size() );
		
		for (auto it = rows.begin(); it != rows.end(); ++it){
			r2.push_back( (*it) );
		}
		
		v = r2.v;
	
	}
	
	inline iterator begin() const{
	
		return iterator(cols, v.data(), v.data()+v.size() );
	}
	
	
	inline iterator beginLast() const{
	
		return iterator(cols, back(), v.data()+v.size() );
	}

	const iterator end;
	
	/*inline void push_back(const CompactTuple& t){
	
		for (int i = 0; i < t.c
	}*/
};

namespace UtilTuple{

	vector<int> without(const vector<int>& v, int x){

		vector<int> ret;
	
		for (int j = 0; j < v.size(); ++j){
	
			if (v[j] != x)
				ret.push_back(v[j]);
		}
		
		return ret;
	}
	
	inline Tuple tuple(const string& s, string split = "_"){
	
		const int n = UtilString::getSplitLength(s, split);
		
		Tuple t;
		
		for (int i = 0; i < n; ++i)
			t.push_back( UtilString::getSplit(s, split, i) );
		
		return t;
	}
	
	inline Tuple tuple(const vector<string>& svec){
	
		Tuple ret;
		
		for (auto it = svec.begin(); it != svec.end(); ++it)
			ret.push_back(*it);
			
		return ret;
	}
	
	inline vector<string> stringvector(const Tuple& t){
	
		vector<string> ret;
		
		for (int j = 0; j < t.size(); ++j)
			ret.push_back(t[j]);
			
		return ret;
	}
	
	inline string getString(const Tuple& t, char sep='_'){
	
		stringstream ss;
	
		for (int j = 0; j < t.size(); ++j){
			if (j > 0)
				ss << sep;
			
			ss << t[j];
		}
		
		return ss.str();
	}


	inline Tuple nonUnique(vector<Tuple> tuples){
	
	
		map<string, int> counts;
		
		for (int j = 0; j < tuples.size(); ++j){
		
			const Tuple& t = tuples.at(j);
		
			for (int k = 0; k < t.size(); ++k)
				counts[ t[k] ]++;
		}
		
		Tuple ret;
		
		
		for (int l = tuples.size(); l >= 2; --l){
			
			for (int j = 0; j < tuples.size(); ++j){
		
				const Tuple& t = tuples.at(j);
		
				for (int k = 0; k < t.size(); ++k){
		
					if (counts[ t[k] ] == l){
				
						ret.push_back(t[k]);
						counts[t[k]] = 0;
					}
				}
			}
		}
		
		
		
		return ret;
		
	}

	inline Tuple select(const Tuple& t, const vector<int>& v){
	
		Tuple ret;
		
		for (int j = 0; j < v.size(); ++j){
		
			if (v.at(j) == -1)
				continue;
		
			assert (v.at(j) >= 0);
			assert (v.at(j) < t.size() );
		
			ret.push_back(t[v.at(j)]);
		}
		
		return ret;
	}
	
	
	inline const string toString(const Tuple& t){
	
		stringstream s;
		
		s << t;
		
		return s.str();
	}



	inline int find(const Tuple& t, const string& s){
		
		for (int i = 0; i < t.size(); i++)
			if (UtilString::equals(s, t[i]))
				return i;
				
		return -1;
	}
	
	
	
	inline int find(const vector<string>& t, const string& s){
		
		for (int i = 0; i < t.size(); i++)
			if (UtilString::equals(s, t[i]))
				return i;
				
		return -1;
	}
	
	
	inline vector<int> find(const Tuple& t, const Tuple& s){
		
		vector<int> ret;
		
		for (int i = 0; i < s.size(); i++){
		
			int f = find(t, s[i]);
			if (f != -1)
				ret.push_back(f);
		}
		
		return ret;
	}
	
	inline Tuple intersection(const Tuple& t1, const Tuple& t2){
	
		Tuple ret;
		
		for (int j = 0; j < t1.size(); ++j){
		
			if (find(t2, t1[j]) != -1)
				ret.push_back(t1[j]);
		}
		
		return ret;
	}
	
	
	inline Tuple minus(const Tuple& t1, const Tuple& t2){
	
		Tuple ret;
		
		for (int j = 0; j < t1.size(); ++j){
		
			if (find(t2, t1[j]) == -1)
				ret.push_back(t1[j]);
		}
		
		//assert (intersection(ret, t2).size() == 0);
		
		return ret;
	}
	
	
	inline Tuple merge(const Tuple& t1, const Tuple& t2){
	
		Tuple ret = t1;
		
		for (int j = 0; j < t2.size(); ++j){
		
			if (find(t1, t2[j]) == -1)
				ret.push_back(t2[j]);
		}
		
		return ret;
	}
	
	inline vector<string> merge(const vector<string>& t1, const vector<string>& t2){
	
		vector<string> ret = t1;
		
		for (int j = 0; j < t2.size(); ++j){
		
			if (find(t1, t2[j]) == -1)
				ret.push_back(t2[j]);
		}
		
		return ret;
	}
	
	
	inline bool containsAll(const Tuple& t1, const Tuple& t2){
	
		return intersection(t1, t2).size() == t2.size();	
	}
	
	inline const string toString(const Tuple& t, const Tuple& a){
	
		hash_set<string> norepeat;
		
		stringstream s;
	
		for (int i = 0; i < a.size(); i++){
		
			if (norepeat.count(a[i]) == 0){

				s << a[i] << "=" << t[i] << " ";
				
				norepeat.insert(a[i]);
			}
		}
		
		return s.str();
	
	}
	
	inline bool equals(const Tuple& a, const Tuple& b){
	
		if (a.size() != b.size() )
			return false;
	
		for (int i = 0; i < a.size(); ++i)
			if (!UtilString::equals(a[i], b[i]))
				return false;
		
		return true;
	}
	
	inline bool equals(const Tuple& a, const vector<int>& v1, const Tuple& b,  const vector<int>& v2 ){
	
		if (v1.size() != v2.size() )
			return false;
		
		for (int i = 0; i < v1.size(); ++i){
		
			if (v1.at(i) == -1 || v2.at(i) == -1)
				continue;
			
			if (!UtilString::equals(a[v1.at(i)], b[v2.at(i)]))
				return false;
		}
		
		return true;
	}
	
	inline shared_ptr< std::unordered_map<string, vector<int> >> getLookupTable(const Tuple& p){
	
		shared_ptr< std::unordered_map<string, vector<int> >> map = std::make_shared< std::unordered_map<string,vector<int>> >(p.size() );
		
		for (int i = 0; i < p.size(); i++){

			const string& s = p[i];

			if (map->count(s) == 0)
				map->insert( {s, vector<int>()} );
			
			(*map)[s].push_back(i);
		}
		
		return map;
	}

	inline vector<vector<int>> getValidityTest(const Tuple& p){
		
		vector<vector<int>> ret;
		
		if (true){
		
			auto map = getLookupTable(p);
			
			for (auto it = map->begin(); it != map->end(); ++it)
				ret.push_back(it->second);
		
			return ret;
		}
	
		for (int i = 0; i < p.size(); i++){
		
			vector<int> inds;
		
			for (int j = 0; j < i; j++){
		
				if (UtilString::equals(p[j], p[i]))
					inds.push_back(j);
			}
			
			ret.push_back(inds);
		}
		
		return ret;
	}

	inline bool isValid(const Tuple& t, const vector<vector<int>>& a){
	
	
		int i = 0;
		for (auto it = a.cbegin(); it != a.cend(); ++it){
		
			const string s = t[i];
		
			for (auto it2 = it->cbegin(); it2 != it->cend(); it2++){
		
				if (!UtilString::equals(s, t[(*it2)] ) )
					return false;
			}	
			i++;
		}
		return true;
	}
	
	inline bool isValid(const Row& t, const vector<vector<int>>& a){
		
		int i = 0;
		for (auto it = a.cbegin(); it != a.cend(); ++it){
		
			for (auto it2 = it->cbegin(); it2 != it->cend(); it2++){
		
				const int j = (*it2);
		
				if (!t.equals(i, j ) )
					return false;
			}
			i++;
		}
		return true;
	}

	inline bool isValid(const Tuple& t, const Tuple& a){
	
	
		for (int i = 0; i < t.size(); i++){
		
			for (int j = 0; j < t.size(); j++){
		
				if (UtilString::equals(a[i], a[j]) && !UtilString::equals(t[i], t[j] ))
					return false;
			}
		}
		
		return true;
	}

	inline const string combine(const Tuple& t, vector<int>& inds){
		
		if (inds.size() == 1)
			return t[inds[0]];
		
		stringstream s;
		
		for (int i = 0; i < inds.size(); i++)
			s << "/-/-/" << t[inds[i]];
			
		return s.str();
	}
	
	inline const string combine(const Row& r, vector<int>& inds){
		
		if (inds.size() == 1)
			return r[inds[0]];
		
		stringstream s;
		
		for (int i = 0; i < inds.size(); i++)
			s << "/-/-/" << r[inds[i]];
			
		return s.str();
	}

	inline vector<int> translate(const Tuple& from, const Tuple& to){
	
		vector<int> ret(from.size() );
		
		for (int i = 0; i < from.size(); ++i){
		
			const string& s = from[i];
			int v = -1;
		
			for (int j = 0; j < to.size(); ++j){
			
				const string& s2 = to[j];
			
				if (UtilString::equals(s, s2)){
				
					v = j;
					break;
				}
			}
			
			ret[i] = v;
		}
		
		return ret;
	}
	
	
	inline vector<int> getMap(const Tuple& from, const Tuple& to){
	
		vector<int> ret;
		
		unordered_map<string, long> map;
		
		for (long j = 0; j < to.size(); j++){
		
			if (map.count(to[j]) == 0){
			
				string s = to[j];
				map.insert( {s, j} );	
			}
		}
		
		for (long j = 0; j < from.size(); j++)
			if (map.count(from[j]) > 0)
				ret.push_back(map[from[j]]);
		
		return ret;
		
	}

	inline vector<int> in(const Tuple& from, const Tuple& to, bool excl = false, bool skipDuplicatesFrom= true, bool skipDuplicatesTo = true){
		
		vector<int> ret;
		
		
		
		for (int i = 0; i < to.size(); ++i){
			
			bool skip = false;
			
			const string& s = to[i];
		
			if (skipDuplicatesTo){
			
				for (int i2 = 0; i2 < ret.size(); ++i2){
				
					const string& s2 = to[ret[i2]];
				
					if (UtilString::equals(s, s2) ){
					
						skip = true;
						break;
					}
				}
			}
			
			if (skip)
				continue;
				
			bool found = false;
			
			for (int j = 0; j < from.size(); ++j){
			
				const string& s2 = from[j];
			
				if (UtilString::equals(s,s2)){
				
					found = true;
				
					if (excl)
						break;
					
					ret.push_back(i);
					
					if (skipDuplicatesFrom)
						break;
				}
			}
			
			if (excl && (!found))
				ret.push_back(i);
		}
		
		//for (int i = 0; i < ret.size(); ++i)
		//	cout << ret[i] << " " << to[ret[i]] << endl;
		
		return ret;
	}
	
};



namespace UtilSortedVector{

	
	template<typename T>
	inline bool isStrictlyIncreasing(const vector<T>& v){
		
		if (v.size() == 0)
			return true;
		
		bool first = true;
		T t0 = v.at(0);
		
		for (auto it = v.begin(); it != v.end(); ++it){
			
			const T& t = (*it);
		
			if (first)
				t0 = t;
			else if (t <= t0)
				return false;
				
			t0 = t;
			first = false;
		}
		
		return true;
	}
	

	template<typename ID, typename T>
	inline void remove(const vector<ID>& sortedIds, const vector<ID>& sortedIdsToRemove, vector<T>& v){

		if (v.size() == 0)
			return;

		assert (isStrictlyIncreasing(sortedIds));
		assert (isStrictlyIncreasing(sortedIdsToRemove));

		if (sortedIds.size() == 0)
			return;
			
		if (sortedIdsToRemove.size() == 0)
			return;

		const long newSize = sortedIds.size()-sortedIdsToRemove.size();
		
		if (newSize == 0){
		
			v.clear();
			return;
		}

		long readPos = 0;
		long writePos = 0;
		
		assert (v.size() == sortedIds.size() );
		
		for (long j = 0; j <= sortedIdsToRemove.size(); j++){
		
			if (j > 0 && j < sortedIdsToRemove.size())
				assert (sortedIdsToRemove.at(j) >= sortedIdsToRemove.at(j-1) );
		
			const long id = j < sortedIdsToRemove.size() ? sortedIdsToRemove.at(j) : 1L << 60;
		
			// skip the remaining ids
			while (readPos < v.size() && sortedIds.at(readPos) < id){
				
				v.at(writePos) = v.at(readPos);
				
				writePos++;
				readPos++;
			}
			
			
			readPos++;
		}
		
		assert (writePos == newSize);
				
		v.resize(newSize);
	}
	
	template<typename ID, typename T>
	inline void removeExcept(const vector<ID>& sortedIds, const vector<ID>& sortedIdsToKeep, vector<T>& v){

		if (sortedIds.size() == 0)
			return;

		long readPos = 0;
		long writePos = 0;
		
		for (long j = 0; j < sortedIdsToKeep.size(); j++){
		
			while (sortedIds[readPos] < sortedIdsToKeep[j])
				readPos++;
			
			assert (sortedIds[readPos] == sortedIdsToKeep[j]);
			
			v[writePos] = v[readPos];
				
			writePos++;
			readPos++;
		}
		
		v.resize(writePos);
	}

};




class WeightingFunction{

public:

	shared_ptr<SQLParser> p;
	vector<long_double> varValues;
	vector<long_double> constValues;
	
	
	
private: 

	long_double MAGIC_VAL = 123.456;

	bool deactivated = true;

	inline void init(){
	
		constValues.reserve(p->getConstants().size());
		
		for (auto it = p->getConstants().begin(); it != p->getConstants().end(); ++it)
			constValues.push_back( UtilHash::to_longdouble(*it) );
		
		varValues.reserve(p->getVariables().size());
		
		for (auto it = p->getVariables().begin(); it != p->getVariables().end(); ++it)
			varValues.push_back( MAGIC_VAL );
	
		const long_double test = p->eval(constValues, varValues);
	
		if (test == UtilExpression::NOT_A_NUMBER){
		
			cout << p->repr() << endl;
		
			assert (false);
		}
		
		for (long j = 0; j < varValues.size(); ++j)
			varValues.at(j) = UtilExpression::NOT_A_NUMBER;
		
		deactivated = test == 1;
		
	}
	
	
	
public:

	inline void disable(bool dis=true){
	
		deactivated = dis;
	}
	
	inline bool isDisabled() const{
	
		return deactivated;
	}

	bool verbose = false;
	
	inline WeightingFunction(){
	
	}
	

	inline void print(){
	
		cout << "DEACTIVATED " << deactivated << endl;
		
		for (auto it = varIndex.begin(); it != varIndex.end(); ++it){
		
			cout << it->first << " " << it->second << endl;
		}
	}
	/*
	inline WeightingFunction(const string& s){
	
		p = shared_ptr<SQLParser>();
		
		p->addBasicOps();
		
		p->parse(s);
		init();
		
	}*/
	
	unordered_map<string, long> varIndex;
	
	
	inline WeightingFunction(shared_ptr<SQLParser> par, const map<string, vector<string>>& aliases){
	
		p = par;
		
		SQLParser& p = (*par);
		
		
		for (auto it = p.expressions.begin(); it != p.expressions.end(); ++it){
		
			const Expression& e = (*it);
		
			const ExpressionOperator& op = p.ops.at(e.op);
			
			if (UtilString::equals(op.name, "rand")){
			
				
				const Expression* order = p.getAncestor(e,"order");
				
				if (order){
				
					// Order By Rand() => weight 1
					if (e.parent == order->id){
					
						Expression& eorder = p.expressions.at(e.parent);
						eorder.mathOp = UtilExpression::ONE;
						eorder.nonMath = false;
						
						//cout << "Order By Rand() => weight 1" << endl;
						
						//assert (false);
						
					} else {
					
						const Expression* log = p.getAncestor(e,"log");
						const Expression* div = p.getAncestor(e,"div");
						
						const Expression* pow = p.getAncestor(e,"pow");
						
						// Order By (a*b*c)/log(Rand() ) => weight a*b*c
						
						if (log != NULL && div != NULL && e.parent == log->id && log->id == div->children.at(1) && div->parent == order->id){
						
							Expression& elog = p.expressions.at(e.parent);
							Expression& ediv = p.expressions.at(elog.parent);
							Expression& eorder = p.expressions.at(ediv.parent);
							
							eorder.mathOp = UtilExpression::IDENTITY0;
							eorder.lenient = true;
							eorder.nonMath = false;
							
							ediv.mathOp = UtilExpression::IDENTITY0;
							ediv.lenient = true;
							
							if ( p.expressions.at( ediv.children.at(0) ).mathOp == UtilExpression::MUL)
								p.expressions.at( ediv.children.at(0) ).lenient = true;
							
							//cout << "Order By (a*b*c)/log(Rand() ) => weight a*b*c" << endl;
							//assert (false);
						}
						
						// Order by Power(Rand(), ...)
						
						if (pow != NULL && e.id == pow->children.at(0) && pow->parent == order->id ){
						
							Expression& epow = p.expressions.at(e.parent);
							Expression& eorder = p.expressions.at(order->id);
							
							eorder.mathOp = UtilExpression::IDENTITY0;
							eorder.lenient = true;
							eorder.nonMath = false;
							
							
							// Order by Power(Rand(), 1 / a*b*c) => weight a*b*c
							if (epow.childrenPtrs.at(1)->mathOp == UtilExpression::DIV){
							
								//cout << "Order by Power(Rand(), 1 / a*b*c) => weight a*b*c" << endl;
							
								Expression& ediv = p.expressions.at(epow.children.at(1));
								
								if (ediv.childrenPtrs.at(0)->mathOp == UtilExpression::MUL)
									p.expressions.at(ediv.children.at(0) ).lenient = true;
								
								if (ediv.childrenPtrs.at(1)->mathOp == UtilExpression::MUL)
									p.expressions.at(ediv.children.at(1) ).lenient = true;
								
								epow.mathOp = UtilExpression::IDENTITY1;
								epow.lenient = true;
								
								ediv.mathOp = UtilExpression::FLIPDIV;
								ediv.lenient = true;
								
								//assert (false);
							
							// Order by Power(Rand(), a*b*c) => weight (1/a)*(1/b)*(1/c)
							} else {
							
								//cout << "Order by Power(Rand(), a*b*c) => weight (1/a)*(1/b)*(1/c)" << endl;
								epow.mathOp = UtilExpression::INVERTED1;
								epow.lenient = true;
								
								//assert (false);
							}
						}
					}
				}
			}
		}
		
		for (auto it = p.expressions.begin(); it != p.expressions.end(); ++it){
		
			Expression& e = (*it);		
			const ExpressionOperator& op = p.ops.at(e.op);
			
			
			if (UtilString::equals(op.name, "var")){
					
				
				const Expression* ptr = &e;
				
				bool unreachable = false;
				
				while (ptr->parent != -1){
				
					unreachable = true;
					
					const Expression* par = &p.expressions.at(ptr->parent);
				
					if ( (par->nonMath))
						break;
				
					if ( (par->mathOp == UtilExpression::IDENTITY0 ||par->mathOp == UtilExpression::INVERTED0) && par->children.at(0) != ptr->id )
						break;
					
					else if ( (par->mathOp == UtilExpression::IDENTITY1 ||par->mathOp == UtilExpression::INVERTED1) && par->children.at(1) != ptr->id )
						break;
					
					else if ( par->mathOp == UtilExpression::NON_MATH)
						break;
						
					else if ( par->mathOp == UtilExpression::ONE)
						break;
						
					ptr = par;
					
					unreachable = false;
				}
				
				if (unreachable)
					e.mathOp = UtilExpression::NON_MATH;
				
				if (e.mathOp == UtilExpression::NON_MATH)
					continue;
				
				if (varIndex.count(e.s) > 0){
				
					e.variableId = varIndex.at(e.s);
					
				} else {
				
					varIndex.insert( {e.s, e.variableId} );
				
					if (aliases.count(e.s) > 0){
				
						const vector<string>& v = aliases.at(e.s);
					
						for (long j = 0; j < v.size(); ++j)
							varIndex.insert( {v.at(j), e.variableId});
					}
				}
			}
		}
		
		init();
	}
	
	
	inline const string toString() const {
		
		if (p)
			return p->toString();
		else
			return "1";
		
	}
	
	
	inline vector<int> getVarIds(const Tuple& t) const {
	
		vector<int> ret;
		
		ret.reserve(t.size() );
		
		for (int j = 0; j < t.size(); ++j){
			if (p)
				ret.push_back( varIndex.count(t[j]) > 0 ? varIndex.at(t[j]) : -1);
			else
				ret.push_back(-1);
		}
		
		/*
		
		vector<long_double> vs;
		
		vs.reserve(varValues.size());
		
		for (long j = 0; j < varValues.size(); ++j)
			vs.push_back(UtilExpression::NOT_A_NUMBER);
		
		
		for (int j = 0; j < ret.size(); ++j)
			if (ret[j] != -1)
				vs.at( ret[j] ) = MAGIC_VAL;
		
		
		const long_double val = p->eval(constValues, vs);
		
		assert (val == p->eval(constValues, vs));
		
		
		for (int j = 0; j < ret.size(); ++j){
			
			if (ret[j] == -1)
				continue;
			
			vs.at( ret[j]) = UtilExpression::NOT_A_NUMBER;
			
			const long_double newval = p->eval(constValues, vs);
			
			if (newval == val){
			
				cout << newval << " == " << val << endl;
				ret[j] = -1;
				
			} else {
			
				cout << newval << " != " << val << endl;
			
				vs.at( ret[j]) = MAGIC_VAL;
				
				assert (p->eval(constValues, vs) == val);
			}
		}*/
		
		return ret;
	}
	
	inline bool usesWeights(const Tuple& t) const {
	
		if (deactivated)
			return false;
			
		vector<int> ids = getVarIds(t);
	
		for (int i = 0; i < t.size(); ++i){
		
			if (ids.at(i) != -1)
				return true;
		}
	
		return false;
	}
	
	
	inline long_double fix(long_double ret){
	
		if (ret < 0)
			return 0;
		
		if ( ((ret-1)+1) != ret)
			return 0;
		
		if ( ((ret-2)+1) == (ret))
			return 0;
			
		if ( ((ret-2)+1) != (ret-1))
			return 0;
		
		return ret;
	} 
	
	template< typename T>
	inline long_double getWeight(const vector<int>& nameIds, const T& r){
	
		if (deactivated)
			return 1;
		
		for (int j = 0; j < nameIds.size(); ++j)
			if (nameIds[j] != -1)
				varValues.at( nameIds[j] ) = UtilHash::to_longdouble(r[j]);
		
		const long_double ret = p->eval(constValues, varValues, verbose);
		
		assert (ret != UtilExpression::NOT_A_NUMBER);
		
		for (int j = 0; j < nameIds.size(); ++j)
			if (nameIds[j] != -1)
				varValues.at( nameIds[j] ) = UtilExpression::NOT_A_NUMBER;
		
		
		
		return fix(ret);
	}
	
	inline long_double getWeight(const vector<int>& nameIds, const Row& r, const vector<int>& sel, bool useSel = true){
	
		if (deactivated)
			return 1;
			
	
		for (int k = 0; k < sel.size(); ++k){
			
			const int j = useSel ? sel[k] : k;
			
			if (nameIds.at(j) != -1)
				varValues.at( nameIds.at(j) ) = UtilHash::to_longdouble(r.get(j), r.getLength(j)-1);
		}
		
		const long_double ret = p ? p->eval(constValues, varValues, verbose) : 1;
		
		assert (ret != UtilExpression::NOT_A_NUMBER);
		
		// RESTORE
		for (int k = 0; k < sel.size(); ++k){
			
			const int j = useSel ? sel[k] : k;
			if (nameIds.at(j) != -1)
				varValues.at( nameIds.at(j) ) = UtilExpression::NOT_A_NUMBER;
		}
		
		
		return fix(ret);
	}
	
	inline long_double getWeight(const vector<int>& nameIds, const Row& r){
	
		return getWeight(nameIds, r, nameIds, false);
	}
	
};

class Weighter{

	WeightingFunction* fptr = NULL;

	const Tuple cols;
	
	const vector<int> nameIds;
	
	shared_ptr<vector<int>> sel;
	
	public:
	
		inline Weighter (WeightingFunction& f, const Tuple& cols) : cols(cols), nameIds(f.getVarIds(cols)){
		
			fptr = &f;
		}
		
		inline Weighter (WeightingFunction& f, const Tuple& cols, const Tuple& selcols) : cols(cols), nameIds(f.getVarIds(cols)){
		
			fptr = &f;
			
			select( UtilTuple::find(cols, selcols) );
		}
		
		inline bool isDisabled() const {
			
			if (nameIds.size() == 0)
				return true;
				
			if (sel)
				if (sel->size() == 0)
					return true;
			
			return fptr->isDisabled();
		}
		
		inline void select(const vector<int>& vsel){
		
			sel = std::make_shared<vector<int>>(vsel.size() );
			
			*sel = vsel;
		}
		
		inline long_double operator()(const Row& r) const{
		
			if (isDisabled() )
				return 1;
				
			if (sel)
				return fptr->getWeight(nameIds, r, *sel, true);
				
			return fptr->getWeight(nameIds, r);
		}
};


namespace JoinColumn{

	const int JOIN = -4;
	const int TARGET_MAIN = -3;
	const int ALL = -2;
	const int NOVEL = -1;
	const int PARENT = 0;
	const int TARGET = 1;
	const int MAIN = 2;
	const int CHILD = 3;
}


class JoinColumns{

typedef long_double Weight;

private:

	bool hasNovel = false;
	bool ready = true;
	
	vector<Tuple> cols;
	vector<vector<vector<int>>> data;
	vector<vector<vector<int>>> data2;
	
	bool allChildrenShareOnlyJoinColumns = true;
	
	vector<vector<int>> sources;
	//vector<vector<int>> sourcesCols;
	
	vector<int> mainVarIds;

	vector<vector<long>> dataHashings;
	
	vector<vector<vector<long>>> dataHashings1;
	vector<vector<vector<long>>> dataHashings2;
	
	shared_ptr< unordered_map<string, long> > hashings;
	
public:
	
	vector<bool> childrenAcyclic;
	vector<int> acyclicChildren;
	vector<int> cyclicChildren;
	
	shared_ptr<WeightingFunction> f;
	
	inline void set(shared_ptr<WeightingFunction> f_){
	
		f = f_;
		
		init();
	}
	
	inline shared_ptr<WeightingFunction> getWeightingFunction() const{
	
		return f;
	}
	
	inline void setHashings(shared_ptr< unordered_map<string, long> > h){
	
		hashings = h;
		
		
		
		if (hasNovel)
			updateHashings();
	}
	
	inline shared_ptr< unordered_map<string, long> > getHashings() const{
	
		return hashings;
	}
	
	vector<long> getHashing(const Tuple& t, const vector<int>& v, bool check=false) const{
	
		vector<long> h;
		
		h.reserve(v.size() );
			
		for (int k = 0; k < v.size(); k++){
			
			if (v.at(k) == -1){
				
				assert (!check);
			
				h.push_back(-1);
				
			} else {
			
				if (hashings){
			
					const string& s = t[v.at(k)];
					h.push_back( hashings->count(s) > 0 ? hashings->at(s) : UtilHash::HASHMAX);
			
				} else {
			
					h.push_back(UtilHash::HASHMAX);
				}
			}
		}	
		
		return h;
	}
	
	vector<long> getHashing(const Tuple& t) const{
	
	
		vector<long> h;
		
		h.reserve(t.size() );
			
		for (int k = 0; k < t.size(); k++){
		
			if (hashings){
			
				const string& s = t[k];
				h.push_back( hashings->count(s) > 0 ? hashings->at(s) : UtilHash::HASHMAX);
			
			} else {
			
				h.push_back(UtilHash::HASHMAX);
			}
		
			
		}
		
		return h;
	}
	
	inline void updateHashings(){
	
		assert (hasNovel);
		
		if (dataHashings.size() != cols.size() ){
			
			dataHashings.clear();
			
			for (int i = 0; i < cols.size(); ++i)
				dataHashings.push_back(vector<long>() );
		}
		
		for (int i = 0; i < cols.size(); ++i)
			dataHashings[i] = getHashing(cols[i]);
		
		for (int i = 0; i < cols.size(); ++i){
			for (int j = 0; j < cols.size(); ++j){
				
				vector<long>& v1 = UtilCollection::tableCell(dataHashings1, i, j);
				v1 = getHashing( getTuple(j), getMap(i,j), true );
				
				vector<long>& v2 = UtilCollection::tableCell(dataHashings2, i, j);
				v2 = getHashing( getTuple(j), getMap2(i,j) );
			}
		}
		
	}
	
	inline const string toString(const Tuple& t) const {
	
		stringstream ss;
		
		for (int i = 0; i < t.size(); ++i)
			ss << " " << t[i];
	
		return ss.str();
	}
	
	inline const string toStringType(int i) const {
	
		
		if (i == (cols.size()+JoinColumn::JOIN) )
			return "JOIN";
		
		if (i == (cols.size()+JoinColumn::TARGET_MAIN) )
			return "TARGET_MAIN";
		
		if (i == (cols.size()+JoinColumn::ALL) )
			return "ALL";
			
		if (i == (cols.size()+JoinColumn::NOVEL) )
			return "NOVEL";
		
		if (i == JoinColumn::PARENT)
			return "PARENT";
			
		if (i == JoinColumn::MAIN)
			return "MAIN";
		
		if (i == JoinColumn::TARGET)
			return "TARGET";
		
		
		if (i >= JoinColumn::CHILD)
			return "CHILD"+to_string(i-JoinColumn::CHILD);
		
		return "?";
	}
	
	inline const string toString(int tabs = 0) const {
	
		stringstream s;
		
		for (int i = 0; i < cols.size(); ++i){
		
			if (i == JoinColumn::TARGET){
			
				for (int ii = 0; ii < tabs; ++ii)
					s << "   ";
			
				s << "JOIN WEIGHTS FOR PARENT" << toString(cols[i]) << endl;
			}
			
			if (i == JoinColumn::MAIN){
			
				for (int ii = 0; ii < tabs; ++ii)
					s << "   ";
			
				s << "TABLE COLUMNS" << toString(cols[i]) << endl;
			}
		
			if (i == JoinColumn::JOIN){
			
				for (int ii = 0; ii < tabs; ++ii)
					s << "   ";
			
				s << "JOIN COLUMNS WITH PARENT" << toString(cols[i]) << endl;
			}
		
			//if (i == cols.size()+ALL ||i == cols.size()+NOVEL)
			//	continue;
		
			//s << i << " " << toString(i) << toString(cols[i]) << endl;
		}
		
		return s.str();
	
	}
	
	inline void ensureOpen(){
		
		if (hasNovel){
		
			for (int i = 1; i <= (-JoinColumn::JOIN); ++i) 
				cols.pop_back();
				
			hasNovel = false;
		}
	}
	
	inline void removeChildren(){
	
		ensureOpen();
		
		while (JoinColumn::CHILD >= cols.size()-1){
		
			cols.pop_back();
		}
	}
	
	inline void set(int j, const Tuple& t){
		
		ensureOpen();
			
		while (j >= cols.size() )
			cols.push_back( Tuple() );	
		
		cols[j].clear();
			
		for (int i = 0; i < t.size(); ++i)	
			cols[j].push_back(t[i]);
			
		//cout << "SET " << toString(j) << " TO " << t << endl;
		
		init();
	}
	
	inline void setTarget(const Tuple& t){
	
		set(JoinColumn::TARGET, t);
		
	}
	
	inline void setParent(const Tuple& t){
	
		set(JoinColumn::PARENT, t);
	}
	
	inline void setMain(const Tuple& t){
	
		set(JoinColumn::MAIN, t);
	}
	
	inline void addChild(const Tuple& t){
	
		ensureOpen();
		
		if (JoinColumn::CHILD >= cols.size())
			set(JoinColumn::CHILD, t);
		else
			set(cols.size(), t);
		
		
	}
	
	
	
	inline long getPos(int j) const{
	
		assert(hasNovel);
	
		return j < 0 ? cols.size()+j : j;
	}
	
	inline const Tuple& getTuple(int j) const{
		
		return cols.at(getPos(j));
	}
	
	
	inline Tuple& get(int j){
		
		return cols.at(getPos(j));
	}
	
	inline const vector<int>& getMap(int i, int j) const{
	
		return data.at( getPos(i) ).at( getPos(j) );
	}
	
	
	inline const vector<int>& getMap2(int i, int j) const{
	
		return data2.at( getPos(i) ).at( getPos(j) );
	}
	
	
	inline const vector<long>& getHashings(int i) const{
	
		return dataHashings.at( getPos(i) );
	}
	
	inline const vector<long>& getHashings(int i, int j) const{
	
		return dataHashings1.at( getPos(i) ).at( getPos(j) );
	}
	
	
	inline const vector<long>& getHashings2(int i, int j) const{
	
		return dataHashings2.at( getPos(i) ).at( getPos(j) );
	}
	
	
	inline const vector<int>& getTranslationTable(int i, int j) const{
	
		return getMap2(i,j);
	}
	
	inline int translate(int i, int k, int j) const{
	
		return getMap2(i,j)[k];
	}
	
	inline const vector<int>& getParentJoinCols() const{
	
		return getMap(JoinColumn::MAIN, JoinColumn::PARENT);
	}
	
	inline const vector<int>& getChildJoinCols(int i) const{
	
		return getMap(JoinColumn::MAIN, JoinColumn::CHILD+i);
	}
	
	inline const vector<int>& getChildJoinColsInMain(int i) const{
	
		return getMap(JoinColumn::CHILD+i, JoinColumn::MAIN);
	}
	
	inline const vector<int>& getNovelInMain() const{
	
		return getMap(JoinColumn::NOVEL, JoinColumn::MAIN);
	}
	
	inline const vector<int>& getTargetInMain() const{
	
		return getMap(JoinColumn::TARGET, JoinColumn::MAIN);
	}
	
	inline const vector<int>& getJoinInMain() const{
	
		return getMap(JoinColumn::JOIN, JoinColumn::MAIN);
	}
	
	
	
	inline const vector<int>& getTargetInChild(int i) const{
	
		return getMap(JoinColumn::TARGET, JoinColumn::CHILD+i);
	}
	
	inline const Tuple& getChildTuple(int j) const{
		
		return getTuple(JoinColumn::CHILD+j);
	}
	
	inline const Tuple& getMainTuple() const {
	
		return getTuple(JoinColumn::MAIN);
	}
	
	inline const Tuple& getNovelTuple() const {
	
		return getTuple(JoinColumn::NOVEL);
	}
	
	inline const Tuple& getAllTuple() const {
	
		return getTuple(JoinColumn::ALL);
	}
	
	inline const Tuple& getJoinTuple() const {
	
		return getTuple(JoinColumn::JOIN);
	}
	
	
	inline const Tuple& getTargetTuple() const {
	
		return getTuple(JoinColumn::TARGET);
	}
	
	inline const Tuple& getParentTuple() const {
	
		return getTuple(JoinColumn::PARENT);
	}
	
	inline int getChildNum() const {
	
		return cols.size()-JoinColumn::CHILD+JoinColumn::JOIN;
	}
	
	inline const vector<int>& getJoinColsWith(int j) const{
	
		return getMap(JoinColumn::MAIN, j);
	}
	
	inline const vector<int>& getJoinColsWithParent() const{
	
		return getJoinColsWith(JoinColumn::PARENT);
	}
	
	inline const vector<int>& getJoinColsWithChild(int j) const{
	
		return getJoinColsWith(JoinColumn::CHILD+j);
	}
	
	template<typename T>
	inline void print(const vector<T>& v) const{
	
		cout << UtilString::toString(v) << endl;
	}
	
	inline void print(const Tuple& t) const{
	
		cout << t << endl;
	}
	
	inline void print(const Tuple& t, const vector<int>& v) const{
	
		cout << "{ " << t << " }[" << UtilString::toString(v) << "] = { ";
		
		for (int k = 0; k < v.size(); ++k)
			cout << t[v[k]] << " ";
			
		cout << "}" << endl;
	}
	
	
	
	
	inline const vector<int>& getSources(int j, int i) const{
	
		/*
		cout << "getMap(j, ALL).at(i) " << getMap(j, ALL).at(i) << endl;
		
		cout << "src0 " << UtilString::toString( sources.at(0)) << endl;
		cout << "src1 " << UtilString::toString( sources.at(1)) << endl;
		cout << "src2 " << UtilString::toString( sources.at(2)) << endl;*/
	
		return sources.at(getMap(j, JoinColumn::ALL).at(i));
	}
	
	/*
	inline const vector<int>& getSourcesCols(int j, int i) const{
	
		return sourcesCols.at(getMap(j, ALL).at(i));
	}*/
	
	inline bool isChildAcyclic(int i) const {
	
		return childrenAcyclic[i];
	}
	
	inline bool areAllChildrenAcyclic() const {
	
		return cyclicChildren.size() == 0;
	}
	
	inline bool areAllChildrenSimple() const {
	
		return areAllChildrenAcyclic() || allChildrenShareOnlyJoinColumns;
	}
	
	inline Weight getNovelWeightInMain(const Row& tmain) const{
		
		if (f)
			return f->getWeight(mainVarIds, tmain, getNovelInMain() );
		else
			return 1;
	}
	
	inline void init(){
		
		data.clear();
		data2.clear();
		
		if (!hasNovel){
		
			for (int i = 1; i <= (-JoinColumn::JOIN); ++i) 
				cols.push_back(Tuple() );
			
			hasNovel = true;
			
		} else {
		
			get(JoinColumn::NOVEL).clear();
			get(JoinColumn::ALL).clear();
			get(JoinColumn::JOIN).clear();
			get(JoinColumn::TARGET_MAIN).clear();
		}
		
		sources.clear();
		//sourcesCols.clear();
		
		childrenAcyclic.clear();
		acyclicChildren.clear();
		cyclicChildren.clear();
		
		const int novel = cols.size()+JoinColumn::NOVEL;
		//const int all = cols.size()+JoinColumn::ALL;
		const int lastChild = cols.size()+JoinColumn::JOIN-1;
		
		std::map<string, vector<int>> msources;
		//std::map<string, vector<int>> msourcesColumns;
		
		const Tuple& parent = getTuple(JoinColumn::PARENT);
		const Tuple& target = getTuple(JoinColumn::TARGET);
		const Tuple& main = getTuple(JoinColumn::MAIN);
		
		std::set<string> tparentset;
		std::set<string> ttargetset;
		std::set<string> tmainset;
		
		for (int i = 0; i < parent.size(); ++i)
			if (tparentset.count(parent[i]) == 0)
				tparentset.insert(parent[i]);
			
		for (int i = 0; i < target.size(); ++i)
			if (ttargetset.count(target[i]) == 0)
				ttargetset.insert(target[i]);
		
		for (int i = 0; i < main.size(); ++i)
			if (tmainset.count(main[i]) == 0)
				tmainset.insert(main[i]);
		
		Tuple& ttarget_main = get(JoinColumn::TARGET_MAIN);
		
		ttarget_main.clear();
		
		Tuple& tjoin = get(JoinColumn::JOIN);
		
		tjoin.clear();
		
		for (int i = 0; i < target.size(); ++i)
			if (tmainset.count(target[i]) > 0)
				ttarget_main.push_back(target[i]);
				
		
		for (int i = 0; i < main.size(); ++i)
			if (ttargetset.count(main[i])+tparentset.count(main[i]) > 1)
				tjoin.push_back(main[i]);
				
		
		for (int j = 0; j <= lastChild; j++){
			
			const Tuple& tj = getTuple(j);
			
			for (int i = 0; i < tj.size(); ++i){
			
				const string& s = tj[i];
			
				if (msources.count(s) == 0){
					
					msources.insert( {s, vector<int>() } );
					//msourcesColumns.insert( std::make_pair(s, vector<int>() ) );
				}
				
				msources[s].push_back(j);
				//msourcesColumns[s].push_back(i);
			}
		}
		
		if (true){
		
			Tuple& tall = get(JoinColumn::ALL);
			std::set<string> tallset;
			
			for (int j = 0; j <= lastChild; j++){
			
				const Tuple& tj = getTuple(j);
			
				for (int i = 0; i < tj.size(); ++i){
				
					const string& s = tj[i];
				
					if ( tallset.count(s) == 0){
					
						tall.push_back(s);
						tallset.insert(s);
					}
				}	
			}
		
			const Tuple& tmain = getTuple(JoinColumn::MAIN);
			Tuple& tnovel = get(JoinColumn::NOVEL);
			
			std::set<string> tnovelset;
		
			for (int i = 0; i < tmain.size(); ++i){
		
				const string& s = tmain[i];
				
				vector<int>& v = msources[s];
				
				//cout << s << " found in " << UtilString::toString(v) << endl;
			
				bool isNovel = true;
			
				for (int j = 0; j < v.size(); j++){
					
					if (v[j] >= JoinColumn::CHILD)
						isNovel = false;
				}
				
				if (isNovel){
				
					//cout << s << " is novel" << endl;
				
					msources[s].push_back(novel);
					//msourcesColumns[s].push_back(tnovel.size() );
					
					tnovel.push_back(s);
					tnovelset.insert(s);
				}
			}
		}
		
		allChildrenShareOnlyJoinColumns = true;
		
		for (int j = JoinColumn::CHILD; j <= lastChild; j++){
		
			const Tuple& tj = getTuple(j);
			
			bool acyclic = true;
		
			for (int l = 0; l < tj.size(); ++l){
			
				const string& s = tj[l];
			
				vector<int>& v = msources[s];
				
				const bool inMain = std::binary_search(v.begin(), v.end(), JoinColumn::MAIN);
				const bool inTarget = std::binary_search(v.begin(), v.end(), JoinColumn::TARGET);
				
				if (inTarget && !inMain){
				
					acyclic = false;
					
					for (int k = JoinColumn::CHILD; k <= lastChild; k++){
						
						if (k == j)
							continue;
						
						if (std::binary_search(v.begin(), v.end(), k)){
							allChildrenShareOnlyJoinColumns = false;
							break;
						}
					}
					
					break;
				}
			}
			
			childrenAcyclic.push_back(acyclic);
			
			if (acyclic)
				acyclicChildren.push_back(j-JoinColumn::CHILD);
			else
				cyclicChildren.push_back(j-JoinColumn::CHILD);
		}
		
		for (int i = 0; i < cols.size(); ++i){
			
			for (int j = 0; j < cols.size(); ++j){
				
				vector<int>& v = UtilCollection::tableCell(data, i, j);
				
				v = UtilTuple::getMap( getTuple(i), getTuple(j) );
				
				vector<int>& v2 = UtilCollection::tableCell(data2, i, j);
				
				v2 = UtilTuple::translate( getTuple(i), getTuple(j) );
			}
		}
		
		const Tuple& tall = getAllTuple();
		
		for (int j = 0; j < tall.size(); ++j){
		
			const string& s = tall[j];
			
			const vector<int> src = msources.at(s);
			//const vector<int> cols = msourcesColumns.at(s);
			
			sources.push_back(src);
			//sourcesCols.push_back(cols);
			
			/*
			cout << "$$" << s << endl;
			
			for (int i = 0; i < src.size(); ++i){
			
				const Tuple& tuple = getTuple(src[i]);
				
				cout << src[i] << " " << tuple << "[" << src[i] << "] = " << tuple[col[i]] << endl;
			}
			
			cout << "/$$" << s << endl;*/
		}
		
		
		if (f){
		
			mainVarIds = f->getVarIds(getMainTuple() );
		}
		
		dataHashings.clear();
		dataHashings1.clear();
		dataHashings2.clear();
		
		updateHashings();
		
		ready = true;
	}
};





namespace WeightedSamplers{

	const int NONE = 0;
	const int FINAL = 1;
	const int POP = 2;
	const int OFFLINE_WR = 3;
	const int ONLINE_WR = 4;
	const int ONLINE_WOR = 5;
	
	
	
	
};

template<typename T>
class WeightedSample{



typedef long_double Weight;

protected:

	// DYNAMIC MANAGEMENT
	
	vector<long> identities; // identities of elements for mark-removal during online phase
	vector<long> toRemove; // identities of elements marked for removal
	long idfactory = 0; // next id
	unique_ptr<RandomNumbers> rnd;
	long currentSize = 0; // size during online phase
	long randomId = 0;
	
	vector<double> randomKeys;
	
	
	vector<T> batchElements;
	vector<Weight> batchWeights;
	vector<int> batchFrequencies;
	
	long batchSize = 0; // size during online phase
	
	Weight beforeBatch = 0;
	
	long maxBatchSize = 0;
	
	int sampler = WeightedSamplers::NONE;
	
	// SAMPLE DATA
public:	
	vector<T> elements;
	const string seed;
	long maxSize; // MAXIMAL SUPPORTED SAMPLE SIZE
	const bool isPop;

protected:
	
	bool weightsCumulative = false; 
	vector<Weight> weights; // probabilistic weight of element
	vector<int> frequencies; // repetitions of elements
	long sampleSize = 0; // final size;
	bool wr = false;
	bool wor = false;
	
	// POPULATION STATISTICS
	
	Weight populationWeight = 0; // probabilistic weight of population
	long populationSize = 0; // number of elements in population
	
	// POPULATION WR SAMPLER
	
	shared_ptr<vector<Weight>> targets;
	
	// ONLINE WOR SAMPLER
	
	long maxRemovalQueue = 0;
	unique_ptr<std::priority_queue< SampleEntry >> heap;
	Weight wc_wi = -1;
	Weight Xw = -1;
	Weight logTw = 1;
	bool exponentialJumps = true;
	
	
	unique_ptr<std::priority_queue< std::pair<Weight,int> >> wrheap;
	vector<Weight> logTws;
	vector<long> wrids;
	
public:

	inline Weight getXw() const {
	
		return Xw;
	}

	class iterator2{

		typename vector<T>::const_iterator elemIt;
		typename vector<T>::const_iterator elemItEnd;
		vector<int>::const_iterator freqIt;
		vector<int>::const_iterator freqItEnd;
		vector<Weight>::const_iterator weightIt;
		vector<Weight>::const_iterator weightItEnd;
		
		long freqsSize = 0;
		Weight lastWeight = 0;
		
		public:
		inline iterator2(const vector<T>& elements, const vector<int>& freqs, const vector<long_double>& weights) :
	
		freqsSize(freqs.size() ),
	
		elemIt(elements.cbegin() ),
		elemItEnd(elements.cend() ),
	
		freqIt(freqs.cbegin() ),
		freqItEnd(freqs.cend() ),
	
		weightIt(weights.cbegin() ),
		weightItEnd(weights.cend() )
	
	
		{
	
		
		
		}
	
	
		inline bool operator!=(const bool& b) const{
	
			if ( !(elemIt != elemItEnd) )
				return false;
	
			if (freqsSize > 0)
			if ( !(freqIt != freqItEnd) )
				return false;
	
			if ( !(weightIt != weightItEnd) )
				return false;
	
			return true;
		}
	
		inline bool operator==(const bool& b) const{
	
			return !( (*this) != b);
		}
	
		inline bool operator()() const{
	
			return true;
		}

		inline const T& operator*() const {
	
			return (*elemIt);
		}
	
		inline long_double getWeight() const {
	
			return (*weightIt)-lastWeight;
		}
	
		inline long_double getCumulativeWeight() const {
	
			return (*weightIt);
		}
	
		inline const long getFrequency() const {
	
			return freqsSize > 0 ? (*freqIt) : 1;
		}
	
		inline iterator2& operator++(){
	
			if (weightIt != weightItEnd)
				lastWeight = (*weightIt);
		
			++elemIt;
			++weightIt;
			
			if (freqIt != freqItEnd)
				++freqIt;
		
			return (*this);
		}
	};
	
	const iterator2 elementsEnd;
	const iterator2 end2;
	
	inline WeightedSample<T>::iterator2 elementsBegin() const{
	
		assert (sampler == WeightedSamplers::FINAL);
		assert (weightsCumulative);
	
		return WeightedSample<T>::iterator2(elements, this->frequencies, this->weights );
	}

	inline WeightedSample<T>::iterator2 begin2() const{
	
		return elementsBegin();
	}
	
	inline void makeWeightsCumulative(bool b){
	
		if (weightsCumulative == b)
			return;
			
		if (b){
		
			for (long j = 1; j < weights.size(); ++j)
				weights.at(j) += weights.at(j-1);
				
		} else {
		
			for (long j = weights.size()-1; j >= 1; --j)
				weights.at(j) -= weights.at(j-1);
		}
		
		weightsCumulative = b;
	}
	
	
	/*
	inline long getIndependentSample(){
		
		requireRandom();
		
		makeWeightsCumulative(true);
		
		const Weight u = rnd->randomDouble() * weights.at(weights.size()-1);
	
		if (u <= (*weights.begin()))
			return 0;
		
		auto it = std::upper_bound(weights.begin(), weights.end(), u);
		
		return it-weights.begin();
	}*/
	
	/*
	inline WeightedSample(const WeightedSample<T>& smp) : isPop(smp.isPop), maxRemovalQueue(smp.maxRemovalQueue), seed(smp.seed), maxSize(smp.maxSize), end2(WeightedSample<T>::iterator2(this->elements, this->frequencies, this->weights )) {
	
		identities = smp.identities;
		toRemove = smp.toRemove;
		weights = smp.weights;
		frequencies = smp.frequencies;
		
		weightsCumulative = smp.weightsCumulative;
		
		idfactory = smp.idfactory;
		populationWeight = smp.populationWeight;
		populationSize = smp.populationSize;
		
		rnd = smp.rnd;
		
		currentSize = smp.currentSize;
		sampleSize = smp.sampleSize;
		
		elements = smp.elements;
		
		wr = smp.wr;
		wor = smp.wor;
		
		assert (frequencies.size() == identities.size() );
	}*/
	
	inline WeightedSample(const string s, long maximalSize=0, bool pop = false) : isPop(pop),seed(s), maxSize(maximalSize), end2(WeightedSample<T>::iterator2(this->elements, this->frequencies, this->weights )), elementsEnd(WeightedSample<T>::iterator2(this->elements, this->frequencies, this->weights )) {
	
		
		if (pop)
			sampler = WeightedSamplers::POP;
	}
	
	
	inline void initOffline( shared_ptr<vector<Weight>> t){
		
		targets = t;
		std::sort(targets->begin(), targets->end(), std::greater<Weight>());
		
		assert (sampler == WeightedSamplers::NONE);
		sampler = WeightedSamplers::OFFLINE_WR;
		
		
	}
	
	inline void initOffline(Weight popWeight, bool withReplacement=true){	
	
		assert (withReplacement);
		
		requireRandom();
	
		shared_ptr<vector<Weight>> vec = std::make_shared<vector<Weight>>();
				
		vec->reserve(maxSize);		
	
		for (long j = 0; j < maxSize; j++)
			vec->push_back( rnd->randomDouble() * popWeight);
				
		initOffline(vec);
		
		rnd.reset();
		
		reserveInternal();
	}
	
	
	inline void init(long size){
	
		maxBatchSize = 0;
		maxRemovalQueue = 0;
		
		cout << "WeightedSample init " << UtilString::valuesToString<int, long, Weight, long>(size, size, size, size) << endl;
		
		
		frequencies.reserve(size);
		identities.reserve(size);
		weights.reserve(size);
		toRemove.reserve(size);
			
		reserveInternal(size);
		
		
		
	}
	
	inline void initPop(){
	
		maxBatchSize = 0;
		maxRemovalQueue = 0;
		
		
			
		
		
		cout << "WeightedSample initPop " << UtilString::valuesToString<int, long, Weight, long>(maxSize, maxSize, maxSize, maxSize) << endl;
		
		frequencies.reserve(maxSize);
		identities.reserve(maxSize);
		weights.reserve(maxSize);
		toRemove.reserve(maxSize);
		
		reserveInternal();
	}
	
	inline void initOnline(bool withReplacement=true, long maxExtra = -1){
	
		if (maxExtra == -1)
			maxExtra = maxSize >> 4;
		
		if (withReplacement){
		
			maxRemovalQueue = maxExtra;
		
			assert (sampler == WeightedSamplers::NONE);
			sampler = WeightedSamplers::ONLINE_WR;
			
			
			cout << "WeightedSample initOnlineWOR heap " << UtilString::valuesToString<std::pair<Weight,int>>(maxSize) << endl;
			
			std::vector<std::pair<Weight,int>> vec;
			vec.reserve(maxSize);
			
			wrheap.reset(new std::priority_queue< std::pair<Weight,int> >(std::less<std::pair<Weight,int>>(), vec) );
			
			const long r = maxSize+maxRemovalQueue;
			
			cout << "WeightedSample initOnlineWR " << UtilString::valuesToString<Weight, long, Weight, long, long>(maxSize, maxSize, r, r, r) << endl;
			
			logTws.reserve(maxSize);
			wrids.reserve(maxSize);
			identities.reserve(maxSize+maxRemovalQueue);
			weights.reserve(maxSize+maxRemovalQueue);
			toRemove.reserve(maxRemovalQueue);
			
			reserveInternal();
			
		} else {
			
			maxRemovalQueue = maxExtra;
		
			assert (sampler == WeightedSamplers::NONE);
			sampler = WeightedSamplers::ONLINE_WOR;
			
			
			cout << "WeightedSample initOnlineWOR heap " << UtilString::valuesToString<SampleEntry>(maxSize) << endl;
			
			
			std::vector<SampleEntry> vec;
			vec.reserve(maxSize);
	
			heap.reset(new std::priority_queue< SampleEntry >(std::less<SampleEntry>(), vec) );
			
			
			
			const long r = maxSize+maxRemovalQueue;
			
			cout << "WeightedSample initOnlineWOR " << UtilString::valuesToString<double, Weight, long, long>(r, r, r, r) << endl;
			
			randomKeys.reserve(maxSize+maxRemovalQueue);
			identities.reserve(maxSize+maxRemovalQueue);
			weights.reserve(maxSize+maxRemovalQueue);
			toRemove.reserve(maxRemovalQueue);
			
			reserveInternal();
		}
		
		
	}
	
	inline bool initOnlineAuto(int sampler, long maxExtra = -1){
		
		cout << "SAMPLER " << sampler << " " << maxSize << endl;
		
		if (sampler == WeightedSamplers::ONLINE_WR){
		
			initOnlineWR(maxExtra);
			
		} else if (sampler == WeightedSamplers::ONLINE_WOR){
		
			initOnlineWOR(maxExtra);
			
		} else {
		
			return false;
		}
		
		return true;
	
	}
	
	inline void initOnlineWR(long maxExtra = -1){
	
		initOnline(true, maxExtra);
	}
	
	inline void initOnlineWOR(long maxExtra = -1){
	
		initOnline(false, maxExtra);
	}
	
	inline void initOfflineWR(long maxExtra = -1){
	
		initOffline(true, maxExtra);
	}
	
	inline void initOfflineWOR(long maxExtra = -1){
	
		initOffline(false, maxExtra);
	}
	
	inline virtual ~WeightedSample(){
	
		
	}
	
	inline long getSizeDistinct() const{
	
		return currentSize;
	}
	
	inline long getSize() const{
	
		return sampleSize;
	}
	
	inline long size() const{
	
		return sampleSize;
	}
	
	inline vector<int>& getFrequencies(){
	
		return frequencies;
	}
	
	inline const vector<Weight>& getSampleWeights() const{
	
		assert (!weightsCumulative);
	
		return weights;
	}
	
	inline const vector<long_double>& getCumulativeSampleWeights() const{
	
		assert (weightsCumulative);
	
		return weights;
	}
	
	inline Weight getPopulationWeight() const{
	
		return populationWeight;
	}
	
	inline void setPopulationWeight(Weight w){
	
		populationWeight = w;
	}
	
	
	
	inline long getPopulationSize() const{
	
		return populationSize;
	}
	
	inline void setPopulationSize(long s){
	
		populationSize = s;
	}
	
	
	inline const string getSeed() const{
	
		return seed;
	}
	
	
	inline void addPopulation(Weight w){
	
		populationWeight += w;
		populationSize++;
	}
	
protected:


	inline virtual void clearToRemove(){
	
		UtilSortedVector::remove(identities, toRemove, elements);
		UtilSortedVector::remove(identities, toRemove, weights);
		
		UtilSortedVector::remove(identities, toRemove, frequencies);
		
		UtilSortedVector::remove(identities, toRemove, randomKeys);
			
		UtilSortedVector::remove(identities, toRemove, identities);
	}

	inline virtual void addInternal(const T& elem){
	
		if (elements.capacity() < identities.capacity())
			elements.reserve(identities.capacity());
		
		elements.push_back(elem);
	}

	inline virtual void reserveInternal(long size =-1){
	
		elements.reserve(size == -1 ? maxSize : size);
	}

	inline virtual void addToBatchInternal(const T& elem){
		
		batchElements.push_back(elem);
	}
	
	inline virtual void compress(){
	
		
	}
	
	inline void requireRandom(){
	
		if (rnd){
		
			
		} else {
		
			rnd.reset(new RandomNumbers(seed+"_"+to_string(++randomId) ) );
		}
	}
	
	inline virtual void shrinkInternal(){
	
		elements.shrink_to_fit();
	}
	
	inline void shrink(){
	
		randomKeys.shrink_to_fit();
		
		frequencies.shrink_to_fit();
		identities.shrink_to_fit();
		toRemove.shrink_to_fit();
		weights.shrink_to_fit();
		
		batchElements.shrink_to_fit();
		batchWeights.shrink_to_fit();
		batchFrequencies.shrink_to_fit();
		
		logTws.shrink_to_fit();
		
	}
	
	

public:

	long queueops = 0;

	inline void finalise(bool withReplacement=false, long m = -1){
	
		assert (sampler != WeightedSamplers::FINAL);
	
		rnd.reset();
		
		clearRemovals(true);
		
			
		if (sampler == WeightedSamplers::ONLINE_WR ||sampler == WeightedSamplers::OFFLINE_WR)
			assert (withReplacement);
		
		
		const long collectedSamps = frequencies.size() == 0 ? weights.size() : UtilCollection::sum(frequencies);
		
		if (!withReplacement){
		
			if (m == -1 || m > collectedSamps)
				m = collectedSamps;
				
		} else {
		
			if (m == -1)
				m = maxSize;
		}
		
		if (m > 0)
		for (int cc = 0; cc < 1; ++cc){
		
			
			if (sampler == WeightedSamplers::POP || sampler == WeightedSamplers::NONE || sampler == WeightedSamplers::OFFLINE_WR){
				
				assert (m == collectedSamps);
				
				if ( sampler == WeightedSamplers::POP)
					assert (UtilCollection::sum(weights) >= populationWeight);
			
				assert (weightsCumulative == false);
				
				break;
			}
			
			if (UtilCollection::sum(weights) >= populationWeight){
			
				frequencies.clear();
				
				cout << "WeightedSampler reserve " << UtilString::valuesToString<int>(weights.size() ) << endl;
				
				frequencies.reserve(weights.size());
	
				for (long j = 0; j < weights.size(); j++)
					frequencies.push_back(1);
			
				requireRandom();
			
				if (withReplacement)
					rnd->getWRfromPOP(weights, frequencies, m);
				else
					rnd->getWORfromPOP(weights, frequencies, m);
			
				break;
			}
			
			if (sampler == WeightedSamplers::ONLINE_WR){
		
				assert (withReplacement);
			
				wrheap.reset();
				rnd.reset();
				
				break;
			}
			
			if (sampler == WeightedSamplers::ONLINE_WOR){
		
				heap.reset();
	
				//if (!isPop)
				//	assert( m <= weights.size() );
			
				assert (!weightsCumulative);
			
				cout << "finalise reserve " << UtilString::valuesToString<int>(weights.size() ) << endl;
			
				frequencies.clear();
				frequencies.reserve(weights.size());
				
				for (long j = 0; j < weights.size(); j++)
					frequencies.push_back(1);
				
				requireRandom();
			
				
				if (withReplacement){
				
					vector<SampleEntry> v;
				
					v.reserve(weights.size() );
				
					if (weights.size() < m)
						m = weights.size();
			
					assert (weights.size() == randomKeys.size() );
			
					for (long j = 0; j < randomKeys.size(); j++)
						v.push_back( SampleEntry(-randomKeys.at(j), j) );
			
					std::sort(v.begin(), v.end() );
			
					vector<long> ordered;
			
					ordered.reserve(v.size() );
			
					for (long j = 0; j < randomKeys.size(); j++)
						ordered.push_back( v.at(j).id );
		
					rnd->getWRfromWOR(populationWeight, weights, frequencies, ordered, m);
				
				} else {
					
					rnd->getWORfromWOR(populationWeight, weights, frequencies, m);
				}
				
				rnd.reset();
			
				assert (weights.size() == identities.size() );
				assert (frequencies.size() == identities.size() );
			
				for (long j = 0; j < weights.size(); j++)
					weights[j] *= frequencies[j];
			
				if (UtilCollection::sum(frequencies) != m){
				
					cout << UtilCollection::sum(frequencies) << " != " << m << endl;
					assert (UtilCollection::sum(frequencies) == m);
				}
				
				break;
			}
		}
		
		//assert (populationWeight > 0);
		
		assert (toRemove.size() == 0);
		
		assert (!weightsCumulative);
		
		if (frequencies.size() > 0)
		for (long j = 0; j < identities.size(); ++j)
			if (frequencies.at(j) == 0)
				remove(identities.at(j) );
		
		clearRemovals(true);
		
		identities.clear();
		
		if (sampler == WeightedSamplers::ONLINE_WOR){
			
			shrink();
			shrinkInternal();
		}
		
		if (frequencies.size() > 0)
			sampleSize = UtilCollection::sum(frequencies);
		else
			sampleSize = maxSize-getInvalidSize();
		
		//compress();
		
		makeWeightsCumulative(true);
		
		wr = withReplacement;
		wor = !withReplacement;
		
		sampler = WeightedSamplers::FINAL;
	}

	inline void shuffle(){
	
		assert (sampler == WeightedSamplers::FINAL);
		
		vector<long> posVec(weights.size());
		
		for (long i = 0; i < posVec.size(); ++i)
			posVec[i] = i;
		
		requireRandom();
		
		rnd->randomPermutation(posVec);
		
		identities.clear();
		identities.reserve(weights.size() );
		
		for (long i = 0; i < identities.size(); ++i)
			identities.push_back(i);
		
		makeWeightsCumulative(false);
		
		UtilCollection::permute(posVec, weights);
		
		makeWeightsCumulative(true);
		
		UtilCollection::permute(posVec, frequencies);
		
		internalPermute(posVec);
	
		rnd.reset();
	}
	
	inline void processBatch(){
	
		if (sampler != WeightedSamplers::ONLINE_WR)
			return;
	
		if (batchSize == 0)
			return;
	
		
		//cout << "batch.start" << endl;
		
		requireRandom();
		
		vector<Weight> newCoords; 
		{
			const Weight batchWeight = populationWeight-beforeBatch;
			const Weight prop1 = batchWeight / populationWeight;
			const Weight prop2 = prop1 * 2 < 1 ? prop1 *2 : prop1;
		
			newCoords.reserve( prop2 * batchSize);
		}
		
		long keep = 0;
		
		//cout << "batch.randomvariates" << endl;
		
		for (long j = 0; j < maxSize; ++j){
		
			const Weight u = rnd->randomDouble() * populationWeight;
			
			if (u < beforeBatch){
				
				keep++;
				
			} else {
				
				if (newCoords.size() == newCoords.capacity() )
					newCoords.reserve(newCoords.size()*2);
				
				newCoords.push_back(u-beforeBatch);
			}
		}
		
		//cout << "batch.simprand" << endl;
		
		rnd->getWRfromWR(frequencies, keep);
		
		//cout << "batch.clear" << endl;
		
		assert (UtilCollection::sum(frequencies) == keep);
		
		for (long j = 0; j < frequencies.size(); ++j)
			if (frequencies.at(j) == 0)
				remove(identities.at(j) );
		
		clearRemovals();
		
		assert (UtilCollection::sum(frequencies) == keep);
		
		//cout << "batch.newcoords" << endl;
		
		std::sort(newCoords.begin(), newCoords.end() );
		
		assert (UtilCollection::sum(batchFrequencies) == batchSize );
		
		rnd->setCoords(batchWeights, batchFrequencies, newCoords);
		
		assert (UtilCollection::sum(batchFrequencies) == newCoords.size() );
		
		//cout << "batch.rows" << endl;
		
		addBatchToSampleInternal();
		
		assert (UtilCollection::sum(frequencies) == maxSize);
		
		batchWeights.clear();
		batchFrequencies.clear();
		batchSize = 0;
		beforeBatch = populationWeight;
		
		//cout << "batch.finished" << endl;
	
	}
	
	inline long addToSample( const T& elem, Weight w=1, long freq = -1, double randomKey = -1){
	
		if (freq == 0)
			return -1;
		
		++idfactory;
		
		
		
	
		const long identity = idfactory;
		
		//cout << "addelem "<<elem << " w " << w << " id " << identity << endl;
		
		identities.push_back(identity);
		addInternal(elem);
		weights.push_back(w);
		
		if (freq != -1){
			frequencies.push_back(freq);
			currentSize += freq;
		} else {
			
			assert (frequencies.size() == 0);
			++currentSize;
		}
		
		if (randomKey != -1){
			randomKeys.push_back(randomKey);
		} else {
			
			assert (randomKeys.size() == 0);
		}
		
		return identity;
	
	}
	
	inline virtual void addBatchToSampleInternal(){
		
		long j = 0;
		for (auto it = batchElements.begin(); it != batchElements.end(); ++it){
		
			addToSample( (*it), batchWeights.at(j), batchFrequencies.at(j));
			++j;
		}
		
		batchElements.clear();
	}
	
	inline void remove(long id){
		
		toRemove.push_back(id);
	
		clearRemovals(false);
		
		--currentSize;
		
		
	}
	
	inline void clearRemovals(bool force = false){
	
		if (force || (toRemove.size() > maxRemovalQueue) ){
		
			std::sort(this->toRemove.begin(), this->toRemove.end() );
			clearToRemove();
			toRemove.clear();
		
		}
	}
	
	inline void add(const T& elem, Weight w=1){
	
		assert (sampler != WeightedSamplers::FINAL);
		
		assert (maxSize > 0);
		
		addPopulation(w);
		
		if (sampler == WeightedSamplers::NONE ||sampler == WeightedSamplers::POP){
			
			addToSample(elem, w);
			return;
		}
		
		if (sampler == WeightedSamplers::ONLINE_WR){
		
			/*addToBatchInternal(elem);
			batchWeights.push_back(w);
			batchFrequencies.push_back(1);
			++batchSize;
			
			if (batchSize > maxBatchSize)
				processBatch();*/
				
			const Weight wi = w;
			
			if (wi == 0)
				return;
			
			assert (wi > 0);
			
			requireRandom();
		
			if (wc_wi == -1)
				wc_wi = 0;
		
			wc_wi += wi;
		
			if ( logTws.size() == 0 ){
			
				for (long k = 0; k < maxSize; ++k){
				
					
					wrids.push_back(addToSample(elem, wi));
					
					const Weight logki = log(rnd->randomDouble() )/wi;
					
					const Weight logTw = logki;
					
					logTws.push_back( logTw );
					
					const Weight logr = log(rnd->randomDouble());
					
					const Weight Xw1 = logr/logTw;
					
					assert (Xw1 > 0);
					
					++queueops;
					
					wrheap->emplace( std::pair<long double, int>(-(wc_wi+Xw1), k) );
					
					
				}
				
				assert (wrheap->size() == maxSize);
				
			} else {
			
				assert (wrheap->size() == maxSize);
				
				
			
				while ( -wrheap->top().first <= wc_wi){
					
					//cout << -wrheap->top().first << " -- " << wc_wi << endl;
					
					const long k = wrheap->top().second;
					
					const long removeId = wrids.at(k);
					
					remove(removeId);
					clearRemovals();
					
					wrids.at(k) = addToSample(elem, wi);
					
					wrheap->pop();
					
					const Weight logtw = logTws.at(k)*wi;
					
					const Weight r2 = 1.0-(1.0-exp(logtw))*rnd->randomDouble();
					
					//assert (r2 >= exp(logtw));
					//assert (r2 <= 1);
					
					const Weight logki = log(r2)/wi;
					
					logTws.at(k) = logki;
					
					const Weight logr = log(rnd->randomDouble());
				
					const Weight Xw1 = logr/logTws.at(k);
					
					assert (Xw1 > 0);
					
					++queueops;
					
					wrheap->emplace( std::pair<Weight, int>(-(wc_wi+Xw1), k) );
				}
				
			}
			
			return;
		}
		
		if (sampler == WeightedSamplers::OFFLINE_WR){
			
			while (targets->size() > 0 && targets->back() <= populationWeight){
			
				addToSample(elem, w);
				targets->pop_back();
			}
			
			return;
		}
		
		//cout << "add*" << endl;
		if (sampler == WeightedSamplers::ONLINE_WOR){
		
			const Weight wi = w;
			
			if (wi == 0)
				return;
				
			Weight ki = 1;
			
			requireRandom();
		
			if ( (currentSize < maxSize) || (!exponentialJumps) ){
			
				ki = log(rnd->randomDouble() )/wi;
			
			} else if (exponentialJumps){
				
				if (wc_wi < 0){
					
					wc_wi = 0;
				
					const Weight r = rnd->randomDouble();
					
					Xw = log(r)/logTw;
				}
				
				assert (Xw >= 0);
				assert (wc_wi >= 0);
				
				const bool cond1 = wc_wi < Xw;
		
				wc_wi += wi;
				
				const bool cond2 = Xw <= wc_wi;
		
				if ( cond1 && cond2 ){
				
					wc_wi = -1;
					Xw = -1;
					
					const Weight tw = exp(logTw*wi);
			
					const Weight min = tw;
					const Weight max = 1;
					
					const Weight r2 = min+rnd->randomDouble()*(max-min);
					ki = log(r2)/wi;
					
				} else {
		
					return;
				}
			}
		
			assert (ki != 1);
			assert (currentSize <= maxSize);
			
			if (currentSize == maxSize){
		
				if (ki < -heap->top().u){
					
					if (exponentialJumps){
					
						cout << "ki " << ki << " currentlySmallestRandomValue " << -heap->top().u << endl;
					}
				
					//assert (!exponentialJumps);
					return;
				}
			
				const long removeId = heap->top().id;
				
				remove(removeId);
				
				clearRemovals();
			
				heap->pop(); // remove item with currently smallest randomValue
			
				assert (heap->size() == currentSize);
			}
		
			assert (currentSize < maxSize);
			const long identity = addToSample(elem, w, -1, ki);
			
			++queueops;
			
			heap->push( SampleEntry(-ki, identity) );
		
			if (exponentialJumps)
				logTw = -heap->top().u;
			
			assert (heap->size() == currentSize);
		
			assert (currentSize <= maxSize);
		
			assert (weightsCumulative == false);
			assert (weights.size() == identities.size() );
			assert( weights.size()-toRemove.size() == currentSize);
		}
		
		//cout << "*add" << endl;
	}
	
	inline bool isWR() const{
	
		return wr;
	}
	
	inline bool isWOR() const{
	
		return wor;
	}
	
	inline bool isPopulation() const{
	
		return isPop;
	}
	
	
	inline void withReplacement(long m = -1){
		
		assert (m <= maxSize);
	
		finalise(true, m);
	}
	
	inline void withoutReplacement(long m = -1){
	
		assert (m <= maxSize);
		
		finalise(false, m);
	}
	
	inline virtual void internalPermute(const vector<long>& v){
	
		UtilCollection::permute(v, elements);
	}
	
	inline virtual long getInvalidSize() const{
	
		return 0;
	}
};




class WeightedRowSample : public WeightedSample<Row>{

typedef long_double Weight;

public:
	const Tuple cols;
	

private:
	Rows rows;
	
	long invalidSize = 0;
	Rows batchRows;
	
public:
	
	class iterator{
	
		Rows::iterator rowsIt;
		bool rowsItEnd;
		vector<int>::const_iterator freqIt;
		vector<int>::const_iterator freqItEnd;
		vector<Weight>::const_iterator weightIt;
		vector<Weight>::const_iterator weightItEnd;
		
		
		Weight lastWeight = 0;
		
		long freqsSize = 0;
		
		public:
		inline iterator(const Rows& rows, const vector<int>& freqs, const vector<long_double>& weights) :
		
		freqsSize(freqs.size()),
		
		rowsIt(rows.begin() ),
		rowsItEnd(rows.end() ),
		
		freqIt(freqs.cbegin() ),
		freqItEnd(freqs.cend() ),
		
		weightIt(weights.cbegin() ),
		weightItEnd(weights.cend() )
		
		
		{
		
			
			
		}
		
		inline iterator(const Rows& rows, const vector<int>& freqs, const vector<long_double>& weights, bool last) :
		
		rowsIt(rows.beginLast() ),
		rowsItEnd(rows.end() ),
		
		freqIt( freqs.begin()+(freqs.size()-1) ),
		freqItEnd(freqs.cend() ),
		
		weightIt(weights.cend()+(freqs.size()-1) ),
		weightItEnd(weights.cend() )
		
		
		{
			
			lastWeight = freqs.size() > 1 ? weights.at(freqs.size()-2) : 0;
			
		}
		
		inline bool operator!=(const bool& b) const{
		
			if ( !(rowsIt != rowsItEnd) )
				return false;
		
			if (freqsSize > 0 && !(freqIt != freqItEnd) )
				return false;
		
			if ( !(weightIt != weightItEnd) )
				return false;
		
			return true;
		}
		
		inline bool operator==(const bool& b) const{
		
			return !( (*this) != b);
		}
		
		inline bool operator()() const{
		
			return true;
		}
	
		inline const Row& operator*() const {
		
			return (*rowsIt);
		}
		
		inline long_double getWeight() const {
		
			return (*weightIt)-lastWeight;
		}
		
		inline long_double getCumulativeWeight() const {
		
			return (*weightIt);
		}
		
		inline const long getFrequency() const {
		
			return freqsSize > 0 ? (*freqIt) : 1;
		}
		
		inline iterator& operator++(){
		
			if (weightIt != weightItEnd)
				lastWeight = (*weightIt);
			
			++rowsIt;
			++weightIt;
			
			if (freqsSize > 0)
				++freqIt;
			
			return (*this);
		}
	};
	
	const iterator end;
	
	
	inline WeightedRowSample(Tuple cols, const string s, long size=0, bool pop = false) : WeightedSample<Row>(s, size, pop), cols(cols), rows(cols.size() ),batchRows(cols.size() ), end(WeightedRowSample::iterator(rows, this->frequencies, this->weights )) {
	
	
		//cout << "WEIGHTED ROW SAMPLE " << cols << " max size " << size << endl;
	}
	
	inline void reserveInternal(long size =-1) override{
	
		if (size == -1)
			size = maxSize;
		
		//cout << "cols.size() " << cols.size() << " size " << size << " maxRemovalQueue " << maxRemovalQueue << endl;
		cout << "reserveInternal " << UtilString::bytesToString(8*cols.size() *size+8*maxRemovalQueue) << endl;
		
		rows.reserve(8*cols.size() *size+8*maxRemovalQueue);
		
		if (false){
		
			cout << "batchrows reserveInternal " << 8*cols.size()*size << endl;
			batchRows.reserve(8*cols.size()*size);	
			
		}
		
		//cout << "reserveInternal " << (16*cols.size() *maxSize+8*maxRemovalQueue ) << " capacity " << rows.capacity() << endl;
	}
	
	
	inline void addBatchToSampleInternal() override{
		
		long j = 0;
		for (auto it = batchRows.begin(); it != batchRows.end(); ++it){
		
			addToSample( (*it), batchWeights.at(j), batchFrequencies.at(j));
			++j;
		}
		
		batchRows.clear();
		
	}
	
	
	inline void shrinkInternal() override{
	
		batchRows.clear();
		batchRows.shrink_to_fit();
	
		rows.shrink_to_fit();
	}
	
	template<typename T>
	inline void getKS(const std::map<Row, long>& index, const Row& last, WeightingFunction& f, T& table, vector<long_double>& less, vector<long_double>& equal) const{
	
		less.clear();
		equal.clear();
		
		const long len = index.size();
		
		less.reserve( len+1);
		equal.reserve( len+1);
		
		const Weighter g(f, cols);
		
		for (long j = 0; j <= len; ++j)
			equal.push_back(0);
		
		for (long j = 0; j <= len; ++j)
			less.push_back(0);
		
		cout << "scan*" << endl;
		for (auto it = table.begin(); it != table.end(); ++it){
			
			long_double weight = g(*it);
			
			const Row& r = (*it);
			
			if (weight > 0){
			
				if (last < r){
				
					less.back() += weight;
					
				} else {
				
					auto it2 = index.lower_bound(r); //std::lower_bound(index.begin(), index.end(), r);
					
					const long ind = it2->second; //it2-index.begin();
					
					if ( r < it2->first ){
					
						less.at(ind) += weight;
						
					} else {
					
						equal.at(ind) += weight;
					}
				}
			}
		}
		
		cout << "*scan" << endl;
		
		long_double wsum = UtilCollection::sum(equal)+UtilCollection::sum(less);
		long_double divmul = 1.0L/wsum;
		
		for (long j = 0; j < equal.size(); ++j)
			equal.at(j) *= divmul;
			
		for (long j = 0; j < less.size(); ++j)
			less.at(j) *= divmul;
			
		for (long j = 1; j < less.size(); ++j)
			less.at(j) += less.at(j-1)+equal.at(j-1); // less.at(2) = less.at(2)+less.at(1)+equal.at(1)
			
		for (long j = 1; j < equal.size(); ++j)
			equal.at(j) += less.at(j); // equal.at(2) = equal.at(2)+less.at(2)+less.at(1)+equal.at(1)
		
	
		//for (long j = 0; j < w.size(); ++j)
		//	cout << w[j] << endl;
	}
	
	inline long_double getMax(long_double max,const vector<long_double>& w1, long j1, const vector<long_double>& w2, long j2) const{
	
		const long_double a = j1 <= 0 ? 0 : j1 >= w1.size() ? 1 : w1.at(j1);
		const long_double b = j2 <= 0 ? 0 : j2 >= w2.size() ? 1 : w2.at(j2);
		
		if ( (b-a) > max)
			return b-a;
		
		if ( (a-b) > max)
			return a-b;
			
		return max;
			
	}
	
	template<typename T>
	inline long_double getKS(WeightingFunction& f, T& table) const{
	
	
		std::map<Row, long> index;
		
		vector<Row> samples;
		
		//index.reserve(sampleSize);
		
		samples.reserve(sampleSize);
		
		for (auto it = begin(); it != end(); ++it){
		
			const Row& r = (*it);
			
			const long freq = it.getFrequency();
		
			for (long k = 0; k < freq; k++)
				samples.push_back(r);
		}
		
		std::sort(samples.begin(), samples.end() );
		
		{
			long ind = 0;
		
			for (auto it = samples.begin(); it != samples.end(); ++it){
		
				const Row& r = (*it);
			
				if (index.count(r) == 0){
				
					index.insert( {r, ind});
					++ind;
				}
			}
		}
		
		/*cout << table.getColNames() << endl;
		
		for (auto it = index.begin(); it != index.end(); ++it){
		
			cout << (*it) << endl;
		}*/
		
		//assert (index.size() == sampleSize);
		
		vector<long_double> smpEqual;
		vector<long_double> smpLess;
		
		vector<long_double> tabLess;
		vector<long_double> tabEqual;
		
		WeightingFunction g;
		
		getKS(index, samples.back(), g, samples, smpLess, smpEqual);
		getKS(index, samples.back(), f, table, tabLess, tabEqual);
		
		long_double max = 0;
		
		for (long j = 0; j <= index.size(); ++j){
			
			// [before: j*2][at: j*2+1]...[at: j*2+1][after: j*2+2]
			
			max = getMax(max, tabLess, j, smpEqual, j-1);
			max = getMax(max, tabLess, j, smpEqual, j);

			max = getMax(max, tabEqual, j-1, smpLess, j);
			max = getMax(max, tabEqual, j, smpLess, j);

		}
		
		return max;
		
	}
	
	inline long getBytes() const {
	
		return rows.getBytes();
	}
	
	inline void internalPermute(const vector<long>& v) override{
	
		rows.permute(v);
	}
	
	
	
	inline long getInvalidSize() const override{
	
		return invalidSize;
	}
	
	inline void setInvalidSize(long w){
	
		invalidSize = w;
	}
	
	
	inline void setMaxSize(long w){
	
		maxSize = w;
	}
	
	
	inline void removeLast(){
	
		assert (wor || wr);
		
		assert (invalidSize == 0);
	
		assert (sampleSize > 0);
		
		--sampleSize;
		--frequencies.back();
	
		if (frequencies.back() > 0)
			return;
		
		rows.pop_back();
		weights.pop_back();
		
		if (frequencies.size() > 0)
			frequencies.pop_back();
		
		if (identities.size() > 0)
			identities.pop_back();
		
		
	}
	
	inline ~WeightedRowSample(){
	
	}
	
	
	inline void push_back_invalid(long_double w = 0){
	
		addPopulation(w);
		
		invalidSize++;
	}
	
	
	inline void push_back(const Row& r1, long_double w = 1){
	
	
		//addPopulation(w);
		
		
		rows.push_back(r1);
		
		if (identities.size() == 0)
			identities.push_back(0);
		else
			identities.push_back(identities.back()+1);
		
		//frequencies.push_back(1);
		weights.push_back(w);
	}
	
	inline void push_back(vector<const Row*>& rs, vector<vector<int>>& sel, long_double w = 1){
	
		rows.push_back(rs, sel);
	
		if (identities.size() == 0)
			identities.push_back(0);
		else
			identities.push_back(identities.back()+1);
		
		frequencies.push_back(1);
		weights.push_back(w);
	}
	
	inline void push_back(const Row& r1, const Row& r2, long_double w = 1){
	
	
		//addPopulation(w);
		
		
		rows.push_back(r1, r2);
		
		if (identities.size() == 0)
			identities.push_back(0);
		else
			identities.push_back(identities.back()+1);
		
		frequencies.push_back(1);
		weights.push_back(w);
	}
	
	inline void push_back(const Row& r1, const vector<int>& v1, const Row& r2, long_double w = 1){
	
	
		//addPopulation(w);
		rows.push_back(r1, v1, r2);
		
		if (identities.size() == 0)
			identities.push_back(0);
		else
			identities.push_back(identities.back()+1);
		
		frequencies.push_back(1);
		weights.push_back(w);
	}
	
	inline void push_back(const Row& r1, const Row& r2, const vector<int>& v2, long_double w = 1){
	
	
		//addPopulation(w);
		rows.push_back(r1, r2, v2);
		
		if (identities.size() == 0)
			identities.push_back(0);
		else
			identities.push_back(identities.back()+1);
		
		frequencies.push_back(1);
		weights.push_back(w);
	}
	
	
	inline void push_back(const Row& r1, const vector<int>& v1, const Row& r2, const vector<int>& v2, long_double w = 1){
	
	
		//addPopulation(w);
		
		
		rows.push_back(r1, v1, r2, v2);
		
		if (identities.size() == 0)
			identities.push_back(0);
		else
			identities.push_back(identities.back()+1);
		
		frequencies.push_back(1);
		weights.push_back(w);
	}
	
	
	
	inline bool compatible(const WeightedRowSample& samp) const {
	
		if (maxSize != samp.maxSize)
			return false;
		
		if (isWOR() != samp.isWOR())
			return false;
		
		if (isWR() != samp.isWR() )
			return false;
			
		if (getPopulationWeight() != samp.getPopulationWeight())
			return false;
			
		if (getPopulationSize() != samp.getPopulationSize())
			return false;
		
		if (!UtilTuple::equals(cols, samp.cols) )
			return false;
			
		return true;
	}
	
	
	
	inline void push_back(const WeightedRowSample& samp){
	
		assert (sampler == WeightedSamplers::FINAL);
		
		assert (samp.sampler == WeightedSamplers::FINAL);
	
		if ( (samp.getSize() == 0) ||(samp.invalidSize == samp.maxSize) ){
		
			invalidSize += samp.invalidSize;
			return;
		}
		
		
		if (!compatible(samp)){
		
			cout << samp << endl;
		
			assert (maxSize == samp.maxSize);
		
			assert (wor == samp.wor);
			assert (wr == samp.wr);
		
			assert( populationWeight == samp.populationWeight);
			assert( populationSize == samp.populationSize);
		
			assert( UtilTuple::equals(cols, samp.cols) );
		}
		
		invalidSize += samp.invalidSize;
		
		//long idmax = identities.back()+1;
		
		
		rows.reserve(rows.getBytes()+samp.rows.getBytes() );
		
		for (auto it = samp.begin(); it != samp.end(); ++it){
		
			//identities.push_back(idmax++);
			
			const Row& r = (*it);
			const long f = it.getFrequency();
			const double long w = it.getCumulativeWeight();
			
			weights.push_back(w );
			frequencies.push_back(it.getFrequency() );
			sampleSize += it.getFrequency();
			
			rows.push_back(r);
		}
	}
	
	
	inline vector<int> getVector(const Tuple& t1, const Tuple& t2, int k){
	
		JoinColumns rel;
		
		Tuple t3 = t1;
		t3.push_back(t2);
		
		rel.setTarget(t3);
		rel.setParent(t1);
		rel.setMain(t2);
		
		//const Tuple& target = rel.getTargetTuple();
		//const Tuple& parent = rel.getParentTuple();
		//const Tuple& main = rel.getMainTuple();
		
		const vector<int>& parentToMain = rel.getTranslationTable(JoinColumn::PARENT, JoinColumn::MAIN);
		const vector<int>& mainToParent = rel.getTranslationTable(JoinColumn::MAIN, JoinColumn::PARENT);
		
		vector<int> ret;
		
		if (k == 0){
		
			for (auto it = parentToMain.begin(); it != parentToMain.end(); ++it)
				if ( (*it) >= 0)
					ret.push_back(*it);
		
		} else if (k == 1){
		
			for (auto it = mainToParent.begin(); it != mainToParent.end(); ++it)
				if ( (*it) >= 0)
					ret.push_back(*it);
		
		} else {
		
			assert (k >= 0);
			assert (k <= 1);
		}
		
		return ret;
	}
	
	inline void compress() override{
		
		if (true)
			return;
		
		assert (weightsCumulative == false);
	
		clearRemovals(true);
			
		Row oldrow;
		
		long newj = 0;
		long oldj = -1;
		
		for (auto it = rows.begin(); it != rows.end(); ++it){
		
			const Row& newrow = (*it);
			
			if (oldj >= 0 && oldrow.equals(newrow)){
			
				// merge new into old
				frequencies.at(oldj) += frequencies.at(newj);				
				weights.at(oldj) += weights.at(newj);
				
				weights.at(newj) = 0;
				frequencies.at(newj) = 0;
				
			} else {
			
				oldrow = newrow;
				oldj = newj;
			}
			
			++newj;
		}
		
		for (long j = 0; j < frequencies.size(); j++){
			if (frequencies.at(j) == 0)
				remove(identities.at(j));
		}
		
		clearRemovals(true);
	}
	/*
	inline void add(const Row& a, const vector<int>& asel, const Row& b, const vector<int>& bsel, long freq = 1, long_double weight = 1){
	
		
		
	}*/
	
	friend inline std::ostream& operator<<(std::ostream &os, const WeightedRowSample& v) { 
    	
    	os << v.toString();
		return os;
	}
	
	
	inline WeightedRowSample::iterator begin() const{
	
		assert (sampler == WeightedSamplers::FINAL);
		assert (weightsCumulative);
	
		return WeightedRowSample::iterator(rows, this->frequencies, this->weights );
	}
	
	inline WeightedRowSample::iterator beginLast() const{
	
		assert (sampler == WeightedSamplers::FINAL);
		assert (weightsCumulative);
		
		return WeightedRowSample::iterator(rows, this->frequencies, this->weights, true );
	}
	
	inline shared_ptr<vector<Row>> getRowVector() const {
		
		shared_ptr<vector<Row>> rows;
		
		auto ret = std::make_shared<vector<Row>>();
		
		ret->reserve(sampleSize);
		
		for (auto it = begin(); it != end(); ++it){
			
			for (int i = (it).getFrequency(); i >= 1; --i)
				ret->push_back(*it);
		}
		
		return ret;
	}
	
	template<typename T>
	inline shared_ptr<T> getTable(string name) const {
		
		shared_ptr<Rows> rows = std::make_shared<Rows>(cols.size() );
		
		long bytes = 0;
		
		long smpsize = 0;
		
		for (auto it = begin(); it != end(); ++it){
			bytes += (*it).getBytes() * (it).getFrequency();
			
			smpsize += (it).getFrequency();
		}
		
		rows->reserve(bytes);
		
		for (auto it = begin(); it != end(); ++it){
			
			for (int i = (it).getFrequency(); i >= 1; --i)
				rows->push_back(*it);
		}
		
		return std::make_shared<T>(name, cols, rows);
	}
	
	
	template<typename T>
	inline shared_ptr<T> getTable(string name, const Tuple& target) const {
		
		shared_ptr<Rows> rows = std::make_shared<Rows>(cols.size() );
		
		long bytes = 0;
		
		long smpsize = 0;
		
		for (auto it = begin(); it != end(); ++it){
			bytes += (*it).getBytes() * (it).getFrequency();
			
			smpsize += (it).getFrequency();
		}
		
		rows->reserve(bytes);
		
		for (auto it = begin(); it != end(); ++it){
			
			for (int i = (it).getFrequency(); i >= 1; --i)
				rows->push_back(*it);
		}
		
		return std::make_shared<T>(name, cols, rows);
	}
	
	inline long_double getSampleWeight() const{
	
		assert (weightsCumulative);
		
		if (weights.size() > 0)
			return weights.back();
			
		return 0;
	}
	
	inline const string json() const {
	
		stringstream ss;
	
		ss << "{";
		
		ss << " \"cols\": \"" << cols << "\"";
		ss << ", ";
		
		ss << " \"popsize\": " << populationSize;
		ss << ", ";
		
		
		ss << " \"smpweight\": " << UtilString::doubleToString(getSampleWeight());
		ss << ", ";
		
		
		double est = 0;
		
		if (sampleSize == 0){
		
			est = ((long) (populationWeight/invalidSize));
			
		} else {
			
			const long_double a = sampleSize;
			const long_double b = sampleSize+invalidSize;
		
			est = (long) (populationWeight/b*a);
		}
		
		ss << " \"maxsmp\": " << maxSize;
		ss << ", ";
		
		ss << " \"invalid\": " << invalidSize;
		ss << ", ";
		
		ss << " \"smpsize\": " << sampleSize;
		ss << ", ";
		
		ss << " \"popweight\": " << UtilString::doubleToString(populationWeight);
		ss << ", ";
		
		ss << " \"estweight\": " << UtilString::doubleToString(est);
		
		ss << "}";
		
		return ss.str();
	
	}
	
	inline const string toString(int j = -1) const {
	
		assert (isPop || (wr || wor));
	
		stringstream ss;
	
		for (auto it = begin(); it != end(); ++it){
		
			if (j == 0)
				break;
			
			ss << j << " " << (*it) << " w " << it.getWeight() << " # " << it.getFrequency() << endl;
			
			--j;
		}
		
		ss << cols << " pop size " << populationSize << " pop weight " << UtilString::doubleToString(populationWeight) << " smp weight " << UtilString::doubleToString(getSampleWeight()); // << " size " << populationSize << endl;
		
		
		ss << " maxSize " << maxSize << " size " << sampleSize << " invalid " << invalidSize << " estimate ";
		
		if (sampleSize == 0){
		
			ss << ">" << ((long) (populationWeight/invalidSize));
			
		} else {
			
			const long_double a = sampleSize;
			const long_double b = sampleSize+invalidSize;
		
			ss << (long) (populationWeight/b*a);
		}
		
		if (wr)
			ss << " WR";
			
		if (wor)
			ss << " WOR";
		
		return ss.str();
	}
	
	inline void push_back(const Tuple& t, long_double w = 1){
	
		
	
		Rows rows(t.size() );
		rows.push_back(t);
		
		Row row;
		rows.back(row);
		
		this->add(row, w);
	}
	
	
	inline void clearToRemove() override {
	
		// CLEAR WEIGHTS
		UtilSortedVector::remove(this->identities, this->toRemove, this->weights);
		UtilSortedVector::remove(this->identities, this->toRemove, this->frequencies);
		UtilSortedVector::remove(this->identities, this->toRemove, this->randomKeys);
		
		rows.remove(this->identities, this->toRemove);
	}
	
	inline void addInternal(const Row& row) override{
		
		this->rows.push_back(row);
	}
	
	inline void addToBatchInternal(const Row& row) override{
		
		this->batchRows.push_back(row);
	}
	
	
	shared_ptr<vector<Row> > rowIndex;
	
	
	inline const Row& getRow(double randomDouble1){
	
		assert (isPop);
	
		if (rowIndex){
		
		
		} else {
		
			rowIndex = std::make_shared<vector<Row>>();
			
			rowIndex->reserve(sampleSize);
		
			for (auto it = begin(); it != end(); ++it){
			
				rowIndex->push_back((*it) );
			}
		}
		
		const long_double sum = weights.back();
		
		// 0.3 0.9 
		
		const long_double u = randomDouble1*sum;
		
		const long ind  = u >= sum ? weights.size()-1 : std::lower_bound(weights.begin(), weights.end(), u )-weights.begin();
		
		assert (u < weights.at(ind));
		
		if (ind > 0)
			assert (weights.at(ind-1) <= u);
		
		assert (ind < weights.size() );
		
		//if (ind >= rowIndex->size() )
		//	return rowIndex->at(rowIndex->size()-1);
		
		assert (ind == 0 ? weights.at(ind) > 0 : weights.at(ind) > weights.at(ind-1));
		
		const Row& ret = rowIndex->at(ind) ;
		
		
		return ret;
		
	}
	
	// treats this sample as full population
	inline shared_ptr<WeightedRowSample> getWithReplacementSampleFromPopulation(long m, const string s ="abc"){
	
		assert (isPop);
		
		if (rowIndex){
		
			
		} else {
			
			rowIndex = std::make_shared<vector<Row>>();
			
			rowIndex->reserve(sampleSize);
		
			for (auto it = begin(); it != end(); ++it){
			
				rowIndex->push_back((*it) );
			}
		}
		
		
		if (rowIndex){
			
			auto ret = std::make_shared<WeightedRowSample>(cols, s, m );
		
			ret->initPop();
		
			ret->populationWeight = populationWeight;
			
			this->requireRandom();
			
			const long_double sum = weights.back();
		
			for (long j = 0; j < m; ++j){
			
				auto it = std::lower_bound(weights.begin(), weights.end(), rnd->randomDouble()*sum );
				
				const long pos = it-weights.begin();
				
				if (pos >= rowIndex->size() )
					ret->add( rowIndex->at(rowIndex->size()-1));
				else
					ret->add( rowIndex->at(pos));
				
			}//ret->add(rowIndex->at(getIndependentSample()));
			
			
			ret->finalise();
			
			rnd.reset();
			
			return ret;
		}
		
		assert (false);
		
	}
};


template<typename T>
inline shared_ptr<WeightedRowSample> withrepl(T& in, const Tuple& t2, long smps, string seed="abc123"){

	shared_ptr<WeightedRowSample> s2 = std::make_shared<WeightedRowSample>(t2, seed, smps );	
	
	/*Weighter weight(f, t2);
	
	double wsum = 0;

	for (auto it = r2.begin(); it != r2.end(); ++it)
		wsum += weight(*it);			
	
	s2->initOffline(wsum, true);
	*/
	
	long wsum = 0;
	
	for (auto it = in.begin(); it != in.end(); ++it)
		wsum += 1;
	
	s2->initOffline(wsum, true);
	//s2->initOnlineWR();

	for (auto it = in.begin(); it != in.end(); ++it)
		s2->add(*it); //, weight(*it)

	s2->withReplacement();
	
	return s2;
}

template<typename T>
inline shared_ptr<WeightedRowSample> withoutrepl(T& in, const Tuple& t2, long smps, string seed="abc123"){

	shared_ptr<WeightedRowSample> s2 = std::make_shared<WeightedRowSample>(t2, seed, smps );	
	
	/*Weighter weight(f, t2);
	
	double wsum = 0;

	for (auto it = r2.begin(); it != r2.end(); ++it)
		wsum += weight(*it);			
	
	s2->initOffline(wsum, true);*/
	
	s2->initOnlineWOR();

	for (auto it = in.begin(); it != in.end(); ++it)
		s2->add(*it); //, weight(*it)

	s2->withoutReplacement();
	
	return s2;
	
}


class TupleStream{

public:

	
	virtual ~TupleStream(){}

	
	inline virtual bool read(Row& p){
	
		return false;
	}
	
	Row dummy;
	
	inline virtual bool skip(long off){
	
		for (long j = 0; j < off; ++j)
			read(dummy);
			
		return true;
	}
	
	inline virtual void close() = 0;

};


/*
class ODBCInterface : public TupleStream{

	SQLHENV env = NULL;
    SQLHDBC dbc = NULL;
    
    //long res;
	
    const static int bufsize = 1; //1000;
    const static int maxlen = 64;
    
    inline SQLCHAR* NTS(const string& s){
	
		return (SQLCHAR*) s.c_str();
	}
    
    struct Column {
 		SQLCHAR data[bufsize][maxlen];
 		SQLLEN orind[bufsize];
	};
    
    vector<Column> columns;
    
public:
	
    inline ODBCInterface( const string& db){
    
        SQLAllocHandle( SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env);
 
        SQLSetEnvAttr( env, SQL_ATTR_ODBC_VERSION, (void*)SQL_OV_ODBC3, 0);
 
		SQLAllocHandle( SQL_HANDLE_DBC, env, &dbc);

		SQLDriverConnect(dbc, NULL, NTS("DSN="+db+";"), SQL_NTS,
                 NULL, 0, NULL, SQL_DRIVER_COMPLETE);
    }
    
    
    inline ODBCInterface( const string& db, const string& query, int cs) : ODBCInterface(db){
    
        startQuery(query, cs);
    }
    
    
    
	inline ~ODBCInterface(){

		close();
	
	}	  
    
    SQLHSTMT    hstmt = NULL;
    
	long row_count = 0;
	
	SQLUSMALLINT row_status[ bufsize ];
 		
 	int cols = 0;
 	int rowpos = -1;
	long rowid = 0;
	
	
    bool closed = true;
    
	
    inline void startQuery(const string& query, int cs){
    
    	assert (closed);
    
    	closed = false;
    	
    	
    	if (hstmt!=SQL_NULL_HSTMT) {
        	SQLFreeHandle(SQL_HANDLE_STMT, hstmt);
			hstmt = SQL_NULL_HSTMT;
		}
		
		cols = cs;
    
    	SQLAllocHandle(SQL_HANDLE_STMT, dbc, &hstmt);
    	
    	SQLExecDirect(hstmt, NTS(query), SQL_NTS);
    	
    	columns = vector<Column>(cols);	
    	
    	if (bufsize > 1){
    		
			SQLSetStmtAttr( hstmt, SQL_ATTR_ROW_ARRAY_SIZE,
				(SQLPOINTER)bufsize, 0 );
			
			for (int i = 1; i <= cols; i++){
			
				SQLBindCol( hstmt,    // Statement handle
				i,                // Column number
				SQL_C_CHAR,       // Bind to a C string
				columns[i-1].data,            // The data to be fetched
				maxlen,        // Maximum length of the data
				columns[i-1].orind );
			}
			
			SQLSetStmtAttr( hstmt, SQL_ATTR_ROWS_FETCHED_PTR, &row_count, 0 );
			SQLSetStmtAttr( hstmt, SQL_ATTR_ROW_STATUS_PTR, row_status, 0 );
		}
		
		row_count = 0;
		rowpos = 0;
		
		hasNextBatch = true;
    }
    
    bool hasNextBatch = true;
    
    inline void readBatch(){
    
    	if (bufsize == 0){
    	
    		row_count = 1;
			rowpos = 0;
				
			hasNextBatch = true;			
			return;
    	
    	}
    	
    	if (bufsize > 1){
    
			SQLRETURN ret = SQLFetchScroll( hstmt, SQL_FETCH_NEXT, 0 );
		
			if (ret == SQL_ERROR || ret == SQL_SUCCESS_WITH_INFO) {  
         		
         		cout << "ERROR FETCH2 !!" << endl;
      		}
		
			hasNextBatch = SQL_SUCCEEDED( ret ) && row_count > 0;
		
			if (hasNextBatch)
				rowpos = 0;
			else
				row_count = 0;
			
		} else {
			
			SQLRETURN ret = SQLFetch(hstmt);
			
			if (ret == SQL_ERROR || ret == SQL_SUCCESS_WITH_INFO) {  
				
         		cout << "ERROR FETCH1 !!" << endl;
      		}
			
			if (ret == SQL_SUCCESS || ret == SQL_SUCCESS_WITH_INFO){
			
				for (int i = 0; i < cols; i++){
					SQLGetData(hstmt, i+1, SQL_C_CHAR, columns[i].data[0], maxlen, &(columns[i].orind[0]));  
				}
				
				row_count = 1;
				rowpos = 0;
				
				hasNextBatch = true;			
				
			} else {
				
				row_count = 0;
				hasNextBatch = false;
			}
			
			//SQLCloseCursor(hstmt);
			
		}
    }
    
    
    inline const string getString(int col, int row) const{
    
    	return string( ((char*)(&columns[col].data[row])) );
    }
    
    
    inline long getInt(const string& cmd){
    
    	startQuery(cmd, 1);
    	
    	Tuple t;
    	
    	read(t);
    	
    	return UtilString::strToInt(t[0]);
    	
    }
    
    
    inline bool read(Tuple& tuple) override{
    		
    	if (rowpos == row_count){
    		
    		if (hasNextBatch)
    			readBatch();	
    		else
    			return false;
    	}
    		
    	if (row_count == 0)
    		return false;
    	
    	
    	tuple.clear();
    	
    	for (int col = 0; col < cols; col++)
    		tuple.push_back( getString(col, rowpos) );
    	
    	//if (rowid % 10000 == 0)
    	//	cout << rowid << endl;
    	
    	rowpos++;
    	rowid++;
    	
    	return true;
    }
    
    inline bool skip(long l) override{
    
    	return false;
    }
    
	inline void close() override{
	
		//if (closed)
		//	return;
		
		//cout << "CLOSE" << endl;
	
		if (hstmt!=SQL_NULL_HSTMT) {
        	SQLFreeHandle(SQL_HANDLE_STMT, hstmt);
			hstmt = SQL_NULL_HSTMT;
		}
	
        if (dbc!=SQL_NULL_HDBC) {
        SQLDisconnect(dbc);
        SQLFreeHandle(SQL_HANDLE_DBC, dbc);
		dbc = SQL_NULL_HDBC;
    	}

    	if (env!=SQL_NULL_HENV) {
			SQLFreeHandle(SQL_HANDLE_ENV, env);
			env = SQL_NULL_HENV;
   		 }
   		
   		closed = true;
	}
};*/


class CSVStream : public TupleStream{

	const Tuple colNames;
	
	bool open = true;
	
	io::LineReader in;
	
	const char sep;
	
	Rows rows;
	
	const long colMask;
	
	const std::bitset<1000> colSel;
	
public:	

	inline void close() override{
	
	
	}
	
	bool addColNames = false;
	
	inline bool skip(long n) override{
	
		if (!open)
			return false;
			
		if (n == 0)
			return true;
	
		for (int i = 0; i < n; i++){
			
			if (!open)
				return false;
		
			char* c = in.next_line();
			
			if (c == 0){
				open = false;
				break;
			}
		}
		
		return true;
	}


	inline bool read(Row& t) override{
		
		if (!open)
			return false;
		
		long colId = 0;
	
		//long col = 1;
		//bool write = col & colMask;
		
		long col = 0;
		
		bool write = colSel.test(col);
		
		const char* c = in.next_line();
		
		//cout << "read row [" << colNames << "]" << endl;
		
		rows.clear();
		
		const bool DEBUG = false;
		
		const int colNameSize = (addColNames || DEBUG) ? colNames.size() : 0; // ; // DEBUG CODE
		
		bool firstInCol = true;  // DEBUG CODE
		
		if (c){
		
			rows.startRow();
			
			for (const char* d = c; true ; d++){
				
				if ( ((*d) == sep) || ((*d) == 0)){
		
					if (write){
					
						if (firstInCol){
						
							if (colId < colNameSize ){
						
								const string& colName = colNames[colId];
						
								for (int j = 0; j < colName.length(); ++j){
							
									rows.push_back_char(colName[j]);
									
								}
								rows.push_back_char('=');
							}
							
							rows.push_back_char('N');
							rows.push_back_char('U');
							rows.push_back_char('L');
							rows.push_back_char('L');
						}
						
						rows.push_back_char(0);
						
						++colId;
						firstInCol = true;
					}
					
					if ((*d) == 0)
						break;
					
					
					
					//col <<= 1;
					//write = (col & colMask) != 0;
					
					++col;
					
					assert (col < 1000);
					
					write = colSel.test(col);
					
				} else {
				
					if (write){
					
						// BEGIN DEBUG CODE
						if (firstInCol){
						
							if (colId < colNameSize ){
						
								const string& colName = colNames[colId];
						
								for (int j = 0; j < colName.length(); ++j){
							
									rows.push_back_char(colName[j]);
									
								}
								rows.push_back_char('=');
							}
							
							firstInCol = false;
						}
						// END DEBUG CODE
						
						rows.push_back_char( (*d) );
						
					}
				}
			}
			
			rows.finishRow();
			
			assert (colId == rows.cols);
			
			rows.back(t);
			
			return true;
		
		} else {
		
			open = false;
			return false;
		}
	}

	static std::bitset<1000> getColSel(const vector<int>& seps){
	
		std::bitset<1000> ret;
	
		for (int i = 0; i < seps.size(); i++)
			ret.set(seps.at(i));
		
		return ret;
	}
	
	static long getColMask(const vector<int>& seps){
	
		long ret = 0;
	
		for (int i = 0; i < seps.size(); i++)
			ret |= 1L << seps[i];
		
		return ret;
	}
	
	inline CSVStream(const string& s, const vector<int>& cols, const Tuple& colnames, char sep=',', bool ignoreFirst = true) : colNames(colnames), in(s), sep(sep), rows(cols.size()), colMask(getColMask(cols)), colSel(getColSel(cols)){
	
		rows.reserve(1000);
		
		if (ignoreFirst){
		
			char* c = in.next_line();
			
			if (c == 0)
				open = false;
		}
		
	}
	
	inline CSVStream(const string& s, const vector<int>& cols, char sep=',', bool ignoreFirst = true) :colNames(),  in(s), sep(sep), rows(cols.size()), colMask(getColMask(cols)), colSel(getColSel(cols)){
	
		rows.reserve(1000);
		
		if (ignoreFirst){
		
			char* c = in.next_line();
			
			if (c == 0)
				open = false;
		}
		
	}

};


class RowTupleStream : public TupleStream{


private:
	
	//int pos = 0;
	
	int cols;
	const char* ptr;
	const char* end;
	
	Row r;	
	
public:	
	inline RowTupleStream(const Rows& rows) : cols(rows.getCols() ), ptr(rows.beginPtr()), end(rows.endPtr()){
		
		
	}
	
	
	inline bool read(Row& p) override{
	
		
		if (ptr == end)
			return false;
	
		ptr = r.set(cols, ptr, end);
		
		p = r;
		
		return true;
	}
	
	
	inline void close() override {
		
	}

};


class SortedSampleTupleStream : public TupleStream{

public:

	SortedSample samp;
	
	shared_ptr<TupleStream> stream;
	
	inline SortedSampleTupleStream( shared_ptr<TupleStream> s, long n, long u, const string& seed): stream(s), samp(n,u, seed){
	
	}

	inline bool read(Row& p) override{
	
		if (samp){
	
			if (stream->skip(samp.getSkip())){
		
				++samp;
		
				return stream->read(p);
			}
		};
		
		return false;
	}
	
	inline void close() override {
		
	}

	
};


// Translate strings to integers
class Translator{

public:
	const static long NOT_MATERIALISED = 1L << 60;
private:
	long hashuniv = UtilHash::HASHMAX;
	
	
	vector<char> keyToStringChars;
	vector<long> keyToStringPositions;
	
	Tuple keyToString;
	hash_map<string, long> stringToKey;
	
	Translator* parent = NULL;

	vector<Translator*> children;

	bool onlyIntegers = true;
	
	
	inline void materialise(){
		
		if (parent == NULL)
			return;
		
		const Translator* t = parent;
		
		while (t->parent != NULL)
			t = t->parent;
		
		assert (keyToString.size() == 0);
		assert (stringToKey.size() == 0);
		
		hashuniv = t->hashuniv;
		
		if (hashuniv == UtilHash::HASHMAX){
		
			keyToStringChars = t->keyToStringChars;
			keyToStringPositions = t->keyToStringPositions;
			stringToKey = t->stringToKey;
		}
		
		parent = NULL;
	}
	
	inline void removeChild(Translator* t){
	
		for (int i = 0; i < children.size(); ++i){
			
			if (children[i] == t){
				
				children[i] = children.back();
				children.pop_back();
			}
		}
	}
	
	inline void addChild(Translator* t){
	
		removeChild(t);
	
		children.push_back(t);
	}
	
	inline void setParent(Translator* t){
	
		assert (t != this);
	
		if (parent != t){
		
			if (parent != NULL)
				parent->removeChild(this);
			
			parent = t;
			parent->addChild(this);
		}
	}

public:

	inline Translator(long m = UtilHash::HASHMAX){
	
		hashuniv = m;
		
		
	}

	inline ~Translator(){
	
		if (parent != NULL){
		
			parent->removeChild(this);
		
			for (auto it = children.begin(); it != children.end(); ++it){
			
				(*it)->parent = parent;
				parent->addChild( (*it) );
			}	
			
			parent = NULL;
			
			
		} else {
		
			for (auto it = children.begin(); it != children.end(); ++it){
		
				(*it)->materialise();
				assert ( (*it)->parent == NULL);
			}
		}
		
		children.clear();
		
		
	}

	inline void set(Translator& t){
	
		keyToStringChars.clear();
		keyToStringPositions.clear();
		stringToKey.clear();
		
		setParent(&t);
		
		hashuniv = t.hashuniv;
		
		//materialise();
	}
	
	inline void print() const{
	
		cout << "print" << endl;
		if (parent != NULL){
		
			cout << "COPY ";
			parent->print();
			return;
		}
			
		for (int i = 0; i < keyToString.size(); ++i){
		
			cout << i << " = " << keyToString[i] << endl;
		}
		
		for (auto it = stringToKey.begin(); it != stringToKey.end(); ++it){
		
			cout << it->first << " = " << it->second << endl;
		}
		
		cout << "/print" << endl;
	}

	inline const string getKeyToString(const long& l) const{
	
		if (parent != NULL)
			return parent->getKeyToString(l);
		
		return string( keyToStringChars.data()+keyToStringPositions.at(l) );
	}

	// get original value for id
	inline const string getValue(long id) const{
	
		if (parent != NULL)
			return parent->getValue(id);
	
		if (id >= 0)
			return to_string(id);
	
		const long stringId = (-id)-1;
	
		if (stringId < keyToString.size() )
			return getKeyToString(stringId);
		else
			return "?"+to_string(id);
	}
	
	inline long getHashUniv() const{
	
		if (parent != NULL)
			return parent->getHashUniv();
	
		return hashuniv;
	}
	
	inline long getKeyConst(const string& ss) const{
	
		if (parent != NULL)
			return parent->getKeyConst(ss);
	
		const long e = UtilString::stringToNonNegativeInteger(ss);
		
		if (e >= 0)
			return hashuniv < UtilHash::HASHMAX ? UtilHash::hash(e, hashuniv): e;
		
		if (ONLYINTEGERS || onlyIntegers)
			return NOT_MATERIALISED;
		
		if (stringToKey.count(ss) == 0)
			return NOT_MATERIALISED;
	
		if (hashuniv == UtilHash::HASHMAX){
		
			if (stringToKey.count(ss) > 0)
				return stringToKey.at(ss);
				
			assert (false);
		}
		
		return UtilHash::hash(ss, hashuniv);
	}
	
	
	
	
	inline long getKey(const string& ss){
	
		if (parent != NULL){
		
			const long ret = parent->getKeyConst(ss);
		
			if (ret != NOT_MATERIALISED)
				return ret;
			
			parent->removeChild(this);
			materialise();
		}
	
		const long e = UtilString::stringToNonNegativeInteger(ss);
		
		if (e >= 0)
			return hashuniv < UtilHash::HASHMAX ? UtilHash::hash(e, hashuniv) : e;
	
		assert (!ONLYINTEGERS);
		onlyIntegers = false;
	
		if (stringToKey.count(ss) == 0){
			
			keyToStringPushBack(ss);
			stringToKey[ss] = -keyToStringPositions.size();
		}
		
		
		
		if (hashuniv == UtilHash::HASHMAX){
		
			if (stringToKey.count(ss) > 0)
				return stringToKey[ss];
				
			assert (false);
		}
		
		return UtilHash::hash(ss, hashuniv);
	}
	
	
	inline long getKeyConst(string& temp, const Row& r, int i) const{
	
		if (parent != NULL)
			return parent->getKeyConst(temp, r,i);
	
		{
			const long key = r.getLong(i);
		
			if (key >= 0)
				return hashuniv < UtilHash::HASHMAX ? UtilHash::hash(key, hashuniv) : key;
		}
		
		if (ONLYINTEGERS || onlyIntegers)
			return NOT_MATERIALISED;
			
		//const string& s = r[i];
		
		temp.clear();
		
		for (const char* c = r.get(i); (*c) != 0; ++c)
			temp.push_back((*c));
	
		if (hashuniv == UtilHash::HASHMAX){
		
			if (stringToKey.count(temp) > 0)
				return stringToKey.at(temp);
			
			return NOT_MATERIALISED;
				
			//assert (false);
		}
		
		return UtilHash::hash(temp, hashuniv);
	}
	
	
	inline void keyToStringPushBack(const string& s){
	
		if (keyToStringChars.capacity()-keyToStringChars.size() < (s.length()+2) ){
		
			cout << "keyToStringReserve " << keyToStringChars.size()*2 << endl;
			keyToStringChars.reserve(keyToStringChars.size()*2);
		}
	
		if (keyToStringPositions.capacity() == keyToStringPositions.size() )
			keyToStringPositions.reserve(keyToStringPositions.size()*2);
	
		keyToStringPositions.push_back(keyToStringChars.size() );
	
		for (int i = 0; i < s.length(); ++i)
			keyToStringChars.push_back(s[i]);
		
		keyToStringChars.push_back(0);		
	}
	
	inline long getKey(string& temp, const Row& r, int i){
	
		if (parent != NULL){
		
			const long ret = parent->getKeyConst(temp, r,i);
		
			if (ret != NOT_MATERIALISED)
				return ret;
			
			parent->removeChild(this);
			materialise();
			
		}
	
		const long key = r.getLong(i);
		
		if (key >= 0)
			return hashuniv < UtilHash::HASHMAX ? UtilHash::hash(key, hashuniv) : key;
	
		assert (!ONLYINTEGERS);
		onlyIntegers = false;
	
		//const string& s = r[i];
		
		temp.clear();
		
		for (const char* c = r.get(i); (*c) != 0; ++c)
			temp.push_back((*c));
	
		if (stringToKey.count(temp) == 0){
			
			keyToStringPushBack(temp);
			stringToKey[temp] = -keyToStringPositions.size();
		}
	
		if (hashuniv == UtilHash::HASHMAX){
		
			if (stringToKey.count(temp) > 0)
				return stringToKey[temp];
				
			assert (false);
		}
		
		return UtilHash::hash(r[i], hashuniv);
	}
	
	/*
	inline long getKey(Translator* old, long key) const{
	
		if (parent != NULL)
			return parent->getKey(old, key);
	
		if (old->getHashUniv() == UtilHash::HASHMAX){
		
			if (hashuniv == UtilHash::HASHMAX)
				return key;
				
			if (key >= 0)
				return UtilHash::hash(key, hashuniv);
				
			return UtilHash::hash( old->getValue(key), hashuniv);
		}
		
		assert (hashuniv <= old->getHashUniv());
		
		return UtilHash::hash(key, hashuniv);
	}*/
	
	// keys from this translator are equal to keys from other translator
	inline bool trivialConversion(const Translator* old){
	
		if ( (old->getHashUniv() != UtilHash::HASHMAX) && (old->getHashUniv() == getHashUniv()) )
			return true;
	
		if (old == this)
			return true;
	
		if (old == parent)
			return true;
	
		if (old->parent == this)
			return true;
			
		if ( (parent != NULL) && (old->parent == this->parent) )
			return true;
		
		
		return UtilTuple::equals(keyToString, old->keyToString);
	}
	
	inline long getKeyConst(const Translator* old, long key) const{
	
		if (parent != NULL)
			return parent->getKey(old, key);
	
		if (old->getHashUniv() == UtilHash::HASHMAX){
		
			/*
			if (hashuniv == UtilHash::HASHMAX)
				return key;*/
			
			if (key >= 0)
				return hashuniv == UtilHash::HASHMAX ? key : UtilHash::hash(key, hashuniv);
			
			const string s = old->getValue(key);
			
			return hashuniv == UtilHash::HASHMAX ? getKeyConst(s) : UtilHash::hash(s, hashuniv);
			
		}
		
		assert (hashuniv <= old->getHashUniv());
		
		return UtilHash::hash(key, hashuniv);
	}
	
	inline long getKey(const Translator* old, long key){
	
		if (parent != NULL)
			return parent->getKey(old, key);
	
		if (old->getHashUniv() == UtilHash::HASHMAX){
			
			if (key >= 0)
				return hashuniv == UtilHash::HASHMAX ? key : UtilHash::hash(key, hashuniv);
			
			const string s = old->getValue(key);
			
			return hashuniv == UtilHash::HASHMAX ? getKey(s) : UtilHash::hash(s, hashuniv);
		}
		
		assert (hashuniv <= old->getHashUniv());
		
		return UtilHash::hash(key, hashuniv);
	}
};


class TableHeader{

	vector<Translator> cols;
	
public:
	
	inline void print(){
	
		cout << "table header " << cols.size() << endl;
	}
	
	inline Translator* getTranslator(int col) {
	
		while (cols.size() <= col){
	
			cols.push_back(Translator() );
		}
	
		return &cols.at(col);
	}
	
	inline const Translator* getTranslatorConst(int col) const {
		
		assert (cols.size() > col);
	
		return &cols.at(col);
	}
	
	inline void setTranslator(int col, Translator& t){
	
		while (cols.size() <= col){
	
			cols.push_back(Translator() );
		}
		
		cols.at(col).set(t);
	}
	
	inline long getKey(Translator* t, int col, long k){
		
		while (cols.size() <= col){
	
			cols.push_back(Translator() );
		}

		return cols.at(col).getKey(t,k);
	}
	
	inline long getKey(int col, const string& s){
		
		while (cols.size() <= col){
	
			cols.push_back(Translator() );
		}
		
		return cols.at(col).getKey(s);
	}
	
	inline long getKey(string& temp, int col, const Row& r, int i){
		
		while (cols.size() <= col){
	
			cols.push_back(Translator() );
		}
		
		return cols.at(col).getKey(temp, r,i);
	}
	
	inline long getKeyConst(string& temp, int col, const Row& r, int i) const{
		
		if (col >= cols.size())
			return Translator::NOT_MATERIALISED;
		
		return cols.at(col).getKeyConst(temp, r,i);
	}
	
	
	inline int getColNum() const {
	
		return cols.size();
	}

};

template<typename Weight>
class RowHistogramNode{

	
	int col = -1;
	
	long_map<long, Weight> counts;
	
	
	long_map<long, unique_ptr<RowHistogramNode<Weight>> > children;
	
	RowHistogramNode<Weight>* parent = NULL;
	
	long size = 0;
	
	string temp;
	
	unique_ptr<std::bitset<100000>> bitset;
	
	bool locked = false;
	
	
public:

	inline void lock(){
	
		locked = true;
		
		assert (children.size() == 0);
	}


	inline Weight getWeightSum() const{
	
		Weight ret = 0;
		
		for (auto it = counts.cbegin(); it != counts.cend(); ++it){
		
			ret += it->second;
		}
		
		for (auto it = children.cbegin(); it != children.cend(); ++it){
		
			ret += it->second->getWeightSum();
		}
		
		return ret;
	}

	TableHeader* tab = NULL;
	
	inline void reserve(long m){
	
		if (locked)
			return;
	
		assert (col != -1);
	
		if (children.size() > 0){
		
			counts.reserve(m);
			
			
			
			//counts.max_load_factor(10.0);
			
			children.reserve( pow(m, 1.0/(col+1) ) );
		} else {
			
			if (false)
			if (m >= 100000){
			
				//cout << "bitset created" << endl;
				
				if (getSize() == 0){
					
					bitset.reset( new std::bitset<100000>() );
					bitset->reset();
				}
			}
			
			counts.reserve( pow(m, 1.0/(col+1) ) );
			
		}	
	}
	
	class iterator{
	
	public:
	
		typedef long_map<long, Weight> IT0;
		typedef typename IT0::const_iterator IT1;
		
		typedef long_map<long, unique_ptr<RowHistogramNode<Weight>> > IT2;
		typedef typename IT2::const_iterator IT3;
	
		vector< IT1 > countit;
		vector< IT1 > countend;
		vector< IT3 > childrenit;
		vector< IT3 > childrenend;
		
		vector<iterator> childit;
		
		bool open = false;
		
		const RowHistogramNode<Weight>* rows;
		
		
		
		inline iterator(){
		
			
		}
		
		inline iterator(const RowHistogramNode<Weight>* rows) : rows(rows){
		
			if (rows == NULL){
			
				open = false;
				return;
			}
		
			if (rows->counts.size() > 0){
			
				countit.push_back( rows->counts.cbegin() );
				countend.push_back( rows->counts.cend() );
				open = true;
			}
			
			if (rows->children.size() > 0){
			
				childrenit.push_back( rows->children.cbegin() );
				childrenend.push_back( rows->children.cend() );
				
				childit.push_back( iterator( childrenit[0]->second.get() ) );
				
				open = childit[0].open;
			}
		}
		
		
		inline Translator* getTranslator(int col) const{
		
			if (col == 0)
				return rows->tab->getTranslator(rows->col);
				
			assert (childit.size() > 0);
			
			return childit[0].getTranslator(col-1);
		}
		
		inline void get(long m, Tuple& t, Weight& f){
		
			t.clear();
		
			for (int i = 0; i < m; i++)
				t.push_back( to_string( getKey(i) ) );
			
			f = getCount();
		}
		
		
		inline long getKey(int col) const {
		
			if (col == 0)
				return countit.size() > 0 ? countit[0]->first : childrenit.size() > 0 ? childrenit[0]->first : 12345;
		
			assert (childit.size() > 0);
		
			return childit[0].getKey(col-1);
		}
		
		inline Weight getCount() const{
		
			if (countit.size() > 0)
				return countit[0]->second;
			
			assert (childit.size() > 0);
		
			return childit[0].getCount();
		}
		
		inline void print(){
			
			if (countit.size() > 0){
			
				cout << countit[0]->first << " " << countit[0]->second;
			}
			
			if (childit.size() > 0){
			
				cout << childrenit[0]->first << " ";
			
				childit[0].print();
			}
		}
		
		inline int operator()() const{
		
			return 123;
		}
		
		
		inline bool operator!=(int x) const{
			
			assert (x == 123);
			
			return open;
		}
		
		inline bool operator==(const bool& b) const{
		
			return !( (*this) != b);
		}
		
		
		inline iterator& operator++(int n){
		
			return ++(*this);
		}
		
		inline iterator& operator++(){
		
			if (countit.size() > 0){
			
				++countit[0];
				
				open = countit[0] != countend[0];
				
				return *this;
			}
		
			if (childit.size() > 0){
			
				++childit[0];
				
				open = childit[0].open;
				
				if (open)
					return *this;
			}
			
			if (childrenit.size() > 0){
			
				++childrenit[0];
				
				open = false;
				
				if (childrenit[0] != childrenend[0]){
					childit[0] = iterator( childrenit[0]->second.get() );
					
					open = childit[0].open;
				}
				
				return *this;
			}
			
			return *this;
		}
		
		
		
	};
	
	inline iterator begin() const{
	
		return iterator(this);
	}
	
	const iterator end;
	
	
	
	inline RowHistogramNode(RowHistogramNode<Weight>* ppp, TableHeader* tab, int col) : col(col){
	
		this->tab = tab;
		setParent(ppp);
		
		temp.reserve(1000);
	}
	
	
	inline void setParent(RowHistogramNode<Weight>* ppp){
	
		assert (!locked);
	
		this->parent = ppp;
		
		if (this->parent){
			this->parent->sizeChange(getSize() );
			this->tab = this->parent->tab;
		}
	}
	
	inline int getColNum() const {
	
		return col;
	}
	
	
	inline Translator* getTranslator(){
	
		return tab->getTranslator(col);
	}
	
	inline const Translator* getTranslatorConst() const{
	
		return tab->getTranslator(col);
	}
	
	inline void merge(const RowHistogramNode<Weight>* rows, bool noAdding=false){
	
		assert (!locked);
	
		//const Weight oldWeightSum = getWeightSum();
	
		if (rows->col != col){
		
			cout << rows->col << " != " << col << endl;
			
			
			cout << endl << "##1" << endl;
			
			auto it = rows->begin();
			
			for (int k = 0; k < 5; k++){
			
				if (it != rows->end() ){
				
					cout << "*" << endl;
					it.print();
					++it;
					
				} else {
				
					break;
				}
			}
			cout << endl << "##2" << endl;	
			
			cout << endl << endl; 
			
			auto it2 = begin();
			
			for (int k = 0; k < 5; k++){
			
				if (it2 != end() ){
				
					cout << "*" << endl;
					it2.print();
					++it2;
					
				} else {
				
					break;
				}
			}
			
			cout << endl << "##3" << endl;
		}
	
		assert (rows->col == col);
		
		
		if (rows->counts.size() > 0){
		
			const long oldsize = counts.size();
			
			Translator* tr1 = getTranslator();
			const Translator* tr2 = rows->getTranslatorConst();
			
			const bool trivial = tr1->trivialConversion(tr2);
			const bool noadd = noAdding;
			
			for (auto it = rows->counts.cbegin(); it != rows->counts.cend(); ++it){
					
				const long newkey = trivial ? it->first : tr1->getKey(tr2, it->first);
				
				if (noadd)
					counts[newkey] = it->second;
				else
					counts[newkey] += it->second;
				
				if (bitset)
					bitset->set(newkey % bitset->size() );
			}
			
			sizeChange( counts.size()-oldsize);
		}
		
		if (rows->children.size() > 0){
			for (auto it = rows->children.cbegin(); it != rows->children.cend(); ++it){
			
				if (children.count(it->first) == 0)
					children[it->first] = make_unique<RowHistogramNode<Weight>>(this, tab, col+1);
				
				children[it->first]->merge(it->second.get(), noAdding );
			}
		}	
	
		//if (!noAdding){
		//	assert (getWeightSum() == (oldWeightSum+rows->getWeightSum() ) );
		//}
	}
	
	inline void switchTranslator(int col1, Translator* tr){
	
		assert (!locked);
	
		//const Weight weightBefore = getWeightSum();
	
		if (col < col1){
		
			if (children.size() > 0){
		
				vector<long> keys;
				
				keys.reserve(children.size() );
		
				for (auto it = children.cbegin(); it != children.cend(); ++it)
					keys.push_back(it->first);
				
				for (long j = 0; j < keys.size(); j++)
					children[keys[j]]->switchTranslator(col1, tr);
					
				return;
			}
			
			//assert (weightBefore == getWeightSum() );
			
			return;
			//assert (false);
		}
		
		assert (col == col1);
		
		Translator* oldtr = tab->getTranslator(col);
		
		assert (oldtr != NULL);
		
		if (counts.size() > 0){
		
			vector<long> keys;
			vector<Weight> vals;
			
			keys.reserve(counts.size() );
			vals.reserve(counts.size() );
		
			for (auto it = counts.cbegin(); it != counts.cend(); ++it){
		
				keys.push_back( tr->getKey(oldtr, it->first) );
				vals.push_back( it->second);
			}
			
			const long oldsize = counts.size();
			
			counts.clear();
			
			if (bitset)
				bitset->reset();
			
			for (int j = 0; j < keys.size(); j++){
				counts[keys[j]] += vals[j];
				
				if (bitset)
					bitset->set(keys[j] % bitset->size());
			}
				
			sizeChange(counts.size()-oldsize);
		}
		
		if (children.size() > 0){
		
			vector<long> keys;
			
			keys.reserve(children.size() );
			
			for (auto it = children.cbegin(); it != children.cend(); ++it)
				keys.push_back(it->first );
			
			vector<unique_ptr<RowHistogramNode<Weight>>> vals;
			
			for (long j = 0; j < keys.size(); j++){
	
				vals.push_back( std::move(children[keys[j]]) );
				sizeChange( -vals.back()->getSize() );
			}
			
			children.clear();
			
			for (int j = 0; j < keys.size(); j++){
				
				const long key = tr->getKey(oldtr, keys.at(j) );
				
				if (children.count(key) > 0){
				
					children[key]->merge(vals.at(j).get() );
					
				} else {
					
					vals.at(j)->setParent(this);
					children[key] = std::move(vals.at(j));
					
				}
			}
		}
		
		//assert (weightBefore == getWeightSum() );
	
	}
	
	inline RowHistogramNode<Weight>* getCreateAlongPath(long key){
	
		if (children.count(key) == 0){
		
			assert (!locked);
		
			children[key] = make_unique<RowHistogramNode<Weight>>(this, tab, col+1);
		}
		
		return children[key].get();
	}
	
	inline RowHistogramNode<Weight>* get(long key) const{
	
		if (children.count(key) > 0)
			return children.at(key).get();
		else
			return NULL;
	}
	
	inline const RowHistogramNode<Weight>* getConst(long key) const{
	
		if (children.count(key) > 0)
			return children.at(key).get();
		else
			return NULL;
	}
	
	inline Weight getCount(long key) const{
	
		if (counts.count(key) == 0)
			return 0;
	
		return counts.at(key);
	}

	inline void sizeChange(long k){
	
		if (k != 0)
			assert (!locked);
	
		size += k;
		
		if (parent)
			parent->sizeChange(k);
	}

	inline void addCount(long key, Weight f){
		
		if (f != 0)
			assert (!locked);
		
		long oldsize = counts.size();
		counts[key] += f;
		
		if (bitset)
			bitset->set(key % bitset->size());
		
		sizeChange(counts.size()-oldsize);
	}
	
	inline RowHistogramNode<Weight>* getCreateAlongPath(Translator* t, long s){
	
		const long k = tab->getKey(t,col,s);
	
		return getCreateAlongPath( k );
	}
	
	inline void addCount(Translator* t, long s, Weight f){
	
		const long k = tab->getKey(t,col,s);
	
		addCount( k, f );
	}
	
	inline RowHistogramNode<Weight>* get(const Row& r, int i){
	
		const long k = tab->getKey(temp, col,r, i);
	
		return get( k );
	}
	
	
	inline RowHistogramNode<Weight>* getConst(string& temp, const Row& r, int i) const{
	
		const long k = tab->getKey(temp, col,r, i);
	
		return get( k );
	}
	
	inline void addCount(const Row& r, int i, Weight f){
	
		const long k = tab->getKey(temp, col,r,i);
		
		//cout << "addCount " << r << " i " << i << " f " << f << endl;
	
		addCount(k, f );
	}
	
	inline long countInternal(long n) const{
	
		if (counts.size() > 0)
			return counts.count(n);
			
		if (children.size() > 0)
			return children.count(n);
		
		return 0;
	}
	inline bool containsConst(string& temp, const Row& r, int i) const{
	
		//cout << "col " << col << " s " << s << endl;
	
		//tab->print();
		
		//cout << "/col " << col << " s " << s << endl;
		
		assert (tab != NULL);
		
		const long k = tab->getKeyConst(temp, col,r, i);
		
		if (k == Translator::NOT_MATERIALISED)
			return false;
			
		if (bitset){
			if ( !bitset->test(k % bitset->size()) )
				return false;
		}
		
		//cout << "k " << k << endl;
	
		const bool ret = countInternal(k) > 0;
		
		//cout << "ret " << ret << endl;
		
		return ret;
	}
	
	inline bool contains(const Row& r, int i){
	
		//cout << "col " << col << " s " << s << endl;
	
		//tab->print();
		
		//cout << "/col " << col << " s " << s << endl;
		
		assert (tab != NULL);
		
		const long k = tab->getKeyConst(temp, col,r, i);
		
		if (k == Translator::NOT_MATERIALISED)
			return false;
			
		if (bitset){
			if ( !bitset->test(k % bitset->size()) )
				return false;
		}
		
		//cout << "k " << k << endl;
	
		const bool ret = countInternal(k) > 0;
		
		//cout << "ret " << ret << endl;
		
		return ret;
	}
	
	inline Weight getCount(const Row& r, int i){
	
		const long k = tab->getKeyConst(temp, col,r,i);
		
		if (k == Translator::NOT_MATERIALISED)
			return 0;
			
		if (bitset){
			if ( !bitset->test(k % bitset->size()) )
				return 0;
		}
	
		return getCount(k);
	}
	
	
	inline Weight getCountConst(string& temp, const Row& r, int i) const{
	
		const long k = tab->getKeyConst(temp, col,r,i);
		
		if (k == Translator::NOT_MATERIALISED)
			return 0;
			
		if (bitset){
			if ( !bitset->test(k % bitset->size()) )
				return 0;
		}
		
		return getCount(k);
	}
	
	inline RowHistogramNode<Weight>* getCreateAlongPath(const Row& r, int i){
	
		const long k = tab->getKey(temp, col,r,i);
		
		return getCreateAlongPath( k);
	}
	
	inline long getSize() const{
	
		return size;
	}
};

template<typename Weight>
class RowHistogram{

	RowHistogramNode<Weight> rows;
	const Tuple cols;
	
	TableHeader tab;
	
	bool locked = false;
	
public:	

	inline void lock(){
	
		locked = true;
		rows.lock();
	}
	
	inline void print(int k = -1) const {
	
		/*
		for (int i = 0; i < cols.size(); i++){
		
			cout << "translator " << cols[i] << endl;
			tab.getTranslator(i)->print();
			cout << "/translator " << cols[i] << endl;
		}*/
	
		long size = 0;
		
		for (auto it = rows.begin(); it != rows.end(); ++it){
		
			if (k != 0){

				cout << "pos=" << size << " ";

				for (int i = 0; i < cols.size(); i++){
			
					const int key = it.getKey(i);
					//
					cout <<cols[i] << "=" <<  tab.getTranslatorConst(i)->getValue(key) << "\t";
				}
				
				cout << "#" << it.getCount() << endl;
			
				--k;
			}
			
			size++;
		}
		
		cout << "total " << size << endl;
	}
	
	inline void check() const {
	
		long size1 = 0;
		for (auto it = rows.begin(); it != rows.end(); ++it){
		
			size1++;
		}
		
		long size2 = rows.getSize();
		
		//assert ( (size1 > 0) == (size2 > 0) );
		assert (size1 == size2);
		
		//cout << "CHECK COMPLETED " << size1 << " == " << size2 << endl;
	}
	
	inline const Tuple& getColNames() const{
	
		return cols;
	}
	
	
	inline RowHistogram(const Tuple& cols) : rows(NULL, &tab, 0), cols(cols){
	
		
	}
	
	inline void reserve(long m){
	
		
	
		rows.reserve(m);
	}
	
	
	inline bool containsConst(string& temp, const Row& t, const vector<int>& v) const{
	
		assert (v.size() > 0);
		
		//cout << "map" << endl;
		
		const RowHistogramNode<Weight>* cr = getConst(temp, t, v, false, true);
		
		if (cr == NULL)
			return false;
		
		//cout << "str " << t<< " " << v.back() << endl;
		
		//const string& s = t[v.back()];
		
		//cout << "ret " << s << " " << cr->getSize() << endl;
		
		const bool ret = cr->containsConst(temp, t, v.back());
		
		//cout << "/ret " << ret << endl;
		
		return ret;
	}
	
	
	inline bool contains(const Row& t, const vector<int>& v){
	
		assert (v.size() > 0);
		
		//cout << "map" << endl;
		
		RowHistogramNode<Weight>* cr = get(t, v, false, true);
		
		if (cr == NULL)
			return false;
		
		//cout << "str " << t<< " " << v.back() << endl;
		
		//const string& s = t[v.back()];
		
		//cout << "ret " << s << " " << cr->getSize() << endl;
		
		const bool ret = cr->contains(t, v.back());
		
		//cout << "/ret " << ret << endl;
		
		return ret;
	}
	
	
	inline Weight getCount(const Row& t, const vector<int>& v) {
	
		assert (v.size() > 0);
		
		return get(t, v, true)->getCount(t, v.back());
	}
	
	
	inline Weight getCountConst(string& temp, const Row& t, const vector<int>& v) const {
	
		assert (v.size() > 0);
		
		return getConst(temp, t, v, false, true)->getCountConst(temp, t, v.back());
	}
	
	inline void addCount(const Row& t, const vector<int>& v, Weight f){
	
		assert (v.size() > 0);
		
		if (f > 0){
		
			assert (!locked);
		
			get(t, v, true, true)->addCount(t, v.back(), f);
		}
	}
	
	inline void addCount(const Row& t, Weight f){
	
		assert (!locked);
	
		vector<int> v;
		
		for (int i = 0; i < t.size(); i++)
			v.push_back(i);
			
		if (f > 0){
			
			assert (!locked);
			addCount(t, v, f);
		}
		
		
	}
	
	
	inline void push_back(const RowHistogram<Weight>& t){
	
		assert (!locked);
	
		assert(UtilTuple::equals(cols, t.cols) );
		
		
		
		rows.merge(&t.rows); // why no adding, true);
	}
	
	inline void push_back_overwriting(const RowHistogram<Weight>& t){
	
		assert (!locked);
	
		assert(UtilTuple::equals(cols, t.cols) );
		
		rows.merge(&t.rows, true);
	}
	
	inline typename RowHistogramNode<Weight>::iterator begin() const{
	
		return rows.begin();
	}
	
	const typename RowHistogramNode<Weight>::iterator end;
	
	
	inline RowHistogramNode<Weight>* get(){
	
		return &rows;
	}
	
	inline const RowHistogramNode<Weight>* getConst() const{
	
		return &rows;
	}
	
	inline const RowHistogramNode<Weight>* getConst(string& temp, const Row& t, const vector<int>& v,bool createAlongPath = false, bool skipLast = false) const {
		
		assert (createAlongPath == false);
		
		const RowHistogramNode<Weight>* ptr = getConst();
		
		const long last = skipLast ? v.size()-2 : v.size()-1;
		
		for (int col = 0; col <= last; col++){
		
			//const string& s = t[v[col]];
		
			if (ptr){
				
				ptr = ptr->getConst(temp, t, v[col]);
			}
		}
		
		return ptr;
	}
	
	inline RowHistogramNode<Weight>* get(const Row& t, const vector<int>& v, bool createAlongPath = false, bool skipLast = false){
		
		RowHistogramNode<Weight>* ptr = get();
		
		const long last = skipLast ? v.size()-2 : v.size()-1;
		
		for (int col = 0; col <= last; col++){
		
			//const string& s = t[v[col]];
		
			if (ptr){
				
				if (createAlongPath)
					ptr = ptr->getCreateAlongPath(t, v[col]);
				else
					ptr = ptr->get(t, v[col]);
			}
		}
		
		return ptr;
	}
	
	
	inline int getColNum(){
	
		return tab.getColNum();
	}
	
	
	inline long getHashing(int col){
	
		return tab.getTranslator(col)->getHashUniv();
	}
	
	inline bool update(unordered_map<string, long>& map){
	
		bool change = false;
	
		for (int i = 0; i < getColNum(); i++){
		
			const string& s = cols[i];
		
			const long global = map.count(s) == 0 ? UtilHash::HASHMAX : map[s];
			
			const long local = getHashing(i);
			
			//bool change1 = false;
			
			if (local < global){
			
				map[s] = local;
				
				cout << "GLOBAL_SET HASHING " << cols[i] << " FROM " << global << " TO " << local << endl;
				
				change = true;
				
			} else if ( global < local){
				
				setHashing(i, global);
			}
		}
		
		if (locked)
			assert (!change);
		
		if (change)
			cout << "UPDATE " << change << endl;
		
		return change;
	}
	
	inline bool usesHashing(){
	
		for (int i = 0; i < getColNum(); i++)
			if (getHashing(i) < UtilHash::HASHMAX )
				return true;
		
		return false;
	}
	
	
	inline void setHashing(int col, long m){
		
		if (locked)
			return;
		
		const long old = getHashing(col);
	
		if (old == m)
			return;
		
		
		assert (!locked);
		assert (m <= old);
		
		cout << "SET HASHING " << col << " TO " << m << endl;
			
		Translator t(m);
			
		get()->switchTranslator(col, &t);
			
		tab.setTranslator(col, t);
	}
	
	
	inline long getSize() const{
	
		return rows.getSize();
	}
	
	inline bool shrinkToSize(const long max){
	
		bool change = false;
		int col = getColNum()-1;
		
		const long min = exp( log(max)/getColNum() );
		
		//cout << "old size " << getSize() << " max " << max << " min " << min<< endl;
		
		int noupdate = 0;
		
		for (int i = 0; (i < 100) && (getSize() > max); i++){
			
			const long m = getHashing(col);
			
			if (m <= min){
			 
				noupdate++;
				
			} else {
			
				const long newm = m > 0 ? UtilMath::maxVal<long>(min, m/4) : min;
			
				cout << "SET HASHING BECAUSE OF OVERLOW ALONG " << getColNames()[col] << " TO " << newm << endl;
				
				setHashing(col, newm );
				
				noupdate = 0;
				change = true;	
			}
			
			if (noupdate > getColNum())
				break;
			
			col = (col-1) % getColNum();
		}
		
		if (locked)
			assert (!change);
		
		return change;
		//cout << "new size " << getSize() << endl;
	}
	
	
	
};


// 
template <typename E, typename F, typename FF>
class MetaIterator{

	public:
	
		vector<F> it;
		vector<MetaIterator> mit;
		
		vector<FF> itend;
		
		vector<E>* sources;
		int m;
		
		bool open = false;
		
		bool skipNull;
	
		int srcSize;
	
		inline MetaIterator(vector<E>& src, bool skipNull = true, int mm = 0) : skipNull(skipNull){
			
			if (skipNull){
				
				while (src[mm] == NULL){
			
					++mm;
				}
			}
			
			srcSize = src.size();
			
			while (srcSize > 0 && (src.at(srcSize-1) == NULL) )
				--srcSize;
			
			m = mm;
			sources = &src;
			open = false;
			
			if (m >= srcSize )
				return;
				
			if (sources->at(m) == NULL){
			
				open = false;
			
				if ( (m+1) < srcSize ){
				
					mit.clear();
					mit.push_back( MetaIterator(*sources, skipNull, m+1) );
			
					open = mit[0].open;
				} 
			}
			
			if (sources->at(m) != NULL){
			
				it.push_back( sources->at(m)->begin() );	
				itend.push_back( sources->at(m)->end() );
				
				if (it[0] != itend[0]){
				
					open = true;
			
					if ( (m+1) < srcSize ){
					
						mit.clear();
						mit.push_back( MetaIterator(*sources, skipNull, m+1) );
				
						open = mit[0].open;
					} 
				}			
			}
			
			
		}
		
		
		inline bool cont(){
		
			return open;
		}
		
		inline F* get(int x){
		
			assert (open);
		
			if (x == 0)
				return &(it[0]);
			
			if (mit.size() > 0)
				return mit[0].get(x-1);
			
			assert (false);
		}
		
		inline F& getIterator(int x){
		
			return *get(x);
		}
		
		inline MetaIterator& operator++(int n){
		
			return ++(*(this));
		}
		
		inline MetaIterator& operator++(){
		
			open = false;
		
			if (mit.size() > 0){
			
				mit[0]++;
				
				if (mit[0].open){
				
					open = true;
					return (*this);
				}
			}
			
			if (it.size() > 0){
			
				it[0]++;
			
				if (it[0] != itend[0]){
			
					open = true;
			
					if ( (m+1) < srcSize ){
				
						mit[0] = MetaIterator(*sources, skipNull, m+1);
						open = mit[0].open;
					}
				
					return (*this);
				}
			}
			
			return (*this);
		}

};


typedef long_double Weight;

class Joiner{

	private:
		shared_ptr<RowHistogram<Weight>> target;

		vector<const RowHistogram<Weight>*> children;

		shared_ptr<JoinColumns> rel;

	public:
		
		
		
		vector<int> semijoins;
		vector<int> fks;
		
		string temp;
		
		bool fk;
		
		inline Joiner(shared_ptr<RowHistogram<Weight>> target_, shared_ptr<JoinColumns> rel_, bool fk = false){
		
			this->fk = fk;
			rel = rel_;
			target = target_;
			
			if (target->update(*(rel->getHashings()) ))
				rel->updateHashings();
			
			//cout << rel->getNovelTuple() << " " << UtilString::toString(novelIds) << endl;
		}
		
		/*
		inline const vector<const RowHistogram<Weight>*>& getChildren() const{
		
			return children;
		}*/
		
		inline const vector<const RowHistogram<Weight>*> getChildren2() const{
		
			vector<const RowHistogram<Weight>*> ret;
			
			for (int j = 0; j < children.size(); ++j){
			
				//assert(semijoins.at(j) == 0);
			
				if (fks.at(j) == 1)
					ret.push_back(NULL);
				else
					ret.push_back(children.at(j) );
					
				
			}	
		
			return ret;
		}
		
		inline void addChild(RowHistogram<Weight>& c, bool semijoin=false, bool fk=false){
		
			const RowHistogram<Weight>& d = c;
		
			semijoins.push_back(semijoin ? 1 : 0);
			fks.push_back(fk ? 1 : 0);
		
			/*
			if (c.update(*(rel->getHashings()) ))
				rel->updateHashings();
			*/
			
			children.push_back( &d );
			
			
			const int childIndex = children.size()-1;
			
			
			if (rel->childrenAcyclic[childIndex])
				return;
			
			cout << "X20" << endl;
			
			const vector<int>& vtarget = rel->getMap(JoinColumn::CHILD+childIndex, JoinColumn::TARGET);
			const vector<int>& vchild = rel->getTargetInChild(childIndex);
			
			
			assert (vtarget.size() == vchild.size() );
			
			for (int ii = 0; ii < vtarget.size(); ++ii){
			
				const int i = vtarget[ii];
				const int j = vchild[ii];
				
				assert (UtilString::equals(rel->getTargetTuple()[i], rel->getChildTuple(childIndex)[j]) );
				
				Translator* tr1 = target->get()->tab->getTranslator(i);
				Translator* tr2 = c.getConst()->tab->getTranslator(j);
			
				tr1->set( (*tr2) );
			}
			
		}
		
		inline void addChild(shared_ptr<RowHistogram<Weight>> c, bool semijoin=false, bool fk=false){
		
			addChild( (*c), semijoin, fk );
		}
		
		
		inline bool joinEmpty(const Row& tmain){
		
			if (fk)
				return false;
		
			for (int j = 0; j < children.size(); j++){
				
				assert ( children.at(j) != NULL);
				
				//cout << rel->getChildTuple(j) << " " << rel->getMainTuple() << endl;
				//cout << UtilString::toString( rel->getChildJoinColsInMain(j) ) << endl;
				
				if (fks.at(j) == 1)
					continue;
				
				if (!children.at(j)->containsConst(temp, tmain, rel->getChildJoinColsInMain(j) ))
					return true;
			}
			
			return false;
		}
		
		inline Weight getJoinSize(const Row& tmain, bool novelweight=true){
		
			if (fk)
				return 1;
			
			//cout << ch.size() << " =?= " << children.size() << endl;
			//assert (ch.size() == children.size() );
		
			if (joinEmpty(tmain)){
			
				//cout << tmain << " " << 0 << endl;
			
				return 0;
			}
		
			const vector<int>& ch = rel->acyclicChildren;
			
			Weight ret = novelweight ? rel->getNovelWeightInMain(tmain) : 1;
			
			//
			
			for (int jj = 0; jj < ch.size(); jj++){
				
				const int j = ch.at(jj);
				
				if (fks.at(j) == 1)
					continue;
				
				const RowHistogram<Weight>* c = children.at(j);
				
				auto d = c->getCountConst(temp, tmain, rel->getChildJoinColsInMain(j) );
				
				if (semijoins.at(j) == 1){
					
					if (d == 0)
						return 0;
				
				} else {
				
					ret *= d;
				}
			}
			
			//cout << tmain << " " << ret << endl;
			
			return ret;
		}
		
		inline void print(){
		
		
		}
		
		inline bool shrinkToSize(long max){
		
			if (target->shrinkToSize(max))
				return true;
			
			return false;
		}
		
		inline void join(const Row& tmain){
			
			if (fk)
				return;
			
			const Weight w = getJoinSize(tmain);
			
			if (w == 0)
				return;
			
			if (rel->areAllChildrenAcyclic() ){
			
				RowHistogramNode<Weight>* t = target->get();
				
				//const vector<int>& target = rel->getTargetInMain();
				
				//assert( target.size() > 0);
				
				const int jmax = rel->getTargetTuple().size()-1;
				
				for (int j = 0; j <= jmax ; j++){
				
					const int k = rel->translate(JoinColumn::TARGET, j, JoinColumn::MAIN);
				
					//assert (UtilString::equals(rel->getTargetTuple()[j], rel->getMainTuple()[k]) );
				
					if (j < jmax)
						t = t->getCreateAlongPath(tmain, k);
					else
						t->addCount(tmain, k, w);
				}
				
				return;
			}
			
			
			if (rel->areAllChildrenSimple() ){
			
				vector<const RowHistogramNode<Weight>*> maps(children.size(), NULL );
			
				const vector<int>& ch = rel->cyclicChildren;
			
				for (int jj = 0; jj < ch.size(); jj++){
				
					const int j = ch.at(jj);
				
					maps.at(j) = children.at(j)->getConst(temp, tmain, rel->getChildJoinColsInMain(j) );
					
					//cout << "##" << j << endl;
				}
			
				const Tuple& ttarget = rel->getTargetTuple();
				
				
				
				for (auto mit = MetaIterator<const RowHistogramNode<Weight>*, typename RowHistogramNode<Weight>::iterator, int>(maps, false); mit.cont(); ++mit){
					
			
					Weight w2 = w;
				
					for (int k = 0; k < maps.size(); k++){
					
						if (maps[k])
							w2 *= mit.get(k)->getCount();
					}
				
					RowHistogramNode<Weight>* map = target->get();
				
					assert (map != NULL);
				
					const int last = ttarget.size()-1;
					
					assert (last == rel->getTargetTuple().size()-1);
				
					for (int k = 0; k <= last; k++){
				
						const vector<int>& sources = rel->getSources(JoinColumn::TARGET,k);
						
						int src = -1;
						
						for (int l = 0; l < sources.size(); l++){
						
							if (sources[l] == JoinColumn::MAIN){
							
								src = JoinColumn::MAIN;
								//col = sourcesCols[l];
								break;
							}
							
							if (sources[l] >= JoinColumn::CHILD){
							
								src = sources[l];
								//col = sourcesCols[l];
								break;
							}
						}
						
						if (src == -1){
						
							cout << rel->getTargetTuple()[k] << endl;
							
							cout << "ALL " << rel->getAllTuple() << endl;
							cout << "MAIN " << rel->getMainTuple() << endl;
							cout << "CHILD " << rel->getChildTuple(0) << endl;
							cout << "CHILD " << rel->getChildTuple(1) << endl;
							cout << UtilString::toString(sources) << endl;
							cout << "SRC " << src << endl;
							//cout << "COL " << col << endl;
						}
						
						const int col = rel->translate(JoinColumn::TARGET, k, src);
						
						//cout << "target " << rel->getTuple(JoinColumn::TARGET) << "[" << k << "]" << endl;
						//cout << "source " << rel->getTuple(src) << "[" << col << "]" << endl;
						
						//assert (UtilString::equals( rel->getTargetTuple()[k], rel->getTuple(src)[col]));
						
						assert (src != -1);
						assert (col != -1);
					
						if (src == JoinColumn::MAIN){
					
							if (k == last)
								map->addCount(tmain, col, w2);
							else
								map = map->getCreateAlongPath(tmain, col);
					
						} else {
					
							const int childIndex = src-JoinColumn::CHILD;
							
							const int col2 = col-rel->getChildJoinColsInMain(childIndex).size();
							
							//assert (UtilString::equals( rel->getTargetTuple()[k], rel->getTuple(src)[col]));
							
							auto m = mit.get(childIndex);
							
							Translator* translator = m->getTranslator(col2);
							
							assert (translator != NULL);
						
							const long key = m->getKey(col2);
							
							if (k == last)
								map->addCount(translator, key, w2);
							else
								map = map->getCreateAlongPath(translator, key);
						}
					
					}
				}
				
				return;
			}
			
			assert (false);
		}
};





class MainMemoryIndex{

	vector<Row> orderedRows;

	// starts and ends of blocks
	
	
	
	RowHashMap<std::pair<const Row*, const Row*>> blocks;
	
	vector<int> proj;
	vector<long> dummy;
	
	vector<int> sel;
	
	vector<long_double> weights;
	
	shared_ptr<Joiner> joiner;
	
	long maxBlockLength = 0;
	long blockNum = 0;
	
public:
	
	inline MainMemoryIndex(vector<int> projection) : proj(projection), dummy(projection.size(), UtilHash::HASHMAX ), sel(projection.size(), -1 ){
		
	}
	
	inline MainMemoryIndex(int projection) : MainMemoryIndex( std::vector<int>(1, projection) ){
		
		
	}
	
	inline void reserve(long m){
	
		orderedRows.reserve(m);
		
		
	}
	

	/*
	inline void startWeights(){
	
		weights.clear();
		weights.reserve(orderedRows.size() );
	}
	
	
	
	inline void addWeight(long_double weight){
	
		weights.push_back(weight);
	}
	
	inline void finaliseWeights(){
	
		
		for (long j = 1; j < weights.size(); ++j)
			weights.at(j) += weights.at(j)-1;
	}*/
	
	inline const Row* begin() const{
		
		return &(*orderedRows.begin());
	}
	
	inline const Row* end() const{
		
		return begin()+orderedRows.size();
	}
	
	inline const Row* get(long ind) const{
	
		assert (ind >= 0);
		assert (ind < orderedRows.size() );
	
		return begin()+ind;
	}
	
	inline long getMaxBlockLength() const {
	
		return maxBlockLength;	
	}
	
	inline void setJoiner(shared_ptr<Joiner> joiner){
	
		this->joiner = joiner;
	}
	
	inline bool hasJoiner() const{
	
		if (joiner)
			return true;
		else 
			return false;
	}
	
	inline const Row* get(const Row* st, const Row* en,  double randomDouble){
	
		assert (hasJoiner() );
		
		assert (st != NULL);
		assert (en != NULL);
	
		const long ind1 = st-begin();
		const long ind2 = en-begin();
		
		assert (ind1 != ind2);
		
		if (joiner){
		
			if (weights.size() == 0){
		
				weights.reserve(orderedRows.size() );
			
				for (long k = 0; k < orderedRows.size(); ++k)
					weights.push_back(-1);
			}
		
			auto beg1 = weights.begin()+ind1;
			auto end1 = weights.begin()+ind2;
		
			if ( (*beg1) == -1 ){
			
				const Row* it = st;

				for (long k = ind1; k < ind2; ++k){
					weights.at(k) = joiner->getJoinSize(*it);
					++it;
				}
				
				if (true){
					
					long sw = ind2-1;
				
					while (sw >= ind1 && weights.at(sw) == 0)
						--sw;
				
				
				
					// sw == last non-zero weight index
				
					Row copy;
				
					for (long k = sw; k >= ind1; --k){
				
						if (weights.at(k) == 0 && (k < sw)){
					
							// swap SW and K
							
							copy = orderedRows.at(sw);
							orderedRows.at(sw) = orderedRows.at(k);
							orderedRows.at(k) = copy;
						
							const long_double copy_w = weights.at(sw);
							weights.at(sw) = weights.at(k);
							weights.at(k) = copy_w;
							
							assert (weights.at(sw) == 0);
							assert (weights.at(k) > 0);
							
							while (sw >= ind1 && weights.at(sw) == 0)
								--sw;
						}
					}
					
					for (long k = sw; k >= ind1; --k){
				
						assert (weights.at(k) > 0);
					}
					
					for (long k = sw+1; k < ind2; ++k){
				
						assert (weights.at(k) == 0);
					}
					
					for (long k = ind1+1; k <= sw; ++k)
						weights.at(k) += weights.at(k-1);
				
					const long_double divmul = 1.0L / weights.at(sw);
			
					for (long k = ind1; k < ind2; ++k)
						weights.at(k) *= divmul;

					weights.at(sw) = 1.0;
					
					for (long k = sw+1; k < ind2; ++k)
						weights.at(k) = 1.0;
					
				} else {
				
					for (long k = ind1+1; k < ind2; ++k)
						weights.at(k) += weights.at(k-1);
					
					assert (it == en);
				
					const long_double divmul = 1.0L / weights.at(ind2-1);
			
					for (long k = ind1; k < ind2; ++k)
						weights.at(k) *= divmul;
				}
			}
			
			assert (randomDouble >= 0.0);
			assert (randomDouble <= 1.0);
			
			auto it = std::lower_bound(beg1, end1, randomDouble);
			
			const long ind = it == end1 ? ind2-1 : it-weights.begin();
			
			assert(ind >= ind1);
			assert(ind < ind2);
			
			if ( (weights.at(ind) == -1) ||((ind > ind1) && (weights.at(ind) <= weights.at(ind-1))) ){
				
				cout << "weights" << endl;
				for (auto it = beg1; it != end1; ++it){
				
					cout << (*it) << endl;
				}
				cout << "/weights" << endl;
			}
			
			assert (weights.at(ind) != -1);
			assert (weights.at(ind) > 0);
			
			if (ind > ind1)
				assert (weights.at(ind) > weights.at(ind-1));
			
			return get( ind );
			
		} else {
		
			return get( ind1+UtilMath::minVal<long>(ind2-1-ind1, floor(randomDouble * (ind2-ind1) )) );
		}
		
		
	}
	
	inline const Row* get(const string& s, double randomDouble){
		
		const Row* st = begin(s);
		const Row* en = end(s);
	
		return get(st, en, randomDouble);
	}
	
	inline const Row* get(const Row& row, int column, double randomDouble){
	
		assert (sel.size() == 1);
	
		sel[0] = column;
	
		const long key = blocks.getLongKey(row, sel, dummy);
		
		if (!blocks.hasKey(row, sel, dummy, key))
			return NULL;
	
		const std::pair<const Row*, const Row*>& p = blocks.get(row, sel, dummy, key);
	
		const Row* st = p.first;
		const Row* en = p.second;
		
		
		if (st == NULL || en == NULL){
		
			cout << row << "[" << column << "] = " << row[column] << endl;
		}
		
		return get(st, en, randomDouble);
	}
	
	inline const Row* begin(const string& s){
	
		return blocks.get(s).first;
	}
	
	
	inline const Row* end(const string& s) {
	
		return blocks.get(s).second;
	}
	
	inline const Row* begin(long s){
	
		return blocks.get(s).first;
	}
	
	inline const Row* end(long s){
	
		return blocks.get(s).second;
	}


	inline const Row* begin(const Row& row, const vector<int>& sel){
	
		if (row.isNull(sel) ||!blocks.hasKey(row, sel, dummy) )
			return NULL;
		
		return blocks.get(row, sel, dummy).first;
	}
	
	inline const Row* end(const Row& row, const vector<int>& sel){
	
		if (row.isNull(sel) ||!blocks.hasKey(row, sel, dummy) )
			return NULL;
		
		return blocks.get(row, sel, dummy).second;
	}

	inline const Row* begin(const Row& row, int j){
		
		assert (sel.size() == 1);
	
		sel[0] = j;
		
		return begin(row, sel);
	}
	
	inline const Row* end(const Row& row, int j){
	
		assert (sel.size() == 1);
			
		sel[0] = j;
		
		return end(row, sel);
	}
	
	
	inline const Row* begin(const Row& row){
		
		return begin(row, proj);
	}
	
	inline const Row* end(const Row& row){
	
		return end(row, proj);
	}
	
	inline void insert(const Row& row){
	
		if (orderedRows.capacity() == orderedRows.size() )
			orderedRows.reserve( orderedRows.size()*1.25);
	
		orderedRows.push_back(row);
	}
	
	inline void finalise(){
	
		for (int k = proj.size()-1; k >= 0; --k){
	
			for (auto it = orderedRows.begin(); it != orderedRows.end(); ++it)
				(*it).setSortKey(proj.at(k));
			
			std::stable_sort( orderedRows.begin(), orderedRows.end() );
		}
		
		const Row* first = NULL;
		const Row* last = NULL;
		
		//starts.reserve(orderedRows.size() );
		//ends.reserve(orderedRows.size() );
		
		BasicAggregator agg;
		
		agg.add(maxBlockLength);
		
		blocks.reserve(orderedRows.size());
		
		for (auto it = orderedRows.begin(); it != orderedRows.end(); ++it){
		
			const Row* row = &(*it);
			
			if (first != NULL){
			
				if (!row->equals(*first, proj) ){
				
					assert (last != NULL);
				
					assert ( first->equals(*last, proj));
				
					if ( !last->isNull(proj)){
						blocks.insert( *last, proj, dummy, std::make_pair(first, last+1) );
						++blockNum;
					}
				
					//ends.insert(*last, proj, dummy, last+1);
					
					
					agg.add( last+1-first );
					maxBlockLength = agg.max();
					
					//cout << "LAST " << (*last) << " ROW " << (*row) << " proj[0] " << proj[0] << endl;
					
					//starts.insert(*row, proj, dummy, row);
					
					first = row;
				}
				
			} else {
			
				first = row;
				//starts.insert(*row, proj, dummy, row);
			}
			
			last = row;
		}
		
		if (last != NULL){
			
			if ( !last->isNull(proj)){
				blocks.insert( *last, proj, dummy, std::make_pair(first, last+1) );
				++blockNum;
			}
			//ends.insert(*last, proj, dummy, last+1);
			
			
			agg.add( last+1-first );
			maxBlockLength = agg.max();
		}
		
		//for (auto it = orderedRows.begin(); it != orderedRows.end(); ++it)
		//	blocks.insert( *it, proj, dummy, std::make_pair(starts.get(*it, proj, dummy), ends.get(*it, proj, dummy) ) );
		
		assert (blocks.getSize() == blockNum);
	}
};


class Profiler{

	public:
		
		BasicAggregator distinctAgg;
		
		string s;
		
		std::uniform_real_distribution<double> distribution;
		std::hash<string> h;
		
		std::bitset<100000> bs;
		
		inline Profiler() : distribution(0.0, 1.0){
		
			s.reserve(100);
			s.clear();
		}
		
		bool open = false;
		
		inline long count() const{
		
			return distinctAgg.count();
		}
		
			
		inline long distinct() const{
		
			return count() == 0 ? 0 : UtilMath::maxVal<long>( distinctMin(), distinctAgg.max()/distinctAgg.min() );
		}
		
		inline long distinctMin() const{
		
			return bs.count();
		}
		
		inline void collect(const char* st, const char* en){
		
		
			for (const char* p = st; p != en; ++p)
				s.push_back(*p);
		}
		
		inline void add(){
			
			size_t hash = h(s);
			
			std::default_random_engine generator (hash);
			
			bs.set( hash % bs.size() );
			
			distinctAgg.add( distribution(generator) );
			
			s.clear();
		}
		
		inline void add(const char* st, const char* en){
		
			collect(st,en);
			add();
		}
		
		inline void collect(const string& st){
		
		
			collect(st.c_str(), st.c_str()+st.length() );
		}
		
		inline void add(const string& st){
		
			collect(st);
			add();
		}
};



class RowComparator : public std::binary_function<const Row, const Row, bool>{

public:

	const vector<int> lex;
	
	inline RowComparator(const vector<int>& v) : lex( UtilTuple::without(v, -1) ){
	
		
	}
	
	
	inline bool operator()(const Row& a, const Row& b) const{
		
		for (int i = 0; i < lex.size(); ++i){
		
			const int j = lex.at(i);
			
			assert (j >= 0 && j < a.size() );
			assert (j >= 0 && j < b.size() );
			
			const char* ptr1 = a.get(j);
			int len1 = a.getLength(j);
			const char* ptr2 = b.get(j);
			int len2 = b.getLength(j);
			
			const int c = a.stringcmp(ptr1,len1, ptr2, len2);
			
			if (c < 0)
				return true;
				
			if (c > 0)
				return false;
		
		}
		
		return false;
	}
	
};


/*
class JoinProc{

	int mode = 0; 
	// mode == 0 for join size
	// mode == 1 for ks-test
	// mode == 2 for sampling
	
	
	long joinsize = 0;
	long double joinweight = 0;
	
	// join row sample
	shared_ptr<Rows> joinRows;	
	shared_ptr<vector<Row>> joinRowsVector;

	Tuple joincols;
	
	vector<Tuple> basecols;
	
	// sorted base tables
	vector<shared_ptr<Rows>> baseRows;
	vector<shared_ptr<Rows::iterator>> baseIters;
	vector<shared_ptr<Weighter>> baseWeighters;
	
	WeightingFunction* f = NULL;
	
	// uniform variates for sample
	vector<double> coords;
	vector<double>::const_iterator coordsIter;

	vector<vector<int>> projsSel;
	
	// join columns for each base table
	vector<vector<int>> projsJoin;
	
	// non-join columns for each base table
	vector<vector<int>> projsNonJoin;
	
	
	vector<RowComparator> cmps;

	inline const Rows& getJoinRows() const {
	
		return *joinRows;
	}
	
	inline void setWeightingFunction(WeightingFunction& ff){
	
		f = &ff;
	}

	inline void addSample(const vector<double>& smp){
		
		assert (mode == 0);	
		mode = 2;
		coords = smp;
	}
	
	inline void init(const Tuple& a, const Tuple& b){
		
		init(vector<Tuple>{Tuple(a), Tuple(b)});
	}
	
	inline void init(const vector<Tuple>& tables){
	
		set<string> setJoinCols;
		
		for (int i = 0; i < tables.size(); ++i){
		
			const Tuple& ti = tables.at(i);
			
			for (int j = 0; j < tables.size(); ++j){
			
				const Tuple& t = tables.at(j);
			
				for (int k = 0; k < ti.size(); ++k){
				
					if ( UtilTuple::find(tj, ti[j]) != -1){
					
					
						if (setJoinCols.count(ti[j]) == 0){
						
							assert (i == 0);	
							setJoinCols.insert(ti[j]);
						}
					}
				}
			}
		}
		
		Tuple joinCols = tables.at(0);
		
		
		for (int i = 0; i < tables.size(); ++i){
		
			const Tuple& ti = tables.at(i);
			
			vector<int> projSel;
			vector<int> projJoin;
			vector<int> projNonJoin;
			
			for (int k = 0; k < ti.size(); ++k){
			
				if (i == 0 || setJoinCols.count(ti[k]) == 0){
				
					projSel.push_back(k);	
					joinCols.push_back(ti[k]);
				}
					
				if (setJoinCols.count(ti[k]) > 0)
					projJoin.push_back(k);
				else
					projNonJoin.push_back(k);
			}
			
			projsSel.push_back(projSel);
			projsJoin.push_back(projJoin);
			projsNonJoin.push_back(projNonJoin);
		}
		
		if (mode == 2){
		
			joinRows.reset(new Rows(joinCols.size() ));
		
			joinRows->reserve(32 * coords.size() );
		}
	}
	
};*/

class JoinProcessor{

public:

	int mode = 0;
	
	const Tuple joinA;
	const Tuple joinB;
	
	double joinsize = 0;
	
	int passes = 1;
	
	Rows sampleRows;
	vector<Row> samples;
	
	
	Tuple smpCols;
	
	shared_ptr<Rows> sampleRows1;
	shared_ptr<Rows> sampleRows2;
	
	shared_ptr<Rows::iterator> it1;
	shared_ptr<Rows::iterator> it2;
	
	vector<double>::const_iterator itw;
	
	shared_ptr<Weighter> weighter1;
	shared_ptr<Weighter> weighter2;
	
	WeightingFunction* f;
	
	vector<double> smpCoords;
	
	inline const Rows& getSample() const {
	
		return sampleRows;
	}
	
	inline JoinProcessor(const Tuple& smp, WeightingFunction& ff, const Tuple& joinA, const Tuple& joinB) :  smpCols(smp), sampleRows(smp.size() ), joinA(joinA), joinB(joinB){
		
		f = &ff;
	}
	
	inline void addSample(const vector<double>& smp){
		
		assert (mode == 0);	
		mode = 2;
		smpCoords = smp;
		
		sampleRows.reserve(32 * smpCoords.size() );
		
		cout << "mode " << mode << endl;
	}
	
	inline void addSample(double joinweight, long size, const string seed = "abc"){
	
		RandomNumbers rnd(seed);
		
		smpCoords.reserve(size);
		
		for (long k = 0; k < size; ++k)
			smpCoords.push_back(joinweight * rnd.randomDouble() );
		
		std::sort(smpCoords.begin(), smpCoords.end() );
		
		addSample(smpCoords);
	}
	
	template<typename T>
	inline void addSample(const T& smp){
	
		assert (mode == 0);	
		mode = 1;
	
		long cnt = 0;
		long bytes = 0;
		
		for (auto it = smp.begin(); it != smp.end(); ++it){
			++cnt;
			//bytes += it.getBytes();
		}
		
		reachedEnd = cnt == 0;
		
		bytes += cnt*32;
		
		sampleRows.reserve(bytes);
		
		for (auto it = smp.begin(); it != smp.end(); ++it)
			sampleRows.push_back(*it);
		
		samples.reserve(sampleRows.getSize() );
		
		for (auto it = sampleRows.begin(); it != sampleRows.end(); ++it)
			samples.push_back(*it);
	}

	JoinColumns rel;
	
	Tuple novCols;
	vector<int> nonjoinTab;
	vector<int> nonjoinExt;
	
	Tuple retCols;
		
	Tuple lex;

	vector<int> tabproj;
		
		
	vector<long> tabhash;
		
	vector<int> extproj;
	vector<long> exthash;
	
	
	vector<int> lex1;
	vector<int> lex2;
	
	vector<int> lex3;
		
	vector<int> smpproj1;
	vector<int> smpproj2;
	
	shared_ptr<Rows> lookup;
	
	shared_ptr<RowComparator> ptrcmp1;
	shared_ptr<RowComparator> ptrcmp2;

	bool reachedEnd = true;
	long ops = 0;
	
	long joincard = 0;
	
	inline void init(const Tuple& tabCols, const Tuple& extCols){
	
		//cout << "lex " << lex << endl;
		joinsize = 0;
		
		rel.setParent(tabCols);
		rel.setMain(extCols);
		
		if (mode != 1)
			smpCols.clear();
		
		{
			for (int j = 0; j < tabCols.size(); ++j){
		
				if (mode != 1)
					smpCols.push_back(tabCols[j]);
		
				if (UtilTuple::find(extCols, tabCols[j]) == -1){
					nonjoinTab.push_back(j);
				}
			}
		
			for (int j = 0; j < extCols.size(); ++j){
		
				if (UtilTuple::find(tabCols, extCols[j]) == -1){
				
					novCols.push_back(extCols[j]);
					nonjoinExt.push_back(j);
				
					if (mode != 1)
						smpCols.push_back(extCols[j]);
				}
			}
		}
		
		
		cout << "tabcols " << tabCols << endl;
		cout << "extcols " << extCols << endl;
		
		cout << "novcols " << novCols << endl;
		
		cout << "smpcols " << smpCols << endl;
	
		
		rel.setTarget(retCols);
		
		
		Tuple uniq;
		
		for (int j = 0; j < tabCols.size(); ++j){
		
			if (UtilTuple::find(joinA, tabCols[j]) == -1)
				continue;
				
			if (UtilTuple::find(joinB, tabCols[j]) == -1)
				continue;
				
			if (UtilTuple::find(uniq, tabCols[j]) != -1)
				continue;
			
			uniq.push_back(tabCols[j]);
			
			const int k = UtilTuple::find(extCols, tabCols[j]);
		
			if (k != -1){
				
				tabproj.push_back(j);
				tabhash.push_back(UtilHash::HASHMAX);
				
				extproj.push_back(k);
				exthash.push_back(UtilHash::HASHMAX);
			}
		}
		
		
		for (int j = 0; j < tabCols.size(); ++j){
		
			const int k = UtilTuple::find(smpCols, tabCols[j]);
		
			assert (k != -1);
		
			if (k != -1)
				smpproj1.push_back(k);
		}
		
		for (int j = 0; j < extCols.size(); ++j){
		
			const int k = UtilTuple::find(smpCols, extCols[j]);
		
			assert (k != -1);
		
			if (k != -1)
				smpproj2.push_back(k);
		}
		
		for (int i = 0; i < tabproj.size(); ++i)
			lex.push_back( tabCols[tabproj.at(i)] );
		
		
		for (int i = 0; i < nonjoinTab.size(); ++i)
			lex.push_back( tabCols[nonjoinTab.at(i)] );
		
		for (int i = 0; i < nonjoinExt.size(); ++i)
			lex.push_back( extCols[nonjoinExt.at(i)] );
		
		
		cout << "lex " << lex << endl;
		
		lookup.reset(new Rows(lex.size() ) );
		
		for (int j = 0; j < lex.size(); ++j){
		
			const int k1 = UtilTuple::find(tabCols, lex[j]);
			const int k2 = UtilTuple::find(extCols, lex[j]);
			const int k3 = UtilTuple::find(smpCols, lex[j]);
			
			if (k1 != -1)
				lex1.push_back(k1);
			
			if (k2 != -1)
				lex2.push_back(k2);
				
			if (k3 != -1)
				lex3.push_back(k3);
		}
		
		if (mode == 1){
			
			sampleRows1 = std::make_shared<Rows>(tabCols.size() );
			sampleRows2 = std::make_shared<Rows>(extCols.size() );
			
			RowComparator cmp3(lex3);
		
			std::sort(samples.begin(), samples.end(), cmp3);
		
			sampleRows1->reserve(sampleRows.getBytes() );
			sampleRows2->reserve(sampleRows.getBytes() );
		
			for (auto it = samples.begin(); it != samples.end(); ++it){
		
				sampleRows1->push_back( *it, smpproj1);
				sampleRows2->push_back( *it, smpproj2);
			}
		
			/*
			for (auto it = sampleRows1->begin(); it != sampleRows1->end(); ++it){
		
				cout << tabCols << " => " << (*it) << endl;
			}
		
			for (auto it = sampleRows2->begin(); it != sampleRows2->end(); ++it){
		
				cout << extCols << " => " << (*it) << endl;
			}*/
		
			it1.reset( new Rows::iterator(sampleRows1->begin() ));
			it2.reset( new Rows::iterator(sampleRows2->begin() ));
		
		
		
			dstat1.reserve(samples.size()+10 );
			dstat2.reserve(samples.size()+10 );
			dstat3.reserve(samples.size()+10 );
		}
		
		if (mode == 2)
			itw = smpCoords.begin();
		
		weighter1.reset(new Weighter(*f, tabCols));
		weighter2.reset(new Weighter(*f, extCols));
		
		weighter2->select(nonjoinExt);
		
		ptrcmp1.reset(new RowComparator(lex1));
		ptrcmp2.reset(new RowComparator(lex2));
	}
	
	inline void add(Rows& rows, const Row& r1, const Row& r2){
	
		//cout << r1 << " " << r2 << endl;
	
		rows.push_back(r1, r2, nonjoinExt);
		
		
	}

	inline long double getWeight1(const Row& r) const{
		
		if (weighter1)
			return (*weighter1)(r);
		else
			return 1;
	}

	inline long double getWeight2(const Row& r) const{
	
		if (weighter2)
			return (*weighter2)(r);
		else
			return 1;
	}
	
	vector<double> dstat1;
	vector<double> dstat2;
	vector<double> dstat3;

	inline double dstat(const string seed = "test123"){
	
		assert (mode == 1);
		
		RandomNumbers rnd(seed);
	
		double ret = 0;
		
		long c = 0;
		
		vector<double> smp;
		
		assert (dstat1.size() > 0);
	
		for (long j = 0; j < dstat1.size(); ++j){
		
			const long cnt = dstat3.at(j);
			
			smp.clear();
			
			smp.push_back(0);
			
			for (int k = 0; k < cnt; ++k)
				smp.push_back( rnd.randomDouble() );
			
			std::sort(smp.begin(), smp.end() );
			
			//cout << "cnt " << cnt << endl;
			
			for (int k = 0; k <= cnt; ++k){
			
				double d = (c+k)*1.0/samples.size() - ( dstat1.at(j)+dstat2.at(j)*smp.at(k) )/joinsize;
			
				if (d < 0.0)
					d = -d;
			
				//cout << "k " << k << " d " << d << endl;
				
				if (d > ret)
					ret = d;
			}
			
			c += cnt;
		}
		
		return ret;
	}
	
	inline bool less(RowComparator& cmp1, const Row& a1, const Row& b1) const{
	
		return cmp1(a1,b1);
	}
	
	inline bool neq(RowComparator& cmp1, const Row& a1, const Row& b1) const{
	
		return less(cmp1, a1, b1) || less(cmp1, b1, a1);
	}
	
	inline bool equal(RowComparator& cmp1, const Row& a1, const Row& b1) const{
	
		return !neq(cmp1, a1, b1);
	}
	
	inline bool less(RowComparator& cmp1, RowComparator& cmp2, const Row& a1, const Row& a2, const Row& b1, const Row& b2) const{
	
		return less(cmp1, a1, b1) ||( equal(cmp1, a1, b1) && less(cmp2, a2, b2) );
	}
	
	inline bool equal(RowComparator& cmp1, RowComparator& cmp2, const Row& a1, const Row& a2, const Row& b1, const Row& b2) const{
	
		return equal(cmp1, a1, b1) && equal(cmp2, a2, b2);
	}
	
	inline bool neq(RowComparator& cmp1, RowComparator& cmp2, const Row& a1, const Row& a2, const Row& b1, const Row& b2) const{
	
		return neq(cmp1, a1, b1) || neq(cmp2, a2, b2);
	}
	
	inline long getJoinSize() const {
	
		return joincard;
	}
	
	inline long getOps() const {
	
		return ops;
	}
	
	inline double getJoinWeight() const {
	
		return joinsize;
	}
	
	inline void proc(const std::vector<Row>::const_iterator& r1, const std::vector<Row>::const_iterator& r1end, const std::vector<Row>::const_iterator& r2, const std::vector<Row>::const_iterator& r2end){
		
		double w1sum = 0;
		double w2sum = 0;
		
		for (auto i1 = r1; i1 != r1end; ++i1)
			w1sum += getWeight1(*i1);
		
		for (auto i2 = r2; i2 != r2end; ++i2)
			w2sum += getWeight2(*i2);
			
		RowComparator& cmp1 = *ptrcmp1;
		RowComparator& cmp2 = *ptrcmp2;
		
		Rows::iterator& iter1 = *it1;
		Rows::iterator& iter2 = *it2;
		
		
		double jj = joinsize;
		double smpw = jj;
		
		joincard += ((long)(r1end-r1))*((long)(r2end-r2));
		joinsize += w1sum*w2sum;
		
		
		if (mode == 2){
			
			if (itw == smpCoords.end() )
				return;
			
			//cout << *itw << endl;
			
			if (*itw > joinsize)
				return;
			
			for (auto i1 = r1; i1 != r1end; ++i1){
		
				const Row& row1 = (*i1);
				
				const double w1 = getWeight1(row1);
		
				if ( *itw > jj+w1*w2sum ){
				
					jj += w2sum * w1;
				
					ops++;
				
					continue;
				}
				
				for (auto i2 = r2; i2 != r2end; ++i2){
				
					const Row& row2 = (*i2);
				
					const double w2 = getWeight2(row2);
					
					jj += w2*w1;
					
					while ( *itw <= jj){
					
						add(sampleRows, row1, row2);
					
						ops++;
						++itw;
					
						if (itw == smpCoords.end() )
							return;
					}
					
					
					ops++;
					
					
				}
			}
		}
		
		if (mode == 1){
		
			for (auto i1 = r1; i1 != r1end; ++i1){
			
			
				if ( (i1+1) != r1end )
					assert ( !less(cmp1, *(i1+1), *i1) );
			
				const Row& row1 = (*i1);
				
				const double w1 = getWeight1(row1);
			
				// result < sample 
				if ( reachedEnd || neq(cmp1, row1, *iter1) ){
			
					//cout << row1 << " " << *iter1 << endl;
			
					assert (reachedEnd ||less(cmp1, row1, *iter1));
			
					jj += w2sum * w1;
				
					smpw = jj;
				
					ops++;
				
					continue;
				}
			
			
			
				//assert( !((*it1) != sampleRows1->end()) || !cmp1(*iter1, row1) );
			
				// result = sample?
		
				for (auto i2 = r2; i2 != r2end; ++i2){
		
					if ( (i2+1) != r2end )
						assert ( !less(cmp2, *(i2+1), *i2) );
		
					const Row& row2 = (*i2);
				
					const double w2 = getWeight2(row2);
				
				
					if ((i2+1) != r2end && equal(cmp2, row2, *(i2+1))){
					
						smpw += w1*w2; 
						ops++;
						continue;
					}
				
					// result != sample
					if (reachedEnd || neq(cmp1, cmp2, row1, row2, *iter1, *iter2)){
				
				
						assert (reachedEnd ||less(cmp1, cmp2, row1, row2, *iter1, *iter2) );
				
						smpw += w1*w2; 
						ops++;
						continue;
					}
				
					long smpcount = 0;

					while( !reachedEnd && equal(cmp1, cmp2, row1, row2, *iter1, *iter2) ){
				
						++smpcount;
						++iter1;
						++iter2;
					
						++ops;
					
						const bool good = iter1 != sampleRows1->end(); 
						reachedEnd = !good;
					}
				
					dstat1.push_back( smpw);
					dstat2.push_back( w1*w2);
					dstat3.push_back( smpcount);
				
					if ( reachedEnd || neq(cmp1, row1, *iter1) )
						break;
				
					++ops;
					smpw += w1*w2;
				}
			
				jj += w2sum * w1;
				ops++;
			
				smpw = jj;
			}
		}
		
		
		 //((long)(r1end-r1))*((long)(r2end-r2));
	}

};


class DatabaseTable{

protected:

	class SampleInfo{
	public:
		const string seed;
		const long size;
	
		inline SampleInfo(const SampleInfo& inf) : size(inf.size), seed(inf.seed){
		
		}
	
		inline SampleInfo(long sampleSize, const string& seed) : size(sampleSize), seed(seed){
		
		}
	
	};
	
	class TableProfile{
	public:
		shared_ptr<DatabaseTable> sample;
	
		long size;
		
		Tuple columns;
		Tuple joinColumns;
		
		vector<long> columnsDistinctValues;

		long distinctJoinRows;
		long distinctRows;
	
		inline const string toString(int j = 5) const{
	
			stringstream ss;
	
			ss << sample->toString(j) << endl;
		
			ss << "columns " << columns << endl;
			ss << "joinColumns " << joinColumns << endl;
		
			ss << "distinct [" << joinColumns <<"] "<< UtilString::toString(columnsDistinctValues) << endl;
		
			ss << "distinctJoinRows " << distinctJoinRows << endl;
			ss << "distinctRows "<< distinctRows << endl;
			ss << "size " << size << endl;
			
			return ss.str();
		}

	};
	
	class TableProfiler{

		public:
			WeightedRowSample sample;
		
			vector<shared_ptr<Profiler>> columnsDistinct;
		
			Profiler allColumnsDistinct;
			Profiler joinColumnsDistinct;
		
			vector<int> join;
		
			inline shared_ptr<TableProfile> get(){
		
				shared_ptr<TableProfile> ret = std::make_shared<TableProfile>();
			
				sample.withoutReplacement();
			
				ret->sample = std::make_shared<DatabaseTable>("", columns);
				
				ret->sample->clear();
				
				ret->sample->reserve(sample.maxSize);
			
				for (auto it = sample.begin(); it != sample.end(); ++it)
					ret->sample->push_back(*it);
			
				ret->distinctJoinRows = joinColumnsDistinct.distinct();
			
				ret->distinctRows = allColumnsDistinct.distinct();
			
				for (long k = 0; k < columnsDistinct.size(); ++k)
					ret->columnsDistinctValues.push_back( columnsDistinct.at(k)->distinct() );
			
				ret->size = allColumnsDistinct.count();
			
				ret->columns = columns;
				ret->joinColumns = joinColumns;
			
			
				return ret;
			}
		
			Tuple columns;
			Tuple joinColumns;
		
			inline TableProfiler(Tuple columns, Tuple joinColumns = Tuple(), int smp=10000) : sample(columns, "abc", smp) {
		
				this->columns = columns;
				this->joinColumns = joinColumns;
		
				for (int j = 0; j < joinColumns.size(); ++j){
				
					join.push_back( UtilTuple::find(columns, joinColumns[j] ) );
				
					columnsDistinct.push_back( std::make_shared<Profiler>() );
				}
			
				sample.initOnlineWOR();
			}
		
		
			template<typename T>
			inline void add(const T& t){
		
				for (auto it = t.begin(); it != t.end(); ++it)
					add( (*it) );
			}
		
		
			template<typename T>
			inline void add(shared_ptr<T> t){
		
				return add(*t);
			}
		
			inline void add(const Row& row){
			
				sample.add(row, 1);
			
				for (int i = 0; i < join.size(); ++i){
			
					const char* p = row.get( join.at(i) );
			
					joinColumnsDistinct.collect(p);
				
					columnsDistinct.at(i)->collect(p);
					columnsDistinct.at(i)->add();
				}
			
				if (join.size() > 0)
					joinColumnsDistinct.add();
			
				for (int i = 0; i < columnsDistinct.size(); ++i){
			
					const char* p = row.get(i);
			
					allColumnsDistinct.collect(p);
				}
			
				allColumnsDistinct.add();
			}

	};
	
	long _size = -1;
	
	shared_ptr<TableProfile> profile;
	
	shared_ptr<Rows> rows;
	
	shared_ptr<DatabaseTable> wrapped;
	
	shared_ptr<SampleInfo> sampleInfo;
	
	vector<SampleInfo> sampleInfoStack;

public:

	inline bool inMemory() const{
	
		if (rows)
			return true;
			
		return false;
	}

	inline const TableProfile& getProfile() const {
	
		return *profile;
	}
	
	inline long estimatedBytes(long m = -1) const {
	
		return 32*colNames.size()*(m == -1 ? getSize() : m);
	}

	inline void disableSampling(){
		
		sampleInfo.reset();
		
		if (sampleInfoStack.size() > 0){
		
			sampleInfo.reset( new SampleInfo(sampleInfoStack.back() ) );
			sampleInfoStack.pop_back();	
		}
		
	}
	
	inline void enableSamplingFrequency(double freq, const string seed="abc"){
	
		if (freq >= 0)
			return enableSampling(getSize()*freq, seed);
	}
	
	
	inline void enableSamplingFrequency(double freq, long max, const string seed="abc"){
	
		if (freq >= 0 && getSize() > max)
			return enableSampling(getSize()*freq, seed);
	}
	
	inline void enableSampling(long sampleSize = -1, const string seed="abc"){
		
		if (sampleSize >= getSize() )
			return;
		
		if (sampleSize >= 0){
		
			if (sampleInfo)
				sampleInfoStack.push_back( SampleInfo(*sampleInfo) );	
		
			sampleInfo.reset(new SampleInfo(sampleSize > getSize() ? getSize(): sampleSize, seed));
		}
	}
	

	class iterator{
	
		private: 
		
			Row current;
			
			Tuple t;
			
			Tuple t2;
			
			bool active = false;
			
	
		public:

			std::shared_ptr<TupleStream> stream;
			
			iterator(){
			
			
			}
			
			iterator(const shared_ptr<TupleStream>& s) : stream(s){
			
				active = stream->read(current);
			}

			inline long operator()() const{
			
				return 123;
			}

			inline const Tuple& getTuple(){
			
				t.clear();
				
				current.write(t);
				
				return t; 
			}
			
			inline const Row& operator*() const{
		
				return current;
			}
			
			inline bool operator ==(long b) const{
	
				assert (b == 123);
	
				return !active;
			}
			
			
			
			inline bool operator !=(long b) const{
	
				assert (b == 123);
				
				return active;
			}
	
			inline iterator& operator +=(const long n){
	
				stream->skip(n-1);
				
				active = stream->read(current);
				
				if (!active)
					stream->close();
				
				return (*this);
			}
			
			inline iterator& operator++(){
	
				return ( (*this) += 1 );
			}
	
			inline iterator& operator++(int n){
			
				return ( (*this) += 1 );
			}
			
		};

	const iterator end;
	
	const string tableName;
	
	Tuple colNames;
	
	vector<shared_ptr<MainMemoryIndex>> indices;
	
	inline DatabaseTable(const string name,const Tuple& colNames) : tableName(name), colNames(colNames){
	
		
	}
	
	
	inline DatabaseTable(const string name,const Tuple& colNames, shared_ptr<Rows> rows) : tableName(name), colNames(colNames){
	
		this->rows = rows;
		
		init();	
	}
	
	
	shared_ptr<BasicAggregator> weightAgg;
	shared_ptr<Weighter> weighter;
	
	
	Weight weightMax = -1;
	Weight weightSum = -1;
	
	inline Weight getMaxWeight(){
	
		return weightMax;
	}
	
	inline Weight getTotalWeight(){
	
		return weightSum;
	}
	
	inline void startWeightScan(){
	
		if (weightAgg){
		
			assert(false);
		} else {
		
			weightAgg = std::make_shared<BasicAggregator>();
		}
	}
	
	inline void weightRegister(const Weight& w){
	
		if (weightAgg)
			weightAgg->add(w);
	}
	
	inline void weightRegister(const Row& r){
	
		if (weighter)
			weightRegister( (*weighter)(r) );
	}
	
	inline void startWeightScan(WeightingFunction& f, Tuple t, Tuple sel){
	
		assert (t.size() != 0);
		
		weighter = std::make_shared<Weighter>(f, t, sel);
		startWeightScan();
	}
	
	inline void endWeightScan(){
	
		if (weightAgg){
	
			weightMax = weightAgg->max();
			weightSum = weightAgg->sum();
	
			weightAgg.reset();
			weighter.reset();
			
			
		
		} else {
		
			assert (false);
		}
		
	}
	
	inline DatabaseTable(const DatabaseTable& tab) : DatabaseTable(tab.getTableName(), tab.getColNames())  {
		
		reserve( tab.getSize() );
		
		assert (inMemory() );
		
		for (auto it = tab.begin(); it != tab.end(); ++it)
			push_back( *it);
		
		rows->shrink_to_fit();
	}
	
	
	inline DatabaseTable(shared_ptr<DatabaseTable> tab) : DatabaseTable(tab->getTableName(), tab->getColNames()) {
	
		tab->init();
		
		reserve( tab->getSize() );
		
		assert (inMemory() );
	
		for (auto it = tab->begin(); it != tab->end(); ++it)
			push_back( *it);
		
		rows->shrink_to_fit();
	}
	
	inline shared_ptr<DatabaseTable> getSample(long samp, const string& seed) {
		
		enableSampling(samp, seed);
		
		auto ret = std::make_shared<DatabaseTable>(*this);
		
		disableSampling();
		
		return ret;
	}
	
	
	virtual inline ~DatabaseTable(){
	
	}
	
	inline void clear(){
	
		profile.reset();
	
		if (rows)
			rows->clear();
		else
			rows.reset( new Rows(colNames.size() ) );
			
		_size = 0;
	}
	
	/*
	inline vector<int> getJoinCols(DatabaseTable& tab){
	
		
	
	}*/
	
	
	template<typename S,typename T>
	inline void join(S& b, T& t){

		DatabaseTable& a = *this;
		
		
		const Tuple tabCols = a.getColNames();
		const Tuple extCols = b.getColNames();
		
		
		const Tuple& joinA = tabCols;
		const Tuple& joinB = extCols;
		
		
		t.init(tabCols, extCols);
		
		auto rows1 = a.getRows();
		auto rows2 = b.getRows();
		
		if (rows1){
		
			cout << "ROWS1" << endl;
			
		} else {
		
			long size1 = 0;
			long bytes1 = 0;
		
			{
				for (auto it = a.begin(); it != a.end(); ++it){
		
					bytes1 += (*it).getBytes();
					size1++;
				}
			}
			
			rows1 = std::make_shared<Rows>(tabCols.size() );
			
			rows1->reserve(bytes1);
			
			for (auto it = a.begin(); it != a.end(); ++it)
				rows1->push_back(*it);
		
		
			rows1->shrink_to_fit();
		}
		
		if (rows2){
			cout << "ROWS2" << endl;
		
		} else {
		
			
			long size2 = 0;
			long bytes2 = 0;
		
			{
				for (auto it = b.begin(); it != b.end(); ++it){
		
					bytes2 += (*it).getBytes();
					size2++;
				}
			}
			
			rows2 = std::make_shared<Rows>(extCols.size() );
			
			rows2->reserve(bytes2);
			
			for (auto it = b.begin(); it != b.end(); ++it)
				rows2->push_back(*it);
				
			rows2->shrink_to_fit();
		}
		
		vector<Row> v1;
		vector<Row> v2;
		
		v1.reserve(rows1->getSize() );
		v2.reserve(rows2->getSize() );
		
		for (auto it = rows1->begin(); it != rows1->end(); ++it)
			v1.push_back(*it);
			
		for (auto it = rows2->begin(); it != rows2->end(); ++it)
			v2.push_back(*it);
		
		RowComparator& cmp1 = *t.ptrcmp1;
		RowComparator& cmp2 = *t.ptrcmp2;
		
		cout << "sort" << endl;
		
		std::sort(v1.begin(), v1.end(), cmp1);
		std::sort(v2.begin(), v2.end(), cmp2);
		
		
		for (auto it = v1.begin(); it != v1.end(); ++it)
			if ( (it+1) != v1.end() )
				assert (!cmp1( *(it+1), *it) );
		
		
		for (auto it = v2.begin(); it != v2.end(); ++it)
			if ( (it+1) != v2.end() )
				assert (!cmp2( *(it+1), *it) );
		
		
		cout << "/sort" << endl;
		
		long pos1 = 0;
		long pos2 = 0;
		
		long pos10 = -1;
		long pos20 = -1;
		
		long joinsize = 0;
		
		while ( (pos10+1) < v1.size() && (pos20+1) < v2.size() ){
		
			pos1 = pos10+1;
			pos2 = pos20+1;
			
			//cout << "pos " << pos1 << " pos " << pos2 << endl;
			
			//cout << v1.at(pos1) << " " << v2.at(pos2) << endl;
			
			int cmp = v1.at(pos1).cmp(t.tabproj, v2.at(pos2), t.extproj );
			
			//cout << v1.at(pos1) << " " << v2.at(pos2) << " cmp " << cmp << endl ;
		
			while (cmp != 0){
		
				if (cmp == +1){
			
					++pos2;
				
					if (pos2 == v2.size())
						break;	
					
				
				} else {
					
					++pos1;
				
					if (pos1 == v1.size())
						break;
				}
			
				Row& r1 = v1.at(pos1);
				Row& r2 = v2.at(pos2);
				
				
				cmp = r1.cmp(t.tabproj, r2, t.extproj);
			}
			
			
			
			pos10 = pos1;
			pos20 = pos2;
			
			while ( (pos10+1) < v1.size() && v1.at(pos10+1).cmp( t.tabproj, v1.at(pos10), t.tabproj) == 0 )
				++pos10;
				
			while ( (pos20+1) < v2.size() && v2.at(pos20+1).cmp( t.extproj, v2.at(pos20), t.extproj) == 0 )
				++pos20;
			
			auto r1 = v1.begin()+pos1;
			auto r1end = v1.begin()+(pos10+1);
			
			auto r2 = v2.begin()+pos2;
			auto r2end = v2.begin()+(pos20+1);
			
			/*
			const Row* r1 = &v1.at(pos1);
			const Row* r1end = (&v1.at(pos10))+1;
			
			const Row* r2 = &v2.at(pos2);
			const Row* r2end = (&v2.at(pos20))+1;*/
			
			t.proc(r1, r1end, r2, r2end);
			/*
			for (long i = pos1; i <= pos10; ++i){
				for (long j = pos2; j <= pos20; ++j){
				
					cout << v1.at(i) << " " << v2.at(j) << endl;
				}
			}*/
			
			
			
			//joinsize += (pos10-pos1+1)*(pos20-pos2+1);
		}
		
		//cout << t.joinsize << endl;
		//cout << "lex " << lex << endl;
		
	}
	
	
	
	inline void push_back(const Row& row){
	
		assert (rows);
	
		profile.reset();
	
		assert (rows);
		
		rows->push_back_parse(row);
		
		++_size;
	}
	
	
	inline void reserveBytes(long bytes){
	
		if (rows){
		
		} else {
		
			clear();
		}
	
		if (rows)
			rows->reserve(bytes);
	}
	
	
	inline void reserve(long m){
	
		reserveBytes( estimatedBytes(m));
	}
	
	virtual inline shared_ptr<DatabaseTable> projection(const string& name, const Tuple& names, const Tuple& aliases, bool ignoreMissingColumns){
	
		assert (false);
	}



	inline bool hasColumnUnion(const Tuple& t) const{
		
		return UtilTuple::find(getColNames(), UtilTuple::getString(t, '_')) != -1;
	}
	


	inline shared_ptr<DatabaseTable> addColumnUnion(const Tuple& t){
	
		
		assert (!hasColumnUnion(t) );
	
		Tuple told = getColNames();
		
		
		vector<Tuple> ts;
		
		for (int i = 0; i < told.size(); ++i)
			ts.push_back( Tuple(told[i]) );
			
		ts.push_back(t);
		
		return projection(ts);
	}

	inline shared_ptr<DatabaseTable> projection(vector<Tuple> cols){
	
		auto rows2 = std::make_shared<Rows>(cols.size() );
		
		vector<vector<int>> vcols;
		
		Tuple tnew;
		
		Tuple told = getColNames();
		
		for (int i = 0; i < cols.size(); ++i){
		
			const Tuple& t = cols.at(i);
			vector<int> v;
			
			v.reserve(t.size() );
		
			for (int j = 0; j < t.size(); ++j){
			
				const int k = UtilTuple::find(told, t[j]);	
				assert (k >= 0);
				
				v.push_back(k);
			}
			
			tnew.push_back(UtilTuple::getString(t, '_') );
			vcols.push_back(v);
		}
		
		if (rows)
			rows2->reserve(rows->getBytes() );
		
			
		for (auto it = begin(); it != end(); ++it){
		
			const Row& t = (*it);
			
			rows2->startRow();
			
			for (int i = 0; i < vcols.size(); ++i){
			
				stringstream ss;
				
				const vector<int>& v = vcols.at(i);
				
				for (int j = 0; j < v.size(); ++j){
				
					if (j > 0){
						rows2->push_back_char('0');	
						rows2->push_back_char('1');	
						rows2->push_back_char('2');	
						rows2->push_back_char('3');	
					}
					
					rows2->push_back_chars( t.get(v.at(j)) );
				}
				
				rows2->push_back_char(0);
			}
			
			rows2->finishRow();
		}
		
		
		
		return std::make_shared<DatabaseTable>(tableName, tnew, rows2);
	}
	
	inline bool hasIndices() const{
	
		return indices.size() > 0;
	}
	
	
	
	inline bool hasIndex(int col) const{
	
		assert (col >= 0);
		assert (col < this->getColNames().size() );
		
		if (!hasIndices() )
			return false;
			
		assert (col < indices.size());
		
		if (indices.at(col))
			return true;
		else
			return false;
	}
	
	inline MainMemoryIndex* getIndex(){
		
		for (int j = 0; j < indices.size(); ++j)
			if (hasIndex(j))
				return getIndex(j);
				
		return NULL;
	}
	
	inline MainMemoryIndex* getIndex(int col){
		
		if (hasIndex(col))
			return indices.at(col).get();
		else
			return NULL;
	}
	
	
	inline void buildIndices(const Tuple& joincols){
	
		assert (inMemory() );
		
		rows->lock();
	
		init();
	
		const Tuple& cols = this->getColNames();
		
		if (indices.size() == cols.size() )
			return;
	
		bool foundOne = false;
	
		for (int i = 0; i < cols.size(); i++){
		
			std::shared_ptr<MainMemoryIndex> index;
		
			if (UtilTuple::find(joincols, cols[i]) != -1 ){
			
				cout << "found " << cols[i] << " from " << cols << " in " << joincols << endl;
				
				index = std::make_shared<MainMemoryIndex>(i);
				index->reserve( getSize() );
				
				foundOne = true;
			} else {
			
			
				cout << "not found " << cols[i] << " from " << cols << " in " << joincols << endl;
			}
			
			indices.push_back( index );
		}
		
		if  (!foundOne){
		
			indices.clear();
			return;
		}
		
		for (auto it = begin(); it != end(); ++it){
		
			const Row& row = (*it);
		
			for (int i = 0; i < cols.size(); i++){
				
				if (indices.at(i))
					indices.at(i)->insert(row);
			}
			
		}
		
		for (int i = 0; i < cols.size(); i++)
			if (indices.at(i))
				indices.at(i)->finalise();
		
	}
	
	inline void buildIndices(){
	
		buildIndices(this->getColNames());
	}
	
	inline shared_ptr<Rows> getRows(){
	
		return rows;
	}
	
	
	virtual inline void init(int buildProfile = 0, Tuple joinColNames = Tuple()){
		
		
			
		if (buildProfile > 0){
		
			
			if (profile){
			
				assert(_size > 0);
			
				return;
			}
		
			TableProfiler tabprof(colNames, joinColNames, buildProfile);
			
			for (auto it = begin(); it != end(); ++it)
				tabprof.add(*it);
				
			profile = tabprof.get();
			
			_size = profile->size;
			
			return;
		}
	
		if (_size != -1)
			return;
		
		if (rows){
		
			_size = rows->size();
			return;
		}
		
		_size = 0;
		
		for (auto it = begin(); it != end(); ++it)
			_size++;
	}
	
	/*
	inline const string toString(int n = 10) const {
	
		stringstream ss;
		
		ss << getTableName() << "{" << getColNames() << "}" << endl;
		
		int kmax = getSize()-1;
		int k = 0;
		for (auto it = begin(); it != end(); ++it){
		
		
			if ( !( k < n || k > (kmax-n) ))
				continue;
		
			ss << (*it) << endl;
		}
		
		return ss.str();
	}*/
	
	friend inline std::ostream& operator<<(std::ostream &os, const DatabaseTable& v) { 
    	
    	assert (false);
    	os << v.toString();
		return os;
	}
	
	inline const string getTableName() const{
	
		return tableName;
	}
	
	inline const Tuple getColNames() const{
	
		return colNames;
	}
	
	inline const string toString(int j = -1) const {
	
		stringstream ss;
		
		ss << getTableName() << "{" << getColNames() << "}" << endl;
	
		for (auto it = begin(); it != end(); ++it){
		
			if (j == 0)
				break;
				
			ss << (*it) << endl;
		
			--j;
		}
	
		return ss.str();
	}
	
	
	inline long commonCols( const DatabaseTable& y) const{
	
		long ret = 0;
		
		for (int i = 0; i < colNames.size(); i++){
		
			for (int j = 0; j < y.colNames.size(); j++){
		
				if (UtilString::equals(colNames[i], y.colNames[j] ))
					ret++;
			}	
		}
		
		return ret;
	}
	
	inline bool hasSize() const {
	
		return _size != -1;
	}
	
	inline virtual long getSize() const {
	
		if (sampleInfo)
			return sampleInfo->size;
		
		assert (hasSize() );
		
		return _size;
	}
	
	
	/*
	inline virtual iterator begin() const{
		
		return iterator(getStream() );
	}*/
	
	
	
	inline virtual shared_ptr<TupleStream> getStream() const{
	
		return std::make_shared<RowTupleStream>(*rows );
	}
	
	
	inline iterator begin() const{
		
		if (sampleInfo){
			
			long siz = UtilMath::minVal<long>(sampleInfo->size, _size);
		
			if (profile){
			
				if (profile->sample->getSize() >= siz)
					return iterator( std::make_shared<SortedSampleTupleStream>( profile->sample->getStream(), siz, profile->sample->getSize(), sampleInfo->seed) );
			}
			
			return iterator( std::make_shared<SortedSampleTupleStream>(getStream(), siz, _size, sampleInfo->seed) );
		}
			
		return iterator(getStream() );
	}
};



class CSVTable : public DatabaseTable{

	const bool ignoreFirst;

	
	
	const vector<int> sel;
	
	const char sep;
	const string file;
	

public:


	static vector<int> getSelection( const Tuple& allCols, const vector<int>& csvColumns, const Tuple& selectedFrom, bool ignoreMissing = true){
		
		vector<int> ret;
		
		
		
		ret.reserve(selectedFrom.size() );
		
		for (int i = 0; i < selectedFrom.size(); i++){	
			
			const long ind = UtilTuple::find(allCols, selectedFrom[i]);
			
			if (ind == -1){
			
				if (ignoreMissing){
					continue;
				} else {
				
					cout << allCols << " " << selectedFrom << endl;
				
					assert (false);
				}
			}
			
			ret.push_back(csvColumns.at(ind) );
		}
		
		return ret;
	}
	
	static vector<int> getSelection( const Tuple& allCols, const Tuple& selectedFrom, bool ignoreMissing = true){
	
		vector<int> v;
		v.reserve(allCols.size());
		
		for (int j = 0; j < allCols.size(); ++j)
			v.push_back(j);
			
		return getSelection(allCols, v, selectedFrom, ignoreMissing);
	}
	
	static const Tuple getAliases( const Tuple& allCols, const Tuple& selectedFrom, const Tuple& aliases){
	
		Tuple ret;
		
		for (int i = 0; i < selectedFrom.size(); i++){	
			
			const long ind = UtilTuple::find(allCols, selectedFrom[i]);
			
			if (ind != -1)
				ret.push_back(aliases[i]);
		}
		return ret;
	}

	inline CSVTable(const string name, const Tuple& aliases, const vector<int>& selection, const string& file, char sep=',', bool ignoreFirst=true) : DatabaseTable(name, aliases), sel(selection), ignoreFirst(ignoreFirst), sep(sep),  file(file){
	
	}

	inline CSVTable(const string name, const Tuple& colNames, const Tuple& cols, const Tuple& aliases, const string& file, char sep=',', bool ignoreFirst=true) : 
	CSVTable( name, aliases, getSelection(colNames, cols, false), file, sep, ignoreFirst ){
	
		
	}
	
	inline shared_ptr<DatabaseTable> projection(const string& name, const Tuple& names, const Tuple& aliases, bool ignoreMissingColumns = true) override{
	
		assert (!hasIndices());
	
		auto ret = std::make_shared<CSVTable>(name, getAliases(colNames, names, aliases), getSelection(colNames, sel, names, ignoreMissingColumns), file, sep, ignoreFirst);
		
		ret->addColNames = addColNames;
		
		ret->_size = _size;
		
		return ret;
	}
	
	/*
	inline shared_ptr<DatabaseTable> projection(const Tuple& names, bool ignoreMissingColumns = true){
	
		return projection(names, names, ignoreMissingColumns);
	}*/
	
	bool addColNames = false;

	inline shared_ptr<TupleStream> getStream() const override{
		
		auto ret = std::make_shared<CSVStream>(file, sel, colNames, sep, ignoreFirst);
		
		ret->addColNames = addColNames;
		
		return ret;
	}

};

/*
class TableSample : public DatabaseTable{

	shared_ptr<DatabaseTable> tab;
	
	const string seed;
	
public:

	inline TableSample(shared_ptr<DatabaseTable> tab, long sampleSize, const string seed="abc") : DatabaseTable(tab->getTableName(), tab->getColNames()), seed(seed), tab(tab){
	
		if (sampleSize > tab->getSize() )
			sampleSize = tab->getSize();
		
		_size = sampleSize;
	}
	
	inline shared_ptr<TupleStream> getStream() const override{
		
		return std::make_shared<SortedSampleTupleStream>(tab->getStream(), _size, tab->getSize(), seed );
	}
};
*/

/*
class MainMemoryTable : public DatabaseTable{

	Rows data;
	

public:
	
	
	inline MainMemoryTable(const string name, const Tuple& colNames) : DatabaseTable(name, colNames), data(colNames.size()){
	
		
	} 
	
	inline MainMemoryTable(const DatabaseTable& tab, long smp, const string& seed) : DatabaseTable(tab.getTableName(), tab.getColNames()), data(tab.getColNames().size()) {
	
		SortedSample samp(smp > tab.getSize() ? tab.getSize() : smp, tab.getSize(), seed);
		
		data.reserve( tab.getSize()*32 * tab.getColNames().size() );
	
		auto it = tab.begin();
		auto itend = tab.end();
		
		if (!(it != itend))
			return;
		
		while (samp){
			
			for (long k = samp.getSkip(); k >= 1; --k){
			
				assert (it != itend );
				++it;
				
				if (!(it != itend))
					return;
			
			}
		
			push_back(*it);
			
			++it;
			++samp;
			
			if (!(it != itend))
				return;
			
		}
	
		
	}
	
	inline MainMemoryTable(const DatabaseTable& tab) : DatabaseTable(tab.getTableName(), tab.getColNames()), data(tab.getColNames().size()) {
	
		data.reserve( tab.getSize()*32 * tab.getColNames().size() );
	
		for (auto it = tab.begin(); it != tab.end(); ++it)
			push_back( *it);
		
	}
	
	inline MainMemoryTable(shared_ptr<DatabaseTable> tab) : DatabaseTable(tab->getTableName(), tab->getColNames()), data(tab->getColNames().size()) {
	
		tab->init();
	
		data.reserve( tab->getSize()*32 * tab->getColNames().size() );
	
		for (auto it = tab->begin(); it != tab->end(); ++it)
			push_back( *it);
		
	}
	
	inline void reserve(long bytes){
	
		data.reserve(bytes);
	}
	
	inline Rows& getRows(){
	
		return data;
	}
	
	inline void clear(){
	
		data.clear();
	}
	
	inline void push_back(const Row& t){
	
		data.push_back(t);
	}
	
	inline void push_back(const Tuple& t){
	
		data.push_back(t);
	}
	
	inline long size() const{
	
		return data.size();
	}
	
	inline long getSize() const override{
	
		if (sampleInfo)
			return sampleInfo->size;
			
		return size();
	}
	
	
	
	inline shared_ptr<TupleStream> getStream() const override{
		
		shared_ptr<TupleStream> ret = std::make_shared<RowTupleStream>(data );
		
		return ret;
	}



};*/



/*
class VectorTable : public DatabaseTable{

	vector<Tuple> data;

public:
	
	
	inline VectorTable(const string name, const Tuple& colNames) : DatabaseTable(name, colNames){
	
		
	}
	
	
	inline VectorTable(shared_ptr<DatabaseTable> tab) : DatabaseTable(tab->getTableName(), tab->getColNames()) {
	
		for (auto it = tab->begin(); it != tab->end(); ++it)
			push_back( it.getTuple() );
	
		
	}
	
	inline const Tuple& operator[](std::size_t idx) const{
		
		return data.at(idx);
	}
	
	inline void push_back(const Tuple& t){
	
		data.push_back(t);
	}
	
	inline long size() const{
	
		return data.size();
	}
	
	inline long getSize() const override{
	
		return size();
	}
	
	inline void clear(){
	
		data.clear();
	}
	
	inline shared_ptr<TupleStream> getStream() const override{
		
		auto ret = std::make_shared<VectorTupleStream>(data);
		
		return ret;
	}

};*/




/*
class ODBCTable : public DatabaseTable{

public:
	
	
	const string conn;
	const string cmd;
	
	long size = -1;
	
	inline ODBCTable(const string tableName, const Tuple& colNames, const string db, const string cmd) : DatabaseTable(tableName,colNames), conn(db), cmd(cmd){
	
		
	}
	
	inline void init() override{
		
		if (size == -1)
			size = ODBCInterface(conn).getInt("SELECT COUNT(*) FROM ("+cmd+") AS COUNTEDTABLE");
	}
	
	inline const string toString() const{
	
		return tableName+" SIZE "+to_string(getSize());
	}
	
	
	inline long getSize() const override{
	
		assert (size != -1);
		
		return size;
	}
	
	inline shared_ptr<TupleStream> getStream() const override{
	
		return std::make_shared<ODBCInterface>(conn, cmd, colNames.size());
	}
};*/


class RowHashMultiset{

	vector<long> linearIndex;
	RowHashMap<long> headCatalog;
	
	unique_ptr<RowHashMap<long>> tailCatalog;
	unique_ptr<GrowingVector> fragmentedIndex;

public:
	inline void reserve(long size){
	
		fragmentedIndex->reserve(size*2);
		tailCatalog->reserve(size);
		headCatalog.reserve(size);
	}
	
	inline RowHashMultiset(){
	
		tailCatalog = make_unique<RowHashMap<long>>();
		fragmentedIndex = make_unique<GrowingVector>();
	}
	
	inline const RowHashMap<long>& getHeadCatalog() const {
	
		return headCatalog;
	}
	
	inline bool insert(const Row& row, const vector<int>& proj, const vector<long>& hash, long entryInd){
	
		if (row.isNull(proj))
			return false;
		
		assert (entryInd >= 0);
		
		assert (tailCatalog);
	
		if (tailCatalog->hasKey(row, proj, hash))
			return false;
		
		const long head = fragmentedIndex->insert();
		
		headCatalog.insert(row, proj, hash, head);
		tailCatalog->insert(row, proj, hash, head);
		
		push_back(row, proj, hash, entryInd);
		
		return true;
	}
	
	inline void push_back(const Row& row, const vector<int>& proj, const vector<long>& hash, long val){
	
		assert(val >= 0);
		assert (tailCatalog->hasKey(row, proj, hash));
		
		long& m = tailCatalog->get(row, proj, hash);
		
		m = fragmentedIndex->push_back(m, val);
	}
	
	inline long get(const Row& row, const vector<int>& proj, const vector<long>& hash){
	
		return headCatalog.get(row, proj, hash);
	}
	
	inline bool hasKey(const Row& row, const vector<int>& proj, const vector<long>& hash){
	
		return headCatalog.hasKey(row, proj, hash);
	}
	
	inline long getEntry(long pos) const {
	
		return linearIndex.at(pos);
	}
	
	inline long getSize(long pos) const {
	
		return linearIndex.at(pos+1);
	}
	
	inline const long* begin(long pos) const{
		
		return linearIndex.data()+(pos+2);
	}
	
	inline const long* last(long pos) const{
		
		return begin(pos)+(getSize(pos)-1);
	}
	
	inline const long* end(long pos) const{
	
		return begin(pos)+getSize(pos);
	}
	
	inline void swap(const long* ptr1, const long* ptr2){
	
		const long v1 = (*ptr1);
		const long v2 = (*ptr2);
		
		linearIndex.at(ptr1-linearIndex.data()) = v2;
		linearIndex.at(ptr2-linearIndex.data()) = v1;
	}
	
	template<typename Key>
	inline void finalise(unordered_map<Key, long>& map, const vector<long_double>* sortKeys = NULL){
	
		vector< std::pair<long_double, long> > sorted;
		
		for (auto it = map.begin(); it != map.end(); ++it){

			const long head1 = fragmentedIndex->next(it->second);
			
			assert (head1 != -1);
			
			it->second = linearIndex.size();
			
			const long entryInd = fragmentedIndex->get(head1);
			
			// entry id
			linearIndex.push_back( entryInd );
			
			const long head = fragmentedIndex->next(head1+1);
			
			assert (head != -1);
			
			long count = 0;
		
			for (long x = fragmentedIndex->next(head); x != -1; x = fragmentedIndex->next(x+1))
				count++;
			
			linearIndex.push_back(count);
			
			if (count == 0)
				continue;
			
			if (sortKeys){
			
				sorted.clear();
				sorted.reserve(count);
				
				for (long x = fragmentedIndex->next(head); x != -1; x = fragmentedIndex->next(x+1)){
				
					const long elemInd = fragmentedIndex->get(x);
				
					sorted.push_back( std::make_pair(sortKeys->at(elemInd), elemInd) );
				}
				
				std::sort(sorted.begin(), sorted.end() );
				
				for (long j = 0; j < sorted.size(); ++j)
					linearIndex.push_back(sorted.at(j).second);
					
				for (long k = it->second+3; k < linearIndex.size(); ++k)
					assert (sortKeys->at(linearIndex.at(k)) >= sortKeys->at(linearIndex.at(k-1)));
					
			} else {
				
				for (long x = fragmentedIndex->next(head); x != -1; x = fragmentedIndex->next(x+1))
					linearIndex.push_back( fragmentedIndex->get(x) );
			}
		}
	}
	
	template<typename Key>
	inline void finalise2(long_map<Key, long>& map, const vector<long_double>* sortKeys = NULL){
	
		vector< std::pair<long_double, long> > sorted;
		
		for (auto it = map.begin(); it != map.end(); ++it){

			const long head1 = fragmentedIndex->next(it->second);
			
			assert (head1 != -1);
			
			it->second = linearIndex.size();
			
			const long entryInd = fragmentedIndex->get(head1);
			
			// entry id
			linearIndex.push_back( entryInd );
			
			const long head = fragmentedIndex->next(head1+1);
			
			assert (head != -1);
			
			long count = 0;
		
			for (long x = fragmentedIndex->next(head); x != -1; x = fragmentedIndex->next(x+1))
				count++;
			
			linearIndex.push_back(count);
			
			if (count == 0)
				continue;
			
			if (sortKeys){
			
				sorted.clear();
				sorted.reserve(count);
				
				for (long x = fragmentedIndex->next(head); x != -1; x = fragmentedIndex->next(x+1)){
				
					const long elemInd = fragmentedIndex->get(x);
				
					sorted.push_back( std::make_pair(sortKeys->at(elemInd), elemInd) );
				}
				
				std::sort(sorted.begin(), sorted.end() );
				
				for (long j = 0; j < sorted.size(); ++j)
					linearIndex.push_back(sorted.at(j).second);
					
				for (long k = it->second+3; k < linearIndex.size(); ++k)
					assert (sortKeys->at(linearIndex.at(k)) >= sortKeys->at(linearIndex.at(k-1)));
					
			} else {
				
				for (long x = fragmentedIndex->next(head); x != -1; x = fragmentedIndex->next(x+1))
					linearIndex.push_back( fragmentedIndex->get(x) );
			}
		}
	}
	
	inline void finalise(const vector<long_double>& sortKeys){
	
		linearIndex.reserve(fragmentedIndex->getSize() );
		
		tailCatalog.reset();
		
		finalise2(headCatalog.longMap, &sortKeys);
		finalise(headCatalog.stringMap, &sortKeys);
		
		fragmentedIndex.reset();
		linearIndex.shrink_to_fit();
	}
	
	inline void finalise(){
	
		linearIndex.reserve(fragmentedIndex->getSize() );
		
		tailCatalog.reset();
		
		finalise2(headCatalog.longMap, NULL);
		finalise(headCatalog.stringMap, NULL);
		
		fragmentedIndex.reset();
		linearIndex.shrink_to_fit();
	}
	
	
	
};

/*
class TableJoin{

	Rows rowsA;
	Rows rowsB;
	
	RowHashMultiset indexa;
	RowHashMultiset indexb;
	
	long getSize;
	long getBytes;
	
	inline TableJoin(){
	
	}

	
	template<typename A, typename B>
	inline void join(A& a, B& b){
	
		long entryId = 0;
	
		for (auto it = a.begin(); it != a.end(); ++it){
		
			if (indexa.insert( (*it), tabCols, tabHash, entryId){
			
			
				++entryId;
			}
		}
		
		for (auto it = b.begin(); it != b.end(); ++it){
		
			if (!a.hasKey( (*it), tabCols, tabHash ))
				continue;
		
			if (indexb.insert( (*it), tabCols, tabHash, entryId){
			
				
				++entryId;
			}
		}
	}
};*/


namespace UtilTables{

	
	
	
	inline shared_ptr<DatabaseTable> getJoin(long& onlyJoinSize, DatabaseTable& a, const Tuple& joinA, DatabaseTable& b, const Tuple& joinB, long size1=-1, long bytes1 = -1, long size2=-1, long bytes2=-1, bool flip=false){
	
		
		if (onlyJoinSize == 0)
			cout << a.getColNames() << " ( size " << size1 << " bytes " << bytes1 << ") is being joined with " << b.getColNames() << " ( size " << size2 << " bytes " << bytes2 << ")" << endl;
		
		const Tuple tabCols = a.getColNames();
		const Tuple extCols = b.getColNames();
		
		JoinColumns rel;
		
		rel.setParent(tabCols);
		rel.setMain(extCols);
		
		Tuple novCols;
		vector<int> novel;
		
		if (!flip){
		
			for (int j = 0; j < extCols.size(); ++j){
			
				if (UtilTuple::find(tabCols, extCols[j]) == -1){
					
					novCols.push_back(extCols[j]);
					novel.push_back(j);
					
				}
			}
		
		} else {
		
			for (int j = 0; j < tabCols.size(); ++j){
			
				if (UtilTuple::find(extCols, tabCols[j]) == -1){
				
					novCols.push_back(tabCols[j]);
					novel.push_back(j);
				}
			}
		}
		/*
		JoinColumns rel3;
		rel3.setMain(!flip ? extCols: tabCols);
		rel3.addChild(!flip ? tabCols: extCols);*/
		
		Tuple retCols;
		
		stringstream name;
		
		if (!flip){
		
			retCols.push_back(tabCols);
			retCols.push_back(novCols);
			
			name << a.tableName << " join " << b.tableName;
			
		} else {
		
			
			retCols.push_back(extCols);
			retCols.push_back(novCols);
			
			name << b.tableName << " join " << a.tableName;
		}
		
		rel.setTarget(retCols);
		
		vector<int> tabproj;
		vector<long> tabhash;
		
		vector<int> extproj;
		vector<long> exthash;
		
		Tuple uniq;
		
		for (int j = 0; j < tabCols.size(); ++j){
		
			if (UtilTuple::find(joinA, tabCols[j]) == -1)
				continue;
				
			if (UtilTuple::find(joinB, tabCols[j]) == -1)
				continue;
				
			if (UtilTuple::find(uniq, tabCols[j]) != -1)
				continue;
			
			uniq.push_back(tabCols[j]);
			
			const int k = UtilTuple::find(extCols, tabCols[j]);
		
			if (k != -1){
				
				tabproj.push_back(j);
				tabhash.push_back(UtilHash::HASHMAX);
				
				extproj.push_back(k);
				exthash.push_back(UtilHash::HASHMAX);
			}
		}
		
		if (extproj.size() == 0){
		
			assert (onlyJoinSize == 1),
			
			onlyJoinSize = a.getSize()*b.getSize();
			
			return nullptr;
		}
		
		assert (extproj.size() == tabproj.size() );
		
		
		for (int j = 0; j < tabproj.size(); ++j){
		
			assert ( UtilString::equals(tabCols[tabproj[j]], extCols[extproj[j]]) );
		}
		
		if (!flip){
						
			for (long j = 0; j < retCols.size(); j++){
			
				const int k = j-tabCols.size();
			
				if (k < 0)
					assert (UtilString::equals(retCols[j], tabCols[j]) );
				else {
					
					/*
				 	cout << "ret[" << retCols << "] tab[" << tabCols << "] ext[" << extCols << "] nov["<<  novCols << "]" << endl;				 	cout << UtilString::toString(tabproj) << endl;
					cout << UtilString::toString(novel) << endl;
					
					cout << retCols[j] << " != " << extCols[novel[k]] << endl;*/
				
					assert (UtilString::equals(retCols[j], extCols[novel[k]]) );
				}
			}			
	
		} else {
			
			for (long j = 0; j < retCols.size(); j++){
			
				const int k = j-extCols.size();
			
				if (k < 0)
					assert (UtilString::equals(retCols[j], extCols[j]) );
				else
					assert (UtilString::equals(retCols[j], tabCols[novel[k]]) );
			}
		}
		
		
		
		if (size1 == -1){
			
			size1 = 0;
			bytes1 = 0;
			size2 = 0;
			bytes2 = 0;
		
			for (auto it = a.begin(); it != a.end(); ++it){
		
				bytes1 += (*it).getBytes();
				size1++;
			}
		
			long bytes2 = 0;
			long size2 = 0;
		
			for (auto it = b.begin(); it != b.end(); ++it){
		
				bytes2 += (*it).getBytes();
				size2++;
			}
			
			if (bytes2 < bytes1 )
				return getJoin(onlyJoinSize, b, joinB, a, joinA, size2, bytes2, size1, bytes1, !flip);
		}
		
		
		const bool hasIndex1 = tabproj.size() == 1 && a.getIndex(tabproj.at(0)) != NULL;
		const bool hasIndex2 = extproj.size() == 1 && b.getIndex(extproj.at(0)) != NULL;
		
		if (hasIndex2 && (!hasIndex1))
			return getJoin(onlyJoinSize, b,joinB, a, joinA, size2, bytes2, size1, bytes1, !flip);
			
		if (hasIndex1){
		
			MainMemoryIndex* ind = a.getIndex(tabproj.at(0) );
		
			long joinSize = 0;
		
			for (auto it = b.begin(); it != b.end(); ++it){
			
				if (ind->begin(*it, extproj.at(0) ) == NULL)
					continue;
			
				joinSize += ind->end(*it, extproj.at(0) )-ind->begin(*it, extproj.at(0) );
			}
			
			if (onlyJoinSize == 1){
			
				onlyJoinSize = joinSize;
				
				return nullptr;
			}
			
			auto ret = std::make_shared<DatabaseTable>(name.str(), retCols);
		
			ret->clear();
			ret->reserve(joinSize);
		
			Rows& retrows = *ret->getRows();
		
			
			for (auto it = b.begin(); it != b.end(); ++it){
				
				const Row& extrow = (*it);
				
				const Row* bg = ind->begin(*it, extproj.at(0));
				const Row* en = ind->end(*it, extproj.at(0));
				
				for (auto it2 = bg; it2 != en; ++it2){
				
					const Row& tabrow = (*it2);
					
					if (!flip){
						
						retrows.push_back(tabrow, extrow, novel);
				
					} else {
				
						retrows.push_back(extrow, tabrow, novel);
					
					}
				}
				
			}
			
			
			
			
			return ret;
		}
		
		Rows rowsdata(a.getColNames().size() );
		vector<Row> rows;
		
		rows.reserve(size1);
		rowsdata.reserve(bytes1);
		
		
		vector<long> index;
		RowHashMap<long> index1;
		
		
		{
			GrowingVector index0;
			index0.reserve(size1*2);
			index1.reserve(size1);
		
			{
				RowHashMap<long> index2;
				index2.reserve(size2);
		
				Row indexrow;
				
				for (auto it = a.begin(); it != a.end(); ++it)
					rowsdata.push_back(*it);	
		
				
				Timer t("growing vector");
				
				if (onlyJoinSize == 0)
					t.start(10, 100000, size1);
				
		
				for (auto it = rowsdata.begin(); it != rowsdata.end(); ++it){
		
					if (onlyJoinSize == 0)
						t.tick();
					
					const Row& indexrow = (*it);
					
					rows.push_back(indexrow);
			
					if (!index2.hasKey(indexrow, tabproj, tabhash)){
				
						const long ind = index0.insert();
				
						index1.insert(indexrow, tabproj, tabhash, ind);
						index2.insert(indexrow, tabproj, tabhash, ind);
					}
					
					long& m = index2.get(indexrow, tabproj, tabhash);
			
					m = index0.push_back(m, rows.size()-1);
				}
				
				if (onlyJoinSize == 0)
					t.end();
			}
		
			index.reserve(size1*2);
		
			long_map<long, long>& map1 = index1.longMap;
			unordered_map<string, long>& map2 = index1.stringMap;
		
			for (auto it = map1.begin(); it != map1.end(); ++it){
	
				const long x1 = it->second;
		
				const long sizeInd = index.size();
		
				it->second = index.size();
		
				index.push_back(0);
				index.push_back(0);
				
				long count = 0;
				long bytes = 0;
	
				for (long x = index0.next(x1); x != -1; x = index0.next(x+1)){
		
					const long ind = index0.get(x);
					index.push_back( ind);
				
					bytes += rows.at(ind).getBytes();
					++count;
				}
			
				index.at(sizeInd) = count;
				index.at(sizeInd+1) = bytes;
			}
			
			for (auto it = map2.begin(); it != map2.end(); ++it){
	
				const long x1 = it->second;
		
				const long sizeInd = index.size();
		
				it->second = index.size();
		
				index.push_back(0);
				index.push_back(0);
		
				long count = 0;
				long bytes = 0;
	
				for (long x = index0.next(x1); x != -1; x = index0.next(x+1)){
		
					const long ind = index0.get(x);
					index.push_back( ind);
				
					bytes += rows.at(ind).getBytes();
					++count;
				}
				
				index.at(sizeInd) = count;
				index.at(sizeInd+1) = bytes;
			}
		}
		
		
		long joinSize = 0;
		long joinBytes = 0;
		
		{
			Timer t("join size");
		
			if (onlyJoinSize == 0)
				t.start(10, 100000, size2);
		
			for (auto it = b.begin(); it != b.end(); ++it){
			
				const Row& extrow = (*it);
				
				if (onlyJoinSize == 0)
					t.tick();
			
				if (!index1.hasKey(extrow, extproj, exthash))
					continue;
				
				const long ind = index1.get(extrow, extproj, exthash);
				const long size = index.at(ind);
				const long bytes = index.at(ind+1);
			
				joinSize += size;
				joinBytes += bytes+extrow.getBytes()*size;
			}
			
			if (onlyJoinSize == 0)
				t.end();
		}
		
		
		if (onlyJoinSize == 1){
		
			onlyJoinSize = joinSize;
			return nullptr;
		}
		
		if (onlyJoinSize == 0)
			cout << "joinSize " << joinSize << " joinBytes " << joinBytes << endl;
		
		if (onlyJoinSize == 0)
			cout << retCols << " joinSize " << joinSize << " joinBytes " << joinBytes << endl;
		
		auto ret = std::make_shared<DatabaseTable>(name.str(), retCols);
		
		ret->clear();
		ret->reserveBytes(joinBytes);
		
		Rows& retrows = *ret->getRows();
		
		{
			Timer t("join");
		
			if (onlyJoinSize == 0)
				t.start(10, 100000, size2);
			
			long k = 0;
		
			for (auto it = b.begin(); it != b.end(); ++it){
			
				t.tick();
			
				const Row& extrow = (*it);
			
				if (!index1.hasKey(extrow, extproj, exthash))
					continue;
				
				const long ind = index1.get(extrow, extproj, exthash);
				const long size = index.at(ind);
				
				const long first = ind+2;
				const long last = ind+1+size;
				
				for (long j = first; j <= last; ++j){
			
					const Row& tabrow = rows.at( index.at(j) );
			
					if (!flip){
						
						retrows.push_back(tabrow, extrow, novel);
				
				
					} else {
				
						retrows.push_back(extrow, tabrow, novel);
						
					}
					
				}
				
				k++;
			}
			
			if (onlyJoinSize == 0)
				t.end();
		}
		
		retrows.shrink_to_fit();
		
		return ret;
		
	}
	
	
	
	inline shared_ptr<DatabaseTable> getJoin( DatabaseTable& a, Tuple joinA, DatabaseTable& b, Tuple joinB, long size1=-1, long bytes1 = -1, long size2=-1, long bytes2=-1, bool flip=false){
	
		long onlyJoinSize = 0;
	
		return getJoin(onlyJoinSize, a,joinA, b,joinB, size1,bytes1,size2,bytes2,flip);
	}
	
	inline shared_ptr<DatabaseTable> getJoin( DatabaseTable& a, DatabaseTable& b, long size1=-1, long bytes1 = -1, long size2=-1, long bytes2=-1, bool flip=false){
	
		return getJoin(a, a.getColNames(), b, b.getColNames(), size1,bytes1,size2,bytes2,flip);
	}

	inline long getJoinSize( DatabaseTable& a, Tuple joinA, DatabaseTable& b, Tuple joinB, long size1=-1, long bytes1 = -1, long size2=-1, long bytes2=-1, bool flip=false){
	
		long onlyJoinSize = 1;
		
		auto x = getJoin(onlyJoinSize, a,joinA, b,joinB, size1,bytes1,size2,bytes2,flip);
		
		return onlyJoinSize;
	}
	
	
	inline long getJoinSize( DatabaseTable& a, DatabaseTable& b, long size1=-1, long bytes1 = -1, long size2=-1, long bytes2=-1, bool flip=false){
	
		return getJoinSize(a, a.getColNames(), b, b.getColNames(), size1,bytes1,size2,bytes2,flip);
	}
	
};








/*
class SQLQuery{

public:

	shared_ptr<ExpressionParser> parsed;
	
	std::map<string, string> tableAliases;
	std::vector<string> selectedTables;
	std::vector<string> selectedTableNames;
	
	std::set<string> selectedColumns;
	std::set<string> joinColumns;
	
	
	vector<string> values;
	std::map<string, long> mapped;
	
	
	std::map<string, string> variableMappings;
	
	inline shared_ptr<ExpressionParser> getParsed() const{
	
		return parsed;
	}
	
	inline bool isSelected(const string& table) const{
		
		return std::find(selectedTables.begin(), selectedTables.end(), table) != selectedTables.end();
	}
	
	inline void addTable(const string& table, const Tuple& colNames){
	
		if (isSelected(table) ){
			
			assert (tableAliases.count(table) > 0);
			
			const string repl = tableAliases.at(table);
			
		
			for (int j = 0; j < colNames.size(); ++j){
			
				const string targ = tableWith(colNames[j]);
				
				
				joinColumns = UtilString::replaceAllIn( joinColumns, targ, repl);
				selectedColumns = UtilString::replaceAllIn( selectedColumns, targ, repl);
				tableAliases = UtilString::replaceAllIn(tableAliases, targ, repl);
				
				mapped = UtilString::replaceAllIn(mapped, targ, repl);
				values = UtilString::replaceAllIn(values, targ, repl);
			}
		}
	}
	
	
	inline const string getAlias(const string& table) const{
	
		assert (tableAliases.count(table) > 0);
	
		return tableAliases.at(table);
	}
	
	
	inline bool isSelected(const string& tableAlias, const string& column) const{
		
		if (selectedColumns.count( combine(tableAlias, column) ) > 0)
			return true;
		
		return false;
	}
	
	inline bool isIncluded(const string& tableAlias, const string& column) const{
		
		if (mapped.count( combine(tableAlias, column) ) > 0)
			return true;
		
		return false;
	}
	
	inline const string onlyTable(const string& s) const{
	
		return UtilString::getSplit(s, ".", 0);
	}
	
	inline const string tableWith(const string& colname) const {
	
		return "?"+colname;
	}
	
	inline const string onlyColumn(const string& s) const{
	
		return UtilString::getSplit(s, ".", 1);
	}
	
	inline const string combine(const string& table , const string& column) const{
	
		return table+"."+column;
	}
	
	inline const string getAlias(const string& tableAlias, const string& column) const{
		
		const string key = combine(tableAlias, column);
		
		assert (mapped.count(key) > 0);
		
		const string key2 = values.at(mapped.at(key));
		
		const string table2 = onlyTable(key2);
		const string column2 = onlyColumn(key2);
		
		const string ret = UtilString::equals(column2, column) ? combine( table2, column2) : column2;
		
		cout << "getAlias " << tableAlias << "." << column << " = " << ret<< endl;
		
		return ret;
	}
	
	inline bool isJoinColumn(const string& tableAlias, const string& column) const{
		
		if (joinColumns.count( combine(tableAlias, column) ) > 0)
			return true;
		
		return false;
	}
	
	inline int get(const Expression& e, int child=-1, int grandchild=-1, int grandgrandchild = -1) const{
	
		if (child == -1)
			return e.mathOp;
		else
			return e.childrenNum <= child ? -1 : get(parsed->expressions.at( e.children[child] ), grandchild, grandgrandchild);
			
	}
	
	inline int& at(Expression& e, int child=-1, int grandchild=-1, int grandgrandchild = -1){
	
		if (child == -1)
			return e.mathOp;
		else
			return at(parsed->expressions.at( e.children[child] ), grandchild, grandgrandchild);
	}
	
	inline Expression& getChild(Expression& e, int child){
	
		return parsed->expressions.at( e.children[child] );
	}
	
	long limit = 0;
	
	
	inline SQLQuery(const string& sql){
	
		parsed = std::make_shared<ExpressionParser>();
		ExpressionParser& p = (*parsed);
		
		
		
		
		p.addBasicOps();
		
		
		
		p.addOp("view", "Create Or Replace View", "", ";", 1, 999);
		
		p.addOp("as", "", "As", "", 2, 2);
		
		p.addOp("limit", "Limit", "", "", 1, 1);
		
		//p.addOp("and", "", "And", "", 2, 999);

		//p.addOp("comma", "", ",", "", 2, 999);
		
		p.addOp("rand", "Rand(", "", ")", 0, 0, UtilExpression::NON_MATH);
		
		p.addOp("weighted", "Order By", "", "", 1, 1, UtilExpression::NON_MATH);
		
		p.addOp("select", "Select *", "", "", 0, 0, UtilExpression::NON_MATH);
		
		p.addOp("select", "As Select", ",", "", 1, 999, UtilExpression::NON_MATH);
		
		p.addOp("select", "Select", ",", "", 1, 999, UtilExpression::NON_MATH);
		p.addOp("from", "From", ",", "", 1, 999, UtilExpression::NON_MATH);
		p.addOp("where", "Where", "And", "", 1, 999, UtilExpression::LENIENT_MUL);
		
		p.addOp("join", "Inner Join", ",", "", 1, 999, UtilExpression::NON_MATH);
		p.addOp("join", "Join", ",", "", 1, 999, UtilExpression::NON_MATH);
		
		
		p.addOp("option", "Option(", ",", ")", 1, 999, UtilExpression::NON_MATH);
		
		p.addOp("option", "Option", ",", "", 1, 999, UtilExpression::NON_MATH);
		
		
		
		p.addOp("weighted", "Weighted By", ",", "", 1, 999, UtilExpression::LENIENT_MUL );
		
		p.addOp("use", "Use Database", "With", ";", 2, 2, UtilExpression::NON_MATH);
		p.addOp("use", "Use", "With", ";", 2, 2, UtilExpression::NON_MATH);
		
		p.addOp("use", "Use Database", "", ";", 1, 1, UtilExpression::NON_MATH);
		p.addOp("use", "Use", "", ";", 1, 1, UtilExpression::NON_MATH);
		
		p.addOp("use", "Use", "", ";", 1, 1, UtilExpression::NON_MATH);
		
		p.addOp("on", "On", ",", "", 1, 999, UtilExpression::NON_MATH);
		*/
		
		//p.addOp("comment", "/*", ",", "*/", 1, 999, UtilExpression::NON_MATH);
		
		/*
		p.addOp("list", "", "", ";", 1, 999, UtilExpression::LENIENT_MUL);
		
		p.addOp("list", "", "", "", 2, 999, UtilExpression::LENIENT_MUL);
		
		
		
		
		
		//p.addOp("list", "", "", "", 2, 999, UtilExpression::LENIENT_MUL);
		
		//p.verbose = true;
	
		const bool check = p.parse(sql);
		
		assert (check);
		
		
		
		
		
		
		//cout << p.toString() << endl;
		
		const Expression& forest = p.expressions.back();
		
		
		
		for (auto it = forest.children.begin(); it != forest.children.end(); ++it){
		
			Expression& e = p.expressions.at(*it);
			
			if (p.isOp(e, "limit")){
			
				if (e.childrenNum == 1){
					
					const string s = e.childrenPtrs[0]->s;
			
					limit = UtilHash::to_longdouble(s.c_str());
				}
			}
			
			if (p.isOp(e, "weighted")){
			
				if (e.children.size() == 1){
					
					at(e) = UtilExpression::ONE;
					
					
					if (get(e, 0) == UtilExpression::DIV && get(e, 0, 1) == UtilExpression::LOG && get(e, 0, 1, 0) == UtilExpression::NON_MATH ){
					
						at(e) = UtilExpression::LENIENT_IDENTITY0;
					
						at(e, 0) = UtilExpression::LENIENT_IDENTITY0;
						
						if (get(e, 0, 0) == UtilExpression::MUL)
							at(e, 0, 0) = UtilExpression::LENIENT_MUL;
						
						
						if (get(e, 0, 0) == UtilExpression::DIV && get(e,0,0,0) == UtilExpression::MUL)
							at(e, 0,0,0) = UtilExpression::LENIENT_MUL;	
						
						if (get(e, 0, 0) == UtilExpression::DIV && get(e,0,0,1) == UtilExpression::MUL)
							at(e, 0,0,1) = UtilExpression::LENIENT_MUL;
													
					}
					
					if (get(e, 0) == UtilExpression::POW && get(e, 0, 0) == UtilExpression::NON_MATH ){
					
						at(e) = UtilExpression::LENIENT_IDENTITY0;
						
						if (get(e, 0, 1) == UtilExpression::DIV){
							
							at(e, 0) = UtilExpression::LENIENT_IDENTITY1;
							at(e, 0, 1) = UtilExpression::FLIPDIV;
							
						} else {
							
							at(e, 0) = UtilExpression::LENIENT_INVERSE1;
						
							if (at(e, 0, 1) == UtilExpression::MUL)
								at(e, 0, 1) = UtilExpression::LENIENT_MUL;
						
						}
					}
					
					
				}
			}
		}
		
		
		for (long j = 0; j < p.expressions.size(); ++j){
		
			Expression& e = p.expressions.at(j);
		
			if ( get(e, 0) == UtilExpression::LOG && get(e, 0, 0) == UtilExpression::EXP){
			
				at(e, 0) = UtilExpression::IDENTITY;
				at(e, 0,0) = UtilExpression::IDENTITY;
			}
			
			if ( get(e, 0) == UtilExpression::EXP && get(e, 0, 0) == UtilExpression::LOG){
			
				at(e, 0) = UtilExpression::IDENTITY;
				at(e, 0,0) = UtilExpression::IDENTITY;
			}
		}
		
		for (auto it = forest.children.begin(); it != forest.children.end(); ++it){
		
			const Expression& root = p.expressions.at(*it);
		
			if (p.isOp(root, "from") ||p.isOp(root, "join")){
			
				for (auto it = root.children.begin(); it != root.children.end(); ++it){
				
					const Expression& selectChild = p.expressions.at(*it);
					
					const string s = UtilString::cleanStringList(selectChild.s, ' ', ' ');
					
					int m = UtilString::getSplitLength(s, " ");
					
					assert (m >= 1);
					assert (m <= 2);
					
					const string key = m == 2 ? UtilString::getSplit(s, " ", 0) : s;
					const string val = m == 2 ? UtilString::getSplit(s, " ", 1) : s;
					
					tableAliases.insert( std::make_pair(key, val));
					
					if (m == 2)
						tableAliases.insert( std::make_pair(val, val));
					
					selectedTables.push_back(key);
					selectedTableNames.push_back(val);
					
				}
			}
			
			if (p.isOp(root, "select")){
			
				const Expression& select = root;
			
				for (auto it = select.children.begin(); it != select.children.end(); ++it){
				
					const Expression& selectChild = p.expressions.at(*it);
					
					string table;
					string column;
					string tablecolumn;
					string alias;
						
					if (p.isOp(selectChild, "as")){
					
						tablecolumn = p.expressions.at(selectChild.children.at(0)).s;
						alias = p.expressions.at(selectChild.children.at(1)).s;
						
					} else {
					
						tablecolumn = selectChild.s;
						alias = selectChild.s;
					}
					
					if (UtilString::getSplitLength(tablecolumn, ".") == 2){
					
						table = UtilString::getSplit(tablecolumn, ".", 0);
						column = UtilString::getSplit(tablecolumn, ".", 1);
					
					} else {
						
						column = tablecolumn;
						table = tableWith(column);
					}
					
					const string key = table+"."+column;
					const string val = table+"."+alias;
					
					
					
					parsed->renameVariable(key, val);
					parsed->renameVariable(column, alias);
					
					if (mapped.count(key) == 0){
					
						values.push_back(val);
						mapped.insert( std::make_pair( key, values.size()-1));
					}
					
					
					if (!UtilString::equals(alias, tablecolumn))
						values.at(mapped.at(key) ) = val;
					
					selectedColumns.insert(key);					
				}
			}
			
			if (p.isOp(root, "where") || p.isOp(root, "on") ){
			
				for (auto it = root.children.begin(); it != root.children.end(); ++it){
			
					const Expression& whereChild = p.expressions.at(*it);
					
					if (p.isOp(whereChild, "equal")){
					
						const string tablecolumn1 = p.expressions.at(whereChild.children.at(0)).s;
						const string tablecolumn2 = p.expressions.at(whereChild.children.at(1)).s;
						
						string table1, table2, column1, column2;
						
						if (UtilString::getSplitLength(tablecolumn1, ".") == 2){
					
							table1 = UtilString::getSplit(tablecolumn1, ".", 0);
							column1 = UtilString::getSplit(tablecolumn1, ".", 1);
					
						} else {
							
							column1 = tablecolumn1;
							table1 = tableWith(column1);
						}
						
						if (UtilString::getSplitLength(tablecolumn2, ".") == 2){
					
							table2 = tableAliases.at(UtilString::getSplit(tablecolumn2, ".", 0));
							column2 = UtilString::getSplit(tablecolumn2, ".", 1);
					
						} else {
						
							column2 = tablecolumn2;
							table2 = tableWith(column2);
						}
						
						const string key1 = combine(table1, column1);
						const string key2 = combine(table2, column2);
						
						joinColumns.insert(key1);
						joinColumns.insert(key2);
						
						const long ind1 = mapped.count(key1) > 0 ? mapped.at(key1) : -1;
						const long ind2 = mapped.count(key2) > 0 ? mapped.at(key2) : -1;
						
						if (ind2 >= 0 && ind1 == -1){
						
							mapped.insert( std::make_pair(key1, ind2) );
							
						} else if (ind1 >= 0 && ind2 == -1){
						
							mapped.insert( std::make_pair(key2, ind1) );
							
						} else if (ind1 >= 0 && ind2 >= 0){
						
							mapped.at(key1) = ind2;
							
						} else {
						
							values.push_back(key1);
							
							mapped.insert( std::make_pair(key1, values.size()-1) );
							mapped.insert( std::make_pair(key2, values.size()-1) );
							
						}	
					}
				}
			}
		}
		
		
	}
	
	
	inline vector<string> getFrom() const{
	
		vector<string> ret;
		
		const ExpressionParser& p = (*parsed);
		
		const Expression& forest = p.expressions.back();
		
		for (auto it = forest.children.begin(); it != forest.children.end(); ++it){
		
			const Expression& root = p.expressions.at(*it);
		
			if (p.isOp(root, "from") ){
			
				for (auto it = root.childrenPtrs.begin(); it != root.childrenPtrs.end(); ++it){
				
					const Expression& from = *(*it);
					
					ret.push_back( from.s);
				}
				
				
			}
		}
		return ret;
	}
	
	inline void print(){
		
		cout << parsed->repr() << endl;
		
		cout << "limit: " << limit << endl;
		
		cout << "join tables: " << endl;
		
		for (auto it = selectedTables.begin(); it != selectedTables.end(); ++it){
		
			const string tableName = (*it);
			
			const string alias = selectedTableNames.at(it-selectedTables.begin());
		
			if (!UtilString::equals(tableName, alias))
				cout << "\t" <<tableName << " renamed to " << alias << endl;
			else
				cout << "\t" <<tableName << endl;
		}
		
		
		for (auto it = mapped.begin(); it != mapped.end(); ++it){
			
			if (it->first[0] == '?')
				cout << "find table " << UtilString::getSplit(it->first, ".", 0) <<" that has column " << UtilString::getSplit(it->first, ".", 1) << endl;
		}
		
		for (auto it = mapped.begin(); it != mapped.end(); ++it){
		
			cout << "rename " << it->first << " to " << values.at(it->second) << endl;
		}
		
		cout << "join along: " << endl;
		
		
		for (auto it = joinColumns.begin(); it != joinColumns.end(); ++it){
		
			cout << "\t" << (*it) << endl;
		}
	}
};
*/


class SQLExpression{
	
	class Column{
	public:
		string table;
		string column;
	
		int scope;

	
		inline Column(const string& table, const string& column, int scope) : table(table), column(column), scope(scope){
		
		}
	
		inline Column(const Column& c) : table(c.table), column(c.column), scope(c.scope){
	
		}
	
		/*
		inline const string getKey() const{
	
			return table+"."+column;
		}*/
		
		inline const string getKey() const{
	
			if (table.length() > 0 && table[0] == '?')
				return column;
				
			return table+"."+column;
		}
	
		
		bool operator<(const Column& col) const {
	
			return scope < col.scope;
		}
	};

	class ColumnGroup{

	public:
		vector<Column> list;
	
		bool join = false; // featured in join conditions
	
		shared_ptr< map<string, ColumnGroup* > > indexPtr;
	
		map<string, ColumnGroup* >& index;
	
		ColumnGroup* owner;
	
		inline ColumnGroup(shared_ptr< map<string, ColumnGroup*> > indexPtr) : indexPtr(indexPtr), index(*indexPtr){
	
			owner = this;
		}
	
		inline void print(){
	
			if (list.size() == 0)
				return;

		
			assert (owner == this);
		
			for (auto it = list.begin(); it != list.end(); ++it){
		
				
		
				cout << it->getKey() << " ";
			}
		
			cout << endl;
		}
	
		inline void add(const Column& col){
	
			if (owner != this)
				return owner->add(col);
				
			const string key = col.getKey();
		
			if (index.count(key) > 0){
			
				ColumnGroup* currentOwner = index[key];
				
				if (currentOwner == this)
					return;
					
				grab (*currentOwner);
				return;
			}
		
			/*
			bool b_add = true;
			
			for (auto it = list.begin(); it != list.end(); ++it){
				if ( UtilString::equals( it->getKey(), key)){
					b_add = false;
					break;
				}
			}
		
		
			if (b_add)
				list.push_back(col);*/
			
			list.push_back(col);
			index[ col.getKey() ] = this;
			
		}
	
		inline void add(const string& view, const string& s, const string& alias, int scope){
	
			if (UtilString::getSplitLength(s, ".") == 2){
		
				add( Column(UtilString::getSplit(s, ".", 0), UtilString::getSplit(s, ".", 1), scope) );
			
				if (view.length() > 0){
					add( Column("?", alias, scope+10) );
					//add( Column("?", alias, scope+11) );
					
				} else{
					add( Column(UtilString::getSplit(s, ".", 0), alias, scope+10) );
				}
			
			} else {
		
				add( Column("?"+s, s, scope) );
				
				if (view.length() > 0){
					add( Column("?", alias, scope+10) );
					
					//add( Column("?", alias, scope+11) );
					
					
				} else {
					add( Column("??"+s, alias, scope+10) );
				}
			}
		
		
		}
	
		inline void addTable(const string& table, const Tuple& t){
	
	
			long siz = t.size();
	
			for (long j = 0; j < siz; ++j){
		
				const string tj = t[j];
		
				for (long k = 0; k < list.size(); ++k){
			
					const string oldKey = list.at(k).getKey();
				
					if ( UtilString::equals(list.at(k).table, "?"+tj)){
				
					
						add( Column(table, tj, list.at(k).scope) );
						break;
					}
				}
			
			}
		}
	
		inline void add(const string& view, const string& s, int scope){
	
			if (UtilString::getSplitLength(s, ".") == 2){
		
				add( Column(UtilString::getSplit(s, ".", 0), UtilString::getSplit(s, ".", 1), scope) );
			
				if (view.length() > 0){
					//add( Column( view, UtilString::getSplit(s, ".", 1), scope+1 ) );
					
					//add( Column( "?", UtilString::getSplit(s, ".", 1), scope+2 ) );
				}
				
			} else {
		
				add( Column("?"+s, s, scope) );
			
				if (view.length() > 0){
					//add( Column( view, s, scope+1 ) );
					
					//add( Column( "?", s, scope+2 ) );
				}
			}
		
		}
	
		inline void sort(){
	
			std::sort(list.begin(), list.end() );
		}
	
		inline void grab(ColumnGroup& group){
		
			assert (owner == this);
		
			for (auto it = group.list.begin(); it != group.list.end(); ++it){		
				
				index[ (*it).getKey() ] = this;
				list.push_back( (*it) );
			}
			
			group.owner = this;
			group.list.clear();
		
		}
	};

	map<string, vector<string>> tableAliases;
	vector<shared_ptr<ColumnGroup>> groups;
	
	std::map<string, Tuple> viewAttributes;
	
	shared_ptr<map<string, ColumnGroup* >> index;

	
public:

	shared_ptr<SQLParser> parsed;
	
	const string raw;

	inline map<string, vector<string>> getAliases(){
	
		map<string, vector<string>> ret;
	
		for (auto it = index->begin(); it != index->end(); ++it){
		
			const ColumnGroup& g = *(*it).second;
		
			if (g.list.size() == 0)
				continue;
			
			vector<string> v;
			
			for (auto it2 = g.list.begin(); it2 != g.list.end(); ++it2)
				v.push_back( (*it2).getKey() );
			
			ret.insert( {(*it).first, v});
		}
		
		return ret;
	}
	
	long limit = 0;

	inline SQLExpression(const string raw) : raw(raw){
	
		parsed = std::make_shared<SQLParser>();
		
		SQLParser& p = (*parsed);
		
		p.parse(raw);
		
		index = std::make_shared< map<string, ColumnGroup* > >();
		
		std::set<string> viewAttributesSet;
		
		for (auto it = p.expressions.begin(); it != p.expressions.end(); ++it){
		
			const Expression& e = (*it);
		
			const ExpressionOperator& op = p.ops.at(e.op);
			
			if (UtilString::equals(op.name, "limit")){
			
				limit = UtilString::stringToLongDouble(e.childrenPtrs.at(0)->s);
				continue;
			}
			
			if (UtilString::equals(op.name, "var")){
			
				const Expression* select = p.getAncestor(e,"select");
				const Expression* from = p.getAncestor(e,"from");
				
				const Expression* join = p.getAncestor(e,"join");
				
				const Expression* where = p.getAncestor(e,"where");
				
				const Expression* on = p.getAncestor(e,"on");
				
				const Expression* weighted = p.getAncestor(e,"weighted");
				const Expression* viewex = p.getAncestor(e,"view");
				
				if (where == NULL)
					where = on;
					
				if (from == NULL)
					from = join;
				
				if (UtilMath::equalsAll( (const Expression*) NULL, select, from, where, weighted) )
					continue;
				
				string view;
				
				int scope = 0;
				
				if (viewex){
					view = viewex->childrenPtrs.at(0)->s;
										
					if (viewAttributes.count(view) == 0)
						viewAttributes.insert( {view, Tuple()} );
				}
				
				if (select){
					
					auto group = std::make_shared<ColumnGroup>(index);
					
					groups.push_back(group);
					
					if (e.parent == select->id){
						
						group->add( view, e.s, scope);
						
						if (view.length() > 0){
							
							const string attr = UtilString::getSplit(e.s, ".", -1);
									
							if (viewAttributesSet.count(view+"."+attr) == 0){
							
								viewAttributes[view].push_back( attr );
								viewAttributesSet.insert(view+"."+attr);
							}
						}
						
					} else {
					
						const Expression* as = p.getAncestor(e, "as");
					
						if (as){
					
							if (e.parent == as->id){
							
								//cout << as->childrenPtrs.at(0)->s << " AS " << as->childrenPtrs.at(1)->s << endl;
							
								group->add( view, as->childrenPtrs.at(0)->s, as->childrenPtrs.at(1)->s, scope);
								
								
								if (view.length() > 0){
								
									if (viewAttributes.count(view) == 0)
										viewAttributes.insert( {view, Tuple() } );
									
									const string attr = UtilString::getSplit(as->childrenPtrs.at(1)->s, ".", -1);
									
									if (viewAttributesSet.count(view+"."+attr) == 0){
									
										viewAttributes[view].push_back( attr );
										viewAttributesSet.insert(view+"."+attr);
									}
								}
							}
						}
				
					}	
				}
			
				if (from){
				
					if (e.parent == from->id){
						
						const string s = UtilString::cleanStringList(e.s,' ', ' ');
					
						const string key = UtilString::getSplitLength(s, " ") == 2 ? UtilString::getSplit(s, " ", 0) : s;
						const string alias = UtilString::getSplitLength(s, " ") == 2 ? UtilString::getSplit(s, " ", 1) : s;
					
						if (tableAliases.count(key) == 0)
							tableAliases.insert( {key, vector<string>() } );
					
						tableAliases[key].push_back(alias);
						
						baseTables.insert(alias);
					}
				
				}
				
				if (where){
				
					const Expression* equal = p.getAncestor(e, "equal");
				
					if (equal){
					
						if (e.parent == equal->id){
						
							auto group = std::make_shared<ColumnGroup>(index);
					
							groups.push_back(group);
							
							group->add( view,  equal->childrenPtrs.at(0)->s, scope);
							group->add( view,  equal->childrenPtrs.at(1)->s, scope);
							
							
						}
					
					}
				
				}
				
				if (weighted){
				
					auto group = std::make_shared<ColumnGroup>(index);
					
					groups.push_back(group);
					group->add( view, e.s, scope);
					
				}
			}
			
		}
		
		for (auto it = viewAttributes.begin(); it != viewAttributes.end(); ++it){
		
			for (auto it2 = groups.begin(); it2 != groups.end(); ++it2){
			
				(*it2)->addTable(it->first, it->second);
			}
		
		}
	}
	
	std::set<string> baseTables;
	
	inline void addTable(const string& tab, const Tuple& t){
	
		//cout << "addTable(" << tab << ", " << t << ")" << endl;
		
		
		
		for (auto it2 = groups.begin(); it2 != groups.end(); ++it2)
			(*it2)->addTable(tab, t);
			
		const vector<string> als = tableAliases.at(tab);
		
		for (long j = 0; j < als.size(); ++j){
		
			for (auto it2 = groups.begin(); it2 != groups.end(); ++it2)
				(*it2)->addTable(als.at(j), t);
		}
		
	}
	
	
	inline void finalise(){
		
		for (auto it = groups.begin(); it != groups.end(); ++it)
			(*it)->sort();
		
		const SQLParser& p = (*parsed);
		
		for (auto it = p.expressions.begin(); it != p.expressions.end(); ++it){
		
			const Expression& e = (*it);
		
			const ExpressionOperator& op = p.ops.at(e.op);
			
			if (UtilString::equals(op.name, "var")){
				
				const Expression* where = p.getAncestor(e,"where");
				
				const Expression* on = p.getAncestor(e,"on");
				
				if (where == NULL)
					where = on;
				
				if (where){
				
					const Expression* equal = p.getAncestor(e, "equal");
					
					if (equal){
					
						if (e.parent == equal->id){
							
							ColumnGroup* g = index->at(e.s);
							
							
							if (where == on){
							
								g->join = true;
								
							} else {
							
								std::set<string> diffBase;
								
								for (auto it = g->list.begin(); it != g->list.end(); ++it){
								
									if (baseTables.count( (*it).table) > 0)
										diffBase.insert((*it).table);
										
								}
								
								g->join = diffBase.size() > 1;
							
							}
							
						}
					}
				
				}
			}
			
		}
	}
	
	inline bool selectedColumn(const string& s) const{
	
		return index->count(s) > 0;
	}
	
	inline bool joinColumn(const string& s) const{
	
		if (!selectedColumn(s))
			return false;
	
		return index->at(s)->join;
	}
	
	inline const string renamedColumn(const string& s) const{
	
		return index->at(s)->list.back().getKey();
	}
	
	inline void print() const{
	
		cout << endl;
		cout << raw << endl << endl;
		
		cout << parsed->repr() << endl;
	
		for (auto it = tableAliases.begin(); it != tableAliases.end(); ++it){
		
			cout << it->first << " => " << UtilString::toString(it->second) << endl;
		}
		
		for (auto it = groups.begin(); it != groups.end(); ++it){
		
			(*it)->print();
		}
		cout << endl;
	}
	

};

#endif