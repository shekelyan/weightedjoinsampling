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

#ifndef SHEKELYAN_UTILCOLLECTION_H
#define SHEKELYAN_UTILCOLLECTION_H

#include <approxdata/utils/utils.h>


class stringvector{

private:
	vector<string> data;

public:
	inline stringvector(){
	
	}
	
	inline void push_back(const string& s){
		
		data.push_back( s);
	}
	
	
	inline stringvector( const stringvector& v){
	
		for (int i = 0; i < v.size(); i++)
			push_back( v[i] );
	}

	
	inline stringvector& operator=(const stringvector& v){
	
		if (this != &v){
		
			clear();
		
			for (int i = 0; i < v.size(); i++)
				this->push_back( v[i] );
		}
		
		return *this;
	}
	
	inline long size() const{
	
		return data.size();
	}

	inline const string operator[](std::size_t idx) const{
		
		return data.at(idx);
	}
	
	inline const string at(std::size_t idx) const{
		
		return data.at(idx);
	}
	
	inline void clear(){
	
		data.clear();
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		for (int i = 0; i < data.size(); i++){
			
			if (i > 0)
				s << ",";
			
			s << data[i];
			
		}
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const stringvector& v){ 
    	
    	os << v.toString();
    	
    	return os;
	}
};

template<typename E>
class constvector{

private:
	//vector<shared_ptr<const E>> data;
	
	unique_ptr<vector<E>> data;
	
public:	
	
	inline constvector(){
	
		clear();
	}
	
	inline constvector( const E& t1){
	
		clear();
		push_back(t1);
	}
	
	inline constvector( const E& t1, const E& t2){
	
		clear();
		push_back(t1);
		push_back(t2);
	}
	
	inline constvector( const E& t1, const E& t2, const E& t3){
	
		clear();
		push_back(t1);
		push_back(t2);
		push_back(t3);
	}
	
	inline constvector( const E& t1, const E& t2, const E& t3, const E& t4){
	
		clear();
		push_back(t1);
		push_back(t2);
		push_back(t3);
		push_back(t4);
	}
	
	inline constvector( const E& t1, const E& t2, const E& t3, const E& t4, const E& t5){
	
		clear();
		push_back(t1);
		push_back(t2);
		push_back(t3);
		push_back(t4);
		push_back(t5);
	}
	
	inline constvector<E>( const constvector<E>& v){
	
		clear();
		for (int i = 0; i < v.size(); i++){
		
			if (i >= v.size())
				break;
		
			push_back( v[i] );
		}
	}
	
	inline constvector<E>& operator=(const constvector<E>& v){
	
		if (this != &v){
		
			clear();
		
			for (int i = 0; i < v.size(); i++)
				push_back( v[i] );
		}
		
		return *this;
	}
	
	inline void push_back(const constvector<E>& v){
		
		for (int i = 0; i < v.size(); i++)
			push_back(v[i]);
	}
	
	inline void push_front(const constvector<E>& v){
	
		const constvector<E> t = (*this);
		
		clear();
		
		push_back(v);
		push_back(t);
	}
	
	inline void push_back(const E& s){
	
		data->push_back(s);
		//shared_ptr<const E> ptr(new E(s));	
		//data.push_back( std::move(ptr) );
	}
	
	inline long size() const{
	
		return data->size();
	}

	inline const E operator[](std::size_t idx) const{
		
		assert (idx >= 0);
		assert (idx < size() );
		
		return data->at(idx); //(*data.at(idx));
	}
	
	inline const E at(std::size_t idx) const{
		
		return data->at(idx); //(*data.at(idx));
	}
	
	inline void clear(){
	
		data.reset(new vector<E>() );
	}
	
	inline const string toString() const { 
    	
    	stringstream os;
    	
    	for (int i = 0; i < size(); i++){
			
			if (i > 0)
				os << ",";
			
			os << (*this)[i];
		}
		
		return os.str();
	}
	
	template <typename T>
	friend inline std::ostream& operator<<(std::ostream &os, const constvector<T>& v) { 
    	
    	for (int i = 0; i < v.size(); i++){
			
			if (i > 0)
				os << ",";
			
			os << v[i];
		}
		
		return os;
	}
};	


class ReservoirSampling{

	private:
		long k;
		
		std::unique_ptr<std::mt19937> g;
		std::seed_seq seed;

	public:
		const long sampleSize;
		
		inline ReservoirSampling(string seedStr, long size) : seed(seedStr.begin(), seedStr.end() ), sampleSize(size){
		
			g.reset(new std::mt19937(seed));
			this->k = 0;
		}
		
		inline void reset(){
	
			k = 0;
		}
		
		inline long pos(){
		
			if (k < sampleSize){
			
				const long j = k;
				++k;
				return j;
			}
			
			std::uniform_int_distribution<long> d(0, k);
			const long j = d( (*g) );
			
			++k;
			return j < sampleSize ? j : -1;
		} 

};

namespace UtilCollection{

	template<typename T>
	inline T& tableCell(vector<vector<T>>& data, int row, int col){
	
		while (data.size() <= row){
	
			vector<T> v;
		
			while (v.size() <= col){
		
				T vv;
				v.push_back(vv);
			}
				
			data.push_back(v);
		}
	
		vector<T>& v = data.at(row);
	
		while (v.size() <= col){
		
			T vv;
			v.push_back(vv);
		}
	
		return v.at(col);
	}


	template<typename T>
	inline const T& tableCell(const vector<vector<T>>& data, int row, int col){
	
		return data.at(row).at(col);
	}

	template<typename E, typename F = std::less<E>>
	inline bool less(E a, E b, F comp=F()){
	
		return comp(a,b);
	}

	template<typename E, typename F = std::less<E>>
	inline bool greater(E a, E b, F comp=F()){
	
		return comp(b,a);
	}

	template<typename E, typename F = std::less<E>>
	inline bool lessOrEqual(E a, E b, F comp=F()){
	
		return !greater(a,b,comp);
	}
	
	template<typename E, typename F = std::less<E>>
	inline bool greaterOrEqual(E a, E b, F comp=F()){
	
		return !less(a,b,comp);
	}
	
	
	
	template<typename E>
	
	inline void remove(vector<E>& v, E elem){
	
	
	}
	
	template<typename E, typename F = std::less<E>>
	inline bool equal(E a, E b, F comp=F()){
	
		return greaterOrEqual(a,b, comp) && greaterOrEqual(b,a, comp);
	}


	template<typename E, typename F = std::less<E>>
	inline long countLessEqualNaive(const vector<E>& v, const E& p, F comp=F() , const long aa = 0, const long bb = -1){
	
		long ret = 0;
		
		for (auto it = v.begin(); it != v.end(); it++)
			if ( !comp( p, (*it) ) )  // !(p < it) <=> p >= it <=> it <= p
				ret++;
		
		return ret;
	}
	
	
	template<typename E, typename F = std::less<E>>
	inline long countLessNaive(const vector<E>& v, const E& p, F comp=F() , const long aa = 0, const long bb = -1){
	
		long ret = 0;
		
		for (auto it = v.begin(); it != v.end(); it++)
			
			if (comp( (*it), p))
				ret++;
		
		return ret;
	}
	
	template<typename E, typename F = std::less<E>>
	inline long countLess(const vector<E>& v, const E& p, F comp=F(), const long aa = 0, const long bb = -1){
		
		const long a = aa;
		const long b = bb == -1 ? v.size() : bb;
		
		assert (b >= a);
		
		if (a == b)
			return a;
		
		auto s = v.begin()+a;
		auto e = v.begin()+b;
		
		return std::distance(v.begin(), std::lower_bound(s, e, p, comp));
	}
	
	template<typename E>
	inline bool equals(const vector<E>& a, const vector<E>& b){
	
		if (a.size() != b.size() )
			return false;
			
		for (int i = 0; i < a.size(); ++i){
		
			if (a[i] != b[i])
				return false;
		}
	
		return true;
	}
	
	
	template<typename E, typename F = std::less<E>>
	inline long countLessEqual(const vector<E>& v, const E& p, F comp=F(), const long aa = 0, const long bb = -1){
		
		const long a = aa;
		const long b = bb == -1 ? v.size() : bb;
		
		assert (b >= a);
		
		if (a == b)
			return a;
		
		auto s = v.begin()+a;
		auto e = v.begin()+b;
		
		return std::distance(s, std::upper_bound(s, e, p, comp));
	}

	template<typename E, typename F = std::less<E>>
	inline long getBucket(const vector<E>& v, const E& p, F comp=F(), const long aa = 0, const long bb = -1){
	
		const long x = countLess(v,p,comp,aa,bb);
		const long y = countLessEqual(v,p,comp,aa,bb);
		
		return y > x ? x : x == 0 ? 0 : x-1;
	}
	
	
	template<typename E, typename F = std::less<E>>
	inline void removeDuplicates(vector<E>& v, F comp=F(), const int maxRepetitions=1){
		
		const long size = v.size();
		const long size2 = v.size()-maxRepetitions;
		
		long current = 0;
		
		for (long j = 0; j < size; j++){
		
			while (j < size2){
				
				bool keepElement = false;
			
				for (long k = j; !(keepElement) && ( k < (j+maxRepetitions) ); k++)
					keepElement = comp(v[k], v[k+1]) || comp(v[k+1], v[k]);
				
				if (keepElement)
					break;
				
				j++;
			}
			
			v[current++] = v[j];
		}
		
		v.resize(current);
	}
	
	template<typename E, typename C, typename F = std::less<E>>
	inline void removeDuplicates(vector<E>& v, vector<C>& weights, F comp=F(), const int maxRepetitions=1){
		
		assert (v.size() == weights.size() );
		
		const long size = v.size();
		const long size2 = v.size()-maxRepetitions;
		
		long current = 0;
		
		for (long j = 0; j < size; j++){
		
			C sum = 0;
			
			while (true){
			
				sum += weights[j];
			
				if (j >= size2)
					break;
				
				bool keepElement = false;
			
				for (long k = j; !(keepElement) && ( k < (j+maxRepetitions) ); k++)
					keepElement = comp(v[k], v[k+1]) || comp(v[k+1], v[k]);
				
				if (keepElement)
					break;
					
				j++;
			}
			
			v[current] = v[j];
			weights[current] = sum;
			
			current++;
		}
		
		v.resize(current);
		weights.resize(current);
	}

/*
	template<typename E, typename F = std::less<E>>
	inline long countLessEqual(const vector<E>& v, const E& p, F comp=F(), const long aa = 0, const long bb = -1){
		
		const long a = aa;
		const long b = bb == -1 ? v.size() : bb;
		
		assert (b >= a);
		
		if (a == b)
			return a;
		
		auto s = v.begin()+a;
		auto t = v.begin()+(b-1);
		auto e = v.begin()+b;
		
		if ( comp(p, (*s) ) )
			return a;
		
		if ( comp((*t), p ) )
			return b-1;
		
		auto it = std::lower_bound(s, e, p, comp);
		
		if (!comp( p, (*it)))
			return it-v.begin();
		
		return (it-v.begin())-1;
	}*/
	
	template<typename T>
	inline void permute( const vector<long>& posVec, vector<T>& v){
	
		assert(posVec.size() == v.size() );
	
		vector<T> v2(v.size() );
		
		for (long j = 0; j < posVec.size(); ++j){
		
			v2[posVec[j]] = v[j];
		}
		
		for (long j = 0; j < v.size(); ++j){
		
			v[j] = v2[j];
		}
	}

	template<typename T>
	inline T sum( const vector<T>& v){
	
		T ret = 0;
		
		for (auto it = v.begin(); it != v.end(); ++it)
			ret += (*it);
			
		return ret;
	}
	
	template<typename T>
	inline T max( const vector<T>& v){
	
		T ret = std::numeric_limits<T>::lowest();
		
		for (auto it = v.begin(); it != v.end(); ++it)
			if ( (*it) > ret)
				ret = (*it);
			
		return ret;
	}
	
	template<typename T>
	inline T min( const vector<T>& v){
	
		T ret = std::numeric_limits<T>::largest();
		
		for (auto it = v.begin(); it != v.end(); ++it)
			if ( (*it) < ret)
				ret = (*it);
			
		return ret;
	}

	template<typename T>
	inline T maxDiff( const vector<T>& v1, const vector<T>& v2){
	
		assert (v1.size() == v2.size() );
		
		T ret = 0;
		
		for (long j = 0; j < v1.size(); ++j){
		
			const T& a = v1[j];
			const T& b = v2[j];
			
			if ((b-a) > ret)
				ret = b-a;
			
			if ((a-b) > ret)
				ret = a-b;
		}
			
		return ret;
	}

	inline shared_ptr< vector<long> > getRandomIndices(long srcSize, long dstSize, string seed){
	
		assert (dstSize < srcSize);
	
		shared_ptr< vector<long> > v(new vector<long>(dstSize, -1) );
	
		ReservoirSampling rs(seed, dstSize);
		
		rs.reset();
		
		for (long k = 0; k < srcSize; k++){
		
			const long j = rs.pos();
			
			if ( (j >= 0) && (j < dstSize) )
				v->at(j) = k;
		}
		
		std::sort(v->begin(), v->end() );
		
		assert (v->at(0) != -1);
		
		return v;
	}
	
	inline void getRandomIndices(vector<long>& v, long srcSize, long dstSize, string seed){
	
		assert (dstSize < srcSize);
		
		for (long k = 0; k < dstSize; k++)
			v.push_back(-1);
	
		ReservoirSampling rs(seed, dstSize);
		
		rs.reset();
		
		for (long k = 0; k < srcSize; k++){
		
			const long j = rs.pos();
			
			if ( (j >= 0) && (j < dstSize) )
				v.at(j) = k;
		}
		
		std::sort(v.begin(), v.end() );
		
		assert (v.at(0) != -1);
	}
	
	template<typename T>
	inline vector<T> nonNegativeOnly(const vector<T>& v){
	
		vector<int> ret;
	
		for (auto it = v.begin(); it != v.end(); ++it){
			
			if ( (*it) >= 0 )
				ret.push_back(*it);
		}
		
		return ret;
	}
	/*
	template<typename E>
	inline long binarySearchGreaterOrEqual(vector<E>* vec, E val){
	
		if (vec->size() == 0)
			return 0L;
	
		auto it = std::lower_bound(vec->begin(), vec->end(), val);
		return it == vec->end() ? vec->size() : it-vec->begin();
	}
	
	template<typename E>
	inline long binarySearchGreater(vector<E>* vec, E val){
	
		if (vec->size() == 0)
			return 0L;
	
		auto it = std::upper_bound(vec->begin(), vec->end(), val);
		return it == vec->end() ? vec->size() : it-vec->begin();
	}
	
	template<typename E>
	inline long binarySearchLessOrEqual(vector<E>* vec, E val){
	
		return binarySearchGreater(vec, val)-1L;
	}
	
	template<typename E>
	inline long binarySearchLess(vector<E>* vec, E val){
	
		return binarySearchGreaterOrEqual(vec, val)-1L;
	}*/
	
};

struct CompObj { 
  
	double key;
    long ind;

	
	bool operator<(const CompObj& other) const{
		return key < other.key;
    }
};

class Sorter{
public:
	vector<CompObj> vec;
	
	long length = 0;
	
	bool sorted = false;
	
	inline void clear(){
		
		sorted = true;
		length = 0;
	}
	
	inline void shrink_to_fit(){
	
		vec.clear();
		vec.shrink_to_fit();
	}
	
	inline void reserve(long size){
	
		vec.reserve(size);
	}
	
	
	inline void add(long e, double key){
		
		if (vec.size() > length){
		
			CompObj& obj = vec[length];
		
			obj.ind = e;
			obj.key = key;
			
		} else {
		
			CompObj obj;
			obj.ind = e;
			obj.key = key;
			vec.push_back(obj);
		}
		
		sorted = false;
		length++;
	}
	
	inline void sort(){
	
		if (sorted)
			return;
	
		std::sort(vec.begin(), vec.begin()+length );
		sorted = true;
	}
	
	inline long getSortedIndex(long ind){
	
		sort();
		
		return vec[ind].ind;
	}
	
	template<typename E>
	inline shared_ptr<vector<E>> getSorted(vector<E>& v){
	
		sort();
		
		shared_ptr<vector<E>> ret(new vector<E>() );
		
		ret->reserve(v.size() );
		
		for (auto it = vec.begin(); it != (vec.begin()+length); it++)
			ret->push_back( v[it->ind] );
		
		return ret;
	}
	
	template<typename E>
	inline void sort(vector<E>& v){
		
		sort();
		
		shared_ptr<vector<E>> sorted = sorted(v);
		v = (*sorted);
	}
};

template<typename E>
class SkipLUT{

private:

	long nonEmpty;

public:
	vector<long> l;
	vector<long> u;
	
	vector<long> rev;
			
	
	inline SkipLUT(vector<E>& v) : l(v.size()), u(v.size() ){
		
		long ind = -1;
		
		for (long j = 0; j < v.size(); j++){
			
			if (v[j] == 0){
			
				l[j] = ind;
				u[j] = ind+1;
				
			} else {
			
				ind++;
				l[j] = ind;
				u[j] = ind;
			}
		}
		
		nonEmpty = ind+1;
		
		rev.reserve(nonEmpty);
		
		for (long j = 0 ; j < v.size(); j++){
			
			if (!isEmpty(j))
				rev.push_back(j);
		}
	}
	
	inline bool isEmpty(long k){
	
		return l[k] != u[k];
	}
	
	inline long getStart(long x){
	
		return u[x];
	}
	
	inline long getSmallIndex(long x){

		if (u[x] != l[x])
			return -1;
	
		return l[x];
	}
	
	inline long getLargeIndex(long x){
	
		return rev[x];
	}
	
	inline long getEnd(long x){
	
		return l[x];
	}
	
	inline long nonEmptyNum() const{
	
		return nonEmpty;	
	}
	

	

};


template<typename T>
class Finder{

private:
	double minVal = std::numeric_limits<double>::max();
	double maxVal = std::numeric_limits<double>::lowest();
	
	T min;	
	T max;

	bool minIsSet = false;
	bool maxIsSet = false;
	
public:	
	inline void feed(const T& object, double val){
	
		if ( (!minIsSet) ||val < minVal){
			min = object;
			minVal = val;
			minIsSet = true;
		}
		
		if ( (!maxIsSet) || val > maxVal){
			max = object;
			maxVal = val;
			maxIsSet = true;
		}
	}
	
	inline bool hasMin(){
	
		return minIsSet;
	}
	
	inline bool hasMax(){
	
		return maxIsSet;
	}
	
	inline void reset(){
		
		minIsSet = false;
		maxIsSet = false;
		
		minVal = std::numeric_limits<double>::max();
		maxVal = std::numeric_limits<double>::lowest();
	}
	
	inline T getMin() const{
	
		assert (minIsSet);
	
		return min;
	}
	
	inline T getMax() const{
	
		assert (maxIsSet);
	
		return max;
	}
};



template<typename E>
class CircularBuffer{

public:
	vector<E> vec;
	
	int frontPos = 0;
	int size = 0;
	
	const int bufferSize;
	
	inline CircularBuffer(int bufferSize) : bufferSize(bufferSize){
	
		
	}
	
	inline bool empty() const{
	
		return size == 0;
	}
	
	inline void push(E e){
	
		if (vec.size() < bufferSize){
		
			vec.push_back(e);
		}
	
		if ((frontPos+size) < vec.size() ){
			
			vec[frontPos+size] = e;
				
		} else {
			
			vec[ (frontPos+size) - vec.size() ] = e;
		}
		
		size++;
		
		if (size > vec.size() )
			cout << "BUFFER OVERFLOW!!!" << endl;
		
		assert (size <= vec.size() );
	}
	
	inline E front() const{
		
		assert (size > 0);
		
		return vec[frontPos];
	}
	
	inline void pop(){
	
		assert (size > 0);
		
		frontPos++;
		size--;
		
		if ( (frontPos == vec.size()) ||size == 0)
			frontPos = 0;
	}
};



class RandomSkipper{

public:

	long current = -1;

	long ind = -1;
	
	long skipIndex = 0;
	
	shared_ptr<vector<long>> indices;

	inline RandomSkipper(long domain, long size, const string seed){
	
		indices = UtilCollection::getRandomIndices( domain, size, seed);
	
		std::sort( indices->begin(), indices->end() );	
		
		reset();
	}

	inline bool skip(){
	
		ind++;
		
		if (ind != current)
			return false;
		
		skipIndex++;
			
		current = -1;
			
		if (skipIndex < indices->size() )
			current = indices->at(skipIndex);
			
		return true;
	}
	
	inline void reset(){
	
		ind = -1;
		skipIndex = 0;
		current = indices->at(skipIndex);
	}

};

#endif