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

#ifndef SHEKELYAN_UTILPOINT_H
#define SHEKELYAN_UTILPOINT_H

#include <approxdata/utils/utils.h>


//typedef vector<double> Point;




class Point{

public:
	
	//double sortKey = 0;
	vector<double> v;
	
	long id = -1;
	
	//byte lexDims = 0;
	
	inline Point(){
	
		
	}
	
	inline Point(int dims) : v(dims, 0){
	
		
	}
	
	inline operator double(){ 
	
		return (*this)[0];
	}
	
	inline Point(int dims, double d, const long id = -1) : v(dims, d), id(id){
	
		setAll(d);
	}
	
	inline void setAll(double d){
	
		for (auto it = v.begin(); it != v.end(); it++)
			(*it) = d;
	}
	
	inline Point(bool b, double x, double y) : v(2, 0){
		
		v[0] = x;
		v[1] = y;
	}
	
	inline double dotProduct(const Point& p) const{
	
		double ret = 0;
		
		for (int i = 0; i < v.size(); i++)
			ret += v[i] * p[i];
		
		return ret;
	}
	
	inline void swap(Point& p){
		
		v.swap(p.v);
		
		const long id2 = id;
		
		id = p.id;
		p.id = id2;
	}
	
	inline Point(Point& p){
	
		v.clear();
		
		for (auto it = p.v.begin(); it != p.v.end(); it++)
			v.push_back( (*it) );
	
		//v = p.v;
		id = p.id;
		//sortKey = p.sortKey;
	}
	
	
	inline double getEnclosedVolume(const Point& max) const{
	
		double ret = 1.0;
	
		for (int i = 0; i < v.size(); i++){
			
			ret *= max[i]-v[i];
			
			if (ret < 0)
				return 0;
		}
		
		return ret;
	}
	
	inline Point(const Point& p){
	
		v.clear();
		
		for (auto it = p.v.begin(); it != p.v.end(); it++)
			v.push_back( (*it) );
		
		//v = p.v;
		id = p.id;
		//sortKey = p.sortKey;
	}
	
	inline const string toString() const{
	
		stringstream s;
	
		s << "[";
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				s << ", ";
				
			const double d = v[i];
			
			s << UtilString::doubleToString(d, 10);
		}
		
		s << "]";
		
		if (id != -1)
			s << "{" << to_string(id) << "}";
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const Point& m) { 
    	
    	return os << m.toString();
	}
	
	
	
	inline Point& operator=(const Point& other){ // copy assignment{
    
    	if (this != &other){ // self-assignment check expected
			
			v = other.v;
			id = other.id;
			//sortKey = other.sortKey;
   		}
   		 
 	   	return *this;
	}
	
	
	inline Point& operator=(const vector<double>& other){ // copy assignment{
    
    	v = other; //sortKey = 0;
   		 
 	   	return *this;
	}
	
	inline bool operator==(const Point& other) const{
    
    	return v == other.v;
	}
	
	inline bool operator!=(const Point& other) const{
    
    	return v != other.v;
	}
	
	inline void clear(){
	
		v.clear();
	}
	
	inline void reserve(int m){
	
		v.reserve(m);
	}
	
	inline const double& operator[](std::size_t idx) const{
		
		return v.at(idx);
	}
	
	inline double& operator[](std::size_t idx){
		
		return v.at(idx);
	}
	
	inline void push_back(double d){
	
		v.push_back(d);
	}
	
	inline double& at(std::size_t idx){
		
		return v.at(idx);
	}
	
	/*
	inline bool operator<(const Point& other) const{ // copy assignment{
    
    	return sortKey < other.sortKey;
	}*/
	
	inline const std::size_t size() const{
	
		return v.size();
	}
	
	inline void set(Point& p){
	
		for (int i = 0; i < getDims(); i++)
			(*this).at(i) = p[i];
			
		id = p.id;
	}
	
	inline void setMin(Point& p){
	
		for (int i = 0; i < getDims(); i++){
			(*this).at(i) = (*this).at(i) <= p.at(i) ? (*this).at(i) : p.at(i);
		}
	}
	
	inline void setMax(Point& p){
	
		for (int i = 0; i < getDims(); i++)
			(*this).at(i) = (*this).at(i) >= p.at(i) ? (*this).at(i) : p.at(i);
	}
	
	inline int getDims() const {
		
		return size();
	}
};

class PointComparator : public std::binary_function<const Point, const Point, bool>{

public:

	bool isSimple;
	bool isLex;
	bool isLexArray;
	
	int dim = -1;
	
	vector<int> lex;
	
	template<bool EQ>
	inline bool f(const Point& a, const Point& b) const{
		
		assert (a.size() == b.size() );
		
		if (isSimple)
			return EQ ? a[dim] <= b[dim] : a[dim] < b[dim];
		
		if (isLex){
		
			if (a[dim] < b[dim])
				return true;
			
			if (a[dim] > b[dim])
				return false;
			
			for (int i = 0; i < dim; i++){
			
				if (a[i] < b[i])
					return true;
			
				if (a[i] > b[i])
					return false;
			}
			
			for (int i = dim+1; i < a.size(); i++){
				
				if (a[i] < b[i])
					return true;
			
				if (a[i] > b[i])
					return false;
			}
			
			if (a.id < b.id)
				return true;
			
			if (a.id > b.id)
				return false;
			
			return EQ;
		} 
		
		if (isLexArray){
		
			for (int ii = 0; ii < lex.size(); ii++){
			
				const int i = lex[ii];
			
				if (a[i] < b[i])
					return true;
					
				if (a[i] > b[i])
					return false;
			}
			
			if (a.id < b.id)
				return true;
			
			if (a.id > b.id)
				return false;
			
			return EQ;
		}
		
		assert(false);
	}
	
public:
	
	inline PointComparator(const PointComparator& pc) : isSimple(pc.isSimple), isLex(pc.isLex), isLexArray(pc.isLexArray), dim(pc.dim), lex(pc.lex){
	
	}
	
	inline PointComparator(PointComparator& pc) : isSimple(pc.isSimple), isLex(pc.isLex), isLexArray(pc.isLexArray), dim(pc.dim), lex(pc.lex){
	
	}
	
	
	inline PointComparator(int dim, bool _isLex=false) : isSimple(!_isLex), isLex(_isLex), isLexArray(false),  dim(dim), lex({}){
	
		
	}

	inline PointComparator(vector<int> _lex) : isSimple(false), isLex(false), isLexArray(true), dim(-1), lex(_lex){
	
		
		
	}
	
	
	inline double interpolation(const Point& p1, const Point& p2, const Point& p3) const{
	
		if (isLex){
		
			assert (p2[dim] >= p1[dim]);
			
			assert (p3[dim] >= p1[dim]);
			assert (p3[dim] <= p2[dim]);
			
			if (p1[dim] == p2[dim]){
			
				for (int i = 0; i < p1.size(); ++i)
					if (p1[i] < p2[i])
						return (p3[i]-p1[i])/(p2[i]-p1[i]);
				
				if (p1.id < p2.id)
					return ((p3.id-p1.id)*1.0)/(p2.id-p1.id);
				
				return 0.5;
			}
			
			if (p3[dim] == p1[dim])
				return 0;
		
			if (p3[dim] == p2[dim])
				return 1;
		
			return (p3[dim]-p1[dim])/(p2[dim]-p1[dim]);
		}
		
		assert(false);
	}
	
	inline PointComparator& operator=(const PointComparator& other){
	
		if (this != &other){ // self-assignment check expected
			
			dim = other.dim;
			lex = other.lex;
			isSimple = other.isSimple;
			isLex = other.isLex;
			isLexArray = other.isLexArray;
   		}
   		 
 	   	return *this;
	}
	
	inline int getMainDim() const{
	
		if (isLexArray)
			return lex[0];
		else
			return dim;
	}
	
	
	inline long countLess(const vector<Point>& points, const Point& p){
	
		if (points.size() == 0)
			return 0;
	
		const int last = points.size()-1;
	
		if ( (*this)(points[last], p))
			return points.size();
		
		return std::upper_bound(points.begin(), points.end(), p, (*this))-points.begin();
	}
	
	template <class E>
	inline void sort(E beg, E end){
	
		//for (auto it = beg; it != end; it++){
		//
		//	cout << (*it) << endl;
		//}
	
		std::sort(beg, end, (*this) );
	}
	
	/*
	inline long countLessOrEqual(const vector<Point>& points, const Point& p){
	
		if (points.size() == 0)
			return 0;
	
		const int last = points.size()-1;
		
		if (comp(points[last], p))
			return points.size();
		
		return (1L+std::lower_bound(points.begin(), points.end(), p, this())-points.begin() );
	}
	
	inline long getBucket(const vector<Point>& points, const Point& p, bool minContained = false, bool maxContained = false, bool beginsAtPoint= true){
	
		long ret = 0;
	
		if (minContained)
			ret = beginsAtPoint ? countLessOrEqual(points, p)-1;
		else
			ret = countLessOrEqual(points, p);
		
		return ret; 
	}*/
	
	inline bool less(const Point& a, const Point& b) const{
		
		return f<false>(a,b);
	}
	
	inline bool lessOrEqual(const Point& a, const Point& b) const{
		
		return f<true>(a,b);
	}
	
	inline bool greater(const Point& a, const Point& b) const{
		
		return !lessOrEqual(a, b);
	}
	
	inline bool greaterOrEqual(const Point& a, const Point& b) const{
		
		return !less(a,b);
	}
	
	inline bool equal(const Point& a, const Point& b) const{
		
		assert (a.size() == b.size() );
		
		if (isSimple)
			return a[dim] == b[dim];
		
		if (isLex){
			
			for (int i = 0; i < a.size(); i++){
			
				if (a[i] != b[i])
					return false;
			}
			
			return true;
		
		} 
		
		if (isLexArray){
		
			for (int ii = 0; ii < lex.size(); ii++){
			
				const int i = lex[ii];
			
				if (a[i] != b[i])
					return false;
			}
			
			return true;
		}
		
		assert(false);
	}
	
	inline bool unequal(const Point& a, const Point& b) const{
	
		return !equal(a,b);
	}
	
	inline bool operator()(const Point& a, const Point& b) const{
		
		return less(a,b);
	}
};


template<int dims>
class ArrayComparator : public std::binary_function<const std::array<double, dims>, const std::array<double, dims>, bool>{

public:

	std::less< const std::array<double, dims> > comp;
	
	inline double interpolation(const std::array<double, dims>& p1, const std::array<double, dims>& p2, const std::array<double, dims>& p3) const{
	
		for (int i = 0; i < dims; ++i)
			if (p1[i] < p2[i])
				return p3[i] <= p1[i] ? 0 : p3[i] >= p2[i] ? 1 : (p3[i]-p1[i])/(p2[i]-p1[i]);
		
		return 0.5;
	}
	
	inline bool operator()(const std::array<double, dims>& a, const std::array<double, dims>& b) const{
		
		return comp(a,b);
	}
	
};

class VectorComparator : public std::binary_function<const vector<double>, const vector<double>, bool>{

protected:

	
	bool isSimple;
	bool isLex;
	bool isLexArray;
	
	int dim;
	
	vector<int> lex;
	
	const int dims;
	
	template<bool EQ>
	inline bool f(const vector<double>& a, const vector<double>& b) const{
		
		//assert (a.size() == b.size() );
		
		if (isSimple)
			return EQ ? a[dim] <= b[dim] : a[dim] < b[dim];
		
		if (isLex){
		
			if (a[dim] < b[dim])
				return true;
			
			if (a[dim] > b[dim])
				return false;
			
			for (int i = 0; i < dim; i++){
			
				if (a[i] < b[i])
					return true;
			
				if (a[i] > b[i])
					return false;
			}
			
			for (int i = dim+1; i < dims; i++){
				
				if (a[i] < b[i])
					return true;
			
				if (a[i] > b[i])
					return false;
			}
			
			return EQ;
		} 
		
		if (isLexArray){
		
			for (int ii = 0; ii < lex.size(); ii++){
			
				const int i = lex[ii];
			
				if (a[i] < b[i])
					return true;
					
				if (a[i] > b[i])
					return false;
			}
			
			return EQ;
		}
		
		assert(false);
	}
	
public:
	
	inline VectorComparator(const VectorComparator& pc) : dims(pc.dims), isSimple(pc.isSimple), isLex(pc.isLex), isLexArray(pc.isLexArray), dim(pc.dim), lex(pc.lex){
	
	}
	
	inline VectorComparator(VectorComparator& pc) : dims(pc.dims), isSimple(pc.isSimple), isLex(pc.isLex), isLexArray(pc.isLexArray), dim(pc.dim), lex(pc.lex){
	
	}
	
	
	inline VectorComparator(int dim, bool _isLex=false) : dims(1), isSimple(!_isLex), isLex(_isLex), isLexArray(false),  dim(dim), lex({}){
	
		
	}

	inline VectorComparator(vector<int> _lex) : dims(lex.size() ), isSimple(false), isLex(false), isLexArray(true), dim(-1), lex(_lex){
	
		
		
	}
	
	inline VectorComparator& operator=(const VectorComparator& other){
	
		if (this != &other){ // self-assignment check expected
			
			lex = other.lex;
			isSimple = other.isSimple;
			isLex = other.isLex;
			isLexArray = other.isLexArray;
   		}
   		 
 	   	return *this;
	}
	
	//
	
	inline VectorComparator(const PointComparator& pc, int dims) : dims(dims), isSimple(pc.isSimple), isLex(pc.isLex), isLexArray(pc.isLexArray), dim(pc.dim), lex(pc.lex){
	
	}
	
	inline VectorComparator(PointComparator& pc, int dims) : dims(dims),isSimple(pc.isSimple), isLex(pc.isLex), isLexArray(pc.isLexArray), dim(pc.dim), lex(pc.lex){
	
	}
	
	inline VectorComparator& operator=(const PointComparator& other){
	
		if (true){ // self-assignment check expected
			
			lex = other.lex;
			isSimple = other.isSimple;
			isLex = other.isLex;
			isLexArray = other.isLexArray;
   		}
   		 
 	   	return *this;
	}
	
	inline int getMainDim() const{
	
		if (isLexArray)
			return lex[0];
		else
			return dim;
	}
	
	
	inline bool less(const vector<double>& a, const vector<double>& b) const{
		
		return f<false>(a,b);
	}
	
	inline bool lessOrEqual(const vector<double>& a, const vector<double>& b) const{
		
		return f<true>(a,b);
	}
	
	inline bool greater(const vector<double>& a, const vector<double>& b) const{
		
		return !lessOrEqual(a, b);
	}
	
	inline bool greaterOrEqual(const vector<double>& a, const vector<double>& b) const{
		
		return !less(a,b);
	}
	
	inline bool equal(const vector<double>& a, const vector<double>& b) const{
		
		assert (a.size() == b.size() );
		
		if (isSimple)
			return a[dim] == b[dim];
		
		if (isLex){
			
			for (int i = 0; i < a.size(); i++){
			
				if (a[i] != b[i])
					return false;
			}
			
			return true;
		
		} 
		
		if (isLexArray){
		
			for (int ii = 0; ii < lex.size(); ii++){
			
				const int i = lex[ii];
			
				if (a[i] != b[i])
					return false;
			}
			
			return true;
		}
		
		assert(false);
	}
	
	inline bool unequal(const vector<double>& a, const vector<double>& b) const{
	
		return !equal(a,b);
	}
	
	inline bool operator()(const vector<double>& a, const vector<double>& b) const{
		
		return less(a,b);
	}
	
};


namespace UtilPoint{

	inline void print(Point& v){
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				cout << ", ";
				
			cout << v.at(i);
		}
		
		cout << endl;
	}
	
	/*
	inline void sort(vector<Point>& pts, int dim){
	
		if (pts.size() == 0)
			return;
	
		for (auto it = pts.begin(); it != pts.end(); it++)
			(*it).sortKey = (*it)[dim];
				
		sort(pts.begin(), pts.end() );
	}
	
	inline void stable_sort(vector<Point>& pts, int dim){
	
		if (pts.size() == 0)
			return;
	
		for (auto it = pts.begin(); it != pts.end(); it++)
			(*it).sortKey = (*it)[dim];
				
		stable_sort(pts.begin(), pts.end() );
	}
	
	inline void sort(const vector<int>& lex, vector<Point>& pts){
	
		for (int ii = lex.size()-1; ii >= 0; ii--)
			stable_sort(pts, lex[ii]);
	}
	
	
	inline void sort_(const vector<int>& lex, vector<Point>& pts){
		
		if (pts.size() <= 1)
			return;
		
		sort(pts, lex[0]);
		
		int prevI = 0;
		
		for (int ii = 1; ii < lex.size(); ii++){
			
			const int i = lex[ii];
			
			auto it = pts.begin();
			
			auto it1 = it;
			auto it2 = it;
			
			double last = (*it)[i];
			
			it++;
			
			auto end = pts.end();
			
			for (; it != end; it++){
			
				const double v = (*it)[prevI];
				
				if (v != last){
				
					if (it2 != end){
					
						bool sortNeeded = false;
						
						double prev = std::numeric_limits<double>::lowest();
						long size = 0;
						
						for (auto it3 = it1; it3 != it2; it3++){
							
							const double d = (*it3)[i];
							
							if ( d < prev)
								sortNeeded = true;
							
							prev = d;
							size++;
						}
						
						if (sortNeeded){
						
							if (false && size == 2){
								
								(*it1).swap( (*it2) );
								
							} else {
							
								for (auto it3 = it1; it3 != it2; it3++)
									(*it3).sortKey = (*it3)[i];
								
								sort(it1, it2);
							}	
						}
					}
					
					it1 = it;
					it2 = end;
					last = v;
					
				} else {
				
					it2 = (it+1);
				}
			}
			
			prevI = i;
		}
	
	}
	
	inline void lexSort(vector<Point>& pts, int dim){
	
		if (pts.size() == 0)
			return;
		
		const int DIMS = pts[0].size();
		
		vector<int> v;
		v.push_back(dim);
		
		for (int i = 0; i < dim; i++)
			v.push_back(i);
			
		for (int i = dim+1; i < DIMS; i++)
			v.push_back(i);
		
		sort(v, pts);
	}
	*/
	inline bool less(const vector<int> lex, const Point& a, const Point& b){
		
		for (int ii = 0; ii < lex.size(); ii++){
			
			const int i = lex[ii];
			
			if (a[i] < b[i])
				return true;
					
			if (a[i] > b[i])
				return false;
		}
			
		return false;
	}
	
	inline bool greater(const vector<int> lex, const Point& a, const Point& b){
		
		for (int ii = 0; ii < lex.size(); ii++){
			
			const int i = lex[ii];
			
			if (a[i] > b[i])
				return true;
					
			if (a[i] < b[i])
				return false;
		}
			
		return false;
	}
	
	inline bool greaterOrEqual(const vector<int> lex, const Point& a, const Point& b){
		
		for (int ii = 0; ii < lex.size(); ii++){
			
			const int i = lex[ii];
			
			if (a[i] > b[i])
				return true;
					
			if (a[i] < b[i])
				return false;
		}
			
		return true;
	}
	
	inline bool equal(const vector<int> lex, const Point& a, const Point& b){
		
		for (int ii = 0; ii < lex.size(); ii++){
			
			const int i = lex[ii];
			
			if (a[i] != b[i])
				return false;
		}
			
		return true;
	}
	
	inline bool lessOrEqual(const vector<int> lex, const Point& a, const Point& b){
		
		for (int ii = 0; ii < lex.size(); ii++){
			
			const int i = lex[ii];
			
			if (a[i] < b[i])
				return true;
					
			if (a[i] > b[i])
				return false;
		}
			
		return true;
	}
	
	
	inline bool isSorted(const vector<int>& lex, vector<Point>& pts){
	
		if (pts.size() <= 1)
			return true;
			
		Point prev(pts[0].size() );
	
		for (int i = 0; i < prev.size(); i++)
			prev[i] = std::numeric_limits<double>::lowest();
		
		for (auto it = pts.begin(); it != pts.end(); it++){
			
			const Point& p = (*it);
			prev = p;
		}
		
		return true;
	
	}
	
	/*
	inline void fastSort(const vector<int>& lex, vector<Point>& pts){
	
		
		if (pts.size() <= 1)
			return;
			
		if (pts.size() <= 100){
		
			sort(lex, pts);
			return;
		}
		 
		
		const int DIMS = pts[0].size();
		
		Point w(DIMS);
		
		const double smallNum = 1e-5;
		const double bigNum = 1;
		const double mul = 1.0/(lex.size()-1);
		
		cout << lex[0] << " " << lex[1] << endl;
		
		for (int i = 0; i < lex.size(); i++){
		
			const double a = i*mul;
			w[ lex[i] ] = smallNum*a+bigNum*(1-a); 
		}
		
		cout << w << endl;
		
		for (auto it = pts.begin(); it != pts.end(); it++){
		
			it->sortKey = it->dotProduct(w);
		}
		
		sort( pts.begin(), pts.end() );
		
		assert ( isSorted(lex, pts));
		
		
		if (isSorted(lex, pts))
			return;
			
		sort(lex, pts);
		
	}*/
	
	/*
	inline void copySortedBetween(const vector<int> lex, vector<Point>& src, vector<Point>& dst, Point& begin, Point& end){
	
		const int i = lex[0];
		
		for (auto it = src.begin(); it != src.end(); it++)
			(*it).sortKey = (*it)[i];
		
		begin.sortKey = begin[i];
		end.sortKey = end[i];
		
		sort( src.begin(), src.end() );
		
		auto it3 = lower_bound(src.begin(), src.end(), begin);
		auto it4 = upper_bound(it3, 		src.end(), end);
		
		const double aval = begin[i];
		const double bval = end[i];
			
		for (auto it = it3; it != it4; it++){
		
			if ( ((*it)[i] > aval) && ((*it)[i] < bval) )
				dst.push_back(*it);
			
			if ( lessOrEqual(lex, begin, *it) && less(lex, *it, end) )
				dst.push_back(*it);
		}
	}*/
	
	/*
	inline void fill(Point& v, int count,...){
	
		va_list args;
		int i;
		
		v.clear();
		va_start(args, count);
		for(i = 0; i < count; i++)
			v.push_back( va_arg(args, double) );
		va_end(args);
	
	}*/
	
	inline bool contains(const Point& amin, const Point& amax, Point& v){

		for (int i = 0; i < amin.size(); i++){

			if ( amin[i] > v[i] || amax[i] < v[i])
				return false;
		}
	
		return true;
	}
	
	inline void copy(const Point& v1, Point& v2){
	
		v2.clear();
		v2.reserve(v1.size());
		
		for (int i = 0; i < v1.size(); i++)
			v2.push_back(v1[i]);
			
		v2.id = v1.id;
	}
	
	inline void min(const Point& v1, Point& v2){
	
		for (int i = 0; i < v1.size(); i++){
		
			if (v1[i] < v2[i] ){
				v2[i] = v1[i];
				v2.id = v1.id;
			}
		}
		
		
	}
	
	inline void max(const Point& v1, Point& v2){
	
		for (int i = 0; i < v1.size(); i++){
		
			if (v1[i] > v2[i] ){
				v2[i] = v1[i];
				v2.id = v1.id;
			}
		}
		
		
	}

};





class Box{

public:

	Point min;
	Point max;
	double count = 0;
	
	inline Box() : min(0), max(0){
		
	}
	
	inline Box(int dims) : min(dims), max(dims){
	
		for (int i = 0; i < dims; i++){
		
			min[i] = 0;
			max[i] = 1;
		}
	}
	
	inline bool isEmpty() const{
	
		return min.size() == 0;
	}
	
	inline void makeEmpty(){
	
		min.v.clear();
		max.v.clear();
	}
	
	inline Box(double x1, double y1, double x2, double y2) : min(x1,y1), max(x2, y2){
		
		
	}
	
	inline Box(const Point& pmin, const Point& pmax, long c=0) : min(pmin.size()), max(pmax.size()), count(c){
	
		min = pmin;
		max = pmax;
	}
	
	
	inline Box(Box& b): min(b.min.size()), max(b.max.size()), count(b.count){
		
		min = b.min;
		max = b.max;
	}
	
	inline Box(const Box& b): min(b.min.size()), max(b.max.size()), count(b.count){
		
		min = b.min;
		max = b.max;
	}
	
	inline void setMin(const Point& v){
	
		min = v;
	}
	
	inline void setMax(const Point& v){
	
		max = v;
	}
	
	inline void setMin(Point& v){
	
		min = v;
	}
	
	inline void setMax(Point& v){
	
		max = v;
	}
	
	inline const Point getCorner(const int dimFlags) const{
			
		const int DIMS = min.getDims();
		
		Point p = min;
		
		for (int i = 0; i < DIMS; i++){
		
			if (dimFlags & (1L << i))
				p[i] = max[i];
		}
	
		return p;
	}
	
	
	inline void intersect(const Box& b){
		
		const int DIMS = min.getDims();
		
		for (int i = 0; i < DIMS; i++){
		
			min[i] = UtilMath::maxVal(min[i], b.min[i]);
			max[i] = UtilMath::minVal(max[i], b.max[i]);
		}
	}
	
	inline bool differenceIsEmptyBox(const Box& b) const{
	
		return UtilMath::isIn( intersection(b), SpatialRelation::EQUAL, SpatialRelation::CONTAINED);
	}
	
	
	inline bool differenceIsNonEmptyBox(const Box& b) const{
		
		int intersectionDims = 0;
		
		const int DIMS = min.getDims();
		
		for (int i = 0; i < DIMS; i++){
				
			const double bmin = UtilMath::maxVal(b.min[i], min[i]);
			const double bmax = UtilMath::minVal(b.max[i], max[i]);
			
			if ( UtilMath::isNotIn(bmin, min[i], max[i]) && UtilMath::isNotIn(bmax, min[i], max[i]))
				return false;
			
			const bool b1 = UtilMath::isOnceIn(bmin, min[i], max[i]);
			const bool b2 = UtilMath::isOnceIn(bmax, min[i], max[i]);
			
			if ( UtilMath::isOnceIn(true, b1, b2) )
				intersectionDims++;
			
			if (intersectionDims > 1)
				return false;
		}
		
		return intersectionDims == 1;
	}
	
	
	// returns false, if set minus would result in non-box shape or empty box
	inline bool minus(const Box& b){
	
		if (!differenceIsNonEmptyBox(b))
			return false;
	
		const int DIMS = min.getDims();
		
		int intersectionDims = 0;
		
		for (int i = 0; i < DIMS; i++){
		
			const double bmin = UtilMath::maxVal(b.min[i], min[i]);
			const double bmax = UtilMath::minVal(b.max[i], max[i]);
				
			const SpatialRelation r1 = UtilMath::intervalIntersection(min[i], max[i], bmin, bmax);
			
			if ( UtilMath::isIn(r1, SpatialRelation::CONTAINS)){
				
				UtilMath::intervalMinus(min[i], max[i], bmin, bmax);
				intersectionDims++;
				continue;
			}
			
			assert (r1 == SpatialRelation::EQUAL);
		}
		
		assert (intersectionDims == 1);
		
		return true;
	}
	
	
	inline virtual Box& operator=(const Box& other){ // copy assignment{
    
    	if (this != &other){ // self-assignment check expected
			
			min = other.min;
			max = other.max;
			count = other.count;
   		}
   		 
 	   	return *this;
	}
	
	inline SpatialRelation intersection(const Box& other) const{
	
		bool contains = true;
		bool contained = true;
		bool equal = true;
	
		for (int i = 0; i < min.size(); i++){
			
			SpatialRelation s = UtilMath::intervalIntersection(min[i], max[i], other.min[i], other.max[i] );
			
			if (s != SpatialRelation::EQUAL){
			
				if (s == SpatialRelation::NO_INTERSECTION)
					return SpatialRelation::NO_INTERSECTION;
			
				equal = false;
				
				if (s != SpatialRelation::CONTAINS)
					contains = false;
			
				if (s != SpatialRelation::CONTAINED)
					contained = false;
			}
		}
		
		if (equal)
			return SpatialRelation::EQUAL;
		
		if (contains)
			return SpatialRelation::CONTAINS;
		
		if (contained)
			return SpatialRelation::CONTAINED;
		
		return SpatialRelation::OVERLAP;
	}
	
	inline bool contains(const Point& p) const{
	
		for (int i = 0; i < p.size(); i++){
		
			if ( min[i] > p[i] || max[i] < p[i])
				return false;
		}
		
		
		
		return true;
	}
	
	inline double volumeDistance(const Point& p) const{
	
		double maxScale = -1;
	
		for (int i = 0; i < p.size(); i++){
		
			const double a = min[i]*0.5;
			const double b = max[i]*0.5;
			
			const double radius = b-a;
			const double dist = UtilMath::absVal(p[i]-(a+b) );
			
			maxScale = UtilMath::maxVal<double>(maxScale, dist/radius);
		}
		
		return maxScale;
	}
	
	inline double getVolume() const{
	
		return min.getEnclosedVolume(max);
		
		/*
	
		double ret = 1.0;
	
		for (int i = 0; i < min.size(); i++){
		
			ret *= max[i]-min[i];
		}
		
		return ret;*/
	}
	
	inline void enclose(const Point& v){
	
		if (min.size() == 0){
		
			UtilPoint::copy(v, min);
			UtilPoint::copy(v, max);
			
		} else {
			
			UtilPoint::min(v, min);
			UtilPoint::max(v, max);
		}
	}
	
	inline void print(long size=-1) const{
		
		if (size==-1)
			cout << "count " << count << " min " << min << " max " << max << endl;
		else
			cout << "count " << count << " sel " << (count*100.0/size) << "%" << " min " << min << " max " << max << endl;
	}


	inline string toString() const{
	
		stringstream s;
		
		s << "box ";
		
		if (count >= 0)
		 	s << "count " << count << " ";
	
		s << "min " << min << " max " << max << endl;
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const Box& m) { 
    	
    	return os << m.toString();
	}
};

/*
class CountBox : public Box{

public:
	long count;
	
	inline CountBox() : Box(), count(0){
	
	}
	
	inline CountBox(const Point& pmin, const Point& pmax, const long c) : Box(pmin, pmax), count(c){
	
	}
	
	inline CountBox(const CountBox& cb) : Box(cb.min, cb.max), count(cb.count){
	
	}
	
	inline CountBox& operator=(const CountBox& other){ // copy assignment{
    
    	if (this != &other){ // self-assignment check expected
			
			min = other.min;
			max = other.max;
			count = other.count;
   		}
   		
 	   	return *this;
	}
	
	inline void print(){
	
		cout << "count " << count << endl;
		UtilPoint::print(min);
		UtilPoint::print(max);	
	}
};*/



class UniformPoints{
public:
	
	std::seed_seq seed;
	int DIMS;
	std::unique_ptr<std::mt19937> g;
	
	inline UniformPoints(int dims, string seedStr="abc") : seed(seedStr.begin(), seedStr.end() ), DIMS(dims)  {
	
		g.reset(new std::mt19937(seed) );
	}
	
	inline void randomPoint(const Box& container, Point& p){
	
		p.clear();
		
		for (int i = 0; i < DIMS; i++){
			
			std::uniform_real_distribution<double> d(container.min[i], container.max[i]);
			p.push_back( d(*g) );
		}
	}
	
	inline void randomBox(const Box& container, Box& created){
		
		Point p1;
		Point p2;
		
		randomPoint(container, p1); 
		randomPoint(container, p2); 
		
		created.enclose(p1);
		created.enclose(p2);
	}
	
	inline void randomBox(const Box& containing, const Box& container, Box& created){
	
		if (containing.min.size() == 0){
			randomBox(container, created);
			return;
		}
	
		Box b1(container.min, containing.min);
		Box b2(containing.max, container.max);
		
		Point p1;
		Point p2;
		
		randomPoint(b1, p1);
		randomPoint(b2, p2);
		
		created.enclose(p1);
		created.enclose(p2);
	}
	
	inline void randomPoint(Box whitelist, Box blacklist, Point& p){
	
		p.clear();
		
		for (int i = 0; i < DIMS; i++){
		
			const double s1 = whitelist.min[i];
			const double w1 = blacklist.min[i]-s1;
			
			const double s2 = blacklist.max[i];
			const double w2 = whitelist.max[i]-s2;
			
			std::uniform_real_distribution<double> d(0, w1+w2);
			
			const double v = d( (*g));
			
			p.push_back( v <= w1 ? v+s1 : v+s2-w1);
			
		}
	}

};

class PointList{

public:	

	
	vector<Point> vec;
	
	inline PointList(){
	
	
	}
	
	inline PointList(const PointList& v){
	
		addAll(v.vec);
	}
	
	inline PointList(const vector<Point>& v){
	
		addAll(v);
	}
	
	inline void reserve(long size){
	
		vec.reserve(size);
	}
	
	inline void addAll(const vector<Point>& v){
	
		std::copy(v.begin(), v.end(), std::back_inserter(vec));
	}
	
	
	inline const Point& operator[](std::size_t idx) const{
		
		return vec[idx];
	}
	
	inline Point& operator[](std::size_t idx){
		
		return vec[idx];
	}
	
	inline Point& front(){
		
		return vec[0];
	}
	
	inline Point& back(){
		
		return vec[vec.size()-1];
	}
	
	inline long size() const{
	
		return vec.size();
	}
	
	inline long getSize() const{
	
		return vec.size();
	}
	
	inline bool operator==(const PointList& l) const{
		
		if (size() != l.size())
			return false;
		
		for (long j = 0; j < size(); j++)
			if ( (*this)[j] != l[j] )
				return false;
		
		return true;
	}
	
	
	inline void push_back(const Point& p){
	
		vec.push_back(p);
	}
	
	inline void add(const Point& p){
	
		vec.push_back(p);
	}
	
	inline void add(double x, double y){
	
		vec.push_back(Point(x,y));
	}
	
	inline void sort(int dim){
	
		std::sort(vec.begin(), vec.end(), PointComparator(dim));
	}
	
	inline void lexSort(int dim){
	
		std::sort(vec.begin(), vec.end(), PointComparator(dim, true));
	}
	
	inline void sort(vector<int> lex){
	
		std::sort(vec.begin(), vec.end(), PointComparator(lex));
	}
	
	inline Box getBoundingBox(){
	
		Box b;
		
		for (auto it = vec.begin(); it != vec.end(); it++)
			b.enclose( (*it) );
			
		return b;
	}
	
	inline string toString() const{
	
		stringstream s;
	
		s << "{";
	
		for (long i = 0; i < vec.size(); i++){
		
			if (i > 0)
				s << ", ";
				
			s << vec[i];
		}
		
		s << "}";
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const PointList& m) { 
    	
    	return os << m.toString();
	}
};

class QueryBox : public Box{

public:
	long flags = 0;	
	
	QueryMode mode = QueryMode::EST;

	const long INF_MIN = 1L << 10;
	const long INF_MAX = 1L << 20;
	
	const long EXC_MIN = 1L << 30;
	const long EXC_MAX = 1L << 40;
	
	inline void setFlags(long f){
	
		this->flags = f;
	}
	
	inline long getFlags() const{
	
		return flags;
	}
	
	inline bool isSet(long f, int i) const{
	
		return (flags & (f << i)) != 0;
	}
	
	inline void set(long f, bool b = true){
	
		flags |= f;
		if (!b)
			flags ^= f;
	}
	
	inline void set(long f, int i, bool b = true){
	
		set(f << i, b);
	}
	
	inline void setUnbounded(int i, bool b1, bool b2){
	
		set(INF_MIN, i, b1);
		set(INF_MAX, i, b2);
	}
	
	inline void setExcluded(int i, bool b1, bool b2){
	
		set(EXC_MIN, i, b1);
		set(EXC_MAX, i, b2);
	}
	
	
	inline long getUnboundedMin() const {
	
		return ((1L << 9)-1) & (flags >> 10);
	}
	
	inline long getUnboundedMax() const {
	
		return ((1L << 9)-1) & (flags >> 20);
	}
	
	inline long getExcludedMin() const {
	
		return ((1L << 9)-1) & (flags >> 30);
	}
	
	inline long getExcludedMax() const {
	
		return ((1L << 9)-1) & (flags >> 40);
	}
	
	inline bool isMinUnbounded(int i) const {
	
		return isSet(INF_MIN, i);
	}
	
	inline bool isMaxUnbounded(int i) const {
	
		return isSet(INF_MAX, i);
	}
	
	inline bool isMinExcluded(int i) const {
	
		return isSet(EXC_MIN, i);
	}
	
	inline bool isMaxExcluded(int i) const {
	
		return isSet(EXC_MAX, i);
	}
	
	
	inline QueryBox(const Point& pmin, const Point& pmax, QueryMode mode) : Box(pmin, pmax), mode(mode){
		
		min.id = -(1L << 40);
		max.id = 1L << 40;
	}
	
	inline QueryBox(const QueryBox& b) : Box(b.min, b.max), mode(b.mode), flags(b.getFlags() ){
	
	}
	
	inline void setMax(int i, const double& pi, bool excl=false){
	
		max[i] = pi;
		
		set(EXC_MAX, i, excl);
		set(INF_MAX, i, false);
	}
	
	inline void setMin(int i, const double& pi, bool excl=false){
		
		min[i] = pi;
		
		set(EXC_MIN, i, excl);
		set(INF_MIN, i, false);
	}
	
	/*
	
	
	inline long getExcludedMin() const {
	
		return ((1L << 10)-1) & (flags >> 30);
	}
	
	inline long getExcludedMax() const {
	
		return ((1L << 10)-1) & (flags >> 40);
	}
	
	inline void setMax(int i, const Point& p, bool excl=false){
	
		max.id = excl ? p.id+1 : p.id;
		max[i] = p[i];
		
		set(EXC_MAX, excl);
		set(INF_MAX, false);
	}
	
	inline void greater(int i, const Point& p){
	
		assert( p[i] >= min[i] );
	
		setMin(i, p, true);
	}
	
	inline void greaterOrEqual(int i, const Point& p){
	
		assert( p[i] >= min[i] );
	
		setMin(i, p, false);
	}
	
	inline void lessOrEqual(int i, const Point& p){
		
		assert( p[i] <= max[i] );
		
		setMax(i, p, false);
	}
	
	
	inline void less(int i, const Point& p){
		
		assert( p[i] <= max[i] );
		
		setMax(i, p, true);
	}*/
	
	
};


#endif