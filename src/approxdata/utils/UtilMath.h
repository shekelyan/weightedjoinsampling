/*
 * ApproxData Library
 * Copyright (c) 2018 Michael Shekelyan <michael.shekelyan@gmail.com>
 */

/*
Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subjecconst T to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions 
of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef SHEKELYAN_UTILMATH_H
#define SHEKELYAN_UTILMATH_H

#include <approxdata/utils/utils.h>



// DEPRECATED
class Stats{

public:
	vector<double> v;

	bool sorted = false;

	const string name;
	const string unit;
	
	inline Stats(const string s="", const string u="") : name(s), unit(u){
	
		
	}

	inline void add(double p){
	
		v.push_back(p);
		sorted = false;
	}
	
	inline void sort(){
	
		if (sorted)
			return;
		
		std::sort(v.begin(), v.end() );
		sorted = true;
	}
	
	inline void reset(){
	
		v.clear();
		sorted = true;
	}
	
	inline double sum() const{
	
		double ret = 0;
	
		for (auto it = v.begin(); it != v.end(); it++)
			ret += (*it);
			
		return ret;
	}
	
	inline double max() const{
	
		assert (v.size() > 0);
	
		if (sorted)
			return v[v.size()-1];
		
		double ret = v[0];
		
		for (auto it = v.begin(); it != v.end(); it++)
			if ( (*it) > ret)
				ret = (*it);
		
		return ret;
	}
	
	inline double min() const{
	
		assert (v.size() > 0);
		
		if (sorted)
			return v[0];
	
		double ret = v[0];
		
		for (auto it = v.begin(); it != v.end(); it++)
			if ( (*it) < ret)
				ret = (*it);
		
		return ret;
	}
	
	inline double avg() const{
		
		return sum()/v.size();
	}
	
	inline long num() const{
		
		return v.size();
	}
	
	inline double median(double q = 0.5) const{
	
		assert (v.size() > 0);
		
		assert (sorted);
		
		return v[ round((v.size()-1)*q) ];
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		s << name;
		s << " \t";
		s << "num " << num();
		s << " \t";
		s << "min " << min() << unit;
		s << " \t";
		s << "avg " << avg() << unit;
		s << " \t";
		s << "max " << max() << unit;
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const Stats& m) { 
    	
    	return os << m.toString();
	}
	
	
	
	
};

// DEPRECATED
class Interval{


public:
	double a;
	double b;
	
	bool aIncluded = true;
	bool bIncluded = true;
	
	inline Interval intersection(Interval& v){
	
		Interval ret;
		
		ret.a = a > v.a ? a : v.a;
		ret.b = b < v.b ? b : v.b;
		
		return ret;
	}
};



class RandomGenerator{

public:

	std::seed_seq seed1;
	std::mt19937 generator;
	
	std::uniform_real_distribution<double> rnd;
	
	inline RandomGenerator(const string seed, double a = 0, double b = 1) : seed1(seed.begin(), seed.end()), generator(seed1), rnd(a,b){
	
	}
	
	template<class E>
	inline void shuffle(E beg, E end){
	
		std::shuffle(beg, end, generator);	
	}
	
	inline double randomDouble(){
	
		return rnd(generator);
	}
	
	inline double randVal(){
		
		return rnd(generator);
	}
	
	inline double laplace(double b){
	
		const double u1 = randomDouble();
		
		const bool sign = u1 >= 0.5;
		
		const double u2 = sign ? u1-0.5 : u1;
		
		return (sign ? +b : -b) * log(1.0-2.0*u2);
	}
	
	inline long randomInteger(long a, long b){
	
		if (a == b)
			return a;
	
		const long ret = round(a+randVal()*(b-a));
		
		return ret < a ? a : ret > b ? b : ret;
	}
	
	template<typename T>
	inline const T& randomElement(const vector<T>& v){
	
		return v.at(randomInteger(0, v.size()-1));
	}
};

class NormGenerator{

public:

	std::seed_seq seed1;
	std::mt19937 generator;
	
	std::normal_distribution<double> rnd;
	
	inline NormGenerator(const string seed, double mean = 0, double std = 1) : seed1(seed.begin(), seed.end()), generator(seed1), rnd(mean,std){
	
	}
	
	inline double randVal(){
	
		return rnd(generator);
	}
};


class BinarySearcher{
	
public:
	long double xmin = std::numeric_limits<long double>::lowest();
	
	long double xminY = std::numeric_limits<long double>::quiet_NaN();
	
	long double xmax = std::numeric_limits<long double>::max();
	
	long double xmaxY = std::numeric_limits<long double>::quiet_NaN();

	bool unreachableMin = false;
	bool unreachableMax = false;

	long double targetY;
	
	long double increasing = true;
	
	long double flip = 1.0;
	
	
	int step = -1;
	
	vector<long double> xs;
	vector<long double> ys;
	
	
	inline void inputGreaterThan(long double x){
	
		xmin = x;
		xminY = std::numeric_limits<long double>::quiet_NaN();
		
		unreachableMin = true;
	}
	
	inline void setOutput(long double y){
	
		targetY = y*flip;
	}
	
	inline void inputLessThan(long double x){
	
		xmax = x;
		xmaxY = std::numeric_limits<long double>::quiet_NaN();
		
		unreachableMax = true;
	}
	
	inline void inputAtLeast(long double x, long double y=std::numeric_limits<long double>::quiet_NaN()){
	
		xmin = x;
		xminY = y;
		
		unreachableMin = false;
	}
	
	inline void inputAtMost(long double x, long double y=std::numeric_limits<long double>::quiet_NaN()){
	
		xmax = x;
		xmaxY = y;
		unreachableMax = false;
	}
	
	inline bool feed(long double x, long double y){
		
		y *= flip;
		
		if (step < 10){
		
			xs.push_back(x);
			ys.push_back(y);
			return false;
		}
		
		if (step == 10){
		
			int inc = 0;
			int dec = 0;
			int eq = 0;
		
			for (int i = 1; i < ys.size(); i++){
			
				if ( ys[i] > ys[i-1])
					inc++;
				
				if ( ys[i] < ys[i-1])
					dec++;
					
				if ( ys[i] == ys[i-1])
					eq++;
			}
			
			if (inc > 0)
				increasing = true;
			else
				increasing = false;
			
			
			if (! ( (increasing && (dec == 0)) || ((!increasing) && (inc == 0)))){
			
				for (int j = 0; j < xs.size(); j++){
				
					cout << "x " << xs[j] << " y " << ys[j] << endl;
				}
				
				cout << "targetY " << targetY << endl;
			}
			
			assert ( (increasing && (dec == 0)) || ((!increasing) && (inc == 0)));
			
			if (!increasing)
				flip = -1;
			
			targetY *= flip;
			xminY *= flip;
			xmaxY *= flip;
			
			for (int i = 0; i < ys.size(); i++)
				ys[i] *= flip;
			
			
			for (int i = 0; i < ys.size(); i++){
				
				if (ys[i] <= targetY && (std::isnan(xminY) || (targetY-ys[i]) < (targetY-xminY)))
					inputAtLeast(xs[i],ys[i]);
				
				if (ys[i] >= targetY && (std::isnan(xmaxY) || (ys[i]-targetY) < (xmaxY-targetY)))
					inputAtMost(xs[i],ys[i]);
			}
			
			
		}
		
		if (step >= 10){
			
			if (y >= targetY)
				inputAtMost(x,y);
			
			if (y <= targetY)
				inputAtLeast(x,y);
		}
		
		return xmin == xmax;
	}
	
	inline long double getInputLB(){
	
		return xmin;
	}

	inline long double getInputUB(){
	
		return xmax;
	}

		
	inline long double getInput(){
	
		step++;
		
		if (step == 0){
		
			if (unreachableMin)
				return xmin*0.99+xmax*0.01;	
		}
		
		if (step == 9){
		
			if (unreachableMax)
				return xmin*0.01+xmax*0.99;	
		}
		
		if (step <= 10){
		
			return xmin+(xmax-xmin)*step/10.0;
		}
		
		assert (!(unreachableMin && unreachableMax));
		
		if (unreachableMin){
		
			return xmin*0.99+xmax*0.01;
		}
		
		if (unreachableMax){
		
			return xmin*0.01+xmax*0.99;
		}
		
		return xmin*0.5+xmax*0.5;
	}

};



class RealFunction{

public:
	inline virtual double operator()(const double x) const = 0;
	inline virtual ~RealFunction(){};
	
	
	inline double inv(double domainMin, double domainMax, double targetY, const bool inc=true, int steps=1000){
		
		const RealFunction& f = (*this);
		
		
		if (false){
		
		
			BinarySearcher b;
			
			b.setOutput(targetY);
			
			b.inputAtLeast(domainMin);
			b.inputAtMost(domainMax);
		
			for (int i = 0; i < steps; i++){
			
				const double x = b.getInput();
				const double y = f(x);
				
				if (b.feed(x,y))
					break;
			}
			
			return b.getInput();
		}
		
		
		//const double a = domainMin*0.9+domainMax*0.1;
		//const double b = domainMin*0.1+domainMax*0.9;
		
		//const bool monotonicallyIncreasing = f(a) < f(b);
		
		
		double xmin = domainMin;
		double xmax = domainMax;
		
		//if (inc)
		//	cout << "f(" << a << ") = " << f(a) << " < " << f(b) << " = f(" << b << ")" << endl;
		//else
		//	cout << "f(" << a << ") = " << f(a) << " >= " << f(b) << " = f(" << b << ")" << endl;
		
		for (int i = 0; i < 1000; i++){
		
			double x = xmin*0.5+xmax*0.5;
			
			const double y = f(x);
			
			//cout << "f(" << x << ") = " << y << endl;
			
			if (y == targetY)
				return x;
				
			if (inc){
			
				if (y > targetY){
				
					xmax = x;
					
				} else {
				
					xmin = x;
				}
				
			} else {
			
				if (y < targetY){
				
					xmax = x;
					
				} else {
				
					xmin = x;
				}
			}
		}
		
		return inc ? xmin : xmax;
	}
};



template<typename X, typename Y>
class BinarySearch{

	
	
public:

	X x1;
	X x2;
	Y y;
	
	X x;
	
	int steps;
	
	X lb;
	X ub;
	
	Y lb_y;
	Y ub_y;
	
	int step = 0;
	
	X lastX;
	
	inline BinarySearch(const X& x1, const X& x2, const Y& y, int steps = 100) : x1(x1), x2(x2), y(y), steps(steps){
		
		x = (x1+x2)/2;
		
		lastX = std::numeric_limits<X>::lowest();
		
		lb = std::numeric_limits<X>::lowest();;
		ub = std::numeric_limits<X>::max();;
		lb_y = std::numeric_limits<Y>::lowest();
		ub_y = std::numeric_limits<Y>::max();
		
	}
	
	inline X toward(const X& x1, const X& x2) const{
	
		if (std::is_same<X, long>::value || std::is_same<X, int>::value){
		
			if (x2 < x1)
				return x1-1;
			else if (x2 > x1)
				return x1+1;
			else
				return x1;
		}
		
		return std::nexttoward(x1, x2);
	}
	
	
	inline bool find(const Y& y3){
		
		++step;
		
		if (x == lastX)
			return false;
			
		lastX = x;
		
		if (lb == ub || step > steps)
			return false;
		
		assert (y3 == y3);
		
		if (y3 > y){ // f(x') > y => x < x'
		
			if (x < ub){
				ub = x;
				ub_y = y3;
			}
			
			x2 = toward(x, x1);
			
			
		} else if (y3 < y){ // f(x') < y => x > x'
		
			if (x > lb){
				lb = x;
				lb_y = y3;
			}
		
			x1 = toward(x, x2);
			
		} else {
		
			lb = x;
			ub = x;
			
			lb_y = y;
			ub_y = y;
			
			x1 = x;
			x2 = x;
		}
		
		x = (x1+x2)/2;
	
		return true;
	}

};

enum class SpatialRelation : int{ NO_INTERSECTION = 1L<<0, CONTAINED=1L<<2, OVERLAP=1L<<3, CONTAINS=1L<<4, EQUAL=1L<<5};




namespace UtilMath{


	template<typename T>
	inline bool isIn(const T t, const T s1){
	
		if (t == s1)
			return true;
			
		return false;	
	}
	
	
	inline long largestDoubleInteger(){
	
		return 1L << 52;
	}
	
	inline long lowestDoubleInteger(){
	
		return -(1L << 52);
	}
	
	template<typename T>
	inline bool equalsAll(const T& t, const T& s1){
	
		return t == s1;	
	}
	
	template<typename T>
	inline bool equalsAll(const T& t, const T& s1, const T& s2){
	
		return t == s2 && equalsAll(t, s1);
	}
	
	template<typename T>
	inline bool equalsAll(const T& t, const T& s1, const T& s2, const T& s3){
	
		return t == s3 && equalsAll(t, s1, s2);
	}
	
	template<typename T>
	inline bool equalsAll(const T& t, const T& s1, const T& s2, const T& s3, const T& s4){
	
		return t == s4 && equalsAll(t, s1, s2, s3);
	}
	
	template<typename T>
	inline bool equalsAll(const T& t, const T& s1, const T& s2, const T& s3, const T& s4, const T& s5){
	
		return t == s5 && equalsAll(t, s1, s2, s3, s4);
	}
	
	template<typename T>
	inline bool isOnceIn(const T t, const T s1){
	
		if (t == s1)
			return true;
			
		return false;	
	}
	
	template<typename T>
	inline bool isOnceIn(const T t, const T s1, const T s2){
	
		int count = 0;
	
		if (t == s1)
			count++;
		
		if (t == s2)
			count++;
			
		return count == 1;	
	}
	
	template<typename T>
	inline bool isOnceIn(const T t, const T s1, const T s2, const T s3){
	
		int count = 0;
	
		if (t == s1)
			count++;
		
		if (t == s2)
			count++;
		
		if (t == s3)
			count++;
		
		return count == 1;	
	}
	
	template<typename T>
	inline bool isOnceIn(const T t, const T s1, const T s2, const T s3, const T s4){
	
		int count = 0;
	
		if (t == s1)
			count++;
		
		if (t == s2)
			count++;
		
		if (t == s3)
			count++;
			
		if (t == s4)
			count++;
		
		return count == 1;	
	}

	template<typename T>
	inline bool isOnceIn(const T t, const T s1, const T s2, const T s3, const T s4, const T s5){
	
		int count = 0;
	
		if (t == s1)
			count++;
		
		if (t == s2)
			count++;
		
		if (t == s3)
			count++;
			
		if (t == s4)
			count++;
			
		if (t == s5)
			count++;
		
		return count == 1;	
	}

	template<typename T>
	inline bool isIn(const T t, const T s1, const T s2){
	
		if (t == s1)
			return true;
		
		if (t == s2)
			return true;
			
		return false;	
	}

	template<typename T>
	inline bool isIn(const T t, const T s1, const T s2, const T s3){
	
		if (t == s1)
			return true;
		
		if (t == s2)
			return true;
	
		if (t == s3)
			return true;
			
		return false;	
	}
	
	template<typename T>
	inline bool isIn(const T t, const T s1, const T s2, const T s3, const T s4){
	
		if (t == s1)
			return true;
		
		if (t == s2)
			return true;
	
		if (t == s3)
			return true;
			
		if (t == s4)
			return true;
			
		return false;	
	}
	
	inline vector<int> shortestPath(vector<vector<int>> adj, int i, int j, int ignoreI = -1, int ignoreJ = -1){
	
		vector<int> pred(adj.size(), -1);
		vector<int> dist(adj.size(), 99999999);
		
		dist[i] = 0;
		
		vector<int> queue;
		vector<int> queue2;
		
		queue.push_back(i);
		
		while (queue.size() > 0){
			
			queue2.clear();
		
			for (int k = queue.size()-1; k >= 0; k--){
			
				const long k1 = queue[k];
			
				if (k1 == -1)
					continue;
				
				vector<int>& v = adj[k1];
				
				for (auto it = v.begin(); it != v.end(); it++){
				
					const long k2 = (*it);
					
					if (k1 == ignoreI && k2 == ignoreJ)
						continue;
					
					if (k1 == ignoreJ && k2 == ignoreI)
						continue;
					
					if (dist[k1]+1 < dist[k2]){
						
						dist[k2] = dist[k1]+1;
						pred[k2] = k1;
						
						queue2.push_back(k2);
					}
				}
			}
			
			queue = queue2;
		}
		
		vector<int> ret;
		
		int cur = j;
		
		while (cur != -1){
		
			ret.insert( ret.begin(), cur);
			
			cur = pred[cur];
		}
		
		return ret;
	}
	
	template<typename T>
	inline bool isIn(const T t, const T s1, const T s2, const T s3, const T s4, const T s5){
	
		if (t == s1)
			return true;
		
		if (t == s2)
			return true;
	
		if (t == s3)
			return true;
			
		if (t == s4)
			return true;
		
		if (t == s5)
			return true;
			
		return false;	
	}

	template<typename T>
	inline bool isNotIn(const T t, const T s1){
	
		return !isIn(t, s1);	
	}

	template<typename T>
	inline bool isNotIn(const T t, const T s1, const T s2){
	
		return !isIn(t, s1, s2);	
	}

	template<typename T>
	inline bool isNotIn(const T t, const T s1, const T s2, const T s3){
	
		return !isIn(t, s1, s2, s3);	
	}

	template<typename T>
	inline bool isNotIn(const T t, const T s1, const T s2, const T s3, const T s4){
	
		return !isIn(t, s1, s2, s3, s4);	
	}

	template<typename T>
	inline bool isNotIn(const T t, const T s1, const T s2, const T s3, const T s4, const T s5){
	
		return !isIn(t, s1, s2, s3, s4, s5);	
	}
	
	
	template <class T>
	inline T log2choose(T n, T k, T limit=0){
	
		if (k > n)
			return 0;
		
		if (k*2 > n)
			k = log2(n-k);
		
		if (k == 0)
			return 1;
	
		T ret = log2(n);
		
		for (int i = 2; i <= k; ++i){
		
			ret += log2(n-i+1);
			ret -= log2(i);
			
			if (limit != 0 && ret > limit)
				return ret;
		}
		
		return ret;
	}
	
	template <class T>
	inline T log10choose(T n, T k, T limit=0){
	
		if (k > n)
			return 0;
		
		if (k*2 > n)
			k = log10(n-k);
		
		if (k == 0)
			return 1;
	
		T ret = log10(n);
		
		for (int i = 2; i <= k; ++i){
		
			ret += log10(n-i+1);
			ret -= log10(i);
			
			if (limit != 0 && ret > limit)
				return ret;
		}
		
		return ret;
	}

	template <class T>
	inline T choose(T n, T k, T limit=0){
	
		if (k > n)
			return 0;
		
		if (k*2 > n)
			k = n-k;
		
		if (k == 0)
			return 1;
	
		T ret = n;
		
		for (int i = 2; i <= k; ++i){
		
			ret *= n-i+1;
			ret /= i;
			
			if (limit != 0 && ret > limit)
				return ret;
		}
		
		return ret;
	}
	
	template <class T>
	inline T multichoose(T n, T k, T limit=0){
	
		return choose(n+k-1, k);
	}

	inline bool isPowerOfTwo(long x){
	
		return x && !(x & (x - 1));
	}
	
	inline int getPowerOfTwo(long x){
    	
    	int ret = -1;
    	
		for ( ; x != 0; x >>= 1)
			ret++;
		
		return ret;
	}
	
	
	inline SpatialRelation intervalIntersection(double amin, double amax, double bmin, double bmax){
	
		if (amin == bmin && amax == bmax)
			return SpatialRelation::EQUAL;
	
		if (amin <= bmin && amax >= bmax)
			return SpatialRelation::CONTAINS;
	
		if (bmin <= amin && bmax >= amax)
			return SpatialRelation::CONTAINED;
			
		if (amin <= bmin && amin >= bmax)
			return SpatialRelation::OVERLAP;
			
		if (amax <= bmin && amax >= bmax)
			return SpatialRelation::OVERLAP;
			
		return SpatialRelation::NO_INTERSECTION;		
	}
	
	
	inline SpatialRelation intervalIntersection(const double amin, const double amax, const double b){
	
		if (amin == b && amax == b)
			return SpatialRelation::EQUAL;
		
		if (amin <= b && b <= amax)
			return SpatialRelation::CONTAINS;
			
		return SpatialRelation::NO_INTERSECTION;		
	}
	
	
	
	inline long gcd(long a, long b){
	
		if ((a <= 1) || (b <= 1))
			return 1;
	
		while (b != 0){
		
			const long b_ = b;
			b = a % b;
			a = b_;
		}
		
		return a;
	}
	
	inline double translate(double v, double min, double max, double newmin, double newmax){
		
		if (max == min)
			return newmin*0.5+newmax*0.5;
		
		if (v <= min)
			return newmin;
		
		if (v >= max)
			return newmax;
		
		return (v-min)/(max-min)*(newmax-newmin)+newmin;
	}

	inline double largestNumberSmallerThan(const double x){
	
		return x-1e-40;
	}
	
	inline double smallestNumberLargerThan(const double x){
		
		return x+1e-40;
	}
	
	
	
	inline void intervalMinus(double& amin, double& amax, const double bmin, const double bmax){
		
		if (bmin <= amin && amin < bmax){
		
			amin = largestNumberSmallerThan(bmax);
			return;
		}
		
		if (bmin < amax && amax <= bmax){
		
			amax = smallestNumberLargerThan(bmin);
		}
	}
	
	
	inline long largestNegativeLong(){
	
		return std::numeric_limits<long>::lowest();
	}
	
	inline long largestLong(){
	
		return std::numeric_limits<long>::max();
	}
	
	
	
	inline double largestNegativeDouble(){
	
		return std::numeric_limits<double>::lowest();
	}
	
	inline double largestDouble(){
	
		return std::numeric_limits<double>::max();
	}

	

	inline double integrate(RealFunction& f, double a, double b, int steps=20){
	
		if (b <= a)
			return 0;
	
		const double step = (b-a)/(steps-1);
		
		assert (step > 0);
	
		const double x1 = a+step;
		const double x2 = b-step;
		
		double sum = 0.5*f(a);
		
		steps *= 2;
		
		for (double x = x1; x <= x2 && steps > 0; x += step){
			sum += f(x);
			
			steps--;
		}
		
		sum += 0.5*f(b);
		
		return step*sum;
	}
	
	template<typename E>
	inline E minVal(E a, E b){
	
		return a < b ? a : b;
	}
	
	template<typename E>
	inline E maxVal(E a, E b){
	
		return a > b ? a : b;
	}
	
	template<typename E>
	inline bool isBetween(const E& v, const E& a, const E& b){
	
		return ( a <= v) && (v <= b);
	}
	
	template<typename E>
	inline E minVal(E a, E b, E c){
	
		return minVal(minVal(a,b),c);
	}
	
	template<typename E>
	inline E maxVal(E a, E b, E c){
	
		return maxVal(maxVal(a,b),c);
	}
	
	template<typename E>
	inline E absVal(E a){
	
		return a < 0 ? -a : a;
	}
	
	
	template<typename E>
	inline E makeBetween(E a, E b, E c){
		
		if (c < a)
			return a;
		
		if (c > b)
			return b;
		
		return c;
		//return maxVal(a, minVal(b, c));
	}
	
	inline double makeMultipleOf(double v, long m){
	
		if (m == 0)
			return v;
	
		return (((long) v)/m)*m;
	}
	
	inline double makeMultipleOfUB(double v, long m){
	
		if (m <= 1)
			return v;
		
		const long x = v;
		const long x2 = (x/m)*m;
		
		if ( x2 == x)
			return x2;
		else
			return x2+m;
	}
	
	inline double truncate(double min, double max, double x){
	
		if (x < min)
			return min;
		
		if (x > max)
			return max;
		
		return x;
	}
	
	inline long powerOfTwoLB(long m){
	
		long ret = (1L << 62);	
		
		if (m <= 0)
			return 0;
		
		while (ret != 0){
		
			if (ret <= m)
				return ret;
				
			ret >>= 1;
		}
		
		return 0;
	}
	
	inline long powerOfTwoUB(long m){
		
		if (m <= 0)
			return 0;
		
		for (long ret = (1L << 62); ret != 0; ret >>= 1){
		
			if (ret == m)
				return ret;
			
			if (ret & m)
				return ret << 1;
		}
		
		return 0;
	}
	
	inline int getPowerOfTwoUB(long x){
    	
    	
		return getPowerOfTwo(powerOfTwoUB(x));
	}
	
	template<typename E>
	inline double sum(const vector<E>& vec){
	
		double ret = 0;
		for (auto it = vec.begin(); it != vec.end(); it++)
			ret += (*it);
		
		return ret;
	}
	
	
	
	inline void optimalRounding2(vector<double>& v){

		const double vsum = sum<double>(v);
		double sumOfFloors = 0;
		
		for (long j = 0; j < v.size(); j++)
			sumOfFloors += floor(v[j]);
		
		const long roundDown = v.size()-(round(vsum)-sumOfFloors);
		
		Sorter s;
		
		for (long j = 0; j < v.size(); j++)
			s.add(j, v[j]-floor(v[j]) );
		
		for (long k = 0; k < roundDown; k++){
			
			const long j = s.getSortedIndex(k);
			v[j] = floor(v[j]);
		}
		
		for (long k = roundDown; k < v.size(); k++){
	
			const long j = s.getSortedIndex(k);
			v[j] = ceil(v[j]);
		}
		
		const double a = round(sum<double>(v));
		const double b = round(vsum);
		
		if (a != b)
			cout << "optimalRounding2 " << a << " != " << b << endl;
		
		//assert ( sum(v) == vsum);
	}
	
	inline double doubleValuesToKb(long n){
	
		return n*(8.0/1024.0);
	}
	
	inline void optimalRounding(vector<double>& v){
	
		const double a = round(sum<double>(v));
	
		vector<double> v2 = v;
		
		for (long j = 0; j < v2.size(); j++)
			if (v2[j] == ((long) v2[j]) )
				v2[j] = 0;
		
		SkipLUT<double> skip(v2);
		
		vector<double> v3(skip.nonEmptyNum() );
		
		for (long j = 0; j < v.size(); j++){
		
			if (!skip.isEmpty(j))
				v3[skip.getSmallIndex(j)] = v[j];
		}
		
		optimalRounding2(v3);
		
		for (long k = 0; k < v3.size(); k++){
		
			v[skip.getLargeIndex(k) ] = v3[k];
		}
		
		const double b = round(sum<double>(v));
		
		if (a != b)
			cout << "optimalRounding " << a << " != " << b << endl;
		
	}
};


namespace UtilHist{

	// input: real "v" between 0 and 1
	// output: integer "discretize(v, b)" between 0 and (b-1)
	
	
	
	// find i in [0, v.size()-2] s.t. p >= v[i] && p < v[i+1]
	
	template<typename E, typename F = std::less<E>>
	inline long getBucket(vector<E> v, E p, F comp=F() ){
		
		auto it = std::upper_bound(v.begin(), v.end(), p, comp);
		
		if (it == v.begin())
			return 0;
		
		if (it == v.end() )
			return v.size()-2;
		
		return (it-v.begin() )-1;
	}
	
	
	template<typename E, typename F = std::less<E>>
	inline long getBucket(vector<E> v, E p, long a, long b, F comp=F() ){
		
		auto s = v.begin()+a;
		auto t = v.begin()+b;
		
		auto it = std::upper_bound(s,t, p, comp);
		
		if (it == s)
			return a;
		
		if (it == t)
			return b-2;
		
		return a+(it-s)-1;
	}
	
	inline long discretizePowerOfTwo(double v, int lod){
		
		if (v <= 0)
			return 0;
		
		if (v >= 1)
			return ((1L << lod)-1);
		
		//return ( (long) (v*(1L << lod)) );
		
		return ( (long) (v*(1L << 62)) ) >> (62-lod);
	}
	
	
	
	// floor( v * buckets )
	
	// [0.0, 0.1[ [0.1, 0.2[ [0.1, 0.3[ [0.1, 0.4[ [0.1, 0.5[ [0.1, 0.6[ [0.1, 0.7[
	
	inline long discretize(double v, long buckets){
	
		assert (buckets > 0);
		
		return v <= 0 ? 0 : v >= 1 ? buckets-1 : (long) (v*buckets);
	}
	
	
	inline double bucketMin(double bucket, double buckets){
		
		return bucket >= buckets ? 1.0 : (0.0+bucket)/buckets;
	}
	
	inline double bucketMax(double bucket, double buckets){
		
		return std::nexttoward( bucketMin(bucket+1, buckets), 0);
	}
	
};

enum class QueryMode : int{ LB=-1, EST=0, UB=1, DEBUG_LB = -2, DEBUG_UB = 2};


class OneDimEstimator{
public:

	virtual ~OneDimEstimator(){};
	inline virtual double intervalCount(double qmin, double qmax, QueryMode mode=QueryMode::EST) const = 0;

protected:
	OneDimEstimator(){}
    OneDimEstimator(const OneDimEstimator&){}
    OneDimEstimator& operator=(const OneDimEstimator&){ return *this; }
};


class Histogram1d{

public:

	virtual ~Histogram1d(){};
	virtual long getBucketNum() const = 0;
	virtual double getBucketMin(long bucket) const = 0;
	
	virtual long getBucket(double d) const = 0;
	virtual double getBucketMax(long bucket) const = 0;
	virtual long& getBucketCount(long bucket) = 0;
	virtual void add(double v) = 0;

protected:
	Histogram1d(){}
    Histogram1d(const Histogram1d&){}
    Histogram1d& operator=(const Histogram1d&){ return *this; }

};

class EquiWidth1d : public Histogram1d, public OneDimEstimator{

private:

	template <typename T>
	inline const string listToString(const vector<T>& v) const{
	
		if (v.size() == 0)
			return "NULL";
	
		if (v.size() > 10){
		
			stringstream s;
	
			for (int i = 0; i < 5; i++){
				
				s << v.at(i);
				s << ", ";
			}
			
			s << "...";
			
			for (int i = 4; i >= 0; i--){
		
				s << ", ";
				s << v.at(v.size()-1-i);
			}
			
			return s.str();
		}
	
		stringstream s;
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				s << ", ";
				
			s << v.at(i);
		}
		
		return s.str();
	}

	inline long getSum(long a, long b) const{
		
		if (a > b)
			return 0;
	
		long ret = 0;
		
		if (!cumulative){
		
			for (long j = a ; j <= b; j++)
				ret += counts[j];
				
			return ret;		
		}
		
		if (a == 0)
			return counts[b];
			
		return counts[b]-counts[a-1];
	}
	

public:
	
	
	bool cumulative = false;
	
	double min;
	double max;
	
	double bucketDiv;
	double normMul;
	
	vector<long> counts;
	
	bool normalize;
	
	int lod = -1;
	
	inline EquiWidth1d(long buckets, double min=0, double max=1) : min(min), max(max), bucketDiv( buckets*1.0/(max-min) ), normMul(1.0/(max-min)), counts(buckets, 0){
	
		if (min != 0 || max != 1)
			normalize = true;	
			
		if (UtilMath::isPowerOfTwo(buckets)){
		
			lod = UtilMath::getPowerOfTwo(buckets);
			
		} else {
			
			lod = -1;
		}
	}
	
	inline EquiWidth1d(const EquiWidth1d& e){
		
		counts.reserve(e.counts.size() );
		
		for (auto it = e.counts.begin(); it != e.counts.end(); it++)
			counts.push_back( (*it) );
		
		cumulative = e.cumulative;
		min = e.min;
		max  = e.max;
		normMul = e.normMul;
		bucketDiv = e.bucketDiv;
		normalize = e.normalize;
		lod = e.lod;
	}
	
	
	~EquiWidth1d() override{}
	
	inline void print(){
	
		for (long j = 0; j < getBucketNum(); j++)
			if (getBucketCount(j) > 0)
				cout << "[" << getBucketMin(j) << "," << getBucketMax(j) << "]:{" << getBucketCount(j) << "}" << endl;
		
		cout << endl;
	}
	
	inline double intervalCount(double qmin, double qmax, QueryMode mode=QueryMode::EST) const override{
		
		const long imin = getBucket(qmin);
		const long imax = getBucket(qmax);
		
		if (imin == imax){
		
			if (mode == QueryMode::LB)
				return 0;
		
			const double c1 = getSum(imin, imin);
			
			if (mode == QueryMode::UB)
				return c1;
		
			return c1*(qmax-qmin)*bucketDiv;
		}
		
		const long c2 = getSum(imin+1, imax-1);
		
		if (mode == QueryMode::LB)
			return c2;
		
		const long c1 = getSum(imin, imin);
		const long c3 = getSum(imax, imax);
		
		if (mode == QueryMode::UB)
			return c1+c2+c3;
		
		const double r1 = UtilMath::absVal<double>( getBucketMax(imin)-qmin )*bucketDiv;
		const double r3 = UtilMath::absVal<double>( qmax-getBucketMin(imax) )*bucketDiv;		
		
		return (r1*c1)+c2+(r3*c3);
	}
	
	
	inline bool makeCumulative(bool b = true){
	
		const bool ret = cumulative;
	
		if (cumulative == b)
			return ret;
		
		if (b){
		
			for (long j = 1; j < counts.size(); j++)
				counts[j] += counts[j-1];
				
		} else {
		
			for (long j = counts.size()-1; j >= 1; j--)
				counts[j] -= counts[j-1];
		}
		
		cumulative = b;
		
		return ret;
	}
	
	
	inline long getBucket(double d) const override{
	
		if (normalize)
			return UtilHist::discretize( (d-min)*normMul, counts.size() );
		
		if (lod > 0)
			UtilHist::discretizePowerOfTwo(d, lod);
		
		return UtilHist::discretize(d, counts.size() );
		
	}
	
	inline void add(double v) override{
	
		assert (!cumulative);
	
		counts[getBucket(v)]++;
	}
	
	inline double getBucketMin(long bucket) const override{
	
		if (normalize)
			return min+( UtilHist::bucketMin(bucket, counts.size()) )*(max-min);
		
		return UtilHist::bucketMin(bucket, counts.size());
	}
	
	inline double getBucketMax(long bucket) const override{
	
		if (normalize)
			return min+( UtilHist::bucketMax(bucket, counts.size()) )*(max-min);
		
		return UtilHist::bucketMax(bucket, counts.size());
	}
	
	inline long getBucketNum() const override{
	
		return counts.size();
	}
	
	inline long getTotalBucketCount() const{
	
		if (cumulative)
			return counts[counts.size()-1];
		
		return UtilMath::sum<long>(counts);
	}
	
	
	inline double getBytes() const{
	
		return counts.size()*8;
	}
	
	inline long& getBucketCount(long bucket) override{
	
		if (cumulative)
			assert (false);
		
		return counts[bucket];
	}
	
	
	
	
	inline long& operator[](size_t bucket){
	
		return getBucketCount(bucket);
	}
	
	inline long& at(size_t bucket){
	
		return getBucketCount(bucket);
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		s << "equiwidth " << getBucketNum() << " buckets between " << min << " and " << max;
		
		if (cumulative)
			s << " cumulative";
		
		s << " count " << getTotalBucketCount();
		
		s << " counts " << listToString<long>(counts);
		
		return s.str();
	}
	

};



#define Calloc2(n, t)   (t *) malloc(n*sizeof(t) )

namespace UtilKS{

	template<typename T>
	inline void Free(T* ptr){
	
		delete[] ptr;
	}
	
	template<typename S, typename T>
	inline double R_pow_di(S s, T t){
	
		return std::pow(s, (double) t);
	}

	static void
	m_multiply(double *A, double *B, double *C, int m)
	{
		/* Auxiliary routine used by K().
		   Matrix multiplication.
		*/
		int i, j, k;
		double s;
		
		cout << 301 << endl;
		
		for(i = 0; i < m; i++)
		for(j = 0; j < m; j++) {
			s = 0.;
			for(k = 0; k < m; k++)
			s+= A[i * m + k] * B[k * m + j];
			C[i * m + j] = s;
		}
		
		cout << 302 << endl;
	}

	static void
	m_power(double *A, int eA, double *V, int *eV, int m, int n)
	{
		/* Auxiliary routine used by K().
		   Matrix power.
		*/
		double *B;
		int eB , i;

		cout << 201 << endl;

		if(n == 1) {
		for(i = 0; i < m * m; i++)
			V[i] = A[i];
		*eV = eA;
		return;
		}
		
		cout << 202 << endl;
		
		m_power(A, eA, V, eV, m, n / 2);
		
		cout << 203 << endl;
		
		B = (double*) Calloc2(m * m, double);
		
		cout << 204 << endl;
		
		m_multiply(V, V, B, m);
		eB = 2 * (*eV);
		
		cout << 205 << endl;
		
		if((n % 2) == 0) {
		for(i = 0; i < m * m; i++)
			V[i] = B[i];
		*eV = eB;
		}
		else {
		
		cout << 206 << endl;
		
		m_multiply(A, B, V, m);
		*eV = eA + eB;
		
		cout << 207 << endl;
		
		}
		
		cout << 208 << endl;
		
		if(V[(m / 2) * m + (m / 2)] > 1e140) {
		for(i = 0; i < m * m; i++)
			V[i] = V[i] * 1e-140;
		*eV += 140;
		}
		
		cout << 209 << endl;
		Free(B);
	}
	
	static double
	K(int n, double d)
	{
		/* Compute Kolmogorov's distribution.
		   Code published in
		 George Marsaglia and Wai Wan Tsang and Jingbo Wang (2003),
		 "Evaluating Kolmogorov's distribution".
		 Journal of Statistical Software, Volume 8, 2003, Issue 18.
		 URL: http://www.jstatsoft.org/v08/i18/.
		*/

	   int k, m, i, j, g, eH, eQ;
	   double h, s, *H, *Q;

		cout << 101 << endl;
	   
		 // The faster right-tail approximation is omitted here.
		  s = d*d*n; 
		  if(s > 7.24 || (s > 3.76 && n > 99)) 
			  return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
	   
	   k = (int) (n * d) + 1;
	   m = 2 * k - 1;
	   h = k - n * d;
	   H = (double*) Calloc2(m * m, double);
	   Q = (double*) Calloc2(m * m, double);
	   
	   cout << 102 << endl;
	   
	   for(i = 0; i < m; i++)
		   for(j = 0; j < m; j++)
		   if(i - j + 1 < 0)
			   H[i * m + j] = 0;
		   else
			   H[i * m + j] = 1;
			   
		cout << 103 << endl;
		
	   for(i = 0; i < m; i++) {
		   H[i * m] -= R_pow_di(h, i + 1);
		   H[(m - 1) * m + i] -= R_pow_di(h, (m - i));
	   }
	   
	   cout << 104 << endl;
	   
	   H[(m - 1) * m] += ((2 * h - 1 > 0) ? R_pow_di(2 * h - 1, m) : 0);
	   for(i = 0; i < m; i++)
		   for(j = 0; j < m; j++)
		   if(i - j + 1 > 0)
			   for(g = 1; g <= i - j + 1; g++)
			   H[i * m + j] /= g;
	   eH = 0;
	   
	   cout << 105 << endl;
	   
	   m_power(H, eH, Q, &eQ, m, n);
	   s = Q[(k - 1) * m + k - 1];
	   
	   cout << 106 << endl;
	   
	   for(i = 1; i <= n; i++) {
		   s = s * i / n;
		   if(s < 1e-140) {
		   s *= 1e140;
		   eQ -= 140;
		   }
	   }
	   cout << 104 << endl;
	   
	   s *= R_pow_di(10.0, eQ);
	   Free(H);
	   Free(Q);
	   return(s);
	}
	
	inline long_double getPValueFast(const long sampleSize, const long_double d){
		
		return 2*exp( (-2)* sampleSize*(d*d));
	}
	
	inline long_double getPValue(const long sampleSize, const long_double d){
	
		return 1.0-K(sampleSize, d);
	}
}


class BasicAggregator{

	long firstMaxIndex = -1;
	long lastMaxIndex = -1;
	long double minValue = std::numeric_limits<long double>::max();
	
	long firstMinIndex = -1;
	long lastMinIndex = -1;
	
	long double maxValue = std::numeric_limits<long double>::lowest();
	
	long double vsum = 0;
	long double vprod = 1;
	
	long vsize = 0;
	
	long vmaxCount = 0;
	long vminCount = 0;

public:
	
	
	inline long count() const{
	
		return vsize;
	}
	
	inline void clear(){
	
		vsum = 0;
		vprod = 1;
		vsize = 0;
		
		vmaxCount = 0;
		vminCount = 0;
		
		firstMaxIndex = -1;
		lastMaxIndex = -1;
		firstMinIndex = -1;
		lastMinIndex = -1;
		minValue = std::numeric_limits<long double>::max();
		maxValue = std::numeric_limits<long double>::lowest();
	}
	
	inline long double min(long double def = -1) const{
	
		return vsize == 0 ? def : minValue;;
	}
	
	
	
	inline long double max(long double def = -1) const{
	
		return vsize == 0 ? def : maxValue;
	}
	
	
	inline long argMinFirst(long def = -1) const{
	
		return vsize == 0 ? def : firstMinIndex;
	}
	
	inline long argMinLast(bool first = true, long def = -1) const{
	
		return vsize == 0 ? def : lastMinIndex;
	}
	
	inline long argMaxFirst(long def = -1) const{
	
		return vsize == 0 ? def : firstMaxIndex;
	}
	
	inline long argMaxLast(bool first = true, long def = -1) const{
	
		return vsize == 0 ? def : lastMaxIndex;
	}
	
	
	
	inline long argMin(bool first = true, long def = -1) const{
	
		return first ? argMinFirst(def) : argMinLast(def);
	}
	
	inline long argMax(bool first = true, long def = -1) const{
	
		return first ? argMaxFirst(def) : argMaxLast(def);
	}
	
	inline long maxCount() const{
	
		return vmaxCount;
	}
	
	inline long minCount() const{
	
		return vminCount;
	}
	
	
	
	inline long double sum(long double def = 0) const{
	
		return vsize == 0 ? def : vsum;
	}
	
	inline long double prod(long double def = 1) const{
	
		return vsize == 0 ? def : vprod;;
	}
	
	inline long double average(long double def = 0) const{
	
		return vsize == 0 ? def : vsum/vsize;
	}
	
	inline void add(long double val){
	
	
		add( vsize, val);
	}
	
	inline void add( long j, long double val){
	
		const bool wasEmpty = vsize == 0;
		++vsize;
	
		if (wasEmpty || val == minValue){
			lastMinIndex = j;
			
			
			++vminCount;
		}
	
		if (wasEmpty || val < minValue){
		
			firstMinIndex = j;
			lastMinIndex = j;
			
			minValue = val;
			
			vminCount = 1;
		}
		
		vsum += val;
		vprod *= val;
		
		if (wasEmpty || val == maxValue){
			lastMaxIndex = j;
			++vmaxCount;	
		}
	
		if (wasEmpty || val > maxValue){
		
			firstMaxIndex = j;
			lastMaxIndex = j;
			maxValue = val;
			
			vmaxCount = 1;
		}
	}
	

};



#endif