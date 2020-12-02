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

#ifndef SHEKELYAN_UTILGRID_H
#define SHEKELYAN_UTILGRID_H

#include <approxdata/utils/utils.h>



class LongCoords: public std::binary_function<const long, const long, bool>{

	private:
	
		const int DIMS;
		
		vector<long> dimMasks;
		
		vector<long> dimAntiMasks;
		
		vector<int> dimOffsets;
		
		long oneCoord;
		

		int sortDim;

	public:
		
		inline LongCoords(int dims) : DIMS(dims), dimMasks(dims), dimAntiMasks(dims), dimOffsets(dims){
		
			assert (dims > 0);
		
			int pos = 0;
			
			const int w = 63/DIMS;
			
			for (int i = 0; i < DIMS; i++){
			
				dimMasks[i] = UtilBits::mask(pos, w);
				dimAntiMasks[i] = UtilBits::mask(0, w*DIMS) ^ UtilBits::mask(pos, w);
				
				
				
				dimOffsets[i] = pos;
 				
 				pos += w;
			}
			
			oneCoord = 0;
			
			for (int i = 0; i < DIMS; i++)
				oneCoord = setDimCoord(oneCoord, i, 1);
		}
		
		inline int getDims() const{
			
			return DIMS;
		}
		
		inline bool operator()(const long& a, const long& b) const{
		
			return getDimCoord(a, sortDim) < getDimCoord(b, sortDim);
		}
		
		
		inline void sort(vector<long>& v, int i){
		
			sortDim = i;
		
			std::sort(v.begin(), v.end(), (*this));
			
		
		}
		
		inline long setAllCoords(long x) const{
		
			if (x < 0)
				return -1;
		
			long ret = 0;
			
			for (int i = 0; i < DIMS; i++)
				ret = setDimCoord(ret, i, x);
				
			return ret;
		}
		
		inline long getMaxima(long c, long x) const{
		
			long ret = 0;
			
			for (int i = 0; i < DIMS; i++)
				ret = setDimCoord(ret, i, UtilMath::maxVal<long>( getDimCoord(c, i),  getDimCoord(x, i)) );
			
			return ret;
		}
		
		/*
		inline long addInSelDims(long c, int dimFlags, long v){
		
			long ret = c;
			
			for (int i = 0; i < DIMS; i++){
			
				if ( (dimFlags & (1 << i)) )
					ret = setDimCoord(ret, i, getDimCoord(ret, i)+v);
			}
			
			return ret;
		}
		
		inline long mulInSelDims(long c, int dimFlags, long v){
		
			long ret = c;
			
			for (int i = 0; i < DIMS; i++){
			
				if ( (dimFlags & (1 << i)) )
					ret = setDimCoord(ret, i, getDimCoord(ret, i)*v);
			}
			
			return ret;
		}*/
		
		inline long divCoords(long c, long v){
			
			for (int i = 0; i < DIMS; i++)
				c = setDimCoord(c, i, getDimCoord(c, i)/v);
			
			return c;
		}
		
		inline long divInSelDims(long c, int dimFlags, long v){
		
			long ret = c;
			
			for (int i = 0; i < DIMS; i++){
			
				if ( (dimFlags & (1 << i)) )
					ret = setDimCoord(ret, i, getDimCoord(ret, i)/v);
			}
			
			return ret;
		}
		
		inline long addInNonSelDims(long c, int dimFlags, long v){
		
			long ret = c;
			
			for (int i = 0; i < DIMS; i++){
			
				if ( !(dimFlags & (1 << i)) )
					ret = setDimCoord(ret, i, getDimCoord(ret, i)+v);
			}
			
			return ret;
		}
		
		inline long mulInNonSelDims(long c, int dimFlags, long v){
		
			long ret = c;
			
			for (int i = 0; i < DIMS; i++){
			
				if ( !(dimFlags & (1 << i)) )
					ret = setDimCoord(ret, i, getDimCoord(ret, i)*v);
			}
			
			return ret;
		}
		
		inline long divInNonSelDims(long c, int dimFlags, long v){
		
			long ret = c;
			
			for (int i = 0; i < DIMS; i++){
			
				if ( !(dimFlags & (1 << i)) )
					ret = setDimCoord(ret, i, getDimCoord(ret, i)/v);
			}
			
			return ret;
		}
		
		inline long getNext(long c, long min, long max) const{
		
			assert (c >= min);
			assert (c <= max);
			assert (max >= min);
		
		
			for (int i = 0; i < DIMS; i++){
			
				if ( (c & dimMasks[i]) < (max & dimMasks[i]) ){
				
					return c + (oneCoord & dimMasks[i]);
					
				} else {
				
					c &= dimAntiMasks[i];
					c |= (min & dimMasks[i]);
				}	
			}
			
			return max+1;
		}
		
		inline bool isBorder(long c, long min, long max) const{
		
			for (int j = 0; j < DIMS; j++){
					
				const long cj = c & dimMasks[j];
				const long minj = min & dimMasks[j];
				const long maxj = max & dimMasks[j];
				
				if ( (cj == minj) || (cj == maxj) )
					return true;
			}
			
			return false;
		}
		
		inline int getBorderCount(long c, long min, long max) const{
		
			int ret = 0;
		
			for (int i = 0; i < DIMS; i++){
					
				const long mi = dimMasks[i];
				const long ci = c & mi;
				const long mini = min & mi;
				const long maxi = max & mi;
				
				if ( (ci == mini) || (ci == maxi) )
					ret++;
			}
			
			return ret;
		}
		
		inline long getBorderFlags(long c, long min, long max) const{
		
			int ret = 0;
		
			for (int j = 0; j < DIMS; j++){
					
				const long cj = c & dimMasks[j];
				const long minj = min & dimMasks[j];
				const long maxj = max & dimMasks[j];
				
				if ( (cj == minj) || (cj == maxj) )
					ret |= 1 << j;
			}
			
			return ret;
		}
		
		inline long getNextBorder(long c, long min, long max) const{
		
			int borders = 0;
			
			for (int i = 0; i < DIMS; i++){
			
				const long mi = dimMasks[i];
				const long ai = dimAntiMasks[i];
				const long ci = c & mi;
				const long mini = min & mi;
				const long maxi = max & mi;
				
				if (ci == maxi){
				
					c &= ai;
					c |= mini;
					
					borders++;
					
					continue;
				}
				
				if (borders == 0)
				for (int j = i+1; j < DIMS; j++){
				
					const long mj = dimMasks[j];
					const long cj = c & mj;
					const long minj = min & mj;
					
					if (cj == minj){
					
						borders++;
						break;
					}
					
					const long maxj = max & mj;
					
					if (cj == maxj){
					
						borders++;
						break;
					}
				}
				
				const bool isOnBorder = borders != 0;
				
				if (isOnBorder)
					return c + (oneCoord & mi);
				else
					return (c & ai) |maxi;
			}
			
			return max+1;
		}
		
		inline long getMinima(long c, long x) const{
		
			long ret = 0;
			
			for (int i = 0; i < DIMS; i++)
				ret = setDimCoord(ret, i, UtilMath::minVal<long>( getDimCoord(c, i),  getDimCoord(x, i)) );
			
			return ret;
		}
		
		inline long decAllCoords(long x) const{
		
			if (x == -1)
				return -1;
		
			for (int i = 0; i < DIMS; i++)
				if ( (x & dimMasks[i]) == 0)
					return -1;
			
			return x-oneCoord;
		}
		
		
		inline long decSelCoords(long x, long selectedDims) const{
		
			if (x == -1)
				return -1;
			
			if (UtilLongSet::isEmpty(selectedDims))
				return x;
		
			long one = oneCoord;
			
			for (int i = 0; i < DIMS; i++){
			
				if (UtilLongSet::contains(selectedDims, i) ){
				
					if ( (x & dimMasks[i]) == 0)
						return -1;
					
				} else {
				
					one &= dimAntiMasks[i];
				}
			}
			
			return x-one;
		}
		
		
		inline long incCoords(long x, long skipDims = 0) const{
		
			if (x == -1)
				return -1;
		
			long one = oneCoord;
			
			for (int i = 0; i < DIMS; i++){
			
				if ( UtilLongSet::contains(skipDims, i) ){
				
					one &= dimAntiMasks[i];
					
				} else {
				
					if ( (x & dimMasks[i]) == dimMasks[i])
						return -1;
				}
			}
			
			return x+one;
		}
		
		inline long decCoords(long x, long skipDims = 0) const{
		
			if (x == -1)
				return -1;
		
			long one = oneCoord;
			
			for (int i = 0; i < DIMS; i++){
			
				if ( UtilLongSet::contains(skipDims, i) ){
				
					one &= dimAntiMasks[i];
					
				} else {
				
					if ( (x & dimMasks[i]) == 0)
						return -1;
				}
			}
			
			return x-one;
		}
		
		inline long incSelCoords(long x, long selectedDims) const{
		
			if (x == -1)
				return -1;
				
			if (UtilLongSet::isEmpty(selectedDims) )
				return x;
		
			long one = oneCoord;
			
			for (int i = 0; i < DIMS; i++){
			
				if ( UtilLongSet::contains(selectedDims, i) ){
				
					if ( (x & dimMasks[i]) == dimMasks[i])
						return -1;
					
				} else {
				
					one &= dimAntiMasks[i];
				}
			}
			
			return x+one;
		}
		
		
		inline long incAllCoords(long x) const{
			
			if (x == -1)
				return -1;
			
			for (int i = 0; i < DIMS; i++)
				if ( (x & dimMasks[i]) == dimMasks[i])
					return -1;
			
			return x+oneCoord;
		}
		
		inline bool contains(long imin, long imax, long x) const{
		
			if (imin == -1)
				return false;
				
			if (imax == -1)
				return false;
		
			for (int i = 0; i < DIMS; i++){
			
				const long m = dimMasks[i];
				
				const long xm = x & m;
				
				if ( xm < (imin & m) )
					return false;
				
				if ( xm > (imax & m) )
					return false;
			}
			
			return true;
		}
		
		inline long getDimCoord(long x, int i) const{
		
			if (x == -1)
				return -1;
		
			return (x & dimMasks[i]) >> dimOffsets[i];
		}
		
		inline long setDimCoord(long x, int i, long v) const{
		
			if ( (x < 0) || (v < 0))
				return -1;
		
			return (x & dimAntiMasks[i]) |(v << dimOffsets[i] );
		
			//return (v << dimOffsets[i]) |UtilBits::setZeroInMask(x, dimMasks[i]);
		}
		
		inline long getSum(long x) const{
		
			long ret = 0;
		
			for (int i = 0; i < DIMS; i++){
			
				ret += getDimCoord(x, i);
			}
			
			return ret;
		
		}
		
		inline long linearDimCoord(long x, int i, long m, long c) const{
		
			return setDimCoord(x, i, getDimCoord(x, i)*m+c);
		}
		
		inline long setDimCoordBetween(long x, int i, long v1, long v2) const{
		
			if ( (x < 0) || (v2 < 0))
				return -1;
			
			const long v = getDimCoord(x, i);
			
			if (v < v1)
				return setDimCoord(x, i, v1);
			
			if (v > v2)
				return setDimCoord(x, i, v2);
			
			return x;
		}
		
		inline long setBetween(long x, long s, long t) const{
		
			if ( (x < 0) || (s < 0) || (t < 0))
				return -1;
			
			for (int i = 0; i < DIMS; i++){
			
				const long m = dimMasks[i];
				
				const long xm = x & m;
				
				if ( xm < (s & m) ){
					
					x &= dimAntiMasks[i];
					x |= (s & dimMasks[i]);
				
				} else if ( xm > (t & m) ){
				
					x &= dimAntiMasks[i];
					x |= (t & dimMasks[i]);
				}
			
			}
			
			return x;
		}
		
		/*
		inline long getArea(long i1, long i2) const{
		
			if (UtilMath::isIn(-1L, i1, i2))
				return -1;
			
			long ret = 1;
		
			for (int i = 0; i < DIMS; i++){
			
				const long a = getDimCoord(i1, i);
				const long b = getDimCoord(i2, i);
			
				if (a > b)
					return 0;
			
				ret *= 1+b-a;
			}
			
			return ret;
		}*/
		
		inline long getArea(long s, long t) const{
		
			if (UtilMath::isIn(-1L, s, t))
				return -1;
		
			long ret = 1;
			
			for (int i = 0; i < DIMS; i++){
				
				if (getDimCoord(t, i) < getDimCoord(s,i) )
					return -1;
			
				ret *= 1+getDimCoord(t, i)-getDimCoord(s, i);
				
				if (ret < 0)
					return -1;
			}
			
			return ret;
		}
		
		inline long getInsideArea(long s, long t) const{
		
			if (UtilMath::isIn(-1L, s, t))
				return -1;
		
			long ret = 1;
		
			for (int i = 0; i < DIMS; i++){
				ret *= -1+getDimCoord(t, i)-getDimCoord(s, i);
				
				if (ret < 0)
					return 0;
			}
			
			return ret;
		}
		
		inline void getCoords(long x, vector<long>& ind) const{
		
			for (int i = 0; i < DIMS; i++)
				ind[i] = getDimCoord(x, i);
		}
		
		
		inline long get2d(long x, long y){
		
			long ret = 0;
			
			ret = setDimCoord(ret, 0, x);
			ret = setDimCoord(ret, 1, y);
		
			return ret;
		}
		
		inline long get3d(long x, long y, long z){
		
			long ret = 0;
			
			ret = setDimCoord(ret, 0, x);
			ret = setDimCoord(ret, 1, y);
			ret = setDimCoord(ret, 2, z);
		
			return ret;
		}
		
		inline long getIndex(long x, long y) const{
		
			return (x << dimOffsets[0]) | (y << dimOffsets[1] ); 
		}
		
		inline long getIndex(long x, long y, long z) const{
		
			return (x << dimOffsets[0]) | (y << dimOffsets[1]) | (z << dimOffsets[2]); 
		}
		
		inline long getIndex(const vector<long>& ind) const{
		
			long ret = 0;
		
			for (int i = 0; i < DIMS; i++){
			
				ret |= ind[i] << dimOffsets[i];
			}	
			
			return ret;
		}

	
		inline const string toString(long x) const{
		
			stringstream s;
			
			for (int i = 0; i < (DIMS-1); i++)
				s << to_string( getDimCoord(x, i) ) << ", ";
			
			s << to_string(getDimCoord(x, DIMS-1));
			
			return s.str();
		}

};

template<bool rowMajor=true, bool bitshift=false>
class MultiDimIndex{

private:

	vector<long> temp;

	vector<long> lens;	
	vector<long> muls;
	vector<int> shifts;	
	
	
	long totalLength = 1;
	
	int log2TotalLength = 0;
	
	inline void prepare(){
		
		temp.clear();
		
		log2TotalLength = 0;
		
		for (int i = 0; i < lens.size(); i++)
			temp.push_back(0);
		
		if (bitshift){
			for (int i = 0; i < lens.size(); i++)
				lens[i] = UtilMath::powerOfTwoUB(lens[i]);
		}
		
		long mul = 1;
		
		if (rowMajor){
			
			for (int i = lens.size()-1; i >= 0; i--){
		
				muls[i] = mul;
				mul *= lens[i];
				
				log2TotalLength += UtilMath::getPowerOfTwo(lens[i]);
			}
			
		} else {
			
			for (int i = 0; i < lens.size(); i++){
		
				muls[i] = mul;
				mul *= lens[i];
				
				log2TotalLength += UtilMath::getPowerOfTwo(lens[i]);
			}
		}
		
		totalLength = mul;
		
		if (bitshift){
			for (int i = 0; i < lens.size(); i++){
				shifts[i] = UtilMath::getPowerOfTwo(muls[i]);
			}
		}
		
		//cout << UtilString::listToString(muls) << endl;
	}

public:

	inline MultiDimIndex(const vector<long> lengths): lens(lengths.size()), muls(lengths.size() ), shifts(lengths.size()){
	
		lens = lengths;
		
		prepare();
	}
	
	inline int getDims() const{
	
		return lens.size();
	}	
	
	inline long getLength(int i) const{
	
		return lens[i];
	}
	
	inline long getLength() const{
	
		return totalLength;
	}
	
	inline long getLog2Length() const{
	
		return log2TotalLength;
	}
	
	inline vector<long> getLengths() const{
	
		return lens;
	}
	
	inline MultiDimIndex(int dims): lens(dims), muls(dims), shifts(dims){
	
		for (int i = 0; i < dims; i++)
			lens[i] = 1L << (62/dims);
		
		prepare();
	}
	
	
	
	inline void setWidth(int i, long k){
	
		lens[i] = k;
		prepare();
	}
	
	
	inline void setAllWidths(long k){
	
		for (int i = 0; i < lens.size(); i++)
			lens[i] = k;
		prepare();
	}
	/*
	inline MultiDimIndex(long x): lens(1), muls(1), shifts(1){
	
		lens[0] = x;
		
		prepare();
	}*/
	
	inline MultiDimIndex(long x, long y): lens(2), muls(2), shifts(2){
	
		lens[0] = x;
		lens[1] = y;
		
		prepare();
	}
	
	inline MultiDimIndex(long x, long y, long z): lens(3), muls(3), shifts(3){
	
		lens[0] = x;
		lens[1] = y;
		lens[2] = z;
		
		prepare();
	}
	
	inline vector<long> begin() const{
	
		vector<long> v = lens;
			
		for (int i = 0; i < lens.size(); ++i)
			v[i] = 0;
			
		return v;
	}
	
	template<typename T>
	inline void begin(vector<T>& v) const{
	
		v.clear();
		increment(v);
	}
	
	template<typename T>
	inline bool endReached(vector<T>& v) const{
	
		
	
		if (v.size() == 0)
			return true;
		
		/*
		assert (v.size() == lens.size() );
		
		const bool order = !rowMajor;
		
		for (int ii = 0; ii < lens.size(); ++ii){
		
			const int i = order ? ii : lens.size()-1-ii;
		
			if (v.at(i) < (lens[i]-1))
				return false;
		}*/
	
		return false;
	}
	
	template<typename T>
	inline void increment(vector<T>& v) const{
	
		if (v.size() == 0){
		
			
			v = lens;
			
			for (int i = 0; i < lens.size(); ++i)
				v.at(i) = 0;
			
			
			
			return;
		}
		
		assert (v.size() == lens.size() );
		
		const bool order = !rowMajor;
		
		for (int ii = 0; ii < lens.size(); ++ii){
		
			const int i = order ? ii : lens.size()-1-ii;
		
			if (v.at(i) < (lens[i]-1)){
			
				++v.at(i);
				return;
			}
			
			v.at(i) = 0;
		}
		
		v.clear();
	
		return;
	}
	
	inline bool contains(long max, long x) const {
		
		if (rowMajor){
		
			if (bitshift){
			
				for (int i = 0; i < lens.size(); i++){
		
					const long maxi = max >> shifts[i];
					const long xi = x >> shifts[i];
					
					if (xi > maxi)
						return false;
					
					max -= maxi << shifts[i];
					x -= xi << shifts[i];
				}
				
			} else {
			
				for (int i = 0; i < lens.size(); i++){
		
					const long maxi = max/muls[i];
					const long xi = x/muls[i];
					
					if (xi > maxi)
						return false;
					
					max -= maxi * muls[i];
					x -= xi * muls[i];
				}
			}
		
		} else {
		
			if (bitshift){
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					const long maxi = max >> shifts[i];
					const long xi = x >> shifts[i];
					
					if (xi > maxi)
						return false;
					
					max -= maxi << shifts[i];
					x -= xi << shifts[i];
				}
				
			} else {
				
				for (int i = lens.size()-1; i >= 0; i--){
		
					const long maxi = max/muls[i];
					const long xi = x/muls[i];
					
					if (xi > maxi)
						return false;
					
					max -= maxi * muls[i];
					x -= xi * muls[i];
				}
			}	
		}
		
		return true;
	}
	
	inline bool contains(long min, long max, long x) const {
		
		if (rowMajor){
		
			if (bitshift){
			
				for (int i = 0; i < lens.size(); i++){
		
					const long mini = min >> shifts[i];
					const long maxi = max >> shifts[i];
					const long xi = x >> shifts[i];
					
					if (xi > maxi || xi < mini)
						return false;
					
					min -= mini << shifts[i];
					max -= maxi << shifts[i];
					x -= xi << shifts[i];
				}
				
			} else {
			
				for (int i = 0; i < lens.size(); i++){
		
					const long mini = min/muls[i];
					const long maxi = max/muls[i];
					const long xi = x/muls[i];
					
					if (xi > maxi || xi < mini)
						return false;
					
					min -= mini * muls[i];
					max -= maxi * muls[i];
					x -= xi * muls[i];
				}
			}
		
		} else {
		
			if (bitshift){
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					const long mini = min >> shifts[i];
					const long maxi = max >> shifts[i];
					const long xi = x >> shifts[i];
					
					if (xi > maxi || xi < mini)
						return false;
					
					min -= mini << shifts[i];
					max -= maxi << shifts[i];
					x -= xi << shifts[i];
				}
				
			} else {
				
				for (int i = lens.size()-1; i >= 0; i--){
		
					const long mini = min/muls[i];
					const long maxi = max/muls[i];
					const long xi = x/muls[i];
					
					if (xi > maxi || xi < mini)
						return false;
					
					min -= mini * muls[i];
					max -= maxi * muls[i];
					x -= xi * muls[i];
				}
			}	
		}
		
		return true;
	}
	
	inline long getIndex(long x) const{
	
		return x;
	}
	
	inline long getIndex(long x, long y) const{
	
		return rowMajor ? (bitshift ? ( y+(x << shifts[0]) ) : (y+x*muls[0]) ) : 
			(bitshift ? ( x+(y << shifts[1]) ) : (x+y*muls[1]) );
	}
	
	inline long getIndex(long x, long y, long z) const{
		
		return rowMajor ? (bitshift ? ( z+(y << shifts[1])+(x << shifts[0]) ) : (z+y*muls[1]+x*muls[0]) ):
			(bitshift ? ( x+(y << shifts[1])+(z << shifts[2]) ) : (x+y*muls[1]+z*muls[2]) );
	}
	
	inline long getIndex(const LongCoords& lc, long coords) const{
	
		long ret = 0;
		
		if (bitshift){
		
			for (int i = 0; i < lens.size(); i++)
				ret += lc.getDimCoord(coords, i) << shifts[i];
		
		} else {
		
			for (int i = 0; i < lens.size(); i++)
				ret += lc.getDimCoord(coords, i)*muls[i];
		}
		
		return ret;
		
	}
	
	inline long getIndex(const vector<long>& coords) const{
		
		long ret = 0;
		
		if (bitshift){
		
			for (int i = 0; i < lens.size(); i++)
				ret += coords[i] << shifts[i];
		
		} else {
		
			for (int i = 0; i < lens.size(); i++)
				ret += coords[i]*muls[i];
		}
		
		return ret;
	}
	
	
	inline long getCoords(const LongCoords& lc, long index) const {
	
		long ret = 0;

		if (rowMajor){
		
			if (bitshift){
			
				for (int i = 0; i < lens.size(); i++){
		
					ret = lc.setDimCoord(ret, i, index >> shifts[i]); 
					index -= lc.getDimCoord(ret, i) << shifts[i];
				}
			} else {
			
				for (int i = 0; i < lens.size(); i++){
		
					ret = lc.setDimCoord(ret, i, index/muls[i]);
					index -= lc.getDimCoord(ret, i)*muls[i];
				}
			}
		
		} else {
		
			if (bitshift){
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					ret = lc.setDimCoord(ret, i, index >> shifts[i]); 
					index -= lc.getDimCoord(ret, i) << shifts[i];
				}
				
			} else {
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					ret = lc.setDimCoord(ret, i, index/muls[i]);
					index -= lc.getDimCoord(ret, i)*muls[i];
				}
			}
		}
		
		return ret;
	}

	inline void getCoords(long index, vector<long>& coords) const{
	
		coords = lens;
		
		if (rowMajor){
		
			if (bitshift){
			
				for (int i = 0; i < lens.size(); i++){
		
					coords[i] = index >> shifts[i];
					index -= coords[i] << shifts[i];
				}
			} else {
			
				for (int i = 0; i < lens.size(); i++){
		
					coords[i] = index/muls[i];
					index -= coords[i]*muls[i];
				}
			}
		
		} else {
		
			if (bitshift){
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					coords[i] = index >> shifts[i];
					index -= coords[i] << shifts[i];
				}
				
			} else {
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					coords[i] = index/muls[i];
					index -= coords[i]*muls[i];
				}
			}
		}
	}
};

template<typename E, bool FAST=false>
class MultiArray{
public:
	MultiDimIndex<true, FAST> m;
	
	vector<E> vec;
	long max;
	
	long wx;
	long wy;
	long wz;
	
	vector<long> lengths;
	
	inline const vector<long> getLengths() const{
	
		return lengths;
	}
	
	inline MultiArray(long wx) : m(wx), wx(wx), wy(1), wz(1){
	
		vec.reserve(m.getLength());
		max = m.getIndex(wx-1);
		
		lengths.push_back(wx);
	}
	
	inline MultiArray(long wx, long wy) : m(wx,wy), wx(wx), wy(wy), wz(1){
		
		vec.reserve(m.getLength());
		max = m.getIndex(wx-1, wy-1);
		
		lengths.push_back(wx);
		lengths.push_back(wy);
	}
	
	inline MultiArray(long wx, long wy, long wz) : m(wx,wy,wz), wx(wx), wy(wy), wz(wz){
		
		vec.reserve(m.getLength());
		max = m.getIndex(wx-1, wy-1, wz-1);
		
		lengths.push_back(wx);
		lengths.push_back(wy);
		lengths.push_back(wz);
	}
	
	
	
	inline long getLengthX() const{
	
		return wx;
	}
	
	inline long getLengthY() const{
	
		return wy;
	}
	
	inline long getLengthZ(){
	
		return wz;
	}
	
	
	inline E& atIndex(long ind) const{
	
		return vec.at(ind);
	}
	
	inline void clear(){
	
		vec.clear();
	}
	
	
	
	inline void fill(E e){
		
		std::fill(vec.begin(), vec.end(), e);
		
		while (vec.size() < m.getLength() )
			vec.push_back(e);
	}
	
	inline long getLength(int i) const{
	
		return m.getLength(i);
	}
	
	inline E& at(long x) {
	
		return vec.at( m.getIndex(x) );
	}
	
	inline E& at(long x, long y) {
		
		return vec.at( m.getIndex(x,y) );
	}
	
	
	inline void set(long x, long y, const E& v) {
		
		vec.at(m.getIndex(x,y)) = v;
	}
	
	
	
	inline E get(long x, long y) const {
		
		return vec.at(m.getIndex(x,y));
	}
	
	
	inline void print() const{
	
		cout << endl;
	
		for (long y = 0; y < getLengthY(); y++){
		
			for (long x = 0; x < getLengthX(); x++){
		
				cout << "\t " << this->get(x,y);
			}
			cout << endl;
		}
	}
	
	inline E& at(long x, long y, long z) {
	
		return vec.at( m.getIndex(x,y,z) );
	}
	
	inline void push_back(E e){
		
		if (FAST){
		
			while (!m.contains(max, vec.size() ))
				vec.push_back(e);
		}
		
		vec.push_back(e);
	}
	
};


template<typename E>
class DenseGridHistogram{
public:
	vector<E> counts;
	MultiDimIndex<> ind;
	
	vector<vector<uint8_t>> countBytes; // compressed representation of counts to avoid wasting bytes
	
	LongCoords lc;
	
	const long cellsPerDim;
	
	bool cumulative = false;
	
	
	inline DenseGridHistogram(int dims, long m, bool alloc=true) : cellsPerDim(m), ind(dims), lc(dims){
	
		assert (dims > 0);
	
		for (int i = 0; i < dims; i++)
			ind.setWidth(i, m);	
		
		if (alloc)
			allocate();
	}
	
	inline DenseGridHistogram(const vector<long>& v, bool alloc=true) : cellsPerDim(-1), ind(v), lc(v.size()){
		
		if (alloc)
			allocate();
	}
	
	inline bool equals(const DenseGridHistogram& h){
	
		const long maxCoord = getMaxCoord();
		
		for (int i = 0; i < ind.getDims(); i++){
			
			if (ind.getLength(i) != h.ind.getLength(i) )	
				return false;
		}
	
		assert (countBytes.size() +h.countBytes.size() == 0);
		
		if (counts.size() != h.counts.size() )
			return false;
		
		if (cumulative != h.cumulative){
		
			for (long l = 0; l <= maxCoord; l = lc.getNext(l, 0, maxCoord) ){
			
				if (sum(l,l) != h.sum(l,l) )
					return false;
			}
			
			return true;
		}
		
		for (long j = 0; j < counts.size(); j++)
			if (counts[j] != h.counts[j])
				return false;
		
		return true;
	}
	
	
	inline long getMinCoord() const{
	
		return 0;
	}
	
	inline long getMaxCoord() const{
	
		long c = 0;
	
		for (int i = 0; i < ind.getDims(); i++)
			c = lc.setDimCoord(c, i, ind.getLength(i)-1);
		
		return c;
	}
	
	inline bool contains(const long& l) const{
	
		return lc.contains(getMinCoord(), getMaxCoord(), l);
	}
	
	// ind.getLength(0) times ind.getLength(1) ... grid is spanned between (0,..,0) and (1,...,1)
	inline long getCoord(const Point& p) const{
	
		long c = 0;
	
		for (int i = 0; i < ind.getDims(); i++)
			c = lc.setDimCoord(c, i, UtilHist::discretize(p[i], ind.getLength(i)));
		
		return c;
	}
	
	inline E debugCount(const Point& qmin, const Point& qmax, const vector<Point>& pts, QueryMode mode=QueryMode::EST, long min = -1, long max = -1, long amin = 0, long amax = 0, long emin = 0, long emax = 0) const{
	
			
		//const long double cellVol = 1.0L/ind.getLength();
		const int DIMS = ind.getDims();
		
		for (int i = 0; i < DIMS; i++){
			assert (qmax[i] >= qmin[i]);
			
			if ( (min != -1)  && (max != -1) )
				assert ( lc.getDimCoord(min, i) <= lc.getDimCoord(max, i) );
			
		}
		
		Box box(qmin, qmax);
		
		if ( (min != -1) && (max != -1) && lc.getArea(min, max) <= 0){
		
			
			for (auto it = pts.begin(); it != pts.end(); it++){
			
				assert (!box.contains(*it) );
			}
			
			return 0;
		}
		
		long ubMin = 0; //getCoord(qmin);
		long ubMax = 0; //getCoord(qmax);
		
		long anchoredDimsMin = amin;
		long anchoredDimsMax = amax;
		
		cout << "qmin " << qmin << endl;
		cout << "qmax " << qmax << endl;
		
		cout << "cellsPerDim " << cellsPerDim << endl;
		
		for (int i = 0; i < DIMS; i++){
			
			const long LEN = ind.getLength(i);
			
			const long mini = (min == -1) ? 0 : lc.getDimCoord(min, i);
			const long maxi = (max == -1) ? LEN-1 : lc.getDimCoord(max, i);
			
			assert (maxi >= mini);
			
			const double bmini = UtilHist::bucketMin( mini, LEN );
			const double bmaxi = UtilHist::bucketMin( maxi+1, LEN );
			
			assert (bmaxi > bmini);
			
			if (qmin[i] >= bmaxi ){
				
				for (auto it = pts.begin(); it != pts.end(); it++){
			
					assert (!box.contains(*it) );
				}
				
				// NO INTERSECTION
				return 0;
			
			} else if ( qmin[i] < bmini ){
			
				ubMin = lc.setDimCoord(ubMin, i, mini);
				anchoredDimsMin = UtilLongSet::add(anchoredDimsMin, i);
				
			} else {
			
				ubMin = lc.setDimCoord(ubMin, i, UtilHist::discretize(qmin[i], LEN ) );
			}
			
			if (qmax[i] < bmini ){
			
				for (auto it = pts.begin(); it != pts.end(); it++){
			
					assert (!box.contains(*it) );
				}
			
				// NO INTERSECTION
				return 0;
			
			} else if ( qmax[i] >= bmaxi ){
			
				ubMax = lc.setDimCoord(ubMax, i, maxi);
				anchoredDimsMax = UtilLongSet::add(anchoredDimsMax, i);
				
			} else {
			
				ubMax = lc.setDimCoord(ubMax, i, UtilHist::discretize(qmax[i], LEN ) );
			}
		}
		
		cout << "ubMin " << lc.toString(ubMin) << endl;
		cout << "ubMax " << lc.toString(ubMax) << endl;
		
		const long lbMin = lc.incCoords( ubMin, anchoredDimsMin);
		const long lbMax = lc.decCoords( ubMax, anchoredDimsMax);
		
		cout << "amin " << anchoredDimsMin << endl;
		cout << "amax " << anchoredDimsMax << endl;
		
		cout << "lbMin " << lc.toString(lbMin) << endl;
		cout << "lbMax " << lc.toString(lbMax) << endl;
		
		cout << "GET AREA LB" << endl;
		
		const long lbCells = UtilMath::maxVal<long>(0, lc.getArea(lbMin, lbMax));
		
		cout << "GET AREA UB" << endl;
		
		const long ubCells = UtilMath::maxVal<long>(0, lc.getArea(ubMin, ubMax));
		
		assert (lbCells >= 0);
		assert (ubCells >= lbCells);
		
		cout << "ubCells " << ubCells << " >= " << lbCells << " lbCells" << endl;
		
		assert (ubCells > 0);
		
		//long double qVol = 1.0;
		//long double ubVol = 1.0;
		
		long double qCells = 1.0; //qVol/cellVol;
		
		for (int i = 0; i < DIMS; i++){
		
			const long mini = (min == -1) ? 0 : lc.getDimCoord(min, i);
			const long maxi = (max == -1) ? ind.getLength(i)-1 : lc.getDimCoord(max, i);
			
			assert (maxi >= mini);
			
			const double bmini = UtilHist::bucketMin( mini, ind.getLength(i) );
			const double bmaxi = UtilHist::bucketMax( maxi, ind.getLength(i) );
			
			
			//const double umini = UtilHist::bucketMin( lc.getDimCoord(ubMin, i), ind.getLength(i) );
			//const double umaxi = UtilHist::bucketMax( lc.getDimCoord(ubMax, i), ind.getLength(i) );
			
			
			assert (bmaxi > bmini);
			
			
			const double a = UtilLongSet::contains(anchoredDimsMin, i) ? bmini : qmin[i];
			const double b = UtilLongSet::contains(anchoredDimsMax, i) ? bmaxi : qmax[i];
			
			
			//const double c = UtilLongSet::contains(dims1, i) ? qmin[i] : bmini;
			//const double d = UtilLongSet::contains(dims2, i) ? qmax[i] : bmaxi;
			
			
			if ( (b > UtilMath::minVal(qmax[i], bmaxi)) || (a < UtilMath::maxVal(qmin[i], bmini)) ){
			
				cout << "a " << a << " b " << b << " qmin " << qmin[i] << " qmax " << qmax[i] << " bmini " << bmini << " bmaxi " << bmaxi << endl;
			}
			
			assert (b <= UtilMath::minVal(qmax[i], bmaxi) );
			assert (a >= UtilMath::maxVal(qmin[i], bmini) );
			
			//assert (a >= umini);
			//assert (b <= umaxi );
			
			//const long double w = b-a;
			
			//assert (w == w);
			
			//if (w < 0)
			//	cout << "a " << a << " b " << b << " qmin " << qmin[i] << " *" << UtilLongSet::contains(dims1, i) << " qmax " << qmax[i] << " *" << UtilLongSet::contains(dims2, i) <<" bmini " << bmini << " bmaxi " << bmaxi << endl;
			
			//assert (w >= 0);
			//assert (w <= 1);
			
			qCells *= (b-a)*ind.getLength(i);
			
			if (qCells < 0)
				qCells = 0;
			
			//ubVol *= (umaxi)-(umini);
			
			//if (qVol < 0)
			//	qVol = 0;
		}
		
		
		
		if (qCells > ubCells)
			qCells = ubCells;
			
		if (qCells < lbCells)
			qCells = lbCells;
		
		//const double qVol = qmin.getEnclosedVolume(qmax);
		
		//assert (qVol == qVol);
		//assert (ubVol == ubVol);
		
		const E ubCount = sum(ubMin, ubMax);
		
		
		
		if (mode == QueryMode::UB)
			return ubCount;
		
		
		//const long double ubVol = ubCells*cellVol;
		
		
		//if (qVol > ubVol)
		//	return ubCount; //cout << "qVol " << qVol << " ubVol " << ubVol << endl;
		
		
		//assert (qVol >= 0);
		//assert (ubVol >= qVol);
		
		if (lbCells <= 0){
		
			if (mode == QueryMode::LB)
				return 0;
			
			return UtilMath::makeBetween<double>(0, ubCount, (qCells/ubCells)*ubCount);
		}
		
		/*
		const long double lbVol = lbCells*cellVol;
		
		assert(lbVol == lbVol);*/
		
		//if (lbVol >= ubVol)
		//	return ubCount;
		
		const E lbCount = sum(lbMin, lbMax);
		
		//if (qVol <= lbVol)
		//	return lbCount;
			
		
		assert (lbCount == lbCount);
		
		if (mode == QueryMode::LB)
			return lbCount;
		
		//assert (ubVol > lbVol);
		
		if (ubCells == lbCells)
			return lbCount;
		
		const double ret = lbCount+(ubCount-lbCount)*((qCells-lbCells)/(ubCells-lbCells) );
		
		assert (ret == ret);
		
		const double ret2 = UtilMath::makeBetween<double>(lbCount, ubCount, ret);
		
		assert (ret2 == ret2);
		
		return ret2;
	}
	
	
	inline std::tuple<long,long,long,long> getGridBox(const QueryBox& box) const{
		
		long ubmin = 0;
		long ubmax = 0;
		
		long lbmin = 0;
		long lbmax = 0;
		
		const int DIMS = ind.getDims();
		
		for (int i = 0; i < DIMS; i++){
		
			assert (box.min[i] <= box.max[i]);
		
			if ( box.isMinUnbounded(i)){
			
				ubmin = lc.setDimCoord(ubmin, i, 0);
				
				if (lbmin != -1)
					lbmin = lc.setDimCoord(lbmin, i, 0);
				
			} else {
			
				const double x = box.isMinExcluded(i) ? std::nexttoward(box.min[i], 1.0) : box.min[i];
				const long y = UtilHist::discretize( x, getCellsPerDim() );
				
				if ( y >= getCellsPerDim())
					return {-1, -1, -1, -1};
				
				ubmin = lc.setDimCoord(ubmin, i, y);
				
				if (lbmin != -1 && (y+1) >= getCellsPerDim()){
				
					lbmin = -1;
					lbmax = -1;
				}
				
				if (lbmin != -1)
					lbmin = lc.setDimCoord(lbmin, i, y+1);
			}
			
			if ( box.isMaxUnbounded(i)){
			
				ubmax = lc.setDimCoord(ubmax, i, getCellsPerDim()-1);
				
				if (lbmax != -1)
					lbmax = lc.setDimCoord(lbmax, i, getCellsPerDim()-1);
				
			} else {
			
				const double x = box.isMaxExcluded(i) ? std::nexttoward(box.max[i], 0.0) : box.max[i];
				const long y = UtilHist::discretize(x, getCellsPerDim() );
				
				if ( y < lc.getDimCoord(ubmin, i) )
					return {-1, -1, -1, -1};
				
				ubmax = lc.setDimCoord(ubmax, i, y);
				
				if (lbmin != -1 && (y-1) < lc.getDimCoord(lbmin, i)){
					
					lbmin = -1;
					lbmax = -1;
				}
				
				if (lbmax != -1)
					lbmax = lc.setDimCoord(lbmax, i, y-1);
			}
		}
		
		return std::make_tuple(ubmin,ubmax, lbmin, lbmax);
	}
	
	inline void getAligned(const QueryBox& b, long& lbMin, long& lbMax, long& ubMin, long& ubMax) const{
	
		const int DIMS = ind.getDims();		
		
		
		ubMin = 0;
		ubMax = 0;
	
		for (int i = 0; i < DIMS; i++){
		
			const long LEN = ind.getLength(i);
			
			ubMin = lc.setDimCoord(ubMin, i, UtilHist::discretize(b.min[i], LEN ) );
			ubMax = lc.setDimCoord(ubMax, i, UtilHist::discretize(b.max[i], LEN ) );
		}
		
		lbMin = lc.incCoords( ubMin, b.getUnboundedMin());
		lbMax = lc.decCoords( ubMax, b.getUnboundedMax());
	}
	
	inline E getCount(const QueryBox& box) const{
	
		const int DIMS = ind.getDims();		
		
		long ubMin = 0; 
		long ubMax = 0;
		
		long lbMin = 0;
		long lbMax = 0;
		
		getAligned(box, lbMin, lbMax, ubMin, ubMax);
		
		const E ubCount = lc.getArea(ubMin, ubMax) >= 0 ? sum(ubMin, ubMax): 0;
		const E lbCount = lc.getArea(lbMin, lbMax) >= 0 ? sum(lbMin, lbMax) : 0;	
		
		if (box.mode == QueryMode::LB)
			return lbCount;
			
		if (box.mode == QueryMode::UB)
			return ubCount;
		
		const long lbCells = UtilMath::maxVal<long>(0, lc.getArea(lbMin, lbMax));
		const long ubCells = UtilMath::maxVal<long>(0, lc.getArea(ubMin, ubMax));
		
		assert (lbCells >= 0);
		assert (ubCells >= lbCells);
		
		assert (ubCells > 0);
		
		long double qCells = 1.0;
		
		for (int i = 0; i < DIMS; i++){
		
			const long mini = 0;
			const long maxi = ind.getLength(i)-1;
			
			assert (maxi >= mini);
			
			const double bmini = UtilHist::bucketMin( mini, ind.getLength(i) );
			const double bmaxi = UtilHist::bucketMax( maxi, ind.getLength(i) );
			
			
			//const double umini = UtilHist::bucketMin( lc.getDimCoord(ubMin, i), ind.getLength(i) );
			//const double umaxi = UtilHist::bucketMax( lc.getDimCoord(ubMax, i), ind.getLength(i) );
			
			assert (bmaxi > bmini);
			
			const double a = UtilLongSet::contains(box.getUnboundedMin(), i) ? bmini : box.min[i];
			const double b = UtilLongSet::contains(box.getUnboundedMax(), i) ? bmaxi : box.max[i];
			
			
			if ( (b > UtilMath::minVal(box.max[i], bmaxi)) || (a < UtilMath::maxVal(box.min[i], bmini)) ){
			
				cout << "a " << a << " b " << b << " qmin " << box.min[i] << " qmax " << box.max[i] << " bmini " << bmini << " bmaxi " << bmaxi << endl;
			}
			
			assert (b <= UtilMath::minVal(box.max[i], bmaxi) );
			assert (a >= UtilMath::maxVal(box.min[i], bmini) );
			
			qCells *= (b-a)*ind.getLength(i);
			
			if (qCells < 0)
				qCells = 0;
			
		}
		
		if (qCells > ubCells)
			qCells = ubCells;
			
		if (qCells < lbCells)
			qCells = lbCells;
		
		if (lbCells <= 0)
			return UtilMath::makeBetween<double>(0, ubCount, (qCells/ubCells)*ubCount);
		
		if (ubCells == lbCells)
			return lbCount;
		
		const double ret = lbCount+(ubCount-lbCount)*((qCells-lbCells)/(ubCells-lbCells) );
		
		assert (ret == ret);
		
		const double ret2 = UtilMath::makeBetween<double>(lbCount, ubCount, ret);
		
		assert (ret2 == ret2);
		
		return ret2;
	}
	
	/*
	inline bool getAligned(const QueryBox& b, long& lbMin, long& lbMax, long& ubMin, long& ubMax) const{
	
		const long emin = b.excludedMin;
		const long emax = b.excludedMax;
		
		const long amin = b.unboundedMin;
		const long amax = b.unboundedMax;
		
		const Point& qmin = b.min;
		const Point& qmax = b.max;
	
		const int DIMS = ind.getDims();
		
		lbMin = -1;
		lbMax = -1;
		ubMin = 0;
		ubMax = 0;
		
		long anchoredDimsMin = amin;
		long anchoredDimsMax = amax;
		
		for (int i = 0; i < DIMS; i++){
			
			const long LEN = ind.getLength(i);
			
			const long mini = 0;
			const long maxi = LEN-1;
			
			assert (maxi >= mini);
			
			const double bmini = UtilHist::bucketMin( mini, LEN );
			const double bmaxi = UtilHist::bucketMin( maxi+1, LEN );
			
			assert (bmaxi > bmini);
			
			if (qmin[i] >= bmaxi ){
				
				// NO INTERSECTION
				return false;
			
			} else if ( qmin[i] < bmini ){
			
				ubMin = lc.setDimCoord(ubMin, i, mini);
				anchoredDimsMin = UtilLongSet::add(anchoredDimsMin, i);
				
			} else {
			
				ubMin = lc.setDimCoord(ubMin, i, UtilHist::discretize(qmin[i], LEN ) );
			}
			
			if (qmax[i] < bmini ){
			
				// NO INTERSECTION
				return 0;
			
			} else if ( qmax[i] >= bmaxi ){
			
				ubMax = lc.setDimCoord(ubMax, i, maxi);
				anchoredDimsMax = UtilLongSet::add(anchoredDimsMax, i);
				
			} else {
			
				ubMax = lc.setDimCoord(ubMax, i, UtilHist::discretize(qmax[i], LEN ) );
			}
		}
		
		if (lc.getArea(ubMin, ubMax) <= 0)
			return false;
		
		lbMin = lc.incCoords( ubMin, anchoredDimsMin);
		lbMax = lc.decCoords( ubMax, anchoredDimsMax);
		
		
		return true;
	}*/
	
	inline bool getAligned(const Point& qmin, const Point& qmax, long& lbMin, long& lbMax, long& ubMin, long& ubMax, long& anchoredDimsMin, long& anchoredDimsMax, long min = -1, long max = -1, long amin = 0, long amax = 0, long emin = 0, long emax = 0) const{
	
	
		const int DIMS = ind.getDims();
		
		ubMin = 0;
		ubMax = 0;
		
		anchoredDimsMin = amin;
		anchoredDimsMax = amax;
		
		for (int i = 0; i < DIMS; i++){
			
			const long LEN = ind.getLength(i);
			
			const long mini = (min == -1) ? 0 : lc.getDimCoord(min, i);
			const long maxi = (max == -1) ? LEN-1 : lc.getDimCoord(max, i);
			
			assert (maxi >= mini);
			
			const double bmini = UtilHist::bucketMin( mini, LEN );
			const double bmaxi = UtilHist::bucketMin( maxi+1, LEN );
			
			assert (bmaxi > bmini);
			
			if (qmin[i] >= bmaxi ){
				
				// NO INTERSECTION
				return false;
			
			} else if ( qmin[i] < bmini ){
			
				ubMin = lc.setDimCoord(ubMin, i, mini);
				anchoredDimsMin = UtilLongSet::add(anchoredDimsMin, i);
				
			} else {
			
				ubMin = lc.setDimCoord(ubMin, i, UtilHist::discretize(qmin[i], LEN ) );
			}
			
			if (qmax[i] < bmini ){
			
				// NO INTERSECTION
				return 0;
			
			} else if ( qmax[i] >= bmaxi ){
			
				ubMax = lc.setDimCoord(ubMax, i, maxi);
				anchoredDimsMax = UtilLongSet::add(anchoredDimsMax, i);
				
			} else {
			
				ubMax = lc.setDimCoord(ubMax, i, UtilHist::discretize(qmax[i], LEN ) );
			}
		}
		
		lbMin = lc.incCoords( ubMin, anchoredDimsMin);
		lbMax = lc.decCoords( ubMax, anchoredDimsMax);
		
		return true;
	}
	
	// Query: qmin, qmax double points, Mode: mode enum, Query limitation: min, max long coords
	inline E getCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST, long min = -1, long max = -1, long amin = 0, long amax = 0, long emin = 0, long emax = 0) const{
	
			
		//const long double cellVol = 1.0L/ind.getLength();
		const int DIMS = ind.getDims();
		
		for (int i = 0; i < DIMS; i++){
			assert (qmax[i] >= qmin[i]);
			
			if ( (min != -1)  && (max != -1) )
				assert ( lc.getDimCoord(min, i) <= lc.getDimCoord(max, i) );
			
		}
		
		if ( (min != -1) && (max != -1) && lc.getArea(min, max) <= 0)
			return 0;
		
		long ubMin = 0; 
		long ubMax = 0;
		
		long lbMin = 0;
		long lbMax = 0;
		
		long anchoredDimsMin = amin;
		long anchoredDimsMax = amax;
		
		if (!getAligned(qmin, qmax, lbMin, lbMax, ubMin, ubMax, anchoredDimsMin, anchoredDimsMax, min, max, amin, amax, emin, emax))
			return 0;
		
		const long lbCells = UtilMath::maxVal<long>(0, lc.getArea(lbMin, lbMax));
		const long ubCells = UtilMath::maxVal<long>(0, lc.getArea(ubMin, ubMax));
		
		assert (lbCells >= 0);
		assert (ubCells >= lbCells);
		
		assert (ubCells > 0);
		
		//long double qVol = 1.0;
		//long double ubVol = 1.0;
		
		long double qCells = 1.0; //qVol/cellVol;
		
		for (int i = 0; i < DIMS; i++){
		
			const long mini = (min == -1) ? 0 : lc.getDimCoord(min, i);
			const long maxi = (max == -1) ? ind.getLength(i)-1 : lc.getDimCoord(max, i);
			
			assert (maxi >= mini);
			
			const double bmini = UtilHist::bucketMin( mini, ind.getLength(i) );
			const double bmaxi = UtilHist::bucketMax( maxi, ind.getLength(i) );
			
			
			//const double umini = UtilHist::bucketMin( lc.getDimCoord(ubMin, i), ind.getLength(i) );
			//const double umaxi = UtilHist::bucketMax( lc.getDimCoord(ubMax, i), ind.getLength(i) );
			
			
			assert (bmaxi > bmini);
			
			
			const double a = UtilLongSet::contains(anchoredDimsMin, i) ? bmini : qmin[i];
			const double b = UtilLongSet::contains(anchoredDimsMax, i) ? bmaxi : qmax[i];
			
			
			//const double c = UtilLongSet::contains(dims1, i) ? qmin[i] : bmini;
			//const double d = UtilLongSet::contains(dims2, i) ? qmax[i] : bmaxi;
			
			
			if ( (b > UtilMath::minVal(qmax[i], bmaxi)) || (a < UtilMath::maxVal(qmin[i], bmini)) ){
			
				cout << "a " << a << " b " << b << " qmin " << qmin[i] << " qmax " << qmax[i] << " bmini " << bmini << " bmaxi " << bmaxi << endl;
			}
			
			assert (b <= UtilMath::minVal(qmax[i], bmaxi) );
			assert (a >= UtilMath::maxVal(qmin[i], bmini) );
			
			//assert (a >= umini);
			//assert (b <= umaxi );
			
			//const long double w = b-a;
			
			//assert (w == w);
			
			//if (w < 0)
			//	cout << "a " << a << " b " << b << " qmin " << qmin[i] << " *" << UtilLongSet::contains(dims1, i) << " qmax " << qmax[i] << " *" << UtilLongSet::contains(dims2, i) <<" bmini " << bmini << " bmaxi " << bmaxi << endl;
			
			//assert (w >= 0);
			//assert (w <= 1);
			
			qCells *= (b-a)*ind.getLength(i);
			
			if (qCells < 0)
				qCells = 0;
			
			//ubVol *= (umaxi)-(umini);
			
			//if (qVol < 0)
			//	qVol = 0;
		}
		
		
		
		if (qCells > ubCells)
			qCells = ubCells;
			
		if (qCells < lbCells)
			qCells = lbCells;
		
		//const double qVol = qmin.getEnclosedVolume(qmax);
		
		//assert (qVol == qVol);
		//assert (ubVol == ubVol);
		
		const E ubCount = sum(ubMin, ubMax);
		
		if (mode == QueryMode::UB)
			return ubCount;
		
		
		//const long double ubVol = ubCells*cellVol;
		
		
		//if (qVol > ubVol)
		//	return ubCount; //cout << "qVol " << qVol << " ubVol " << ubVol << endl;
		
		
		//assert (qVol >= 0);
		//assert (ubVol >= qVol);
		
		if (lbCells <= 0){
		
			if (mode == QueryMode::LB)
				return 0;
			
			return UtilMath::makeBetween<double>(0, ubCount, (qCells/ubCells)*ubCount);
		}
		
		/*
		const long double lbVol = lbCells*cellVol;
		
		assert(lbVol == lbVol);*/
		
		//if (lbVol >= ubVol)
		//	return ubCount;
		
		const E lbCount = sum(lbMin, lbMax);
		
		//if (qVol <= lbVol)
		//	return lbCount;
			
		
		assert (lbCount == lbCount);
		
		if (mode == QueryMode::LB)
			return lbCount;
		
		//assert (ubVol > lbVol);
		
		if (ubCells == lbCells)
			return lbCount;
		
		const double ret = lbCount+(ubCount-lbCount)*((qCells-lbCells)/(ubCells-lbCells) );
		
		assert (ret == ret);
		
		const double ret2 = UtilMath::makeBetween<double>(lbCount, ubCount, ret);
		
		assert (ret2 == ret2);
		
		return ret2;
	}
	
	/*
	inline E getCount(const QueryBox& b) const{
	
		const long emin = b.excludedMin;
		const long emax = b.excludedMax;
		
		const long amin = b.unboundedMin;
		const long amax = b.unboundedMax;
		
		const Point& qmin = b.min;
		const Point& qmax = b.max;
	
		return getCount(qmin, qmax, b.mode, -1, -1, amin, amax, emin, emax);	
	}*/
	
	
	
	inline long getIndex(const Point& p) const{
	
		return ind.getIndex(lc, getCoord(p) );
	}
	
	inline void count(const Point& p) {
	
		counts.at(getIndex(p))++;
	}
	
	inline void countSub(const Point& p) {
	
		counts.at(getIndex(p))--;
	}
	
	
	
	inline void setAllCounts(long n) {
	
		for (long j = 0; j < counts.size(); j++)
			counts.at(j) = n;
	}
	
	inline bool allCountsEqualTo(long n) {
	
		for (long j = 0; j < counts.size(); j++)
			if (counts.at(j) != n)
				return false;
				
		return true;
	}
	
	inline long allCountsEqualTo() {
	
		long ret = -1;
		
		for (long j = 0; j < counts.size(); j++){
		
			long m = counts.at(j);
		
			if (ret != -1 && m != ret)
				return -1;
				
			ret = m;
		}
				
		return ret;
	}
	
	inline E getCountAtIndex(long index) const{
	
		if (countBytes.size() > 0){
		
			long ret = 0;
		
			int shift = 0;	
		
			for (int i = 0; i < 8; i++){
			
				if (countBytes.at(i).size() > 0)
					ret += ((long) countBytes.at(i).at(index)) << shift;
				
				shift += 8;
			}
			
			assert (ret >= 0);
			
			return ret;
		}
	
		//assert (counts.at(index) >= 0);
	
		return counts.at( index );
	}
	
	inline double distMax(const DenseGridHistogram<E>& g, double w1, double w2) const{
	
		assert (cumulative);
	
		assert (cellsPerDim == g.cellsPerDim);
		
		assert (ind.getDims() == g.ind.getDims());
		
		assert (ind.getLength() == g.ind.getLength());
	
		assert (w1 >= 0);
		assert (w2 >= 0);
	
		double ret = 0;
	
		for (long j = 0; j < ind.getLength(); j++){
			
			const double d = UtilMath::absVal<double>( getCountAtIndex(j)*w1- g.getCountAtIndex(j)*w2);
			
			if (d > ret)
				ret = d;
		}
		
		return ret;
	}
	
	inline void compress(){
		
		if (std::is_same<float, E>::value)
			return;
		
		if (std::is_same<double, E>::value)
			return;
		
		countBytes.clear();
	
		for (int i = 0; i < 8; i++){
		
			vector<uint8_t> v;
			countBytes.push_back(v);
		}
		
		const long byteSel = (1L << 8)-1;
		
		for (long j = 0; j < counts.size(); j++){
		
			long count = counts.at(j);
			
			for (int i = 0; i < 8; i++){
			
				if (count == 0)
					continue;
				
				const uint8_t b = count & byteSel;
				
				if (b != 0){
				
					if (countBytes.at(i).size() == 0){
					
						countBytes.reserve(counts.size() );
					
						for (long jj = 0; jj < counts.size(); jj++)
							countBytes.at(i).push_back(0);
					}
					
					countBytes.at(i).at(j) = b;
				}
				
				count >>= 8;
			}
		}
		
		counts.clear();
		counts.shrink_to_fit();
	}
	
	inline long getSize() const{
	
		return ind.getLength();
	}
	
	inline double getBytes() const{
	
		if (countBytes.size() > 0){
		
			double ret = 0;
		
			for (int i = 0; i < 8; i++)
				ret += countBytes.at(i).size();
		
			return ret;
		}
		
		return counts.size()*8;
	}
	
	inline long getCellsPerDim() const{
	
		return cellsPerDim;
	}
	
	inline void allocate(){
	
		countBytes.clear();
		counts.clear();
		counts.reserve(ind.getLength() );
		
	
		for (long j = 0; j < ind.getLength(); j++)
			counts.push_back(0);
			
		//cout << "counts.size() " << counts.size() << endl;
			
		cumulative = false;
	}
	
	inline long getIndex(const vector<long>& coords) const{
	
		return ind.getIndex(coords);
	}
	
	inline long getCoord(const vector<long>& coords) const{
	
		return lc.getIndex(coords);
	}
	
	inline E getCount(long l1) const{
	
		if (l1 == -1)
			return 0;
	
		const long index = ind.getIndex(lc, l1);
	
		return getCountAtIndex(index);
	}
	
	
	inline long getFirstNonEmpty(long s, long t) const{
	
	
	
		for (long l = s; l <= t; l = lc.getNext(l, s, t) ){
		
			if ( sum(l,l) > 0){
			
				assert (lc.contains(s, t, l));
			
				return l;
			}
		}
		
		assert (sum(s, t) > 0);
		
		assert (false);
	
	}
	
	
	inline long getRandomNonEmpty(RandomGenerator& gen, long s, long t, vector<long>& v, long dstSize) const{
	
		//const long srcSize = sum(s,t);
		
		v.clear();
		
		for (long k = 0; k < dstSize; k++)
			v.push_back(-1);
	
		ReservoirSampling rs(to_string(gen.randVal() ), dstSize);
		
		rs.reset();
		
		for (long l = s; l <= t; l = lc.getNext(l, s, t) ){
		
			const long ss = sum(l,l);
			
			for (long jj = 0; jj < ss; jj++){
			
				const long j = rs.pos();
			
				if ( (j >= 0) && (j < dstSize) )
					v.at(j) = l;
			}
		}
		
		std::sort(v.begin(), v.end() );
		
		assert (v.at(0) != -1);
	
	}
	
	inline long getRandom(RandomGenerator& gen, long s, long t) const{
	
	
		ReservoirSampling rs(to_string(gen.randVal() ), 1);
		
		rs.reset();
		
		long ret = s;
		
		for (long l = s; l <= t; l = lc.getNext(l, s, t) ){
		
			const long j = rs.pos();
			
			if (j == 0)
				ret = l;
			
		}
	
		return ret;
	}
	
	inline long getRandomNonEmpty(RandomGenerator& gen, long s_, long t_) const{
	
		if (false){
		
			const int DIMS = ind.getDims();
			
			long s = s_;
			long t = t_;
			
			for (int i = 0; i < DIMS; i++){
		
				cout << "getRandomNonEmpty " << i << endl;
		
				const E ss = sum(s,t);
				
				assert (ss > 0);
				
				const double dd = gen.randomDouble()*ss;
				
				assert (dd >= 0);
				assert (dd <= ss);
				
				const long a1 = lc.getDimCoord(s, i);
				
				long a = a1;
				long b = lc.getDimCoord(t, i);
				
				assert (b >= a);
			
				while ( b > a){
				
					long mid = (a+b) >> 1;
					
					cout << "a " << a << " mid " << mid << " b " << b << endl;
					
					assert (mid >= a);
					assert (mid <= b);
					
					//cout << "dd " << dd << " a " << a << " mid " << mid << " b " << b << endl;
					
					//cout << "sum( s, lc.setDimCoord(t, i, mid-1) " << sum( s, lc.setDimCoord(t, i, mid-1) ) << endl;
					//cout << "sum( s, lc.setDimCoord(t, i, mid) ) " << sum( s, lc.setDimCoord(t, i, mid) ) << endl;
					
					
					cout << "s " << lc.toString(s) << endl;
					cout << "t " << lc.toString(t) << endl;
					cout << "t' " << lc.toString(lc.setDimCoord(t, i, mid)) << endl;
					cout << "t''" << lc.toString(lc.setDimCoord(t, i, mid-1)) << endl;
					
					
					
					const E bb = sum( s, lc.setDimCoord(t, i, mid) );
					const E aa = mid == a1 ? 0 : sum( s, lc.setDimCoord(t, i, mid-1) ); 
					
					
					cout << "aa " << aa << " bb " << bb << endl;
					
					if ( aa < dd ){
				
						a = mid+1;
						continue;
					}
					
					if ( bb > dd ){
						
						b = mid-1;
						continue;
					}
					
					
					bool br = false;
					
					for (long j = a; j <= b; j++){
					
						const E bb = sum( s, lc.setDimCoord(t, i, j) );
						const E aa = j == a1 ? 0 : sum( s, lc.setDimCoord(t, i, j-1) ); 
						
						if ( (aa >= dd) && (bb <= dd) && (bb > aa)){
						
							a = j;
							b = j;
							br = true;
							break;
						}
					}
				
					assert (br);
					
					break;
				}
				
				assert (a == b);
				
				s = lc.setDimCoord(s, i, a);
				t = lc.setDimCoord(t, i, a);
				
				
				assert (sum(s,t) != 0);
				
				cout << "s " << lc.toString(s) << " t " << lc.toString(t) << endl;
				cout << "/getRandomNonEmpty " << i << endl;
			}
			
			assert (lc.contains(s_, t_, s) );
			assert (lc.contains(s_, t_, t) );
				
			cout << "s " << lc.toString(s) << " t " << lc.toString(t) << endl;
			
			assert ( s == t);
			
			assert (sum(s,t) != 0);
			
			assert (sum(s,t) > 0);
			
			
			return s;
		}
		
		E ss = 0;
		
		
		long s = s_;
		long t = t_;
		
		const double dd = gen.randomDouble()*sum(s,t);
		
		for (long l = s; l <= t; l = lc.getNext(l, s, t) ){
		
			ss += sum(l,l);
			
			
		
			if (ss > dd){
			
				assert (lc.contains(s, t, l));
			
				return l;
			}
		}
		
		assert (sum(s, t) > 0);
		
		assert (false);
	
	}
	
	inline long getMax(long s, long t) const{
		
		long max = s;
		
		E maxVal = sum(s,s);
		
		for (long l = s; l <= t; l = lc.getNext(l, s, t) ){
		
			const E sum1 = sum(l,l);
		
			if (sum1 > maxVal){
				max = l;
				maxVal = sum1;
			}
		}
		
		return max;
	}
	
	inline void addCount(const long l1, const E v){
		
		assert (!cumulative);
		
		counts.at( ind.getIndex(lc, l1) ) += v;
	}
	
	inline void mulCount(const long l1, const E v){
		
		assert (!cumulative);
		
		counts.at( ind.getIndex(lc, l1) ) *= v;
	}
	
	inline void subtractCount(const long l1, const E v){
		
		assert (!cumulative);
		
		counts.at( ind.getIndex(lc, l1) ) -= v;
	}
	
	inline void setCount(const long l1, const E v){
		
		assert (!cumulative);
		
		counts.at(ind.getIndex(lc, l1)) = v;
	}
	
	inline long nextCoord(const long c, int skipDims=0){
	
		const int dims = ind.getDims();
		
		long ret = c;
		
		for (int ii = 0; ii < dims; ii++){
					
			if (skipDims & (1 << ii) )
				continue;
						
			const long nextVal = lc.getDimCoord(c, ii)+1;
			
			if (nextVal == ind.getLength(ii) ){
			
				ret = lc.setDimCoord(ret, ii, 0);
				continue;
			}
			
			return lc.setDimCoord(ret, ii, nextVal);
		}
		
		return -1;
	}
	
	inline long prevCoord(const long c, int skipDims=0){
	
		const int dims = ind.getDims();
		
		long ret = c;
		
		for (int ii = (dims-1); ii >= 0; ii--){
					
			if (skipDims & (1 << ii) )
				continue;
						
			const long prevVal = lc.getDimCoord(c, ii)-1;
			
			if (prevVal < 0 ){
			
				ret = lc.setDimCoord(ret, ii, ind.getLength(ii)-1);
				continue;
			}
			
			return lc.setDimCoord(ret, ii, prevVal);
		}
		
		return -1;
	}
	
	// make differentially private
	inline void addLaplacianNoise(double lambda){
	
		assert (!cumulative);
	
		RandomGenerator rndGen("");
	
  	 	for (long j = 0; j < ind.getLength(); j++)
  	 		counts.at(j) += rndGen.laplace(lambda);
	}
	
	
	inline void setCumulative(bool b){
	
		if (cumulative == b)
			return;
		
		const int dims = ind.getDims();
		
		if (b){
		
			for (int i = 0; i < dims; i++){
		
				const long m = ind.getLength(i)-1;
				const int dimSel = 1 << i;
			
				for (long j = 1; j <= m; j++){
				
					for (long c = lc.setDimCoord(0, i, j); c != -1; c = nextCoord(c, dimSel))
						counts.at( ind.getIndex(lc, c) ) += getCount(lc.setDimCoord(c, i, j-1));
				}
			}
			
		} else {
		
			for (int i = (dims-1); i >= 0; i--){
		
				const long m = ind.getLength(i)-1;
				const int dimSel = 1 << i;
			
				for (long j = m; j >= 1; j--){
				
					for (long c = lc.setDimCoord(0, i, j); c != -1; c = nextCoord(c, dimSel))
						counts.at( ind.getIndex(lc, c) ) -= getCount(lc.setDimCoord(c, i, j-1));
				}
			}
		}
		
		cumulative = b;	
	}
	
	inline void finalise(){
	
		setCumulative(true);
		
		/*
		for (int i = 0; i < ind.getDims(); i++){
		
			const long m = ind.getLength(i)-1;
		
			for (long j = 0; j < ind.getLength(); j++){
			
				const long l1 = ind.getCoords(lc, j);
				
				//assert (j == ind.getIndex(lc, c) );
				
				const long n = lc.getDimCoord(l1, i);
				
				if (n < m){
					
					const long l2 = lc.setDimCoord(l1, i, n+1);
					addCount(l2, getCount(l1) );
				}
			}
		}*/
		
		compress();
		
		
		
	}
	
	inline void fillDebug(){
	
		for (long j = 0; j < ind.getLength(); j++){
		
			if (j % 2 == 0)
				counts.at(j) = j;
			else
				counts.at(j) = -j;
		} 
	}
	
	inline void print(){
		
		for (long y = 0; y < ind.getLength(1); y++){
		
			for (long x = 0; x < ind.getLength(0); x++){
			
				cout << counts.at(ind.getIndex(x,y) ) << "\t";
			}
			
			cout << endl;
		}
		
	}



	// Sum between (sx, sy) and (tx,ty)
	
	// A_1 = getCount(tx, sy-1)
	// A_2 = getCount(sx-1, ty)
	
	// A_1 cap A_2 = getCount(sx-1, sy-1)
	
	// Union( 
	

	// s and t in long coords
	inline E sum(const long s, const long t) const{
	
	
		if (!cumulative){
		
			assert (s <= t);
	
		
			E ret = 0;
		
			for (long l = s; l <= t; l = lc.getNext(l, s,t) )
				ret += getCount(l);
			
			return ret;
		}
	
		if (s <= 0)
			return getCount(t);
		
		if (t <= 0)
			return 0;
		
		
		const int dimMask = UtilBits::mask(ind.getDims() );
		
		// $dimMask$ is the set of all dimensions
		
		// go through all non-empty subsets $f$ of all dimensions
		
		
		long anti = 0;
		
		for (int f = 1; f <= dimMask; f++){
			
			long x = t;
			
			for (int i = 0; i < ind.getDims(); i++){
				
				if (f & (1 << i))
					x = lc.setDimCoord(x, i, lc.getDimCoord(s, i)-1 );
			}
			
			if (x == -1)
				continue;
			
			const int n = __builtin_popcount(f);
			
			const bool neg = (n-1) & 1;
			
			if (neg)
				anti -= getCount(x);
			else
				anti += getCount(x);
		}
		
		//assert (anti >= 0);
		
		//assert (anti <= getCount(t) );
		long ret = getCount(t)-anti;
		
		if (ret < 0)
			return 0;
		
		//assert (ret >= 0);
		
		return ret;
	}
	
	inline E sum() const{
	
		return sum(getMinCoord(), getMaxCoord() );
	}
};

typedef DenseGridHistogram<long> DenseGrid ;

class MultipleChoiceKnapsack{

private:
	
	
	const long knapsacks;
	const long items;
	
	MultiArray<double> profits;
	MultiArray<double> weights;

public:

	inline MultipleChoiceKnapsack(long knapsacks, long items) : knapsacks(knapsacks), items(items), profits(knapsacks, items),  weights(knapsacks, items){
	
		profits.fill( std::numeric_limits<double>::lowest() );
		weights.fill( std::numeric_limits<double>::max() );
	}
	
	
	
	inline void add(int knapsack, int item, double profit, double weight){
	
		profits.at(knapsack, item) = profit;
		weights.at(knapsack, item) = weight;
	}
	
	inline vector<int> maximizeProfit(double weightConstraint){
	
		if (knapsacks <= 5 && items <= 100)
			return maximizeProfitBruteForce(weightConstraint);
	
		assert (false);
	}
	
	inline double getProfit(vector<int> v){
	
		double ret = 0;
	
		for (int i = 0; i < v.size(); i++)
			ret += profits.at(i, v[i]);
			
		return ret;
	}
	
	inline double getWeight(vector<int> v){
	
		double ret = 0;
	
		for (int i = 0; i < v.size(); i++)
			ret += weights.at(i, v[i]);
			
		return ret;
	}
	
	inline vector<int> maximizeProfitBruteForce(double weightConstraint){
	
		vector<int> ret(knapsacks, -1);
		
		double maxProfit = std::numeric_limits<double>::lowest();
		
		for (int j1 = 0; j1 < items; j1++){
	
			const double p1 = profits.at(0, j1);
			const double w1 = weights.at(0,j1);
		
			if (w1 > weightConstraint)
				continue;

			if (knapsacks == 1){
					
				if (p1 > maxProfit){
			
					ret.at(0) = j1;
					maxProfit = p1;
				}
				continue;
			}

			for (int j2 = 0; j2 < items; j2++){
		
				const double p2 = p1+profits.at(1, j2);
				const double w2 = w1+weights.at(1, j2);
			
				if (w2 > weightConstraint)
					continue;
				
				if (knapsacks == 2){
					
					if (p2 > maxProfit){
				
						ret.at(0) = j1;
						ret.at(1) = j2;
						maxProfit = p2;
					}
					continue;
				}
			
				for (int j3 = 0; j3 < items; j3++){
	
					const double p3 = p2+profits.at(2, j3);
					const double w3 = w2+weights.at(2, j3);
				
					if (w3 > weightConstraint)
						continue;
					
					if (knapsacks == 3){
					
						if (p3 > maxProfit){
					
							ret.at(0) = j1;
							ret.at(1) = j2;
							ret.at(2) = j3;
							maxProfit = p3;
						}
						continue;
					}
				
					for (int j4 = 0; j4 < items; j4++){
				
						const double p4 = p3+profits.at(3, j4);
						const double w4 = w3+weights.at(3, j4);
					
						if (w4 > weightConstraint)
							continue;
					
						if (knapsacks == 4){
							
							if (p4 > maxProfit){
					
								ret.at(0) = j1;
								ret.at(1) = j2;
								ret.at(2) = j3;
								ret.at(3) = j4;
								maxProfit = p4;
							}
							continue;
						}
						
						for (int j5 = 0; j5 < items; j5++){
				
							const double p5 = p4+profits.at(4, j5);
							const double w5 = w4+weights.at(4, j5);
					
							if (w5 > weightConstraint)
								continue;
							
							if (knapsacks == 5){
							
								if (p5 > maxProfit){
					
									ret.at(0) = j1;
									ret.at(1) = j2;
									ret.at(2) = j3;
									ret.at(3) = j4;
									ret.at(4) = j5;
									maxProfit = p5;
								}
								continue;
							}
						}//5
					}//4
				}//3	
			}//2
		}//1
		
		return ret;
	}
};


class BipartiteMatcher{ 
	
	

	const int INF_DISTANCE = 99999999;
	const int NIL = 0;
	
	
	
	
	vector<int> v;
	
	//int queuepos = -1;
	
	//vector<bool> Q;
	
	std::random_device r;
		
	CircularBuffer<int> Q;
	vector<bool> Qbool;
	
	std::seed_seq seed;
	
	std::mt19937 eng;
	
	MultiArray<bool> g;
	
	int m, n; 
	
	vector<int> pairU;
	vector<int> pairV;
	vector<int> dist;
	
	
public: 
	inline BipartiteMatcher(int mm, int nn) : Q(mm+nn+1),  Qbool(mm+nn+1), seed{r(), r(), r(), r(), r(), r(), r(), r()}, eng(seed), g(mm+1, nn+1), m(mm), n(nn), pairU(mm+1), pairV(nn+1), dist(mm+1){
	
		
	}
	
	inline bool bfs(){
	
		
		
		for (int u = 1; u <= m; u++){
		
			if (pairU[u] == NIL){
				
				dist[u] = 0; 
				
				if (!Qbool[u])
					Q.push(u); 
				
				/*if (!Q[u]){
					Q[u] = true;
					queuesize++;
				}
				
				if (u < front)
					front = u;*/
				
			} else {
			
				dist[u] = INF_DISTANCE; 
			}
		}
		
		dist[NIL] = INF_DISTANCE;
		
		//long count = 0;
		
		while (!Q.empty()){
			
			const int u = Q.front(); 
			Qbool[u] = false;
			Q.pop(); 
			
			if (dist[u] < dist[NIL]){
				 
				for (int v = 1; v <= n ; v++){ 
				
					if (!g.get(u,v))
						continue;

					if (dist[pairV[v]] == INF_DISTANCE){
						dist[pairV[v]] = dist[u]+1; 
						
						const int q = pairV[v];
						
						if (!Qbool[q])
							Q.push(q);
					} 
				}
			}
		}
		
		return dist[NIL] != INF_DISTANCE;
	}
	
	bool dfs(int u){
	
		if (u == NIL)
			return true;
	
		for (int v = 1; v <= n ; v++){ 
			
			if (!g.get(u,v))
				continue;
				
			if (dist[pairV[v]] == (dist[u]+1) ){ 
				
				if (dfs(pairV[v]) == true){
					 
					pairV[v] = u; 
					pairU[u] = v; 
					return true;
				} 
			}
		}
		
		dist[u] = INF_DISTANCE; 
		return false;
	}
	
	inline void getMatching(const MultiArray<bool>& in, MultiArray<bool>& out, bool random=true){
	
		vector<long> px;
		vector<long> py;
		
		for (long i = 0; i <= m; i++){
			px.push_back(i);
		}
		
		for (long i = 0; i <= n; i++){
			py.push_back(i);
		}
		
		if (random){
			std::shuffle(px.begin()+1, px.end(),eng );
			std::shuffle(py.begin()+1, py.end(),eng );
		}
		
		g.fill(false);
		
		for (long u = 1; u <= m; u++){
			for (long v = 1; v <= n; v++){
			
				if (in.get(px.at(u)-1, py.at(v)-1) ){
				
					g.set(u,v, true);
				}
			}
		}
		
		for (int v = 0; v <= n; v++)
			pairV.at(v) = NIL;

		for (int u = 0; u <= m; u++){
			pairU.at(u) = NIL; 
			dist.at(u) = INF_DISTANCE;
		}
			
		out.fill(false);
		
		while ( bfs() == true )
			for (int u = 1; u <= m; u++)
				if (pairU[u] == NIL)
					dfs(u);
		
		for (int u = 1; u <= m; u++)
			out.set(px.at(u)-1, py.at(pairU[u])-1, true);
	}
	
	inline void getMatching(MultiArray<bool>& out, bool random=true){
	
		return getMatching(out, out, random);
	}
}; 


#endif