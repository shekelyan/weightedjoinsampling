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

#ifndef SHEKELYAN_UTILBITS_H
#define SHEKELYAN_UTILBITS_H

#include <approxdata/utils/utils.h>


typedef unsigned char byte;


namespace UtilBits{

	inline int countBits(long x){
		
		return __builtin_popcountl(x);
	}
	
	
	
	
	inline long mask(int bits){
		
		return (1L << bits)-1;
	}
	
	inline long mask(int pos, int bits){
	
		return ((1L << bits)-1) << pos;
	}
	
	inline int getNumOfVariableBytes(long x){
	
		const byte sel = 0x7F;
		
		int bytes = 1;
		
		while (x > sel){
			x >>= 7;
			bytes++;
		}
		
		return bytes;
	}
	
	
	inline long setZeroInMask(long x, long mask){
	
		return (x | mask)^mask;
	}
	
	/*
	inline int getPosOfBitFlag(long x){
	
		int field_set(uint64_t input) {
  	  	uint64_t field = (x * 0x20406080a0c0e1ULL) >> 60;
    	return (field >> 60) & 15;
	}*/
	
	/*
	inline int getLoneBitPos(long x){
    	
    	int ret = 0;
		while (x >>= 1) ++ret;
		
		return ret;
	}*/
	
	inline long spreadBitsOverMask(long x, long mask){
		
		long ret = 0;
		long sel = 1;
	
		for (long m = mask; m != 0; m >>= 1, sel <<= 1){
	
			if (m & 1L){
			
				if (x & 1L)
					ret |= sel;
			
				x >>= 1;
			}
		}
		
		return ret;
	}
	
	inline long collectBitsFromMask(long x, long mask){
	
		long ret = 0;
		long sel = 1;
	
		for (long m = mask; m != 0; m >>= 1, x >>= 1){
			if (m & 1L){
			
				if (x & 1L)
					ret |= sel;
			
				sel <<= 1L;
			}
		}
		return ret;
	}
	
	/*
	inline long collectBits(long x, int pos, int bits){
	
		return collectBitsFromMask(x, pos, bits);
	}*/
	
	inline long incrementInMask(long x, long mask){
	
		long ret = x;
		
		for (long m = mask, sel = 1 ; m != 0; m >>= 1, sel <<= 1){
		
			if ( (mask & sel) != 0){
			
				if ( (x & sel) == 0){
				
					ret |= sel;
					return ret;
				}
				
				ret ^= sel;
			}
		}
		
		return -1;
	}
	
	inline long decrementInMask(long x, long mask){
	
		long ret = x;
		
		for (long m = mask, sel = 1 ; m != 0; m >>= 1, sel <<= 1){
		
			if ( (mask & sel) != 0){
			
				if ( (x & sel) != 0 ){
			
					ret ^= sel;
					return ret;
				}
				
				ret |= sel;
			}
		}
		
		return -1;
	}
};


namespace UtilLongSet{


	inline long emptySet(){
	
		return 0;
	}
	
	inline bool isEmpty(long l){
	
		return l == 0;
	}
	
	inline long getIndex( long l){
	
		for (int i = 0; i < 63; i++){
		
			if (l == (1L << i))
				return i;
		}
		
		return -1;
	}

	inline bool contains(long l, int i){
	
		return (l & (1L << i)) != 0;
	}
	
	inline long add(long l, int i){
	
		//assert (!contains(l, i));
		
		return l | (1L << i);
	}
	
	inline long remove(long l, int i){
	
		//assert (contains(l, i));
	
		return l ^ (1L << i);
	}
	
	inline long addSet(long l, long l2){
	
		return l | l2;
	}
	
	inline long removeSet(long l, long l2){
	
		return (l | l2) ^ l2;
	}
	
	
	inline long sequence(int i1, int i2){
	
		return UtilBits::mask(i1, (i2-i1)+1 );
	}
	
	inline long addSequence(long l, int i1, int i2){
	
		return l | sequence(i1, i2);
	}
	
	inline const string toString(long l){
	
		stringstream s;
		
		bool first = true;
		
		s << "{";
		
		for (int i = 0; i < 63; i++){
		
			if (contains(l, i)){
			
				if (!first)
					s << ", ";
					
				s << i;
				first = false;
			}
		}
		
		s << "}";
		
		return s.str();
	}
};

class BitSet{

private:

	vector<unsigned long> words;
	const long size;
	
public:
	inline BitSet(long size) : words(1+(size >> 6), 0), size(size){
	
		
	}
	
	inline void set(long k){
	
		if (k < 0 || k >= size)
			throw Exception("");
	
		words[k >> 6] |= 1L << (k-((k >> 6) << 6));
	}
	
	inline bool get(long k){
	
		if (k < 0 || k >= size)
			throw Exception("");
		
		return words[k >> 6] & ( 1L << (k-((k >> 6) << 6)) );
	}
};

class ByteReader{

public:
	inline virtual int readByte() = 0;
	virtual ~ByteReader(){};	
};


class ByteWriter{

public:
	inline virtual void writeByte(int b) = 0;
	virtual ~ByteWriter(){};
};






class IntegerInputStream{

public:
	inline virtual long read() = 0;
	virtual ~IntegerInputStream(){};

};

class IntegerOutputStream{

public:
	inline virtual void write(long x) = 0;
	virtual ~IntegerOutputStream(){};

};


class VariableByteVector : public IntegerOutputStream, public IntegerInputStream{

protected:

	int cursor = 0;
	
	vector<byte> bytes;
	
	const byte cont = 0x80;
	const byte sel = 0x7F;
	
public:
	
	~VariableByteVector() override {}
	
	inline void startRead(){
	
		cursor = 0;
	}
	
	inline long getSize(){
	
		return bytes.size();
	}
	
	inline long read() override{
	
		long ret = 0;
	
		for (int j = 0; j < 64; j+=7){
		
			const int b = bytes[cursor++];
			
			ret |= (b & sel) << j;
			
			if (! (b & cont) )
				return ret;
		}
		
		return ret;
	}
	
	inline void write(long v) override{
		
		if (v == 0){
			
			bytes.push_back(v);
			return;
		}
		
		while (v > 0){
		
			if (v > sel)
				bytes.push_back( (v & sel)|cont);
			else
				bytes.push_back( (v & sel) );
			
			v >>= 7;
		}
	}
};


#endif