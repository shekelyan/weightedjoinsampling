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

#ifndef SHEKELYAN_UTILHASH_H
#define SHEKELYAN_UTILHASH_H

#include <approxdata/utils/utils.h>


	
namespace UtilHash{

	const long HASHMAX = std::numeric_limits<long>::max();
	
	std::hash<std::string> str_hash;
	
	//std::hash<const char*> char_hash;
	
	inline long hash(const long s, const long hashuniv){

		assert (hashuniv != 0);
		assert (hashuniv > 0);
		assert (hashuniv < HASHMAX);
		
		assert (s >= 0);
		
		//RandomNumbers rnd(to_string(s));
		
		//return rnd.randomInteger(0, hashuniv);
		
		return s % hashuniv;
	}
	
	inline long hash32_32(uint64_t key, long seed){
	
		uint32_t h;
		
		MurmurHash3_x86_32(&key, 4, seed, &h);
		
		return h;
	}
	
	inline long hash64_32(uint64_t key, long seed){
	
		uint32_t h;
		
		MurmurHash3_x86_32(&key, 8, seed, &h);
		
		return h;
	}
	
	
	
	
	inline long hash(const string& s, const long hashuniv){
		
		assert (hashuniv != 0);
		assert (hashuniv > 0);
		assert (hashuniv < HASHMAX);
		
		/*
		const long ret = str_hash(s) % hashuniv;
		RandomNumbers rnd(to_string(ret));
		return rnd.randomInteger(0, hashuniv);
		*/
		
		return str_hash(s) % hashuniv;
	}
	
	inline long_double to_longdouble(const char* s, int length, bool b = true){

		long_double val = UtilString::stringToLongDouble(s, length);
	
		const long_double NAN_RESERVED = std::numeric_limits<long_double>::max();
		
		if (val != NAN_RESERVED)	
			return val;
	
		val = UtilTime::daysSince1900(s, length);
		
		if (!b && val == -1)
			return NAN_RESERVED;
			
		if (val == -1)
			return str_hash(string(s, s+length) ) % (1L << 32);
			
		return val;
	}
	
	inline long_double to_longdouble(const string& s, bool b = true){
	
		long_double val = UtilString::stringToLongDouble(s);
	
		const long_double NAN_RESERVED = std::numeric_limits<long_double>::max();
		const char* c = s.c_str();
		
		if (UtilString::endOfNumber(c) != c+s.length() )
			val = NAN_RESERVED;
		
		if (val != NAN_RESERVED)	
			return val;
		
		val = UtilTime::daysSince1900(s);
		
		if (!b && val == -1)
			return NAN_RESERVED;
		
		if (val == -1)	
			return str_hash(s) % (1L << 32);;
			
		return val;
	}
}

#endif