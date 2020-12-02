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

#ifndef SHEKELYAN_UTILSTRING_H
#define SHEKELYAN_UTILSTRING_H

#include <approxdata/utils/utils.h>
#include <regex>
#include <external/MurmurHash3.cpp>

#include <set>

namespace UtilString{


	inline const string replaced(const string& s, char pattern, char repl, char beg='\0', char end='\0'){

		stringstream ss;

		int open = 0;

		for (int i = 0; i < s.length(); i++){
	
			if (s[i] == beg)
				open++;
		
			if (s[i] == end)
				open--;
		
			if (open == 0 && s[i] == pattern){
			
				ss << repl;
			
			} else {
			
				ss << s[i];
			}
		}
	
		return ss.str();
	}
	
	
	
	inline const string toUpperCase(const string& s){

		string ret = s;

		for (auto & c: ret)
			c = std::toupper(c);
		
		return ret;
	}
	
	inline const string toLowerCase(const string& s){

		string ret = s;

		for (auto & c: ret)
			c = std::tolower(c);
		
		return ret;
	}
	
	
	// contract white space
	inline const string cleanStringList(const string& s, char sep=',', char newsep=','){
	
	
		stringstream ss;
		
		
		bool open = false;
		
		for (int j = 0; j < s.length(); ++j){
			
			if (s[j] > 0 && (s[j] < '!' || s[j] == sep)){
			
				if (j > 0)
					open = true;
				continue;
			} 
			
			if (open){
				
				if (newsep != 0)
					ss << newsep;
				open = false;
			}	
			
			ss << s[j];
		}
		
		return ss.str();
		
	}

	inline const string replaced(const string& s, const string& pattern, const string& repl, char beg='\0', char end='\0'){

		stringstream ss;

		int open = 0;

		for (int i = 0; i < s.length(); i++){
	
			if (s[i] == beg)
				open++;
		
			if (s[i] == end)
				open--;
		
		
	
			if (open > 0 || i >= s.length()-pattern.length() ){
				ss << s[i];
			
			} else {
		
				int count = 0;
			
				for (int j = pattern.length()-1; j >= 0; j--){
				
					if ( (i+j) >= s.length())
						break;
					
					if (s[i+j] == pattern[j])
						count++;
				}
			
				if (count == pattern.length() ){
			
					ss << repl;
					i += pattern.length()-1;
				} else {
			
					ss << s[i];
				}
			}
		
		}
	
		return ss.str();
	}

	inline const string trim(const string& s, char c){
		
			
		int a = 0;
		int b = s.length()-1;
	
		if (b < 0)
			return "";
	
		while (s[a] == c){
			a++;
		
			if (a >= s.length() )
				return "";
		}
		
		while (s[b] == c){
			b--;
	
			if (b < 0)
				return "";	
		}
		
		return s.substr(a,b-a+1);
	}	

	inline long stringToNonNegativeInteger(const string& s, int off = 0, int end = -1){
		
		long ret = 0;
  		long b = 1;
  		
  		const int n = end >= 0 ? end : s.length();
  		
  		if (n == 0)
  			return 1L << 62;
  		
  		//if (s[0] == '0')
  		//	return n == 1 ? 0 : -1;
  		
  		for (int i = n-1; i >= off; --i){
    		
    		if (s[i] > '9')
    			return -1;
    		else if (s[i] < '0')
    			return -1;
    		else
    			ret += ((int)(s[i]-'0'))*b;
    			
    		b *= 10;
  		}
	
		return ret;
	}
	
	
	inline int getLength(const char* s){
	
		int ret = 0;
	
		for ( ; (*s) != 0; ++s)
			ret++;
			
		return ret;
	}
	
	// map integer strings to integers and one-char strings to 2^62-char
	inline static long stringToNonNegativeInteger(const char* s){
	
		const char c0 = (*s);
	
		if (c0 == 0)
  			return -1;
  		
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
  				return -1;
  			
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
	
	inline long stringToNonNegativeInteger(const char* s, int length, int off = 0, int end = -1){
		
		assert (getLength(s) == length);
		
		long ret = 0;
  		long b = 1;


  		const int n = end >= 0 ? end : length;
  		
  		if (n == 0)
  			return 1L << 62;
  		
  		if (s[0] == '0')
  			return n == 1 ? 0 : -1;
  		
  		for (int i = n-1; i >= off; --i){
    		
    		if (s[i] > '9')
    			return -1;
    		else if (s[i] < '0')
    			return -1;
    		else
    			ret += ((int)(s[i]-'0'))*b;
    			
    		b *= 10;
  		}
	
		return ret;
	}
	
	/*
	inline const char* endOfNumber(const char* b){
		
		const char* s = b;
		const bool minus = (*s) == '-';
		
		if (minus)
			++s;

		if ( (*s) < '0' || (*s) > '9')
			return b;
		
		const bool pre = (*s) == '0';
		++s;
		
		while ( (*s) >= '0' && (*s) <= '9'){
			
			if (pre)
				return b;

			++s;
		}
		
		if ( (*s) == 0)
			return s;
		
		if ((*s) == '.'){
		
			++s;
		
			if ( (*s) < '0' || (*s) > '9')
				return b;

			++s;
		
			while ( (*s) >= '0' && (*s) <= '9')
				++s;
			
			if ( (*s) == 0)
				return s;
		}
		
		if ( (*s) != 'e' && (*s) != 'E')
			return b;
		++s;
		
		const bool expminus = (*s) == '-';
		
		if (expminus ||(*s) == '+')
			++s;
		
		if ( (*s) < '0' || (*s) > '9')
			return b;
		
		const bool exp = (*s) == '0';
		++s;
		
		while ( (*s) >= '0' && (*s) <= '9'){
		
			if (exp)
				return b;
			
			++s;
		}
		
		if ( (*s) == 0)
			return s;
		
		return b;
	}*/
	
	inline long_double stringToLongDouble(const char* s, int length){
		
		const char* e = s+length;
		
		if (s == e)
			return std::numeric_limits<long_double>::max();
		
		const bool minus = (*s) == '-';
		
		if (minus){
			++s;
			
			if (s == e)
				return std::numeric_limits<long_double>::max();	
		}
		
		
		unsigned long pre = 0;
		
		while ( (*s) >= '0' && (*s) <= '9'){
		
			pre = pre*10 + ((*s)-'0');
			++s;
			
			if (s == e)
				return minus ? -pre : pre;
				
			
		}
		
		unsigned long post = 0;
		unsigned long div = 1;
		
		//
		if ((*s) == '.'){
		
			++s;
			
			if (s == e)
				return std::numeric_limits<long_double>::max();
			
			while ( (*s) >= '0' && (*s) <= '9'){
			
				div *= 10;
				//
				post = post*10 + ((*s)-'0');
				++s;
				
				if (s == e)
					return minus ? -(pre+((long double)post)/((long double) div) ) : (pre+((long double)post)/((long double) div) );
			}
		} else {
		
			if (pre == 0)
				return std::numeric_limits<long_double>::max();
		}
		
		if ( (*s) != 'e' && (*s) != 'E')
			return std::numeric_limits<long_double>::max();
		++s;
		//
		if (s == e)
			return std::numeric_limits<long_double>::max();
		//
		const bool expminus = (*s) == '-';
		
		if (expminus ||(*s) == '+'){
			++s;
			
			if (s == e)
				return std::numeric_limits<long_double>::max();
		}
		//
		
		if ( (*s) < '0' || (*s) > '9')
			return std::numeric_limits<long_double>::max();
		
		
		unsigned long exp = 0;
		
		while ( (*s) >= '0' && (*s) <= '9'){
			
			exp = exp*10 + ((*s)-'0');
			++s;
			
			if (s == e){
	
				const long signedexp = expminus ? -exp : exp;
				const long mantissa = minus ? -(pre+((long double)post)/div) : (pre+((long double)post)/div);		
	
				return mantissa * pow(10, signedexp);
			}	
			
			if (exp == 0)
				return std::numeric_limits<long_double>::max();
		}
		
		return std::numeric_limits<long_double>::max();
		
		/*
		
		int minus = -1;
		int dot = length;
		
		if (s[0] == '-')
			minus = 0;
		
		for (int i = 0; i < length; ++i){
		
			if (i > 0 && s[i] == '-'){
				
				return std::numeric_limits<long_double>::max();
			}
			
			if (s[i] == '.'){
				assert (dot == length);
				dot = i;
			}
		}
		
		
		long pre = dot == (minus+1) ? 0 : stringToNonNegativeInteger(s, length, minus+1, dot);
		
		if (pre == -1)
			return std::numeric_limits<long_double>::max();
		
		if (dot == length)
			return pre;
		
		long after = dot == length ? 0 : stringToNonNegativeInteger(s, length, dot+1, length);
		
		if (after == -1)
			return std::numeric_limits<long_double>::max();
		
		
		long div = 1;
		
		for (int i = length-dot; i > 1; --i)
			div *= 10;
			
			
		if ( minus >= 0)
			pre = -pre;
			
		//cout << pre << "." << after << " div " << div << endl;
		
		if ( minus >= 0)
			after = -after;
		
		return pre + ((long_double)after) / div;*/
	}
	
	
	
	inline long_double stringToLongDouble(const string& s){
		
		return stringToLongDouble(s.c_str(), s.length() );
		/*
		int minus = -1;
		int dot = s.length();
		
		if (s[0] == '-')
			minus = 0;
		
		for (int i = 0; i < s.length(); ++i){
		
			if (i > 0 && s[i] == '-'){
				
				return std::numeric_limits<long_double>::max();
			}
			
			if (s[i] == '.'){
				assert (dot == s.length());
				dot = i;
			}
		}
		
		
		long pre = dot == (minus+1) ? 0 : stringToNonNegativeInteger(s, minus+1, dot);
		
		if (pre == -1)
			return std::numeric_limits<long_double>::max();
		
		long after = dot == s.length() ? 0 : stringToNonNegativeInteger(s, dot+1, s.length() );
		
		if (after == -1)
			return std::numeric_limits<long_double>::max();
		
		if (after == 0)
			return pre;
		
		long div = 1;
		
		for (int i = s.length()-dot; i > 1; --i)
			div *= 10;
			
			
		if ( minus >= 0)
			pre = -pre;
			
		//cout << pre << "." << after << " div " << div << endl;
		
		if ( minus >= 0)
			after = -after;
		
		return pre + ((long_double)after) / div;*/
	}

	template<typename T>
	inline const string toString( const vector<T>& v){
	
		stringstream ss;
		
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				ss << ",";
		
			ss << v[i];
		}
		
		return ss.str();
	}
	
	
	
	
	inline const string toString( const vector<int>& v, const vector<string>& s){
	
		stringstream ss;
		
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				ss << ",";
		
			ss << s.at(v[i]);
		}
		
		return ss.str();
	}
	
	
	template<typename T>
	inline const string toString( constvector<T>& v){
	
		stringstream ss;
		
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				ss << ",";
		
			ss << v[i];
		}
		
		return ss.str();
	}

	inline const string twoDigit(long l){
	
		if (l >= 10)
			return to_string(l);
		else
			return "0"+to_string(l);
	}

	inline const string binaryToString(long l, int bits){
	
		stringstream s;
		
		for (long j = bits-1; j >= 0; j--){
		
			if (l & (1L << j)){
				
				s << "1";
				
				
			} else {
				s << "0";
			}
		}
		
		return s.str();
	}

	inline const string binaryToString(long l){
	
		stringstream s;
		
		bool zeros = false;
		
		for (long j = 60; j >= 0; j--){
		
			if (l & (1L << j)){
				
				s << "1";
				zeros = true;
				
			} else if (zeros){
				s << "0";
			}
		}
		
		if (!zeros)
			return "0";
		
		return s.str();
	}
	
	inline const string replace(const string& s, const string& regex, const string& t){
		
		std::regex r(regex);
			
		return std::regex_replace(s, r, t);
	}
	
	
	inline const string integerToShortString(long n){
		
		if ( (n % 1000000000) == 0)
			return to_string(n/1000000)+"g";
	
		if ( (n % 1000000) == 0)
			return to_string(n/1000000)+"m";
			
		if ( (n % 1000) == 0)
			return to_string(n/1000)+"k";
			
		return to_string(n);
	}
	
	inline const string bitString(long l){
	
		char* c = new char[64];
	
		for (int i = 0; i < 64; i++){
		
			c[63-i] = (l & (1L << i)) != 0 ? '1' : '0';
		}
		
		for (int i = 0; i < 64; i++){
		
			if (c[i] == '1')
				break;
				
			c[i] = '0';
		}
		
		string s(c);
		
		return s;
	
	}
	
	
	inline const char* endOf(const char* c, const string& s){
	
		assert (c != NULL);
		assert (s.length() > 0);
	
		const char* ret = c;
		
		for (int j = 0; j < s.length(); j++){
		
			if ( (*ret) != s[j])
				return c;
	
			++ret;
			
			if (j == s.length()-1)
				return ret;
		}
		
		assert (false);
	}
	
	
	inline const char* endOfNumber(const char* c){
	
		const char* ret = c;
		
		
		while ( UtilMath::isBetween(*ret, '0', '9') )
			++ret;
		
		if (ret == c)
			return c;
		
		if ( (*ret) == '.'){
		
			++ret;
			const char* beg = ret;
		
			while ( UtilMath::isBetween(*ret, '0', '9') )
				++ret;
			
			if (ret == beg)
				return c;
		}
		
		if (ret == c)
			return c;
			
		if ( UtilMath::isIn( *ret, 'e', 'E') ){
			++ret;
			
			if (UtilMath::isIn( *ret, '+', '-'))
				++ret;
			
			const char* beg = ret;
			
			while ( UtilMath::isBetween(*ret, '0', '9') )
				++ret;
				
			if (beg == ret)
				return c;
			
			return ret;
		}
		
		if (ret == c)
			return c;
			
		return ret;
	}
	
	inline const char* endOfRegion(const char* c, const string& begin, const string& end){
	
		assert (c != NULL);
		
		const char* ret = endOf(c, begin);
		
		if (ret == c)
			return c;
		
		for (const char* c2 = ret; (*c2) != 0; ++c2){
		
			ret = endOf(c2, end);
			
			if (ret != c2)
				return ret;
		}
		
		return c;
	}
	
	inline int getSplitLength(const string& stringbuf, const string& splitString){
	
		const char* beg = stringbuf.c_str();
		const char* end = beg+stringbuf.length();
		
		int counter = 0;
		
		for (const char* ptr = beg; ptr != end; ){
			
			const char* p = endOf(ptr, splitString);
			
			if (p != ptr){
			
				++counter;
				ptr = p;
			} else {
			
				++ptr;
			}
		}
		
		return counter+1;
		
	}
	
	
	inline const string getSplit(const string& stringbuf, const string& splitString, int n){
	
		if (n < 0){
		
			n += getSplitLength(stringbuf, splitString);
		}
		
		const char* beg = stringbuf.c_str();
		const char* end = beg+stringbuf.length();
		
		int counter = 0;
		
		const char* start = beg;
		
		for (const char* ptr = beg; ptr != end; ){
			
			const char* p = endOf(ptr, splitString);
			
			if (p != ptr){
			
				if (counter == n)
					return string(start, ptr);
				
				start = p;
				
				++counter;
				ptr = p;
				
			} else {
			
				++ptr;
			}
		}
		
		assert (counter == n);
		
		return string(start, end);
	}
	

	inline bool split(const string& str, const string& beginStr, string& a, string& b){
							
		std::size_t x = str.find_first_of(beginStr);
		
		if (x == std::string::npos)
			return false;
		
		a = str.substr(0,x);
		b = str.substr(x+1, str.size()-x-1);
		
		return true;
	}
	
	inline bool equals(const string& a, const string& b){
		
		return a.compare(b) == 0;
	}
	
	
	inline bool equalsAny(const string& s1, const string& s2){
	
		return equals(s1, s2);
	}
	
	inline bool equalsAny(const string& s1, const string& s2, const string& s3){
	
		return equalsAny(s1, s2) || equals(s1, s3);
	}
	
	inline bool equalsAny(const string& s1, const string& s2, const string& s3, const string& s4){
	
		return equalsAny(s1, s2, s3) || equals(s1, s4);
	}
	
	inline bool equalsAny(const string& s1, const string& s2, const string& s3, const string& s4, const string& s5){
	
		return equalsAny(s1, s2, s3, s4) || equals(s1, s5);
	}	
	
	
	inline const string beforeSplit(const string& str, const string& split){
							
		std::size_t x = str.find_first_of(split);
		
		if (x == std::string::npos)
			return "";
		
		return str.substr(0,x);
	}
	
	//inline bool startsWith(const string& str, const string& split){
	//						
	//	return str.find_first_of(split) == 0;
	//}
	
	inline const string afterSplit(const string& str, const string& split){
							
		std::size_t x = str.find_first_of(split);
		
		if (x == std::string::npos)
			return str;
		
		return str.substr(x+split.length());
	}
	
	
	inline const string trim(const string& s){
	
		string str = s;
		
		{
		size_t endpos = str.find_last_not_of(" \t");
		size_t startpos = str.find_first_not_of(" \t");
		if( std::string::npos != endpos )
		{
			str = str.substr( 0, endpos+1 );
			str = str.substr( startpos );
		}
		else {
			str.erase(std::remove(std::begin(str), std::end(str), ' '), std::end(str));
		}
		}
		// trim leading spaces
		{
		size_t startpos = str.find_first_not_of(" \t");
		if( string::npos != startpos )
		{
			str = str.substr( startpos );
		}
		}
		
		return str;
	
	}
	
	
	
	inline const string afterLastSplit(const string& str, const string& split){
							
		std::size_t x = str.find_last_of(split);
		
		if (x == std::string::npos)
			return str;
		
		return str.substr(x+split.length());
	}
	
	
	
	inline const string getUnixFolder(const string& path){
	
		int lastSlash = -1;
		
		for (int i = 0; i < path.length(); i++)
			if (path[i] == '/')
				lastSlash = i;
	
		//if (lastSlash == -1)
		//	return "";
			
		return path.substr(0, lastSlash+1);
	}
	
	inline const string getUnixFile(const string& path){
	
		int lastSlash = -1;
		
		for (int i = 0; i < path.length(); i++)
			if (path[i] == '/')
				lastSlash = i;
			
		return path.substr(lastSlash+1);
	}
	
	inline const string getUnixFileName(const string& path){
	
		int lastDot = -1;
		
		const string file = getUnixFile(path);
		
		for (int i = 0; i < file.length(); i++){	
			if (file[i] == '.')
				lastDot = i;
		}
		
		return file.substr(0, lastDot);
	}
	
	inline const string getFileEnding(const string& path){
	
		int lastDot = -1;
		
		const string file = getUnixFile(path);
		
		for (int i = 0; i < file.length(); i++)
			if (file[i] == '.')
				lastDot = i;
		
		if (lastDot == -1)
			return "";
		
		return file.substr(lastDot);
	}
	
	inline const string getFile(const string& path){
	
		const string f = "/"+path+"$";
		
		try {
		  std::regex re("[\\/]([^\\/]+)[$]");
		  std::smatch match;
		  if (std::regex_search(f, match, re) && match.size() > 1) {
		  
			return match.str(1);
		  } else {
			return "$";
		  } 
		} catch (std::regex_error& e) {
		  
		  	return "error";
		}
	}
	
	inline const string getFileName(const string& path){
	
		const string f = "/"+path+"$";
		
		try {
		  std::regex re("[\\/]([^\\/]+)[.][^.]+[$]");
		  std::smatch match;
		  if (std::regex_search(f, match, re) && match.size() > 1) {
		  
			return match.str(1);
		  } else {
		  	
			return "$";
		  } 
		} catch (std::regex_error& e) {
		  
		  	return "error";
		}
	}
	
	inline double strToDouble(const string& valstr){
	
		try{
		
			std::string::size_type sz;
			 
			const double d = std::stod(valstr, &sz);
			 
			 if (sz == valstr.size() )
				return d;
			
		} catch (...){
		
			
		}
		
		return std::numeric_limits<double>::infinity();
	}
	
	inline int strToInt(const string& valstr){
	
		try{
		
			std::string::size_type sz;
			 
			const int d = std::stoi(valstr, &sz);
			 
			 if (sz == valstr.size() )
				return d;
			
		} catch (...){
		
			
		}
		
		return std::numeric_limits<int>::min();
	}
	
	template <typename E>
	inline void getVal(const string& s, E& val){
	
		istringstream numStream(s);
			
		if (!(numStream >> val))
			assert(false);
	}
	
	
	template <typename E>
	inline bool splitGetVal(const string& str, const string& beginStr, E& a, E& b){
							
		std::size_t x = str.find_first_of(beginStr);
		
		if (x == std::string::npos)
			return false;
		
		string as = str.substr(0,x);
		string bs = str.substr(x+1, str.size()-x-1);
		
		getVal(as,a); 
		getVal(bs,b);
		
		return true;
	}
	
	
	inline bool split(const string& str,  const string& beginStr, const string endStr, string& a, string& b, string& c){
							
		std::size_t x = str.find_first_of(beginStr);
		
		if (x == std::string::npos)
			return false;
		
		std::size_t y = str.find_first_of(endStr);
		
		if (y == std::string::npos)
			return false;
		
		a = str.substr(0,x);
		b = str.substr(x+1, y-x-1);
		c = str.substr(y+1,str.size()-y-1);
		
		return true;
	}
	
	
	inline const string digitToStr(int i){
	
		if (i == 0)
			return "zero";
		
		if (i == 1)
			return "one";
	
		if (i == 2)
			return "two";
	
		if (i == 3)
			return "three";
			
		if (i == 4)
			return "four";
		
		if (i == 5)
			return "five";
	
		if (i == 6)
			return "six";
	
		if (i == 7)
			return "seven";
			
		if (i == 8)
			return "eight";
		
		if (i == 9)
			return "nine";
	
		if (i == 10)
			return "ten";
	
		if (i == 11)
			return "eleven";
			
		if (i == 12)
			return "twelve";
			
		if (i == 13)
			return "thirteen";
			
		if (i == 14)
			return "fourteen";
		
		if (i == 15)
			return "fifteen";
		
		if (i == 16)
			return "sixteen";
		
		if (i == 17)
			return "seventeen";
			
		if (i == 18)
			return "eighteen";
			
		if (i == 19)
			return "nineteen";
		
		if (i == 20)
			return "twenty";
		
		return "overtwenty";
	}
	
	
	inline bool endsWith(const string& value, const string& ending){
		
    	if (ending.size() > value.size())
    		return false;
    	
  		for (int k = 0; k < ending.size(); ++k)
  			if (value[value.length()-ending.length()+k] != ending[k])
  				return false;
      	
      	return true;
    	//return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
    }
    
    inline vector<string> getSplits(const string& s, char sep = ' ', bool ignoreWhiteSpaces = true){
    
    	char lim = ignoreWhiteSpaces ? ' ' : 0;
    
		vector<string> vec;
		shared_ptr<stringstream> ss;
	
		for (int j = 0; j < s.length(); ++j){
	
			if (s[j] != sep && s[j] > lim){
			
				if (ss){
			
				} else {
			
					ss = std::make_shared<stringstream>();
				}
			
				(*ss) << s[j];
			
			} else {
		
				if (ss){
					vec.push_back( ss->str() );
					ss.reset();
				}
			}
		}
	
		if (ss){
			vec.push_back( ss->str() );
			ss.reset();
		}
		
		return vec;
	}
    
    inline bool startsWith(const string& value, const string& ending){
		
    	if (ending.size() > value.size())
    		return false;
    	
  		for (int k = 0; k < ending.size(); ++k)
  			if (value[k] != ending[k])
  				return false;
      	
      	return true;
    	//return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
    }
	
	
	/*
	inline bool startsWith(const string& value, const string& prefix){
	
		if (prefix.length() >= value.length() )
			return false;
			
		for (int j = 0; j < prefix.length(); j++){
		
			if (value[j] != prefix[j])
				return false;
		}
		
		return true;
	
		//return value.substr(0, prefix.size()).compare(prefix) == 0;
    }*/
	
	
	template <typename E>
	inline string getDigitRepresentation(E n, int radix, int digitNum){
	
		stringstream s;
		
		for (int i = digitNum-1; i >= 0; --i){
		
			long coef = pow( radix, i);
		
			int digit = n /coef;
			
			s << digit;
			
			n -= digit*coef;
		}
		
		return s.str();	
	}
	
	template <typename T>
	inline const string iterableToString(const T& v){
	
		if (v.size() == 0)
			return "NULL";
	
		if (v.size() > 10){
		
			stringstream s;
	
			for (int i = 0; i < 5; i++){
				
				s << v[i];
				s << ", ";
			}
			
			s << "...";
			
			for (int i = 4; i >= 0; i--){
		
				s << ", ";
				s << v[v.size()-1-i];
			}
			
			return s.str();
		}
	
		stringstream s;
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				s << ", ";
				
			s << v[i];
		}
		
		return s.str();
	}
	
	template <typename T>
	inline const string listToString(const vector<T>& v){
	
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
	
	inline const string repeat(const string& s, int n){
	
		stringstream ss;
		
		for (int i = 0; i < n; i++){
		
			ss << s;
		}
		
		return ss.str();
	}
	
	template <typename T>
	inline const string horizontalLatexTable(const vector<T>& v){
	
		stringstream ss;
		
		ss << "\\begin{tabular}{|";
		
		for (int i = 0; i < v.size(); i++)
			ss << "c";
		
		ss << "|}\n";
		
		ss << "\\hline" << endl;
		
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				ss << " & ";
			
			ss << v.at(i);
		}
		ss << "\\\\" << endl;
		ss << "\\hline" << endl;
		
		ss << "\\end{tabular}" << endl;
	
		return ss.str();
	}
	
	/*
	inline vector<string>* newStringVector(string v1, string v2){
	
		vector<string>* ret = new vector<string>();
		
		ret->push_back(v1);
		ret->push_back(v2);
		
		return ret;
	}*/
	
	template <typename T>
	inline const string latexTable(const vector<vector<T>>& vs, const vector<string>& rownames, const string betweenRows=""){
	
		stringstream ss;
		
		ss << "\\begin{tabular}{|";
		
		ss << "l|";
		
		for (int i = 0; i < vs.at(0).size(); i++)
			ss << "c|";
		
		ss << "}\n";
		
		ss << "\\hline" << endl;
		
		for (int j = 0; j < vs.size(); j++){
		
			ss << rownames.at(j) << " & ";
			
			vector<T>& v = vs.at(j);
			
			for (int i = 0; i < v.size(); i++){
			
				if (i > 0)
					ss << " & ";
			
				ss << v.at(i);
			}
			ss << "\\\\" << endl;
			ss << "\\hline" << endl;
			ss << betweenRows << endl;
		}
		
		ss << "\\end{tabular}" << endl;
	
		return ss.str();
	}
	
	template <typename T>
	inline const string listToStringJoin(const vector<T>& v, const string join=","){
	
		stringstream s;
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				s << join;
			
			s << v.at(i);
		}
		
		return s.str();
	}
	
	template <typename T>
	inline const string arrayToStringJoin(int size, T* v, const string join=","){
	
		if (v == NULL)
			return "NULL";
		
		stringstream s;
	
		for (int i = 0; i < size; i++){
		
			if (i > 0)
				s << join;
			
			s << v[i];
		}
		
		return s.str();
	}
	
	/*
	*/
	
	inline bool contains(const string& s1, const string& s2){
	
		return s1.find(s2) != string::npos;
	}
	
	
	inline bool containsAll(const string& s1, const string& s2){
	
		return contains(s1, s2);
	}
	
	inline bool containsAll(const string& s1, const string& s2, const string& s3){
	
		return contains(s1, s2) && contains(s1, s3);
	}
	
	inline bool containsAll(const string& s1, const string& s2, const string& s3, const string& s4){
	
		return containsAll(s1, s2, s3) && contains(s1, s4);
	}
	
	inline bool containsAny(const string& s1, const string& s2){
	
		return contains(s1, s2);
	}
	
	inline bool containsAny(const string& s1, const string& s2, const string& s3){
	
		return contains(s1, s2) || contains(s1, s3);
	}
	
	inline bool containsAny(const string& s1, const string& s2, const string& s3, const string& s4){
	
		return containsAny(s1, s2, s3) || contains(s1, s4);
	}
	
	inline bool containsAny(const string& s1, const string& s2, const string& s3, const string& s4, const string& s5){
	
		return containsAny(s1, s2, s3, s4) || contains(s1, s4);
	}
	
	/*
	inline bool containsAny(const string& s1, const vector<const string>& s2){
	
		for (auto it = s2.begin(); it != s2.end(); it++){
			
			if (contains(s1, (*it)))
				return true;
		}
		
		return false;
	}*/
	
	/*
	inline bool containsAll(const string s1, const vector<const string>& s2){
	
		for (auto it = s2.begin(); it != s2.end(); ++it){
		
			if (!contains(s1, (*it)))
				return false;
		}
		
		return true;
	}*/
	
	template <typename T>
	inline void printList(const vector<T>& v){
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				cout << ", ";
				
			cout << v.at(i);
		}
		
		cout << endl;
	}
	
	
	template <typename T>
	inline void print(const vector<T>& v){
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				cout << ", ";
				
			cout << v.at(i);
		}
		
		cout << endl;
	}
	/*
	inline void printList(Point& v){
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				cout << ", ";
				
			cout << v.at(i);
		}
		
		cout << endl;
	}*/
	
	
	template <typename T>
	inline void printList(shared_ptr<vector<T>> v){
	
		printList( (*v) );
	}
	
	
	
	inline void printException(std::exception& exception){
	
		cout << "exception " << exception.what() << endl;
		cerr << "exception " << exception.what() << endl;
	}
	template <typename T>
	inline void printArray(int len, T* v){
	
		for (int i = 0; i < len; i++){
		
			if (i > 0)
				cout << ", ";
			
			cout << v[i];
		}
		
		cout << endl;
	}
	
	
	inline void splitFirst(const string& source, char split, string& first, string& rest){
	
		stringstream lineStream(source);
		
		string s;
			
		first = "";
		rest = "";
			
		if (!getline(lineStream, s, split))
			return;
		
		first = s;
		
		std::stringstream ss;	
		
		for (bool first = true; getline(lineStream, s, split); first = false){
		
			if(!first)
				ss << split;
			ss << s;
		}
		
		rest = ss.str();
	}
	
	inline const string removeAll(const string str, char remove) {
		
		stringstream s;
		
		for (int i = 0; i < str.length(); i++){
		
			if (str[i] != remove)
				s << str[i];
		}
		
		return s.str();
	}
	
	template<typename T>
	inline const string replaceAll(const string& str, const string& from, const T& to) {
		
		stringstream ss;
		
		for (int i = 0; i < str.length(); ++i){
		
			int matches = 0;
			
			for (int k = 0; k < from.length(); k++){
			
				if ((i+k) < str.length() && str[i+k] == from[k])
					matches++;
			}
			
			if (matches == from.size() ){
			
				ss << to;
				i += from.length()-1;
				
			} else {
			
				ss << str[i];
			}	
		}
		
		return ss.str();
	}
	
	inline std::set<string> replaceAllIn(const std::set<string>& set1, const string& targ, const string& repl){
	
		std::set<string> ret;
	
		for (auto it = set1.begin(); it != set1.end(); ++it)
			ret.insert( UtilString::replaced(*it, targ, repl) );
	
		return ret;
	}
	
	inline vector<string> replaceAllIn(const vector<string>& vec, const string& targ, const string& repl){
	
		vector<string> ret;
	
		for (auto it = vec.begin(); it != vec.end(); ++it)
			ret.push_back( UtilString::replaced(*it, targ, repl) );
	
		return ret;
	}
	
	template<typename T>
	inline std::map<string,T> replaceAllIn(const std::map<string,T>& map, const string& targ, const string& repl) {
	
		std::map<string,T> ret;
	
		for (auto it = map.begin(); it != map.end(); ++it)
			ret.insert( std::make_pair( UtilString::replaced(it->first, targ, repl), it->second) );
	
		return ret;
	}
	/*
	inline const string replaceAll(const string& str, const string& from, const string& to) {
		
		stringstream ss;
		
		for (int i = 0; i < str.length(); ++i){
		
			int matches = 0;
			
			for (int k = 0; k < from.length(); k++){
			
				if ((i+k) < str.length() && str[i+k] == from[k])
					matches++;
			}
			
			if (matches == from.size() ){
			
				ss << to;
				i += from.length()-1;
				
			} else {
			
				ss << str[i];
			}	
		}
		
		return ss.str();
		
		
		string s = str;
		
		size_t start_pos = 0;
		while((start_pos = str.find(from, start_pos)) != std::string::npos) {
			s.replace(start_pos, from.length(), to);
			start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
		}
		
		const string ret = s;
		
		return ret;
	}*/

	inline const string removeLast(const string s, int m){
	
		if (m >= s.length())
			return s;
	
		return s.substr(0, s.length()-m);
	}

	inline const string numToStr(double d){
	
		if (d < 0)
			return "-"+numToStr(-d);
	
		if (d != d)
			return " ";
	
		if (d == ((double)((int)d)))
			return to_string((int) d);
		
		return replaceAll(to_string(d), ".", ",");
	}
	
	
	
	inline const string numToStr(double value, int decimals){
		
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(decimals) << value;
		std::string s = ss.str();
		if(decimals > 0 && s[s.find_last_not_of('0')] == '.') {
			s.erase(s.size() - decimals + 1);
		}
		return s;
	}
	
	inline const string doubleToString(double value, int decimals=-1){
	
		if (decimals == -1){
		
			const string t = "Q";
		
			const string s = to_string(value)+t;
			
			std::regex r("([.]?[0]+)"+t);
			
			const string s2 = std::regex_replace(s, r, t);
			
			return s2.substr(0, s2.length()-t.length() );
		}
	
		
		if (value == floor(value))
			return to_string( ((long) value));
			
		std::ostringstream ss;
	
		if (decimals == -1){
		
			ss << std::fixed << std::setprecision(100) << value;
			
		} else{
			ss << numToStr(value, decimals);
			
		}
		
		return ss.str();
	}
	
	
	
	
	inline void print(const string str, const vector<int>& coords){
	
		cout << str;
		
		for (int i = 0; i < coords.size(); i++)
			cout << " " << coords.at(i);
			
		cout << endl << std::flush;
	}
	
	
	
	inline void printErr(const string str, const vector<int>& coords){
	
		std::cerr << str;
	
		for (int i = 0; i < coords.size(); i++)
			std::cerr << " " << coords.at(i);
		
		std::cerr << endl << std::flush;
	}
	
	inline void print(const string str, const vector<double>& coords){
	
		cout << str;
	
		for (int i = 0; i < coords.size(); i++)
			cout << " " << coords.at(i);
		
		cout << endl << std::flush;
	}
	
	inline void print(std::string str, double d){
	
		cout << str << " " << d << endl << std::flush;
	}
	
	inline void printErr(std::string str, int dims, vector<double>& coords){
	
		std::cerr << str;
	
		for (int i = 0; i < coords.size(); i++)
			std::cerr << " " << coords.at(i);
			
		std::cerr << endl << std::flush;
	}
	
	inline const string msToShortString(double ms){
	
		if (ms < 0.1)
			return "<0.1 ms";
	
		const double SEC = 1000;	
		const double MIN = 60*SEC;
		const double HRS = 60*MIN;
		const double DAY = 24*HRS;
	
		if ( (ms/DAY) >= 0.5)
			return UtilString::doubleToString(ms/DAY, 1)+" days";
		
		if ( (ms/HRS) >= 0.5)
			return UtilString::doubleToString(ms/HRS, 1)+" hs";
		
		if ( (ms/MIN) >= 0.5)
			return UtilString::doubleToString(ms/MIN, 1)+" mins";
		
		if ( (ms/SEC) >= 0.5)
			return UtilString::doubleToString(ms/SEC, 1)+" secs";
			
		return UtilString::doubleToString(ms/SEC, 1)+" ms";	
	}
	
	
	inline const string msToString(double ms){
	
		if (ms <= 0)
			return "0 ms";
	
		stringstream s;
	
		bool first = true;
		
		for (int i = 0; i <= 4; i++){
		
			long mul = 0;
			
			if (i == 0)
				mul = 1000*60*60*24;
			if (i == 1)
				mul = 1000*60*60;
			if (i == 2)
				mul = 1000*60;
			if (i == 3)
				mul = 1000;
			if (i == 4)
				mul = 0;
			
			double val = UtilMath::makeMultipleOf(ms, mul);
			
			if (i == 3 && (val > 0 || (!first) ) ){
			
				if (!first)
					s << " ";
				
				s << UtilString::doubleToString(ms/mul, 3);
				
				s << " ";
				
				if (i == 0)
					s << "d";
				else if (i == 1)
					s << "hs";
				else if (i == 2)
					s << "mins";
				else if (i == 3)
					s << "secs";
				else if (i == 4)
					s << "ms";
				
				break;
			}
			
			if (val > 0){
			
				ms -= val;
			
				if (!first)
					s << " ";	
			
				first = false;
			
				s << UtilString::doubleToString(mul == 0 ? val : val/mul, 3);
				s << " ";
			
				if (i == 0)
					s << "d";
				else if (i == 1)
					s << "hs";
				else if (i == 2)
					s << "mins";
				else if (i == 3)
					s << "secs";
				else if (i == 4)
					s << "ms";
			}
		}
		
		return s.str();	
	}
	
	
	inline void writeToFile(const string file, const string text){
	
		ofstream out(file);
		
		out << text;
		
		out.close();
	}
	
	inline const string bytesToString(double bytes){
	
	
		if (true){
		
			const double d = 0.1;
		
			if (bytes < d*pow(1024,1)*10 )
				return numToStr(bytes,2)+" B";
		
			if (bytes < d*pow(1024,2)*10)
				return numToStr(bytes/1024.0,2)+" KB";
		
			if (bytes < d*pow(1024,3)*10)
				return numToStr(bytes/1024.0/1024.0,2)+" MB";
		
			if (bytes < d*pow(1024,4)*10)
				return numToStr(bytes/1024.0/1024.0/1024.0,2)+" GB";
			
			if (bytes < d*pow(1024,5)*10)
				return numToStr(bytes/1024.0/1024.0/1024.0/1024.0,2)+" TB";
			
			if (bytes < d*pow(1024,6)*10)
				return numToStr(bytes/1024.0/1024.0/1024.0/1024.0/1024.0,2)+" PB";
			
			
			return numToStr(bytes/1024.0/1024.0/1024.0/1024/1024/1024,2)+" EB";
		
		}
	
		if (bytes <= 0)
			return "0 bytes";
	
		stringstream s;
	
		bool first = true;
		
		for (int i = 0; i <= 4; i++){
		
			long mul = pow(1024, 4-i);
			
			double val = UtilMath::makeMultipleOf(bytes, mul);
			
			if (i == 3 && (val > 0 || (!first) ) ){
			
				if (!first)
					s << " ";
				
				s << UtilString::doubleToString(bytes/mul, 3);
				
				s << " ";
				
				if (i == 0)
					s << "TB";
				else if (i == 1)
					s << "GB";
				else if (i == 2)
					s << "MB";
				else if (i == 3)
					s << "KB";
				else if (i == 4)
					s << "B";
				
				break;
			}
			
			if (val > 0){
			
				bytes -= val;
			
				if (!first)
					s << " ";	
			
				first = false;
			
				s << UtilString::doubleToString(mul == 0 ? val : val/mul, 3);
				s << " ";
			
				if (i == 0)
					s << "TB";
				else if (i == 1)
					s << "GB";
				else if (i == 2)
					s << "MB";
				else if (i == 3)
					s << "KB";
				else if (i == 4)
					s << "B";
			}
		}
		
		return s.str();	
	}
	
	
	
	template<typename T1>
	inline const string valuesToString(long vals){
	
		return bytesToString( vals*sizeof(T1) );
	}
	
	template<typename T1, typename T2>
	inline const string valuesToString(long vals1, long vals2){
	
		return bytesToString( vals1*sizeof(T1) + vals2*sizeof(T2) );
	}
	
	template<typename T1, typename T2, typename T3>
	inline const string valuesToString(long vals1, long vals2, long vals3){
	
		return bytesToString( vals1*sizeof(T1) + vals2*sizeof(T2)+vals3*sizeof(T3) );
	}
	
	template<typename T1, typename T2, typename T3, typename T4>
	inline const string valuesToString(long vals1, long vals2, long vals3, long vals4){
	
		return bytesToString( vals1*sizeof(T1) + vals2*sizeof(T2)+vals3*sizeof(T3) +vals4*sizeof(T4) );
	}
	
	template<typename T1, typename T2, typename T3, typename T4, typename T5>
	inline const string valuesToString(long vals1, long vals2, long vals3, long vals4, long vals5){
	
		return bytesToString( vals1*sizeof(T1) + vals2*sizeof(T2)+vals3*sizeof(T3) +vals4*sizeof(T4)+vals5*sizeof(T5) );
	}
	
	
	
	inline const string secsToString(long secs){
	
		return msToString(secs*1000);
	}
};



/*
class TextTable{

	public:

	vector<vector<const string>> rows;
	
	vector<int> columnMaxLengths;
	
	int colnum;
	const string sep = "  ";
	
	inline TextTable(int colnum) : columnMaxLengths(colnum, 0), colnum(colnum){
	
		
	}
	
	inline void addRow(vector<const string> cols){
	
		for (int i = 0; i < colnum; i++)
			columnMaxLengths.at(i) = UtilMath::maxVal<int>(columnMaxLengths.at(i), cols.at(i).size() );
			
		rows.push_back(cols);	
	}
	
	inline const string toString(){
	
		std::stringstream sbuf;
		
		for (int i = 0; i < rows.size(); i++){
		
			vector<const string>& row = rows.at(i);
		
			for (int j = 0; j < colnum; j++){
		
				const string s = row.at(j);
				
				sbuf << s;
				
				for (int k = s.size(); k < columnMaxLengths.at(j); k++)
					sbuf << ' ';
					
				sbuf << sep;
			}
			
			sbuf << '\n';
		}
		
		return sbuf.str();
	}
	
	
	
	
	

};*/


struct Expression{

	vector<int> children;

	string s;
	
	int variableId;
	int constantId;
	int id;
	int op;
	
	int parent = -1;
	
	vector<const Expression*> childrenPtrs;
	int childrenNum;
	bool nonMath;
	bool lenient = false;
	int mathOp;
};

/*
enum OperatorType{

	IDENTITY, LENIENT_MUL, NON_MATH, NEG, SUM, SUB, MUL, DIV, MOD, BETWEEN, ABS, FLOOR, CEIL, ROUND, LOG, EXP, POW, SQRT, AND, OR, LT, LEQ, EQ, GT, GEQ


};*/

namespace UtilExpression{
	
	
	const int IDENTITY = 0;
	
	
	
	const int LENIENT_MUL = 1;
	const int NON_MATH = 2;
	
	
	const int NEG = 3;
	const int SUM = 4;
	const int SUB = 5;
	const int MUL = 6;
	const int DIV = 7;
	
	/*
	const long ACOS = 'a'*'c';
	const long ATAN = 'a'*'t';
	const long ASIN = 'a'*'s';
	const long SIN = 's';
	const long COS = 'c';
	const long TAN = 't';
	
	const long SIN = 's';
	const long COS = 'c';
	const long TAN = 't';*/
	
	const int MOD = 8;
	
	const int BETWEEN = 9;
	
	const int ABS = 10;
	const int FLOOR = 11;
	const int CEIL = 12;
	const int ROUND = 13;
	const int LOG = 14;
	const int EXP = 15;
	const int POW = 16;
	const int SQRT = 17;
	
	const int AND = 18;
	const int OR = 19;
	
	const int LT = 20;
	const int LEQ = 21;
	const int EQ = 22;
	const int GT = 23;
	const int GEQ = 24;
	
	const int ONE = 25;
	const int IDENTITY0 = 26;
	const int IDENTITY1 = 27;
	
	const int FLIPDIV = 28;
	
	const int NEQ = 29;
	
	const int INVERTED0 = 30;
	const int INVERTED1 = 31;
	
	const long_double NOT_A_NUMBER = std::numeric_limits<long_double>::max();
	
	inline long_double sum(const vector<long_double>& v){
	
		long_double ret = 0;
			
		for (int j = 0; j < v.size(); ++j)
			ret += v[j];
				
		return ret;
	}
	
	inline long_double prod(const vector<long_double>& v){
	
		long_double ret = 1;
			
		for (int j = 0; j < v.size(); ++j)
			ret *= v[j];
				
		return ret;
	}
	
	inline long_double lenientProd(const vector<long_double>& v){
	
		long_double ret = 1;
			
		for (int j = 0; j < v.size(); ++j)
			if ( v[j] != NOT_A_NUMBER)
				ret *= v[j];
			
		return ret;
	}
	
	
	
	
	
	
	
	
	
	inline long_double eval(bool lenient, const int& op, const long_double& v0){
	
		
		if (v0 == NOT_A_NUMBER)
			return lenient ? 1 : NOT_A_NUMBER;
		
		switch (op){
		
			case NON_MATH:
				return NOT_A_NUMBER;
			
			case ONE:
				return 1;
				
			case IDENTITY0:
				return v0;
			
			case IDENTITY:
				return v0;
				
			case AND:
				return (v0 != 0);
			
			case OR:
				return (v0 != 0);
			
			case ABS:
				return abs(v0);
		
			case NEG:
				return -v0;
				
			case FLOOR:
				return floor(v0);
		
			case CEIL:
				return ceil(v0);
		
			case SQRT:
				return sqrt(v0);
			
			case EXP:
				return exp(v0);
		
			case LOG:
				return log(v0);
		}
		
		cout << "OP not found " << op << endl;
		assert (false);
	}
	
	inline long_double eval(bool lenient, const int& op, const long_double& v0, const long_double& v1){
	
		
		const bool n0 = v0 == NOT_A_NUMBER;
		const bool n1 = v1 == NOT_A_NUMBER;
	
		if(lenient){
		
			if (n0 && n1)
				return 1;
	
			switch (op){
	
				case IDENTITY0:
					return !(n0) ? v0 : 1;
	
				case IDENTITY1:
					return !(n1) ? v1 : 1;
			
				case INVERTED0:
					return (!n0) ? 1/v0 : 1;
				
				case INVERTED1:
					return (!n1) ? 1/v1 : 1;
				
				case DIV:
					return (!n1) ? 1/v1 : (!n0) ? v0 : v0/v1;
				
				case FLIPDIV:
					return (!n1) ? v1 : !(n0) ? 1/v0 : v1/v0;
				
				case SUB:
					return (!n1) ? -v1 : !(n0) ? v0 : v0-v1;
			}
			
			cout << "OP NOT FOUND: " << op << endl;
			
			assert (false);
		
		} else {
			
			if (n0 || n1)
				return NOT_A_NUMBER;
	
		
			switch (op){
	
				case IDENTITY0:
					return v0;

				case IDENTITY1:
					return v1;
		
				case INVERTED0:
					return 1/v0;
		
				case INVERTED1:
					return 1/v1;
	
				case SUB:
					return v0-v1;
		
				case DIV:
					return v0/v1;
	
				case FLIPDIV:
					return v1/v0;
		
				case POW:
					return pow(v0, v1);
	
				case MOD:	
					return ((long) v0) % ((long)v1);
	
				case LOG:
					return log(v0)/log(v1);
		
				case EQ:
					return v0 == v1;
		
				case LT:
					return v0 < v1;
	
				case GT:
					return v0 > v1;
		
				case LEQ:
					return v0 <= v1;
			
				case NEQ:
					return v0 != v1;
		
				case GEQ:
					return v0 >= v1;
	
				case AND:
					return (v0 != 0) && (v1 != 0);
		
				case OR:
					return (v0 != 0) || (v1 != 0);
			}
	
		}
		
		cout << "lenient " << lenient << " mathop " << op << endl;
		assert (false);
	}
	
	inline long_double eval(bool lenient, const int& op, const long_double& v0, const long_double& v1, const long_double& v2){
	
		switch (op){
		
			case BETWEEN:
			
				if ( (v0 == NOT_A_NUMBER) ||(v1 == NOT_A_NUMBER) || (v2 == NOT_A_NUMBER))
					return NOT_A_NUMBER;
			
				return (v1 <= v0) && (v0 <= v2);
		}
		
		assert (false);
	}
	
};

struct ExpressionOperator{

	string name;
	int id;
	string prefix;
	string midA;
	string midB;
	string suffix;
	int min;
	int max;
	
	bool active = true;
	
	int mathOp = 0;
	
	bool lenient = false;
	
};

class ExpressionParser{

public:
	vector<ExpressionOperator> ops;
	
	string stringbuf;
	
	unordered_map<string, long> catalog;
	
	vector<Expression> expressions;
	
	vector<string> expressionNames;
	
	vector<string> constants;
	
	unordered_map<string, long> constantsIndex;
	
	
	vector<string> variables;
	
	unordered_map<string, long> variablesIndex;
	
	std::set<string> protectedStrings;
	
	inline const vector<string>& getVariables() const{
	
		return variables;
	}
	
	inline const vector<string>& getConstants() const{
	
		return constants;
	}
	
	
	inline void renameVariable(const string& var, const string& alias, long expressionid = -1){
	
		if (expressionid == -1){
		
			if (variablesIndex.count(var) == 0)
				return;
				
			if (variablesIndex.count(alias) == 0){
				
				variablesIndex.insert( std::make_pair(alias, variablesIndex.at(var) ));
				return;
			}
		}
		
		if (expressions.size() == 0)
			return;
		
		if (expressionid == -1)
			expressionid = expressions.size()-1;
		
			
		const vector<int>& v = getChildren(expressionid);
		
		if (v.size() == 0){
		
			Expression& e = expressions.at(expressionid);
			
			
			if (e.variableId != -1 && UtilString::equals(e.s, var))
				e.variableId = variablesIndex.at(alias);
				
			variablesIndex.at(var) = variablesIndex.at(alias);
		}
		
		for (auto it = v.begin(); it != v.end(); ++it)
			renameVariable(var, alias, (*it));
		
	}
	
	inline const vector<int>& getChildren(int id) const {
	
		return expressions.at(id).children;
	}
	
	inline bool isOp(const Expression& e, const string& q) const{
	
		return UtilString::equals( ops.at(e.op).name, q);
	}
	
	inline int getConstantId(const string& s) const{

		if (constantsIndex.count(s) > 0)
			return constantsIndex.at(s);
			
		return -1;
	}
	
	inline int getLastExpressionId() const{
	
		return expressions.size()-1;
	}
	
	inline int getVariableId(const string& s) const{

		if (variablesIndex.count(s) > 0)
			return variablesIndex.at(s);
			
		return -1;
	}
	
	/*
	inline long_double eval2(const vector<long_double>& constVals, const vector<long_double>& variableVals, const Expression& e, bool verbose=false) const{
	
		const long_double ret = eval2(constVals, variableVals, e, verbose);
		
		if (verbose){
		
			const ExpressionOperator& op = ops.at(e.op);
			
			//cout << repr(e) << endl;
			cout << "eval("<< op.name << "(" << e.s << ") "<< e.mathOp << (e.lenient ? "L" : "S") << ") == " <<ret << endl;
		}
		
		return ret;
	}
	*/
	inline long_double eval(const vector<long_double>& constVals, const vector<long_double>& variableVals, const Expression& e, const bool verbose=false) const{
		
		
		const long_double def = e.lenient ? 1 : UtilExpression::NOT_A_NUMBER;
		const long_double defsum = e.lenient ? 0 : UtilExpression::NOT_A_NUMBER;
		
		switch (e.mathOp){
		
			case UtilExpression::SUM:
				
				{
				switch (e.childrenNum){
				
					
				
					case 0:
						return defsum;
						
					case 1:
						{
						const Expression& e0 = *(e.childrenPtrs[0]);
						return e0.nonMath ? defsum: eval( constVals, variableVals , e0, verbose );
						}
				
					default:
						
						{
						long_double ret = 0;
		
						for (int j = 0; j < e.childrenNum; ++j){
	
							const Expression& e0 = *(e.childrenPtrs[j]);
							const long_double val = e0.nonMath ? UtilExpression::NOT_A_NUMBER : eval(  constVals, variableVals , e0, verbose );
				
							if (val != UtilExpression::NOT_A_NUMBER)
								ret += val;
							else if (!e.lenient)
								return UtilExpression::NOT_A_NUMBER;
										
						}
		
						return ret;
						}
						
				};
				}
		
			case UtilExpression::MUL:
				
				{
				//if (verbose)
				//	cout << "MUL" << e.childrenNum << endl;
				
				switch (e.childrenNum){
				
					case 0:
						return def;
						
					case 1:
						{
						const Expression& e0 = *(e.childrenPtrs[0]);
						return e0.nonMath ? def : eval(  constVals, variableVals , e0, verbose );
						}
				
					default:
					
						{
						long_double ret = 1;
		
						for (int j = 0; j < e.childrenNum; ++j){
	
							const Expression& ex = *(e.childrenPtrs[j]);
							
							if (verbose)
								cout << "mul" << e.lenient << ".eval[" << j << "/" << e.childrenNum << "]( " << ex.mathOp << " " << ops.at(ex.op).name << ex.lenient << ") = ?" << endl;
							
							const long_double val = ex.nonMath ? UtilExpression::NOT_A_NUMBER : eval(  constVals, variableVals , ex, verbose );
				
							//if (verbose)
							
							if (verbose)
								cout << "/mul" << e.lenient << ".eval[" << j << "/" << e.childrenNum << "]( " << ex.mathOp << " " << ops.at(ex.op).name << ex.lenient << ") = " << val << endl;
				
							if (val != UtilExpression::NOT_A_NUMBER)
								ret *= val;
							else if (!e.lenient)
								return UtilExpression::NOT_A_NUMBER;
										
						}
		
						return ret;
						}
						
				};
				}
		
			default:
			
				{
				switch (e.childrenNum){
				
					case 0: 
					
						if (e.nonMath)
							return e.lenient ? 1 : UtilExpression::NOT_A_NUMBER;
					
						if (e.constantId != -1)
							return constVals.at(e.constantId);
		
						if (e.variableId != -1){
						
							if (verbose)
								cout << "VAR " << e.s << " = " << variableVals.at(e.variableId) << endl;
						
							return variableVals.at(e.variableId);
						}
				
						return UtilExpression::NOT_A_NUMBER;
				
					case 1:
						{
						const Expression& e0 = *e.childrenPtrs[0];
			
						const long_double v0 = e0.nonMath ? def : eval( constVals, variableVals , e0, verbose);
			
						return UtilExpression::eval(e.lenient, e.mathOp, v0);
						}
				
					case 2:
					
						{
						const Expression& e0 = *e.childrenPtrs[0];
						const Expression& e1 = *e.childrenPtrs[1];
			
						const long_double v0 = e0.nonMath ? def : eval( constVals, variableVals , e0, verbose);
						const long_double v1 = e1.nonMath ? def : eval( constVals, variableVals , e1, verbose);
		
						return UtilExpression::eval(e.lenient, e.mathOp, v0, v1);
						}
						
					case 3:
					
						{
						const Expression& e0 = *e.childrenPtrs[0];
						const Expression& e1 = *e.childrenPtrs[1];
						const Expression& e2 = *e.childrenPtrs[2];
		
						const long_double v0 = e0.nonMath ? def : eval( constVals, variableVals , e0, verbose);
						const long_double v1 = e1.nonMath ? def : eval( constVals, variableVals , e1, verbose);
						const long_double v2 = e2.nonMath ? def : eval( constVals, variableVals , e2, verbose);
			
						return UtilExpression::eval(e.lenient, e.mathOp, v0, v1, v2);	
						}
				
					default:
					
						return def;
				}
				}
		};
	}
	
	
	inline long_double eval(const vector<long_double>& constVals, const vector<long_double>& variableVals, const bool verb=false) const{
		
		if (expressions.size() == 0)
			return 1;
		
		if (expressions.back().nonMath)
			return UtilExpression::NOT_A_NUMBER;
		
		return eval(constVals, variableVals, expressions.back(), verb );
	}
	
	inline void addOp(string name, string prefix, string midA, string midB, string suffix, int min, int max, int mathop){
	
		ExpressionOperator op;
		
		op.name = name;
		op.id = ops.size();
		
		op.lenient = mathop == UtilExpression::LENIENT_MUL;
		
		op.mathOp = mathop == UtilExpression::LENIENT_MUL ? UtilExpression::MUL : mathop;
		
		if (UtilString::equals(name, "var")){
		
			op.prefix = prefix;
			op.midA = midA;
			op.midB = midB;
			op.suffix = suffix;
		
		} else {
		
			op.prefix = UtilString::toLowerCase(prefix);
			op.midA = UtilString::toLowerCase(midA);
			op.midB = UtilString::toLowerCase(midB);
			op.suffix = UtilString::toLowerCase(suffix);
		}
		
		
		op.min = min;
		op.max = max;
		
		ops.push_back(op);
		
		if (min >= 0)
		for (int capcase = -1; capcase <= 1; ++capcase){
		
			const string prefix2 = capcase == -1 ? UtilString::toLowerCase(prefix) : capcase == 0 ? prefix : UtilString::toUpperCase(prefix) ;
			const string midA2 = capcase == -1 ? UtilString::toLowerCase(midA) : capcase == 0 ? midA : UtilString::toUpperCase(midA) ;
			const string midB2 = capcase == -1 ? UtilString::toLowerCase(midB) : capcase == 0 ? midB : UtilString::toUpperCase(midB) ;
			const string suffix2 = capcase == -1 ? UtilString::toLowerCase(suffix) : capcase == 0 ? suffix : UtilString::toUpperCase(suffix) ;
		
			if ((prefix2.length()) > 0 && (protectedStrings.count(prefix2) == 0))
				protectedStrings.insert(prefix2);
			
			if ((midA2.length()) > 0 && (protectedStrings.count(midA2) == 0))
				protectedStrings.insert(midA2);
			
			if ((midB2.length()) > 0 && (protectedStrings.count(midB2) == 0))
				protectedStrings.insert(midB2);
			
			if ((suffix2.length()) > 0 && (protectedStrings.count(suffix2) == 0))
				protectedStrings.insert(suffix2);
		
		}
		
		/*
		if (rec)
		for (int capcase = 0; capcase <= 1; ++capcase){
		
			for (int pre = 0; pre <= 1; ++pre){
		
				const string prefix2 = pre == 0 ? prefix : capcase == 0 ? UtilString::toLowerCase(prefix) : UtilString::toUpperCase(prefix) ;
				
				if (pre == 1 && UtilString::equals(prefix, prefix2))
					continue;
		
				for (int mi = 0; mi <= 1; ++mi){
		
					const string mid2 = mi == 0 ? mid : capcase == 0 ? UtilString::toLowerCase(mid) : UtilString::toUpperCase(mid) ;
				
					if (mi == 1 && UtilString::equals(mid, mid2))
						continue;
		
					for (int suf = 0; suf <= 1; ++suf){
		
						const string suffix2 = suf == 0 ? suffix : capcase == 0 ? UtilString::toLowerCase(suffix) : UtilString::toUpperCase(suffix) ;
						
						if (suf == 1 && UtilString::equals(suffix, suffix2))
							continue;
		
						addOp(name, prefix2, mid2, suffix2, min, max, false);
						
						//cout << prefix2 << " " << mid2 << " " << suffix2 << endl;
					}
				}
		
			}
		}*/
		
	}
	
	
	inline void addOp(string name, string prefix="", string midA="",string suffix="", int min=-1, int max=-1, long mathop = 0){
	
		addOp(name, prefix, midA, midA, suffix, min, max, mathop);
	}
	
	inline void addElementaryOps(){
	
		ExpressionParser& p = (*this);
	
		p.addOp("str", "'", "", "", "'", -1, -1, UtilExpression::IDENTITY);
		
		p.addOp("str", "\"", "", "", "\"", -1, -1,  UtilExpression::IDENTITY);
		p.addOp("var", "", "azAZ__", "azAZ09..__", "", -1, -1,  UtilExpression::IDENTITY);
		p.addOp("num");
		
	}
	
	inline void addBasicOps(){
	
		ExpressionParser& p = (*this);
	
		p.addOp("str", "'", "", "", "'", -1, -1, UtilExpression::IDENTITY);
		
		p.addOp("str", "\"", "", "", "\"", -1, -1 , UtilExpression::IDENTITY);
		p.addOp("var", "", "azAZ__[[", "azAZ09..__]]", "", -1, -1, UtilExpression::IDENTITY);
		p.addOp("num");
	
		p.addOp("log", "log(", "", ")", 1, 2, UtilExpression::LOG);
		
		//p.addOp("ind", "ind(", "", ")", 1, 1);
		
		/*
		p.addOp("sin", "sin(", "", ")", 1, 1, UtilExpression::SIN);
		p.addOp("cos", "cos(", "", ")", 1, 1, UtilExpression::COS);
		p.addOp("tan", "tan(", "", ")", 1, 1, UtilExpression::TAN);
		
		p.addOp("sinh", "sinh(", "", ")", 1, 1, UtilExpression::SINH);
		p.addOp("cosh", "cosh(", "", ")", 1, 1, UtilExpression::COSH);
		p.addOp("tanh", "tanh(", "", ")", 1, 1, UtilExpression::TANH);*/
		
		
		p.addOp("between", "", "Between", "And", "", 3, 3, UtilExpression::BETWEEN);
		
		
		//p.addOp("date", "", "", "", 1, 1, UtilExpression::EXP);
		
		p.addOp("abs", "Abs(", "", ")", 1, 1, UtilExpression::ABS);
		
		p.addOp("floor", "Floor(", "", ")", 1, 1, UtilExpression::FLOOR);
		p.addOp("ceil", "Ceil(", "", ")", 1, 1, UtilExpression::CEIL);
		
		p.addOp("exp", "Exp(", "", ")", 1, 1, UtilExpression::EXP);
		
		p.addOp("pow", "Power(", ",", ")", 2, 2, UtilExpression::POW);
		
		p.addOp("pow", "Pow(", ",", ")", 2, 2, UtilExpression::POW);
		p.addOp("pow", "", "^", "", 2, 2, UtilExpression::POW);
		p.addOp("pow", "", "**", "", 2, 2, UtilExpression::POW);
		p.addOp("div", "", "/", "", 2, 2, UtilExpression::DIV);
		
		
		p.addOp("mul", "", "*", "", 2, 999, UtilExpression::MUL);
		p.addOp("add", "", "+", "", 2, 999, UtilExpression::SUM);
		p.addOp("sub", "", "-", "", 2, 2, UtilExpression::SUB);
		p.addOp("neg", "-", "", "", 1, 1, UtilExpression::NEG);
		
		p.addOp("mod", "", "%", "", 2, 2, UtilExpression::MOD);
		
		p.addOp("par", "(", "", ")", 1, 1, UtilExpression::LENIENT_MUL);
		
		//p.addOp("dot", "", ".", "", 2, 2);
		
		p.addOp("notequal", "", "<>", "", 2, 2, UtilExpression::NEQ);
		
		p.addOp("notequal", "", "!=", "", 2, 2, UtilExpression::NEQ);
		
		p.addOp("equal", "", "==", "", 2, 2, UtilExpression::EQ);
		
		p.addOp("equal", "", "=", "", 2, 2, UtilExpression::EQ);
		
		p.addOp("less", "", "<", "", 2, 2, UtilExpression::LT);
		p.addOp("greater", "", ">", "", 2, 2, UtilExpression::GT);
		
		p.addOp("leq", "", "<=", "", 2, 2, UtilExpression::LEQ);
		p.addOp("geq", "", ">=", "", 2, 2, UtilExpression::GEQ);
	
	}
	
	
	inline static const string protectedBegin(){
	
		return "[__";
	}
	
	inline static const string protectedEnd(){
	
		return "__]";
	}
	
	
	inline static bool valid(const char& c, const string& mask){
	
		assert ( (mask.length() & 1) == 0);
	
		for (int j = 0; j < mask.size(); j += 2){
		
			if ( (mask[j] <= c ) && (c <= mask[j+1]))
				return true; 
		}
		
		return false;
	}
	
	inline static const char* endOfName(const char* c, const string& mask1, const string& mask2){
	
		const char* ret = c;
	
		if (valid(*ret, mask1))
			++ret;
			/*
		if (UtilMath::isBetween(*ret, '_', 'z') || UtilMath::isBetween(*ret, 'A', 'Z')) 
			++ret;*/
			
		if (ret == c)
			return c;
		
		/*
		while ( UtilMath::isIn( true, UtilMath::isBetween(*ret, '0', '9'), UtilMath::isBetween(*ret, '_', 'z'), UtilMath::isBetween(*ret, 'A', 'Z')))
			++ret;*/
			
		while (valid(*ret, mask2))
			++ret;
			
		if (ret == c)
			return c;
			
		if ( (*ret) == '(' )
			return c;
			
		return ret;
	}
	
	
	inline static const char* endOfProtected(const char* c){
	
		assert (c != NULL);
		
		const string s1 = protectedBegin();
		const string s2 = protectedEnd();
	
		return UtilString::endOfRegion(c, s1, s2);
	}
	
	inline const char* endOfSkipped(const char* c){
	
		
		const char* ret = NULL;
		const char* nov = c;
		
		while (nov != ret ){
			
			ret = nov;
			
			nov = endOfProtected(ret);
			
			
			if (nov == ret){
				
				for (auto it = protectedStrings.begin(); it != protectedStrings.end(); ++it){
				
					const string& s = (*it);
					
					assert (s.length() > 0);
				
					nov = UtilString::endOf(nov, s);
					
					if (nov != ret)
						break;
				}
			}
		}
		return ret;
	}
	
	inline static const string replace(const char* beg, const char* end, const char* a, const char* b, const string& repl){
	
		stringstream ss;
		
		
		for (const char* ptr = beg; ptr != a; ++ptr)
			ss << (*ptr);
			
		ss << repl;
		
		for (const char* ptr = b; ptr != end; ++ptr)
			ss << (*ptr);
		
		return ss.str();
	}
	
	
	inline const string addExpression(const string& type, const ExpressionOperator& op, const vector<int>& children, const char* beg=NULL, const char* end=NULL){
	
		long id = expressions.size();
		string name = protectedBegin()+type+to_string(id)+protectedEnd();
		
		Expression ex;
		
		ex.variableId = -1;
		ex.constantId = -1;
		
		
		
		if (UtilString::equalsAny(op.name, "num", "str")){
		
			const string con = string(beg,end);
			
			if (constantsIndex.count(con) > 0){
			
				
			
				ex.constantId = constantsIndex.at(con);
				
				//cout << "found constant " << con << " at " << ex.constantId << endl;
			
			} else {
			
				ex.constantId = constants.size();
				constantsIndex.insert( std::make_pair(con, ex.constantId) );
				constants.push_back( con );
				
				//cout << "created constant " << con << " at " << ex.constantId << endl;
			}
		}
		
		if (UtilString::equals(op.name, "var")){
		
			const string var = string(beg,end);
			
			if (variablesIndex.count(var) > 0){
			
				ex.variableId = variablesIndex.at(var);
			
			} else {
			
				ex.variableId = variables.size();
				variablesIndex.insert( std::make_pair(var, ex.variableId) );
				variables.push_back( var );
			}
		}
		
		
		ex.id = id;
		
		if (beg != NULL && end != NULL)
			ex.s = string(beg,end);
			
		ex.op = op.id;
		ex.children = children;
		
		
		
		ex.mathOp = op.mathOp;
		ex.lenient = op.lenient;
		
		for (long j = 0; j < children.size(); ++j)
			expressions.at(ex.children[j]).parent = id;
		
		
		
		ex.childrenNum = ex.children.size();
		
		ex.nonMath = ex.mathOp == UtilExpression::NON_MATH;
		
		//ex.temp.reserve(children.size() );
		
		
		expressions.push_back(ex);
		
		
		
		
		expressionNames.push_back(name);
		
		
		catalog[name] = id;
		
		return name;
	}
	
	inline const Expression* getAncestor(const Expression& ex, const string& opName) const{
	
		const Expression* ptr = &ex;
		
		while (ptr->parent != -1){
		
			ptr = &expressions.at(ptr->parent);
			
			if (UtilString::equals(ops.at(ptr->op).name, opName))
				return ptr;
			
		}
		
		return NULL;
	
	}
	
	/*
	inline vector<Expression*> getAncestorPtrs(const Expression& ex){
		
		vector<Expression*> ret;
		ret.reserve(10);
		const Expression* ptr = &ex;
		
		while (ptr->parent != -1){
		
			ptr = &expression.at(ptr->parent);
			
			ret.push_back(ptr);
		}
		
		return ret;
	}*/
	
	inline bool parse(const ExpressionOperator& op){
	
		const char* beg = stringbuf.c_str();
		const char* end = beg+stringbuf.length();
		
		vector<int> children;
		
		
		
		if (UtilString::equals(op.name, "var")){
		
		
		
			for (const char* ptr = UtilString::endOfNumber(endOfSkipped(beg)); ptr != end; ptr = UtilString::endOfNumber(endOfSkipped(ptr+1)) ){
			
				
				const char* x = endOfName(ptr, op.midA, op.midB);
				
				
				if (x != ptr){
				
					
					const char* y = NULL;
				
					while (x != y){
					
						y = x;
						
						//cout << string(ptr, y) << endl;
						
						if (endOfSkipped(y) != y)
							break;
					
						while ( y != end && (*y) == ' ')
							++y;
						
						if (y == end)
							break;
						
						if (endOfSkipped(y) != y)
							break;
						
						const char* e = endOfName(y, op.midA, op.midB);
						
						if (e != y){
						
							x = e;
							
							if (x == end)
								break;	
						}
						
					}
					
					stringbuf = replace(beg, end, ptr, x, addExpression("var", op, children, ptr, x) );
					return true;
				}
			}
			
			return false;
		}
		
		
		if (UtilString::equals(op.name, "str")){
		
			for (const char* ptr = endOfSkipped(beg); ptr != end; ptr = endOfSkipped(ptr+1)){
		
				const char* x = UtilString::endOfRegion(ptr, op.prefix, op.suffix);
			
				if (x != ptr){
			
					stringbuf = replace(beg, end, ptr, x, addExpression("str", op, children, ptr+op.prefix.size(), x-op.suffix.size() ) );
					return true;
				}
			}
			return false;
		}
		
		if (UtilString::equals(op.name, "num")){
		
			for (const char* ptr = endOfSkipped(beg); ptr != end; ptr = endOfSkipped(ptr+1)){
			
				const char* x = UtilString::endOfNumber(ptr);
				
				if (x != ptr){
				
					stringbuf = replace(beg, end, ptr, x, addExpression("num", op,children, ptr, x) );
					return true;
				}
			}
			
			return false;
		}
		
		for (const char* ptr = beg; ptr != end; ++ptr){
		
			const char* p = ptr;
			
			children.clear();
			
			if ( op.prefix.length() > 0){
			
				const char* e = UtilString::endOf(p, op.prefix);
				
				if (e == p)
					continue;
				
				p = e;
			}
			
			int count = 0;
			
			while (count < op.max){
			
				
				{
					const char* e = endOfProtected(p);
			
					if (e == p)
						break;
					
					string s(p, e);
				
					children.push_back( catalog.at(s) );
					
					p = e;
				}
				
				++count;
				
				if (count >= op.max)
					break;
				
				const string& mid = count == 1 ? op.midA : op.midB;
				
				if (mid.length() > 0){
					const char* e = UtilString::endOf(p, mid);
					
					if (e == p)
						break;
					
					p = e;
				}
				
			}
			
			if (count < op.min)
				continue;
			
			if ( op.suffix.length() > 0){
			
				const char* e = UtilString::endOf(p, op.suffix);
				
				if (e == p)
					continue;
				
				p = e;
			}
			
			if (UtilString::equals(op.name, "par")){
			
				assert (children.size() == 1);
			
				stringbuf = replace(beg, end, ptr, p, expressionNames.at(children.at(0)) );
				return true;
			}
			
			
			stringbuf = replace(beg, end, ptr, p, addExpression(op.name, op, children) );
			return true;
		}
		
		
		return false;
	
	}
	
	
	inline const string repr() const{
		
		return repr(expressions.back() );
		
	}
	
	inline const string repr(const Expression& ex, int lvl = 0) const{
	
		stringstream ss;
		
		if (ex.children.size() > 1){
		
		ss << "$NL$";
		
		for (int j = 0; j < lvl; ++j)
			ss << "$WS$";
		}
	
		const ExpressionOperator& op = ops.at(ex.op);
		
		bool par = false;
		
		if (UtilString::equalsAny(op.name, "sum","mul", "sub", "div") && op.prefix.length()+op.suffix.length() == 0)
			par = true;
		
		if (UtilString::equalsAny(op.name, "neg"))
			par = true;
		
		if (ex.children.size() == 0 && UtilString::equalsAny(op.name, "str", "var", "num") ){
		
			if (UtilString::equals(op.name, "str")){
				ss << " ";
				ss << UtilString::toUpperCase(op.prefix);

			} else {
			
				ss << " ";
			}
		
			
			ss << ex.s;
			
			
			if (UtilString::equals(op.name, "str")){
			
				ss << op.suffix;
				ss << " ";
			} else {
				
				ss << " ";
			}
			
			return ss.str();
		}
		
		//if (op.prefix.length()*op.suffix.length() == 0)
		if (par)
			ss << /*op.name <<*/ "(" ;
		else
			ss << " ";
		
		ss << UtilString::toUpperCase(op.prefix);
		ss << " ";
		
		long count = 0;
		
		for (auto it = ex.children.begin(); it != ex.children.end(); ++it){
		
			const Expression& e = expressions.at( *it );
			
			ss << " ";
			ss << repr(e, lvl+1);
			ss << " ";
			
			if ( (it+1) != ex.children.end()){
			
				const string& mid = count == 0 ? op.midA : op.midB;
				
				ss << " ";
				ss << UtilString::toUpperCase(mid);
				ss << " ";
			}
			count++;
		}
		
		if (par)
			ss << ")" /*<< op.name*/;
		else
			ss << " ";
		ss << UtilString::toUpperCase(op.suffix);
		ss << " ";
		
		if (lvl == 0){
		
			string ret = ss.str();
		
			ret.reserve(ret.length()*2);
		
			ret = UtilString::cleanStringList(ret, ' ', ' ');
			ret = UtilString::replaced(ret, " , ", ", ");
			ret = UtilString::replaced(ret, " ; ", "; ");
			ret = UtilString::replaced(ret, "( ", "(");
			ret = UtilString::replaced(ret, " )", ")");
		
			long j = 0;
		
			while (j < ret.length() && ret[j] < '!')
				++j;
		
			ret = ret.substr(j);
		
			j = ret.length()-1;
		
			while (j >= 0 && ret[j] < '!')
				--j;
		
			ret = ret.substr(0, j+1);
			
			ret = UtilString::replaced(ret, "$NL$", "\n");
			ret = UtilString::replaced(ret, "$WS$", "  ");
			
			return ret;
		
		}
		
		
		return ss.str();
		//if (op.prefix.length()*op.suffix.length() == 0)
		//	ss << ")";
		
		
		
	}
	
	inline const string toString() const{
	
		return toString(expressions.back() );
		
	}
	
	inline const string toString(const Expression& ex) const{
	
		if (true)
			return repr(ex);
		
		stringstream ss;
	
		const ExpressionOperator& op = ops.at(ex.op);
		
		ss << "{" << ex.mathOp << "}";
	
		if (ex.children.size() == 0 && UtilString::equalsAny(op.name, "str", "var", "num") ){
		
			if (UtilString::equals(op.name, "str"))
				ss << op.prefix;
		
			ss << ex.s;
			
			if (UtilString::equals(op.name, "str"))
				ss << op.suffix;
			
			return ss.str();
		}
		
		
		
		if (op.prefix.length()*op.suffix.length() == 0)
			ss << "(";
		
		ss << op.prefix << " ";
		
		long count = 0;
		
		for (auto it = ex.children.begin(); it != ex.children.end(); ++it){
		
			const Expression& e = expressions.at( *it );
			
			ss << toString(e);
			
			
			
			if ( (it+1) != ex.children.end()){
			
				const string& mid = count == 0 ? op.midA : op.midB;
				
				if (ex.lenient)
					ss << " " << mid << "' ";
				else
					ss << " " << mid << " ";
			}
			count++;
		}
		
		ss << " " << op.suffix;
		
		if (op.prefix.length()*op.suffix.length() == 0)
			ss << ")";
		
		return ss.str();
	}
	
	bool verbose = false;
	inline bool parse(const string& input){
	
		vector<std::pair<int, const char*>> toSort;
		
		
		expressions.reserve(100);
		
		for (auto it = protectedStrings.begin(); it != protectedStrings.end(); ++it)
			toSort.push_back( std::make_pair( -(*it).length(), (*it).c_str() ) );
		
		
		std::sort(toSort.begin(), toSort.end() );
		
		vector<string> whiteList;
		
		whiteList.reserve(toSort.size() );
		
		for (auto it = toSort.begin(); it != toSort.end(); ++it){
		
			whiteList.push_back( string(it->second) );
		}
		
		
	
		bool changed = true;
		
		stringbuf = input;
		
		//stringbuf.reserve(10000);
	
		bool cleanWhiteSpaces = false;
	
		while (changed){
		
			if (verbose)
				cout << stringbuf << endl;
		
			changed = false;
			
			for (auto it = ops.begin(); it != ops.end(); ++it){
			
				const ExpressionOperator& const_op = (*it);
				
				if (!const_op.active)
					continue;
			
				if (!cleanWhiteSpaces && const_op.min >= 0)
					continue;

				if (cleanWhiteSpaces && const_op.min < 0)
					continue;
				
				if (parse(const_op)){
				
					changed = true;
					break;
					
				} else {
					
					/*
					if (const_op.min >= 0){
					
						ExpressionOperator& op = ops.at(const_op.id);
					
						if (op.active && ( op.prefix.length() >= 0 ) )
							op.active = stringbuf.find(op.prefix) != std::string::npos;
					
						if (op.active && ( op.min >= 2 && op.midA.length() > 0 ) )
							op.active = stringbuf.find(op.midA) != std::string::npos;
					
						if (op.active && ( op.min >= 3 && op.midB.length() > 0 ) )
							op.active = stringbuf.find(op.midB) != std::string::npos;
						
						if (op.active && ( op.suffix.length() > 0 ) )
							op.active = stringbuf.find(op.suffix) != std::string::npos;
						
					}*/
					
				}	
			}
		
			if (!changed && !cleanWhiteSpaces){
				
				vector<string> protec;
				
				{
				bool chang = true;
				
				while (chang){
				
					chang = false;
					
					const char* beg = stringbuf.c_str();
					const char* end = beg+stringbuf.length();
					
					for (auto it = whiteList.begin(); it != whiteList.end(); ++it){
					
						const string& its = (*it);
						
						bool br = false;
						
						for (const char* ptr = beg; ptr != end; ++ptr){
						
							const char* p = UtilString::endOf(ptr, its);
							
							if (ptr != p){
							
								br = true;
								
								protec.push_back( string( ptr, p) );
								
								long j = protec.size()-1;
								
								stringbuf = replace(beg, end, ptr, p, protectedBegin()+"temp"+to_string(j+expressions.size())+protectedEnd() );
								
								chang = true;
								break;
							}
						}
						
						if (br)
							break;
					}
				
				}
				}
				
				stringbuf = UtilString::cleanStringList(stringbuf, ' ', 0);
				
				{
				
				for (long j = 0; j < protec.size(); ++j){
				
					const string torepl = protectedBegin()+"temp"+to_string(j+expressions.size())+protectedEnd();
				
					const char* beg = stringbuf.c_str();
					const char* end = beg+stringbuf.length();
					
					for (const char* ptr = beg; ptr != end; ++ptr){
						
						const char* p = UtilString::endOf(ptr, torepl);
						
						if (ptr != p){
						
							protec.push_back( string( ptr, p) );
							
							stringbuf = replace(beg, end, ptr, p, protec.at(j) );
							
							break;
						}
					}
				}
				}
				
				stringbuf = UtilString::toLowerCase(stringbuf);
				
				cleanWhiteSpaces = true;
				changed =true;
			}
		}
		
		expressions.shrink_to_fit();
		
		for (auto it = expressions.begin(); it != expressions.end(); ++it){
		
			Expression& ex = (*it);
			
			ex.childrenPtrs.clear();
			
			for (long j = 0; j < ex.children.size(); ++j)
				ex.childrenPtrs.push_back( &expressions.at(ex.children.at(j) ) );
		}
		
		const bool ret =  UtilString::equals(stringbuf, protectedBegin()+ops.at(expressions.back().op).name+to_string(expressions.back().id)+protectedEnd()) ;
		
		if (!ret){
		
			cout << stringbuf << endl;
		}
		
		assert ( ret);
		
		return ret;
		
		//cout << stringbuf;
		
	}
	
	inline void print() const{
	
		if (expressions.size() > 0){
			cout << expressions.at(expressions.size()-1).children.size() << endl;
		
			cout << expressionNames.back() << " = "<< toString(expressions.back() ) << endl;
		}
		
		for (auto it = variables.begin(); it != variables.end(); ++it)
			cout << (*it) << endl;
	}
};


class SQLParser : public ExpressionParser{

public:
	inline SQLParser() : ExpressionParser(){
		
		ExpressionParser& p = (*this);
	
		p.addBasicOps();
		
		p.addOp("select", "Select *", "", "", 0, 0, UtilExpression::LENIENT_MUL);
		
		
		p.addOp("as", "", "As", "", 2, 2, UtilExpression::NON_MATH);
		
		
		p.addOp("select", "As Select", ",", "", 1, 999, UtilExpression::NON_MATH);
		
		p.addOp("limit", "Limit", "", "", 1, 1, UtilExpression::NON_MATH);
		
		//p.addOp("and", "", "And", "", 2, 999);

		//p.addOp("comma", "", ",", "", 2, 999);
		
		p.addOp("rand", "Rand(", "", ")", 0, 0, UtilExpression::NON_MATH);
		
		p.addOp("order", "Order By", "", "", 1, 1, UtilExpression::NON_MATH);
		
		
		
		
		
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
		
		p.addOp("view", "Create Or Replace View", "", ";", 1, 999, UtilExpression::LENIENT_MUL);
		
		
		p.addOp("list", "", "", ";", 1, 999, UtilExpression::LENIENT_MUL);
		
		p.addOp("comment", "/*", ",", "*/", 1, 999, UtilExpression::NON_MATH);
		
		
		p.addOp("list", "", "", "", 2, 999, UtilExpression::LENIENT_MUL);		
	}
};



#endif