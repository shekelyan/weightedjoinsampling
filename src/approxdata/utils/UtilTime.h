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

#ifndef SHEKELYAN_UTILTIME_H
#define SHEKELYAN_UTILTIME_H

#include <approxdata/utils/utils.h>


namespace UtilTime{

	inline void readDate(const char* s, long length, int& year, int& month, int& day){
	
		if (length == 10 && ( (s[4] == '-' && s[7] == '-') ||(s[4] == '/' && s[7] == '/') ) ){
			
			year = UtilString::stringToNonNegativeInteger(s, length, 0, 4);
			month = UtilString::stringToNonNegativeInteger(s, length, 5, 7);
			day = UtilString::stringToNonNegativeInteger(s, length, 8, 10);
		
		} else if(length == 10 && (s[2] == '.' && s[5] == '.') ){
				
			
			day = UtilString::stringToNonNegativeInteger(s, length, 0, 2);
			month = UtilString::stringToNonNegativeInteger(s, length, 3, 5);
			year = UtilString::stringToNonNegativeInteger(s, length, 6, 10);
			
			
		} else {
	
			year = -1,
			month = -1;
			day = -1;
		
		}
	}
	
	inline void readDate(const string& s, int& year, int& month, int& day){
	
		return readDate(s.c_str(), s.length(), year, month, day);
	}

	inline int leapYearsSince1900(int y){
		
		assert (y >= 1900);
		assert (y <= 2100);
	
		if (y == 2100)
			return 49;
		
		return (floor(0.25*y)*4-1904)/4+1;
	}
	
	
	inline int daysSinceBeginningOfYear(int y, int m, int d){
	
		assert (m >= 1);
		assert (m <= 12);
	
		if (m == 1)
			return d;
			
		if (m == 2)
			return 31+d;
			
		const bool isLeapYear = (y > 1900 && y < 2100) && y % 4 == 0;
		
		assert (y >= 1900);
		assert (y <= 2100);
		
		const long february = isLeapYear ? -1 : -2;
		
		if (m >= 8)
			return 31*4+30*3+february + 31*((m-7)/2)+30*((m-8)/2)+d;
		
		return 31*((m)/2)+30*((m-1)/2)+february+d;
	}
	
	
	inline int daysSince1900(int y, int m, int d){
	
		if (y == 1900)
			return daysSinceBeginningOfYear(y, m,d);
	
		const int leapyears = leapYearsSince1900(y-1);
		const int years = y-1900;
		
		return (years-leapyears)*365+leapyears*366+daysSinceBeginningOfYear(y, m,d);
	}


	inline int daysBetween(int y1, int m1, int d1, int y2, int m2, int d2){
		
		/*
		struct tm t1 = {0,0,0,d1,m1-1,y1-1900};
    	struct tm t2 = {0,0,0,d2,m2-1,y2-1900};
    	time_t x = mktime(&t1);
    	time_t y = mktime(&t2);
    	//cout <<  << endl;
		
		return difftime(y,x)/3600/24;*/
		
		return daysSince1900(y2,m2,d2)-daysSince1900(y1,m1,d1);
	}
	
	inline long daysSince1900(const char* s, int length){

		int y = -1;
		int m = -1;
		int d = -1;
		
		readDate(s, length, y,m,d);
		
		if (y == -1)
			return -1;
		
		if (m == -1)
			return -1;
		
		if (d == -1)
			return -1;
		
		return daysSince1900(y,m,d);
	}
	
	inline long daysSince1900(const string& s){

		int y = -1;
		int m = -1;
		int d = -1;
		
		readDate(s, y,m,d);
		
		if (y == -1)
			return -1;
		
		if (m == -1)
			return -1;
		
		if (d == -1)
			return -1;
		
		return daysSince1900(y,m,d);
	}

	inline long daysTill2100(const string& s){

		int y = -1;
		int m = -1;
		int d = -1;
		
		readDate(s, y,m,d);
		
		assert (y != -1);
		assert (m != -1);
		assert (d != -1);
		
		return daysBetween(y,m,d, 2100,1,1);
	}
	
	inline void daysSince1900test(){
	
		int last = 0;
		
		vector<int> months = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
		
		for (int year = 1900; year <= 2100; year++){
		
			months[2] = 28;
		
			if (year > 1900 && year < 2100 && (year % 4 == 0)){
			
				months[2] = 29;
			}
			
			for (int m = 1; m <= 12; m++){
		
				for (int d = 1; d <= months[m]; d++){
			
					int now = daysSince1900(year, m, d);
				
					cout << year << "-" << m << "-" << d << ": " << now << endl; 
					assert (now-last == 1);
				
					last = now;
				}
			}
		}
	}

	inline void dateTimeExtract(long t, long& y, long& m, long& d, long& h, long& min, long& sec){
	
		struct tm * now = localtime( & t );
		
		y = now->tm_year+1900;
		m = now->tm_mon+1;
		d = now->tm_mday;
		h = now->tm_hour;
		min = now->tm_min;
		sec = now->tm_sec;
	}
	
	inline long dateTimeNow(){
	
		return time(0);
	}
	
	inline const string dateTimeString(long l, bool printOrig = false){
	
		long y, m, d, h, min, sec;
		
		dateTimeExtract(l, y, m, d, h, min, sec);
		
		return "["+UtilString::twoDigit(h)+":"+UtilString::twoDigit(min)+"]["+UtilString::twoDigit(sec)+"s]";
	
		//return to_string(y)+"-"+twoDigit(m)+"-"+twoDigit(d)+"_"+twoDigit(h)+":"+twoDigit(min)+"::"+twoDigit(sec)+(printOrig ? (" ["+to_string(l)+"]") : "");
	}
	
	inline void dateTimePrint(bool printOrig = true){
	
		cout << dateTimeString( dateTimeNow(), printOrig) << endl;
	}
	
	inline string dateTimeNowString(bool printOrig = true){
	
		return dateTimeString( dateTimeNow(), printOrig);
	}

};


class Timeinterval{

public:
	const double ms;

	inline Timeinterval(double ms) : ms(ms){
	
	
	}
	
	inline double getMilliseconds() const{
	
		return ms;
	}
	
	inline double getSeconds() const{
	
		return ms/1000.0;
	}
	
	inline double getMinutes() const{
	
		return ms/(1000.0*60.0);
	}
	
	inline double getHours() const{
	
		return ms/(1000.0*60.0*60.0);
	}
	
	inline const string toString() const{
	
		return UtilString::msToString(ms);
	}

	friend inline std::ostream& operator<<(std::ostream &os, const Timeinterval& m) { 
    	
    	return os << m.toString();
	}
};

class Timestamp{

public:
	
	
	const string name;
	const std::chrono::time_point<std::chrono::system_clock> t;
	
	inline Timestamp(const string name = "") : name(name), t(std::chrono::system_clock::now()){
	
		
	}
	
	
	inline Timestamp(long n, const string name = "") : name(name), t(std::chrono::system_clock::from_time_t( (std::time_t) n)){
	
		
	}
	
	inline Timestamp(const Timestamp& ts) : name(ts.name), t(ts.t){
	
		
	}
	
	inline bool operator<(const Timestamp& other) const{ // copy assignment{
    
    	return t < other.t;
	}
	
	inline double msPassed() const{
	
		const std::chrono::time_point<std::chrono::system_clock> tnow = std::chrono::system_clock::now();
	
		return std::chrono::duration_cast<std::chrono::milliseconds>(t-tnow).count();
	}
	
	inline double msPassed(const Timestamp& ts) const{
	
		return std::chrono::duration_cast<std::chrono::milliseconds>(t-ts.t).count();
	}
	
	inline const Timeinterval since(const Timestamp& ts) const{
	
		return Timeinterval( msPassed(ts));
	}
	
	
	inline const string stringTimePassed(const Timestamp& ts) const{
	
		return since(ts).toString();
	}
	
	inline void printTimePassed(const Timestamp& ts) const{
	
		cout << since(ts) << endl;
	}
	
	inline const string toString() const{
	
		std::time_t tt = std::chrono::system_clock::to_time_t(t);
	
		if (name.length() == 0)
			return UtilString::removeLast(std::ctime(&tt), 1);
	
		return name+" "+std::ctime(&tt);
	}
	
	inline const long toLong() const{
	
		return (long) std::chrono::system_clock::to_time_t(t);
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const Timestamp& m) { 
    	
    	return os << m.toString();
	}
};


class Timer{

private:
	string s;

	long tickCounter = 0;
	
	clock_t timeCounter = 0;
	
	//long printCounter = 0;
	clock_t firstPrint = 0;
	clock_t lastPrint = 0;
	
	double secsToWait = 0;
	
	long ticksToWait = 0;
	
	long tickPrecision = 1000;
	
	long ticks = 0;
	long maxTicks = 0;
	
	long timeRes = 0;
	
public:
	
	inline void setMs(long ms){
	
		secsToWait = ms/1000.0;
	}
	
	inline Timer(string str) : s(str){
	
	
	}
	
	inline void setSecs(long secs){
	
		secsToWait = secs;
	}
	
	inline void setTicks(long ticks){
	
		this->ticksToWait = ticks;	
	}
	
	inline void setTimeResolution(long res){
	
		timeRes = res;
	}
	
	inline void msMode(){
	
		setTimeResolution(1);
	}
	
	inline void secMode(){
	
		setTimeResolution(1000);
	}
	
	
	inline void setTickPrecision(long prec){
	
		this->tickPrecision = prec;
	}
	
	long mt1 = 0;
	long mt2 = 0;
	long mt3 = 0;
	
	inline Timer& max(long maxTicks = 0){
	
		mt3 = maxTicks;
		
		return *this;
	}
	
	inline Timer& each(long ticks = 0){
	
		mt2 = ticks;
		
		return *this;
	}
	
	inline Timer& secs(long secs = 0){
	
		mt1 = secs;
		
		return *this;
	}
	
	
	inline void start(long secs=0, long ticks=0, long maxTicks = 0){
	
		if (secs == 0)
			secs = mt1;
	
		if (ticks == 0)
			ticks = mt2;
			
		if (maxTicks == 0)
			maxTicks = mt3;
		
		mt1 = 0;
		mt2 = 0;
		mt3 = 0;
	
		this->maxTicks = maxTicks;
	
		setSecs(secs);
		setTicks(ticks);
		setTickPrecision(ticks);
		
		this->ticks = 0;
		
		firstPrint = clock();
		lastPrint = clock();
		
		print(-1);
		
					
		tickCounter = ticksToWait;
		timeCounter = firstPrint+((long) (CLOCKS_PER_SEC*secsToWait));
	}
	
	inline void end(){
	
		print(1);
	}
	
	
	inline double msPassed(clock_t later, clock_t earlier){
	
		return ((long)(later-earlier))/(CLOCKS_PER_SEC/1000.0);
	}
	
	
	inline string msToString(double ms){
	
		return UtilString::msToString( UtilMath::makeMultipleOf(ms, timeRes) );
	}
	
	inline void print(int mode=0){
	
		clock_t now = clock();
		
		if (mode == -1){
		
			cout << "[BEG]\t ";
		} else if (mode == 0){
		
			cout << "[MID]\t ";
			
		} else if (mode == 1){
		
			cout << "[END]\t ";
		}
		
		cout << s << "\t ";
		
		cout << UtilTime::dateTimeString(UtilTime::dateTimeNow() ) << " ";
		
		if (mode >= 0){
		
		const double msFirst = msPassed(now, firstPrint);
		//const double msLast = msPassed(now, lastPrint);
		
			if (mode == 0){
				cout << "\t step " << ticks;
		
				if (maxTicks > 0)
					cout << " (" << UtilString::doubleToString(ticks*100.0/maxTicks,1) << "%)";
		
				//cout << " [" << msToString( msLast ) << "]";
				
			}
		
		cout << "\t total [" << msToString( msFirst ) << "]";
		
		if (mode == 0)
		if (ticks > 0 && ticks < maxTicks){
		
			cout << "\t est [" << msToString( msFirst/ticks*maxTicks ) << "]";
		}
		
		cout << " " << s << ".secs = " << to_string(msFirst/1000.0);
		}
		
		
		
		cout << endl;
		
		lastPrint = now;
	}

	inline bool tick(){
	
		tickCounter--;
		
		ticks++;
		
		if (tickCounter <= 0){
		
			if ((tickCounter % tickPrecision) == 0){
		
				clock_t now = clock();
			
				if (now >= timeCounter){
					
					print(0);
					
					tickCounter = ticksToWait;
					timeCounter = now+((long) (CLOCKS_PER_SEC*secsToWait));
					
					return true;
				}
			}
			
		}
		
		return false;
	}
};



#endif