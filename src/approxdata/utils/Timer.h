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

#ifndef HEADER_TIMER_H
#define HEADER_TIMER_H

#include <approxdata/utils/utils.h>

class Timer{

private:
	const string s;

	long tickCounter = 0;
	
	clock_t timeCounter = 0;
	
	long printCounter = 0;
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
	
	inline Timer(const string& str) : s(str){
	
	
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
	
	inline Timer& ticks(long ticks = 0){
	
		mt2 = ticks;
		
		return *this;
	}
	
	inline Timer& secs(long secs = 0){
	
		mt1 = setSecs(secs);
		
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
	
	
		if (ticks <= 0)
			ticks = 1;
		
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
		
		cout << UtilString::dateTimeString(UtilTime::dateTimeNow() ) << " ";
		
		if (mode >= 0){
		
			const double msFirst = msPassed(now, firstPrint);
			const double msLast = msPassed(now, lastPrint);
		
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
		}
		
		
		
		cout << " " << s << ".secs = " << msFirst/1000.0;
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
