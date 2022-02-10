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

#ifndef HEADER_DATASUMMARY_H
#define HEADER_DATASUMMARY_H

#include <approxdata/datasummaries/datasummaries.h>









/*
class Params{

public:

	map<string, double> values;
	
	const string s;
	
	inline Params(const Params& p) : s(p.s){
	
		values = p.values;
	}
	
	
	
	~Params(){}
	
	inline Params() : s(""){
		
		
	}
	
	inline Params(const string s) : s(s){
	
		std::stringstream ss(s);
    	std::string item;
    	std::string item2;
    	
    	while (std::getline(ss, item, '-')) {
        	
        	std::stringstream ss2(item);
        	
        	int k = 0;
        	
        	string cmd = "undefined";
        	double val = 1.0;
        	
        	while (std::getline(ss2, item2, ' ')) {
        	
        		if (k == 0){
        		
        			cmd = item2;
        		}
        		
        		if (k == 1){
        		
        			string clean = item2;
        			
        			clean.erase(std::remove(clean.begin(), clean.end(), ' '), clean.end());
        			clean.erase(std::remove(clean.begin(), clean.end(), '\t'), clean.end());
					
        			try{
		
						std::string::size_type sz;
			 
						const double d = std::stod(clean, &sz);
			 
						 if (sz == clean.size() )
							val = d;
			
					} catch (...){
		
						
					}
        		}
        		
        		k++;
        	}
        	
        	if (UtilString::equals(cmd, "undefined"))
        		;
        	else
        		values[cmd] = val;
    	}
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		for (auto it = values.begin(); it != values.end(); it++){
		
			s << "-" << it->first << " " << it->second << " ";
		}
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const Params& m) { 
    	
    	return os << m.toString();
	}
	
	inline void print() const{
	
		cout << toString() << endl;
	}
	
	inline void set(const string s, double v){
	
		values[s] = v;
	}
	
	inline double get(const string s) const{
	
		if (values.count(s) == 0)
			return 0.0;
		
		return values.at(s);
	}
	
};*/


class SummaryObserver{


	map<string, double> values;
	
	map<string, shared_ptr<Timestamp> > begins;
	map<string, shared_ptr<Timestamp> > ends;
	
public:	

	inline SummaryObserver(){
	
	}

	inline SummaryObserver(const SummaryObserver& obv){
	
		values = obv.values;
		begins = obv.begins;
		ends = obv.ends;
	}
	
	inline long count(const string& s) const{
	
		return values.count(s);
	}	
	
	
	
	inline double& operator[](const string& s){
		
		return values[s];
	}
	
	
	inline void begin(const string& s){
	
		
		shared_ptr<Timestamp> ptr(new Timestamp() );
		begins[s] = ptr;
	}
	
	inline void end(const string& s){
	
		shared_ptr<Timestamp> ptr(new Timestamp() );
		
		assert (begins.count(s) > 0);
		
		ends[s] = ptr;
		
		const string s2 = s+" minutes";
		
		values[s2] = ends[s]->since( (*begins[s]) ).getMinutes();
	}
	
	inline void fromJson(const json& j){
	
		if (j.count("values") > 0){
		
			for (auto& it : j["values"].items() )
				values[it.key()] = it.value();
		}
	}
	
	inline const string toJson(){
	
		vector<string> a1;
		vector<string> v1;
		
		a1.push_back("values");
		
		v1.push_back(UtilJson::toJson(values));
		
		a1.push_back("timepoints");
		
		vector<string> a2;
		vector<string> v2;
		
		if (begins.size() > 0){
			
			for (auto& it : begins){
			
				const string key = it.first;
				Timestamp& t1 = *(begins[key]);
				Timestamp& t2 = *(ends[key]);
			
				vector<string> a3;
				vector<string> v3;
				
				a3.push_back("begin");
				v3.push_back( UtilJson::toJson(t1) );
				
				a3.push_back("end");
				v3.push_back( UtilJson::toJson(t2) );
				
				
				a2.push_back(key);
				v2.push_back(UtilJson::toJson(a3, v3));
			}
		}
		
		v1.push_back(UtilJson::toJson(a2, v2) );
	
		return UtilJson::toJson(a1, v1);
	}
	
};

class DataSummary{



public:

	SummaryObserver observer;
	
	vector<Point>* debugPoints = NULL; 
	
	virtual ~DataSummary(){};
	
	inline virtual double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const = 0;
	
	inline void activateDebug(vector<Point>* v){
	
		debugPoints = v;
	}
	
	inline virtual void drawRandomPoint(RandomGenerator& rg, Point& p) const{}
	
	inline virtual double getBytes() const = 0;
	
	inline virtual void debug(){}
	
	inline virtual const string toString() const = 0;
	
	//inline virtual void summarize(DataSet& data, Params& params) = 0;
	
protected:
    DataSummary(){}
    DataSummary(const DataSummary&){}
    DataSummary& operator=(const DataSummary&){ return *this; }
    
    
};



class UnifOneDimEstimator : public OneDimEstimator{
public:
	inline double intervalCount(double qmin, double qmax, QueryMode mode=QueryMode::EST) const override{
	
		return qmax-qmin;
	}
	
	~UnifOneDimEstimator() override{}
};



#endif
