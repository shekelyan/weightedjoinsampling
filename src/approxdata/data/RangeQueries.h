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

#ifndef HEADER_RANGEQUERIES_H
#define HEADER_RANGEQUERIES_H

#include <approxdata/data/data.h>

class RangeQueries{

public:
	vector<Box> ranges;
	
	vector<string> sources;
	
	shared_ptr<const string> name;
	
	inline RangeQueries(){
	
		
	}
	
	inline RangeQueries(const string file){
		
		add(file);
	}
	
	inline void setName(const string s){
	
		name = std::make_shared<const string>(s);
	}
	
	
	inline void addFromString(const string s){
	
		stringstream ss(s);
		
		nlohmann::json j;
		
		ss >> j;
		
		const long size = j["queries"].size();
		
		if (j.count("name") > 0){
		
			const string n = j["name"];
			setName(n);
		}
		
		for (long k = 0; k < size; k++){
		
			long count = j["queries"][k]["count"];
			
			Point min(j["queries"][k]["min"].size() );
			
			for (int i = 0; i < min.size(); i++)
				min[i] = j["queries"][k]["min"][i];
				
			Point max(j["queries"][k]["max"].size() );
			
			for (int i = 0; i < max.size(); i++)
				max[i] = j["queries"][k]["max"][i];
			
			if (j["queries"][k].count("source") > 0){
			
				string src = j["queries"][k]["source"];
				
				sources.push_back(src);
			}
			
			Box b(min, max, count);
			
			ranges.push_back(b);
			
		}
	}
	
	
	inline void discretize(){
	
		const string s = toString();
		
		ranges.clear();
		
		addFromString(s);
	}
	
	inline const string getName(){
	
		if (name)
			return (*name);
			
		return "ranges";
	}
	
	// read binary
	inline void addLegacy(const string file){
	
		BinaryReader r(file);
	
		const long size = r.readLong();
		const int DIMS = r.readInt();
		
		assert (size > 0);
		assert (size < 1000000000);
		
		for (long j = 0; j < size; j++){
			
			Point qmin(DIMS);
			Point qmax(DIMS);
			
			for (int i = 0; i < DIMS; i++)
				qmin[i] = r.readDouble();
			
			for (int i = 0; i < DIMS; i++)
				qmax[i] = r.readDouble();
			
			Box b(qmin, qmax, r.readDouble() );
			
			ranges.push_back(b);
			
			sources.push_back(file);
		}
		
		r.close();
	}
	
	inline void print() const{
	
		cout << "number of ranges " << ranges.size() << endl;
		
		long k = 0;
		for (auto it = ranges.begin(); it != ranges.end(); it++){
		
			if (k < 5 || k > (ranges.size()-10)){
			
				it->print();
			}
			k++;
		}
	}
	
	
	inline void add(ifstream& in, const string source){
	
		nlohmann::json j;
		
		in >> j;
		
		const long size = j["queries"].size();
		
		if (j.count("name") > 0){
		
			const string n = j["name"];
			setName(n);
		}
		
		for (long k = 0; k < size; k++){
		
			long count = j["queries"][k]["count"];
			
			Point min(j["queries"][k]["min"].size() );
			
			for (int i = 0; i < min.size(); i++)
				min[i] = j["queries"][k]["min"][i];
				
			Point max(j["queries"][k]["max"].size() );
			
			for (int i = 0; i < max.size(); i++)
				max[i] = j["queries"][k]["max"][i];
				
			Box b(min, max, count);
			
			if (j["queries"][k].count("source") > 0){
			
				string src = j["queries"][k]["source"];
				
				sources.push_back(src);
			} else {
			
				sources.push_back(source);
			}
			
			
			ranges.push_back(b);
			
		}
	
	}
	
	// read JSON file
	inline void add(const string file){
	
		ifstream in(file);
		
		add(in, file);
		
		in.close();
	}
	
	inline const string toString() const {
	
		stringstream out;
	
		out << "{" << endl;
		
		if (name)
			out << "\"name\":\""+(*name)+"\"," << endl << endl;
		
		out << "\"queries\":[" << endl << endl;
		
		for (long j = 0; j < ranges.size(); j++){
		
			if (j > 0)
				out << ",";
		
			out << "{" << endl;
			
			
			Point pmin = ranges[j].min;
			Point pmax = ranges[j].max;
			
			pmin.id = -1;
			pmax.id = -1;

			if (sources.size() > 0)
				out << "\"source\":\"" << sources.at(j) << "\"," << endl;

			out << "\"count\":" << UtilString::doubleToString(ranges[j].count) << "," << endl;
			out << "\"min\": " << pmin << "," << endl;
			out << "\"max\": " << pmax << "" << endl;
			
			
			out << "}";
		}
		
		out << endl << endl << "]" << endl;
		out << "}" << endl;
		
		return out.str();
	}
	
	inline void write(ofstream& out) const{
		
		out << toString();
	}
	
	// write JSON file
	inline void write(string file) const{
	
		ofstream out(file);
		
		write(out);
		
		out.close();
	}
	
	
	inline void print(long size) const{
	
		
		long k = 0;
		for (auto it = ranges.begin(); it != ranges.end(); it++){
		
			if (k < 5 || k > (ranges.size()-10)){
			
				it->print(size);
			}
			k++;
		}
	}
	
	
	const Box& operator[](std::size_t idx) const{
		
		return ranges[idx];
	}
	
	
	inline void add(const Box& cb){
	
		ranges.push_back(cb);
	}
	
	
	inline void add(shared_ptr<vector<Box>> b){
	
		add( (*b) );
	}
	
	inline void setBox(long k, Box& b){
	
		ranges[k].setMin(b.min);
		ranges[k].setMax(b.max);
	}
	
	inline void add(vector<Box>& b){
	
		for (auto it = b.begin(); it != b.end(); it++){
		
			add( (*it) );
		}
	}
	
	inline void add(vector<string>& b){
	
		for (auto it = b.begin(); it != b.end(); it++){
		
			sources.push_back( (*it) );
		}
	}
	
	
	inline void add(RangeQueries& rq){
	
		add(rq.ranges);
		add(rq.sources);
	}
	
	inline void add(shared_ptr<RangeQueries> rq){
	
		add( (*rq) );
	}
	
	
	inline long size() const{
	
		return ranges.size();
	}
	
	inline long getSize() const{
	
		return ranges.size();
	}
	
	inline void setCount(long ind, long count){
	
		ranges[ind].count = count;
	}
	
	
	inline shared_ptr<RangeQueries> getWithCounts(double min, double max) const{
	
		shared_ptr<RangeQueries> ret(new RangeQueries() );
		
		for (auto it = ranges.begin(); it != ranges.end(); it++)
			if ( it->count >= min && it->count <= max)
				ret->add( (*it) );
		
		return ret;
	}
	
	
	inline shared_ptr<RangeQueries> getRandomSubset(long s, const string seed="123") const{
		
		if (s >= size() ){
		
			shared_ptr<RangeQueries> ret(new RangeQueries() );
			
			for (auto it = ranges.begin(); it != ranges.end(); it++)
				ret->add( (*it) );
			
			return ret;
		}
		
		vector<long> inds(s);
		
		ReservoirSampling rs(seed, s);
		
		for (long pos = 0; pos < size(); pos++){
		
			const long j = rs.pos();
			
			if ( (j >= 0) && (j < s))
				inds[j] = pos;
		}
		
		shared_ptr<RangeQueries> ret(new RangeQueries() );
		
		for (auto it = inds.begin(); it != inds.end(); it++)
			ret->add( ranges[(*it)] );
		
		return ret;
	}
};






#endif
