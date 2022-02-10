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

#ifndef HEADER_UTILJSON_H
#define HEADER_UTILJSON_H


#include <approxdata/external/json.h>

using json = nlohmann::json;



namespace UtilJson{

	inline json parse(const string& s){
	
		json ret = json::parse(s.begin(), s.end());
		
		return ret;
	
	}

	inline void readFile(json& j, const string &file){
		
		ifstream in(file);
		
		assert (in.is_open() );
		
		in >> j;
		
		in.close();
	}
	
	inline vector<double> getArray(const json& j){
	
		vector<double> ret(j.size() );
		
		for (int i = 0; i < ret.size(); i++)
			ret[i] = j[i];
		
		return ret;
	}
	
	inline Point getPoint(const json& j){
		
		Point p(j.size() );
		
		for (int i = 0; i < p.size(); i++)
			p[i] = j[i];
		
		return p;
	}

	inline Box getBox(const json& j){
			
		const Box b(getPoint(j["min"]), getPoint(j["max"]), j["count"]);
			
		return b;
	}
	
		
	
	inline Timestamp getTimestamp(const json& j){
	
		const long unixtime = j["unixtime"];
			
		Timestamp t(unixtime);
			
		return t;
	}
	
	
	inline const string toJson(const string& s){
	
		return "\""+s+"\"";
	}
	
	inline const string toJson(double d){
	
		return UtilString::doubleToString(d);
	}
	
	inline const string toJson(const Point& p){
	
		stringstream out;
		
		out << p;	
	
		return out.str();
	}
	
	inline const string toJson(const vector<string>& values){
	
		stringstream out;
	
		out << "[ ";
		
		for (long j = 0; j < values.size(); j++){
		
			if (j > 0)
				out << ", ";
		
			out << values[j];
		}
		
		out << " ]" << endl;
		
		return out.str();
	}
	
	inline const string toJson(const vector<string>& attributes, const vector<string>& values){
	
		stringstream out;
	
		out << "{ ";
		
		for (long j = 0; j < values.size(); j++){
		
			if (j > 0)
				out << ", ";
		
			out << toJson(attributes[j]) << ": " << values[j];
		}
		
		out << " }" << endl;
		
		return out.str();
	}
	
	inline const string toJson(const Box& b){
	
		vector<string> attributes;
		vector<string> values;
		
		attributes.push_back("count");
		values.push_back(toJson(b.count));
	
		attributes.push_back("min");
		values.push_back(toJson(b.min));
		
		attributes.push_back("max");
		values.push_back(toJson(b.max));
		
		return toJson(attributes, values);
	}
	
	
	inline const string toJson(const map<string, double>& m){
	
		vector<string> attributes;
		vector<string> values;
	
		for (auto& it : m){
		
			attributes.push_back(it.first);
			values.push_back(toJson(it.second));
		}
		
		return toJson(attributes, values);
	}
	
	inline const string toJson(const Timestamp& t){
	
		vector<string> attributes;
		vector<string> values;
		
		attributes.push_back("date");
		values.push_back(toJson(t.toString()));
	
		attributes.push_back("unixtime");
		values.push_back(toJson(t.toLong() ));
		
		return toJson(attributes, values);
	}
	
	inline const string paramsToJsonStringOld(const string s0){


		const string s1 = s0+" --";
		
		std::regex r1("[-][-]([^\"'\\s-]+)[\\s]+[\"']([^\"']+)[\"'][\\s]+[-][-]");
		
		const string s2 = std::regex_replace(s1, r1, "--\"$1\":[\"$2\"] --");
		
		std::regex r2("[-][-]([^\"'\\s-]+)[\\s]+([^-]+)\\s+[-][-]");
		
		const string s3 = std::regex_replace(s2, r2, "--\"$1\":$2 --");
		
		std::regex r3("[-][-]");
		
		const string s4 = std::regex_replace(s3, r3, ",");
		
		cout << "s1 " << s1 << endl;
		cout << "s2 " << s2 << endl;
		cout << "s3 " << s3 << endl;
		
		return "{"+s4.substr(1, s4.length()-2)+"}";
	}
	
	inline const string paramsToJsonString(const string s0){


		


		const string s1 = s0+" ";
		
		cout << s1 << endl;
		
		std::regex r1("[=]");
		const string s2 = std::regex_replace(s1, r1, " ");
		
		cout << s2 << endl;
		
		std::regex r10("[-][-]([^#\"'\\s-]+)[\\s]+");
		const string s20 = std::regex_replace(s2, r10, "-$1 ");
		
		std::regex r2("[-]([^#\"'\\s-]+)[\\s]+[\"']([^\"']+)[\"'][\\s]+");
		
		cout << s20 << endl;
		
		const string s3 = std::regex_replace(s20, r2, "#\"$1\":[\"$2\"] ");
		
		std::regex r3("[-]([^#\\s]+)[\\s]+([0-9]+)[\\s]+");
		
		cout << s3 << endl;
		
		const string s4 = std::regex_replace(s3, r3, "#\"$1\":$2 ");
		
		std::regex r4("[-]([^#\"'\\s-]+)[\\s]+([^-\\s]+)[\\s]+");
		
		cout << s4 << endl;
		
		const string s5 = std::regex_replace(s4, r4, "#\"$1\":\"$2\" ");
		
		cout << s5 << endl;
		
		std::regex r5("[#]");
		
		const string s6 = std::regex_replace(s5, r5, ",");
		
		return "{"+s6.substr(1, s6.length())+"}";
	}
};



inline const string translateUnits(const string s_, const string p1, double mul){

	std::stringstream ret;

	std::smatch m;
	std::regex e("([-+]?[0-9]*[.]?[0-9]+(?:[eE][-+]?[0-9]+)?)"+p1);

	string s = s_;

	bool notFound = true;
	
	while (std::regex_search (s,m,e)) {

		nlohmann::json j = nlohmann::json::parse("["+m[1].str()+"]");
		
		const double d = j[0];
		
		ret << m.prefix().str();
		ret << UtilString::doubleToString(d*mul);
		
		s = m.suffix().str();
		
		notFound = false;
	}
	
	ret << s;
	
	return ret.str();
}


inline const string translateUnits(const string s0){

	
	const string s1 = translateUnits(s0, "\\s*[kK][bB]", 1);
	const string s2 = translateUnits(s1, "\\s*[mM][bB]", 1000);
	const string s3 = translateUnits(s2, "\\s*[gG][bB]", 1000*1000);
	
	const string s5 = translateUnits(s3, "\\s*million", 1000*1000);
	const string s6 = translateUnits(s5, "\\s*billion", 1000*1000*1000);
	
	const string s4 = translateUnits(s6, "\\s*[%]", 0.01);
	
	cout << s0 << " translated to " << s4 << endl;
	
	return s4;
}

inline constvector<string> variants(const string& s_){

	std::smatch m;
	
	string s = s_;
	
	const string sep = "|";
	
	std::regex e("[{]([^"+sep+"{}]+(?:["+sep+"][^"+sep+"{}]+)*)[}]");
	
	constvector<string> ret;
	
	ret.push_back("");
	
	while (std::regex_search (s,m,e)) {
		
		{
		constvector<string> ret2;
	  	for (int i = 0; i < ret.size(); i++)
			ret2.push_back( ret[i]+m.prefix().str() );
			
		ret.clear();
			
		for (int i = 0; i < ret2.size(); i++)
			ret.push_back( ret2[i] );
		}
			
		string ss = m[1].str();
		std::regex ee("([^"+sep+"]+)");
		std::smatch mm;
		
		constvector<string> ret2;
		
		while (std::regex_search (ss,mm,ee)) {

				
			for (int i = 0; i < ret.size(); i++)
				ret2.push_back( ret[i]+mm[1].str() );
			
			ss = mm.suffix().str();
		}

		{
			ret.clear();
			
			for (int i = 0; i < ret2.size(); i++)
				ret.push_back( ret2[i] );
		}

		//cout << m.prefix().str() << " " << m[1] << " " << m.suffix().str() << endl;
		
		//cout << m[1] << endl;
		
    	s = m.suffix().str();
  	}
  	
  	constvector<string> ret2;
  	for (int i = 0; i < ret.size(); i++)
		ret2.push_back( translateUnits(ret[i]+s));
	 
	return ret2;
}

/*
inline vector<string> variants2(const vector<string>& v){
	
	vector<string> ret;
	
	for (auto it = v.begin(); it != v.end(); it++){
	
		vector<const string> s = variants( (*it) );
	
		for (auto it2 = s.begin(); it2 != s.end(); it2++){
		
			ret.push_back( (*it2) );
		}
	}
	
	return ret;
}*/

inline vector<string> variants2(const vector<string>& v){
	
	vector<string> ret;
	
	for (auto it = v.begin(); it != v.end(); it++){
	
		constvector<string> s = variants( (*it) );
	
		for (int i = 0; i < s.size(); i++){
		
			ret.push_back( s[i] );
		}
	}
	
	return ret;
}




inline const nlohmann::json parseParamsToJson(const string& s_){



	const string s = s_+" -";

	std::regex r0("[-]([^\\s-]+)\\s+[-]");
	
	std::regex r1("[-]([a-z]+)[ ]+([-+]?[0-9]*[.]?[0-9]+(?:[eE][-+]?[0-9]+)?)");
	std::regex r2("[-]([^\\s-]+)\\s+([^\\s-]+)");
	
	//cout << s << endl;
	
	const string s0 = std::regex_replace(s, r0, "\"$1\":1, -");
	
	//cout << s0 << endl;
	
	const string s1 = std::regex_replace(s0, r1, "\"$1\":$2,");
	
	//cout << s1 << endl;
	
	const string s2 = std::regex_replace(s1, r2, "\"$1\":\"$2\",");
	
	//cout << s2 << endl;
	
	const string s3 = "{"+s2.substr(0,s2.length()-3)+"}";
	
	//cout << s3 << endl;
	
	cout << "parseParamsToJson " << s3 << endl;

	
	return nlohmann::json::parse(s3);
}

class Params{

public:

	map<string, double> values;
	
	const string name;
	const string s;
	
	
	
	inline Params(const Params& p) : s(p.s){
	
		values = p.values;
	}
	
	
	
	~Params(){}
	
	inline Params() : s(""){
		
		
	}
	
	
	
	
	inline Params(bool b, const nlohmann::json j) : name(j.value("name", "")), s(""){
		
		for (auto it = j.begin(); it != j.end(); it++){
		
			if (!UtilString::equals( it.key(), "name"))
				values[it.key()] = it.value();
		}
	}
	
	
	inline Params(const string& s) : Params(true, parseParamsToJson("-name "+s) ){
		
		
	}
	
	
	inline const string toString() const{
	
		stringstream s;
		
		if (name.length() > 0)
			s << name << " ";
		
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
	
	inline void set(const string& s, double v){
	
		values[s] = v;
	}
	
	inline double get(const string& s) const{
	
		if (values.count(s) == 0)
			return 0.0;
		
		return values.at(s);
	}
	
};


#endif
