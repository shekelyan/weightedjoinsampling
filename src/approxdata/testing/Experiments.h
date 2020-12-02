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

/*40979 ?        R    49226:56 bin/approx --data (file data/gaia_5_7_47_50.json) --ranges (file data/queryboxes/gaia4d_{0.01|0.1|1|5|10}.json) --summary 
56755 ?        S      0:00 make vldbgaia4dslow
56776 ?        S      0:00 /bin/sh -c bin/approx "--data (file data/gaia_5_7_47_50.json) --ranges (file data/queryboxes/gaia4d_u{0.01|0.1|1|5|10}.json)
56779 ?        R    37310:43 bin/approx --data (file data/gaia_5_7_47_50.json) --ranges (file data/queryboxes/gaia4d_u{0.01|0.1|1|5|10}.json) --summary
64421 ?        S      0:00 sshd: shekelyan@notty
*/

#ifndef SHEKELYAN_EXPERIMENTS_H
#define SHEKELYAN_EXPERIMENTS_H

#include <approxdata/testing/testing.h>

// import PLY: nohup sed -e 's/[ ]\+/,/g' test.txt > test.csv &

// nohup sed -e 's/[ ]\+/,/g' sg27_station8_intensity_rgb.txt > pointcloud.csv


/*
class ExpRow{

	constvector<string> cells;
};


class ExpTable{

	vector<ExpRow> rows;


	

};*/

namespace Plots{
	
	
	const string EQUIDEPTH = "equidepth";
	const string EQUIDEPTH_STYLE = "color=black,mark=x";
	
	const string DYADIC = "dyadichist";
	const string DYADIC_STYLE = "color=blue,every mark/.append style={solid,fill=blue},mark=diamond*";
	
	const string DIGITHIST = "digithist";
	const string DIGITHIST_STYLE = "color=violet,every mark/.append style={solid,fill=violet},mark=otimes*";
	
	const string SLICEHIST = "slicehist";
	const string SLICEHIST_STYLE = "color=red,every mark/.append style={solid,fill=white},mark=square*";
	
	
	const string RANDOMSAMPLING = "sampling";
	const string RANDOMSAMPLING_STYLE = "color=gray,mark=*";
	
	
	const string INDSLICEHIST = "indslicehist";
	const string INDSLICEHIST_STYLE = "dashed, color=red,every mark/.append style={solid},color=red, mark=square*";
	
	const string RANKHIST = "rankhist";
	const string RANKHIST_STYLE = "color=red,every mark/.append style={solid,fill=red},mark=square*";
};



/*
class CSVImporter{
public:
	vector<double> p;
	
	long size = 0;
	
	shared_ptr<DatWriter> out;	
	
	char sep = ',';
	
	vector<int> sel;
	
	const string path;
	
	int dims = 0;
	
	size_t sz;
	
	inline CSVImporter(const string& path) : path(path){
	
		
	}
	
	inline void close(){
	
		out->close();
		
		cout << "written data file " << path << ".dat" << " with " << size << " points" << endl;
		
		DataDescription desc(dims);
		
		desc.setDataFile(path+".dat");
		
		desc.setSize(size);
		
		Point domainMin(out->min.size() );
		Point domainMax(out->max.size() );
		
		for (int i = 0; i < out->min.size(); i++)
			domainMin[i] = out->min[i];
		
		for (int i = 0; i < out->max.size(); i++)
			domainMax[i] = out->max[i];
		
		Box box(domainMin, domainMax);
		
		desc.setBoundingBox(box);
		
		UtilString::writeToFile( path+".json", desc.toJson() );
		
		cout << "written metadata file " << path << ".json" << endl;
	}
	
	inline bool notSelected(int i) const{
	
		return i >= sel.size() || sel[i] == 0;
	}
	
	
	inline void selectDim(int i){
	
		while (i >= sel.size() )
			sel.push_back(0);
	
		if (sel[i] == 0){
		
			dims++;
			sel[i] = 1;
		}
	}
	
	
	inline void selectDims(const vector<int>& v){
	
		for (int i = 0; i < v.size(); i++)
			selectDim(v[i]);
	}
	
	inline const string readGzip(const string& filename) const {

		std::ifstream stream(filename, std::ios_base::in | std::ios_base::binary);
		
		if (!stream.is_open()){
			throw std::runtime_error("could not open: '" + filename + "'");
		}
		std::string str_compressed((std::istreambuf_iterator<char>(stream.rdbuf())),
								   std::istreambuf_iterator<char>());
		stream.close();

		std::size_t limit = 500 * 1024 * 1024; 
		// file should be about 500 mb uncompressed
		gzip::Decompressor decomp(limit);
		std::string output;
		decomp.decompress(output, str_compressed.data(), str_compressed.size());
	
		return output;
	}
	
	inline void readLine(const string& line){
	
		p.clear();
				
		stringstream symbs(line);

		bool error = false;

		int i = 0;
		
		//cout << "T" ;

		for (std::string symb; std::getline(symbs, symb, sep); i++) {
	
			if (notSelected(i))
				continue;
		
			try {
				double d = std::stod(symb,&sz);
			
				if (d < -10e+12)
					error = true;
			
				if (d > 10e+12)
					error = true;
			
				p.push_back(d);
				
				//cout << symb << ", ";
			
			} catch (std::invalid_argument& e){
			
				p.push_back(NAN);
			
				error = true;
				
				//cout << symb << ", ";
			}
			
			
		}
		
		
		cout << endl;
		
		for (int j = 0; j < p.size(); j++)
			cout << p[j] << ", ";
				
		cout << endl;
		
		if (error)
			cout << line << endl;
		
		
		if (out->writeDirect(p))
			size++;
		
		if (size % 1000000 == 0)
			cout << "written " << size << " points so far ..." << endl;
	}
	
	inline void readFile(const string& file){
	
		if (out){
		
			// nothing to be done
		
		} else {
		
			out.reset(new DatWriter(path+".dat", dims));
		}
		
		if (UtilString::endsWith(file, ".csv.gz")){
		
			const string s = readGzip(file);
			
			stringstream f(s);
			
			for (std::string line; std::getline(f, line); )
				readLine(line);
			
			return;
		}
		
		if (UtilString::endsWith(file, ".csv")){
		
			ifstream f;
			
			f.open(file);
    		
    		assert (f.is_open() );
			
			for (std::string line; std::getline(f, line); )
				readLine(line);
				
			f.close();
			
			return;
		}
		
	}

	
};*/


template<typename E=double>
class SimpleAggregator{

	E min;
	E max;
	E sum;

	long count = 0;
	
	vector<E> vs;
	bool vsSorted = true;

public:
	inline void feed(E e){
	
		if (count == 0|| e < min)
			min = e;	
		
		if (count == 0|| e > max)
			max = e;	
		
		if (count == 0)
			sum = e;
		else
			sum += e;
		
		count++;
		
		vs.push_back(e);
		vsSorted = false;
	}
	
	inline void reset(){
	
		count = 0;
		
		vs.clear();
		vsSorted = true;
	}
	
	
	inline long getCount() const{
	
		return count;
	}
	
	inline E getMin() const{
	
		assert (getCount() > 0);
	
		return min;
	}
	
	inline E getMax() const{
	
		assert (getCount() > 0);
	
		return max;
	}
	
	inline E getSum() const{
	
		assert (getCount() > 0);
	
		return sum;
	}
	
	inline E getAvg() const{
	
		assert (getCount() > 0);
	
		return sum/getCount();
	}
	
	inline void finalise(){
	
		if (!vsSorted){
			
			std::sort(vs.begin(), vs.end() );
			vsSorted = true;
		}
	}
	
	inline E get(int aggregationMode) const{
	
		if (aggregationMode > 1){
		
			
			const long j = (aggregationMode/100.0)*(vs.size()-1);
			
			return vs[j];
			
			/*
			E sum = 0;
			E count = 0;
			
			for (auto it = vs.begin()+j; it != vs.end(); it++){
			
				sum += (*it);
				count++;
			}
			
			return sum/count;*/
		}
	
		if (aggregationMode < 0)
			return getMin();
			
		if (aggregationMode > 0)
			return getMax();
			
		return getAvg();
	}
	
};

inline const string getPgfLatexString(const string s){
	
	stringstream ss;
	
	ss << "\\documentclass{minimal}" << endl;
	
	ss << "\\usepackage{tikz}" << endl;
	
	ss << "\\usepackage{pgfplots}" << endl;
	
	ss << "\\begin{document}" << endl;
	
	ss << s;
	
	
	ss << "\\end{document}" << endl;

	return ss.str();
}


inline int createImage(const DataSet& data, const string file, long xres=3000, long yres=2000, bool flipAxes = false){

	return 0;
	/*
	DataImage img(xres,yres);
	array<float,3> black = {0,0,0};
		
	Point p(2);
	
	double w = 1.0/(data.getSize());

	RandomGenerator rg("abc",0,1);
	
	const double a = 1.0/100;
	
	for (auto it = data.begin(); it != data.end(); it++){
		{		
			
			if (flipAxes){
			
				p[0] = (*it)[1];
				p[1] = (*it)[0];
			} else {
			
				p[0] = (*it)[0];
				p[1] = (*it)[1];
			}
			
			img.drawPoint(p, black, w);
		}	
	}
	
	img.save(file);*/
}

class Plot{

	const string title;
	const string x;
	const string y;

	unordered_map<string, vector<Point>> map;
	
	constvector<string> datasets;
	constvector<string> styles;
	
	PointComparator pc;
	
public:	

	double xmin = numeric_limits<double>::max();
	double xmax = numeric_limits<double>::lowest();
	double ymin = numeric_limits<double>::max();
	double ymax = numeric_limits<double>::lowest();

	vector<Point> allPoints;

	int num = 0;
	
	inline Plot(const string title="title", const string x="x", const string y="y") : title(title), x(x), y(y), pc(0, true){
	
	}
	
	inline void add(const string s, const string style){
	
		datasets.push_back(s);
		styles.push_back(style);
		
		vector<Point> v;
		map[s] = v;
	}
	
	
	inline void sort(){
	
		for (auto it = map.begin(); it != map.end(); it++){
		
			vector<Point>& v = it->second;
					
			std::sort( v.begin(), v.end(), pc);
		}
	}
	
	inline void aggregate(int aggregationMode){
	
		sort();
		
		
	
		for (auto it = map.begin(); it != map.end(); it++){
		
			vector<Point>& v = it->second;
			vector<Point> v2;
			
			
			
			
			double currentX = -1;
			
			SimpleAggregator<double> ag;
			
			for (auto it = v.begin(); it != v.end(); it++){
			
				double x = (*it)[0];
				double y = (*it)[1];
				
				if (x != currentX){
				
					if (ag.getCount() > 0)
						v2.push_back( Point(true, currentX, ag.get(aggregationMode) ) );
					
					currentX = x;
					ag.reset();
				}
				
				ag.feed(y);
			}
			
			ag.finalise();
			
			if (ag.getCount() > 0)
				v2.push_back( Point(true, currentX, ag.get(aggregationMode) ) );
			
			it->second = v2;
		}
	}
	
	inline void add(string s, const Point& p){
	
		num++;
		
		if (p[0] < xmin)
			xmin = p[0];
				
		if (p[0] > xmax)
			xmax = p[0];
		
		if (p[1] < ymin)
			ymin = p[1];	
		
		if (p[1] > ymax)
			ymax = p[1];
		
	
		allPoints.push_back(p);
	
		if (map.count(s) == 0){
		
			vector<Point> v;
		
			map[s] = v;
			
			datasets.push_back(s);
			styles.push_back("");
		}
		
		map[s].push_back(p);
		
		//std::sort( map[s].begin(), map[s].end(), pc);
	}
	
	inline const string getPgfString(const string s = ""){
	
		stringstream ss;
		
		{
		
			ss << endl;
			ss << endl;
			
			ss << "%%%%%%%%%%%%%%%" << x << " -> " << y << " %%%%%%%%%%%%%%%" << endl;
				
			ss << endl << endl;
				
			ss << "\\begin{tikzpicture}" << endl;
			ss << "\\begin{axis}[" << endl;
			
			ss << "title={" << title << "}," << endl;
			ss << "xlabel={" << x << "},"  << endl;
			ss << "ylabel={" << y << "}" << endl;
			ss << s << endl;
			ss << "]" << endl;
			
			
			for (int i = 0; i < datasets.size(); i++){
			
				const vector<Point>& vec = map[datasets[i]];
				
				if (vec.size() == 0)
					continue;
			
				ss << endl;
				ss << "%%% BEGIN " << datasets[i] << endl;
			
				ss << "\\addplot[" << styles[i] << "] coordinates{" << endl;
			
				
			
				for (int j = 0; j < vec.size(); j++){
				
					const Point& p = vec[j];
					ss << "(" << p[0] << ",\t" << p[1] << ")" << endl;
				}
				
				ss << "};" << endl;
				
				ss << "\\addlegendentry{" << datasets[i] << "}" << endl;
				
				ss << "%%% END " << datasets[i] << endl;
				ss << endl;
			}
			
			ss << "\\end{axis}" << endl;
			ss << "\\end{tikzpicture}" << endl;
		}
		
		return ss.str();
	}
	
	
	inline const string getLatexString(){
	
		return getPgfLatexString(getPgfString() );
		
	}
};


namespace UtilErr{


	const int MAX_QERR = 100;


}

class SummaryResult{
public:

	string name = "";
	double lowerbound = -1;
	double upperbound = -1;
	double estimate = -1;
	Box query;
	shared_ptr<Timestamp> timestamp;
	long datasize = -1;
	int datadims = -1;
	
	string source;
	
	bool useCenterEstimates = false;
	
	bool visible = true;
	
	inline void fromJson(const json& j){
	
		lowerbound = j["bounds"][0];
		upperbound = j["bounds"][1];
		estimate = j["estimate"];
		
		if (j.count("source") > 0)
			source = j["source"];
		else
			source = "unknown";
		
		query = UtilJson::getBox( j["query"] );
		timestamp.reset( new Timestamp(UtilJson::getTimestamp( j["timestamp"] ))  );
	}
	
	inline const string toJson(){
	
		vector<string> attributes;
		vector<string> values;
		
		attributes.push_back("estimate");
		
		values.push_back(UtilJson::toJson(estimate) );
		
		attributes.push_back("source");
		
		values.push_back(UtilJson::toJson(source) );
		
		attributes.push_back("bounds");
		
		Point p(2);
		p[0] = lowerbound;
		p[1] = upperbound;
		
		values.push_back(UtilJson::toJson(p) );
		
		attributes.push_back("timestamp");
		
		values.push_back(UtilJson::toJson( (*timestamp) ) );
		
		attributes.push_back("query");
		
		values.push_back(UtilJson::toJson(query) );
		
		return UtilJson::toJson(attributes, values);
	}
	
	inline double get(const string& s){
	
		
		const double T = query.count;
		const double D = datasize;
		const double L = lowerbound;
		const double U = upperbound;
		
		
		
		if (UtilString::containsAny(s, "volgroup")){
		
			const double V = query.getVolume();
			
			return pow(100, round(log(V)/log(100)))*100;
		}
		
		
		if (UtilString::containsAny(s, "densgroup")){
		
			const double V = query.getVolume();
			
			const double ret = (T/D)/V;
			
			return pow(100, round(log(ret)/log(100)));
		}
		
		if (UtilString::containsAny(s, "selgroup")){
		
			const double ret = T/D;
			
			if (ret >= 0.075)
				return 0.1 * 100;
		
			if (ret >= 0.025)
				return 0.05 * 100;
		
			if (ret >= 0.005)
				return 0.01 * 100;
		
			if (ret >= 0.0005)
				return 0.001 * 100;
			
			if (ret >= 0.00005)
				return 0.0001 * 100;
			
			return 0.00001 * 100;
		}
		
		if (UtilString::containsAny(s, "sel")){
		
			const double ret = T/D;
			
			if (true) return ret*100;
			
		}
		
		const bool useCenterEstimates = UtilString::containsAny(s, "center");
		
		const double E = useCenterEstimates ? L*0.5+U*0.5 : estimate;
		
		const bool BIAS = UtilString::containsAny(s, "bias");
		const bool ERROR = UtilString::containsAny(s, "error", "err", "deviation" );
		const bool NORM = UtilString::containsAny(s, "normalised", "normalized", "norm" );
		const bool ABS = UtilString::containsAny(s, "absolute", "abs" );
		const bool BOUNDS = UtilString::containsAny(s, "wid", "width", "bounds");
		const bool REL = UtilString::containsAny(s, "rel", "relative", "mul");
		
		const bool VIOL = UtilString::containsAny(s, "viol");
		
		const bool Q = UtilString::containsAny(s, "q", "qerr", "qwid");
		
		if (T < L)
			cout << name << " T " << T/D << " < " << L/D << " L"<< endl;
			
		if (T > U)
			cout << name << " T " << T/D << " > " << L/D << " U" << endl;
		
		if (VIOL)
			return T > U ? (U-T)/D : T < L ? (L-T)/D : 0;
		
		if (BIAS){
		
			if (NORM)
				return 100*(E-T)/D;
			
			if (ABS)
				return E-T;
		}
		
		if (ERROR){
		
			if (ABS)
				return abs(E-T);
			
			if (NORM)
				return 100*abs(E-T)/D;
			
			if (REL)
				return E > T ? 100*E/T : 100*T/E;
			
			if (Q)
				return UtilMath::minVal<double>(UtilErr::MAX_QERR, E == 0 ? UtilErr::MAX_QERR : UtilMath::maxVal<double>( E/T, T/E));  //T == 0 && E != T ? D/T : E == 0 ? D/T : UtilMath::maxVal<double>( E/T : T/E );
		}
		
		if (BOUNDS){
			
			if (ABS)
				return U-L;
			
			if (NORM)
				return 100*(U-L)/D;
				
			if (REL)
				return 100*U/L;
				
			if (Q)
				return UtilMath::minVal<double>(UtilErr::MAX_QERR, (L == 0 || E == 0) ? UtilErr::MAX_QERR : UtilMath::maxVal<double>( U/L, UtilMath::maxVal<double>( E/T, T/E) )); //UtilMath::maxVal<double>( U/E, E/ ( L < 1 ? 1 : L) );
		}
		
		return -1;
	}
};

/*

	avg|max 

*/

class SummaryResults{
public:
	const string summaryName;

	const string dataName;

	int datadims = 0;
	long datasize = 0;
	
	vector<shared_ptr<SummaryResult>> results;
	
	
	double bytes = 0;
	
	shared_ptr<Timestamp> constructionStart;
	shared_ptr<Timestamp> constructionEnd;
	
	SummaryObserver obs;
	
	bool visible = true;
	
	vector<long> mergeInds;
	
	inline SummaryResults(const SummaryObserver& obv, const string summaryName, const string dataName) : obs(obv), summaryName(summaryName), dataName(dataName){
	
		
	}
	
	inline void merge(SummaryResults& r){
	
		mergeInds.push_back(results.size() );
	
		for (auto it = r.results.begin(); it != r.results.end(); it++){
		
			results.push_back( (*it) );		
		}
	}
	
	inline void filterSelectivity(double selMin, double selMax){
	
		for (auto it = results.begin(); it != results.end(); it++){
		
			const double sel = (*it)->get("sel")/100;
		
			(*it)->visible = (selMin <= sel) && (sel <= selMax);
		}
	}
	
	inline void filter(bool b){
	
		visible = b;
	
		for (auto it = results.begin(); it != results.end(); it++){
		
			(*it)->visible = b;
		}
	}
	
	inline void selectRandom(long m){
	
		std::random_shuffle(results.begin(), results.end() );
	
		long k = 0;
	
		for (auto it = results.begin(); it != results.end(); it++){
			
			(*it)->visible = k < m;
			k++;
		}
	}
	
	inline void fromJson(const json& j){
	
		//summaryName = j["summary"];
		datadims = j["datadims"];
		datasize = j["datasize"];
		bytes = j["bytes"];
		
		constructionStart.reset( new Timestamp(UtilJson::getTimestamp(j["constructionStart"])) );
		constructionEnd.reset( new Timestamp(UtilJson::getTimestamp(j["constructionEnd"])) );
		
		obs.fromJson(j["observer"]);
		
		for (int i = 0; i < j["results"].size(); i++){
		
			shared_ptr<SummaryResult> r(new SummaryResult() );
			
			r->name = dataName+"_"+summaryName;
			r->datadims = datadims;
			r->datasize = datasize;
			r->fromJson(j["results"][i]);
			  
			results.push_back( std::move(r) );
		}
	}
	
	inline const string toJson(){
	
		vector<string> attributes;
		vector<string> values;
		
		
		attributes.push_back("summary");
		values.push_back(UtilJson::toJson(summaryName) );
			
		attributes.push_back("observer");
		values.push_back(obs.toJson() );
			
		attributes.push_back("summarysize");
		values.push_back(UtilJson::toJson(UtilString::bytesToString(bytes) ));
		
		attributes.push_back("constructionTime");
		values.push_back(UtilJson::toJson( constructionEnd->since( (*constructionStart) ).toString() ) );
		
		attributes.push_back("queryMs");
		values.push_back(UtilJson::toJson(UtilJson::toJson( get("avg query milliseconds") )+"ms") );
		
		attributes.push_back("epsilonLowerBound");
		values.push_back(UtilJson::toJson(UtilJson::toJson( get("max norm abs wid") )+"%") );
		
		attributes.push_back("error bias");
		values.push_back(UtilJson::toJson(UtilJson::toJson( get("bias") )) );
		
		attributes.push_back("avgError");
		values.push_back(UtilJson::toJson(UtilJson::toJson( get("avg norm abs err") )+"%") );
		
		attributes.push_back("avgWidth");
		values.push_back(UtilJson::toJson(UtilJson::toJson( get("avg norm abs wid") )+"%") );
		
		
		attributes.push_back("dataname");
		values.push_back(UtilJson::toJson(dataName) );
		
		attributes.push_back("datasize");
		values.push_back(UtilJson::toJson(datasize) );
		
		attributes.push_back("datadims");
		values.push_back(UtilJson::toJson(datadims) );
		
		attributes.push_back("bytes");
		values.push_back(UtilJson::toJson(bytes) );
		
	
		
		attributes.push_back("constructionStart");
		values.push_back(UtilJson::toJson( (*constructionStart) ) );
		
		attributes.push_back("constructionEnd");
		values.push_back(UtilJson::toJson( (*constructionEnd) ) );
		
		attributes.push_back("results");
		
		vector<string> res;
		
		for (auto it = results.begin(); it != results.end(); it++)
			res.push_back( (*it)->toJson() );
		
		values.push_back( UtilJson::toJson(res) );
		
		return UtilJson::toJson(attributes, values);	
	}
	
	inline void reserve(long queryNum){
		
		results.reserve(queryNum);
	}
	
	inline vector<double> getVector(const string s){
	
		vector<double> ret;
	
		for (auto it = results.begin(); it != results.end(); it++){
			
			if (!(*it)->visible)
				continue;
		
			
			const double d = (*it)->get(s);
			
			if (d == -1) 
				ret.push_back( get(s) );
			else
				ret.push_back(d);
		}
	
		return ret;
	}
	
	inline double get(const string& s){
		
		if (UtilString::equalsAny(s, "bytes"))
			return bytes;
	
		if (UtilString::equalsAny(s, "kilobytes", "kb", "kbs"))
			return bytes/1024.0;
	
		if (UtilString::equalsAny(s, "megabytes", "mb", "mbs", "summary size"))
			return bytes/(1024.0*1024.0);
		
		if (UtilString::equalsAny(s, "gigabytes", "gb", "gbs"))
			return bytes/(1024.0*1024.0*1024.0);
	
		if (UtilString::equalsAny(s, "terabytes", "tb", "tbs"))
			return bytes/(1024.0*1024.0*1024.0*1024.0);
		
		
		if (UtilString::equalsAny(s, "construction hours")){
		
			if (obs.count("constructionhours") > 0)
				return constructionEnd->since( (*constructionStart) ).getHours()+obs["constructionhours"];
			else 
				return constructionEnd->since( (*constructionStart) ).getHours();
		}
		
		if (UtilString::equalsAny(s, "construction mins", "construction time")){
		
			if (obs.count("constructionhours") > 0)
				return constructionEnd->since( (*constructionStart) ).getMinutes()+obs["constructionhours"]*60;
			else 
				return constructionEnd->since( (*constructionStart) ).getMinutes();
		
		}
		
		if (UtilString::equalsAny(s, "quantile ops"))
			return obs["quantileInserts"]+obs["quantileLookups"];
			
		if (UtilString::equalsAny(s, "eps", "epsilon", "guarantee"))
			return (obs.count("eps") > 0 ? obs["eps"] : 1)*100;
		
		
		if (UtilString::equalsAny(s, "dimensions", "dimensionality", "dims"))
			return datadims;
			
			
		if (UtilString::equalsAny(s, "data size"))
			return datasize;
	
		if (UtilString::startsWith(s, "p95 ")){
			
			SimpleAggregator<double> ag;
			
			const string s2 = UtilString::afterSplit(s, "p95 ");
			
			for (auto it = results.begin(); it != results.end(); it++){
			
				if (!(*it)->visible)
					continue;
			
				ag.feed( (*it)->get(s2) );
			}
			
			ag.finalise();
			
			return ag.get(95);
		}
		
		if (UtilString::startsWith(s, "p75 ")){
			
			SimpleAggregator<double> ag;
			
			const string s2 = UtilString::afterSplit(s, "p75 ");
			
			for (auto it = results.begin(); it != results.end(); it++){
			
				if (!(*it)->visible)
					continue;
			
				ag.feed( (*it)->get(s2) );
			}
			
			ag.finalise();
			
			return ag.get(75);
		}
	
		if (UtilString::startsWith(s, "avg ")){
		
			const string s2 = UtilString::afterSplit(s, "avg ");
			
			auto end1 = mergeInds.size() > 0 ? results.begin()+mergeInds.at(0) : results.end();
			
			if (UtilString::equalsAny(s2, "query ms", "query time")){
			
				shared_ptr<Timestamp> min;
				shared_ptr<Timestamp> max;
				
				auto it = results.begin();
				
				while (it != end1 && !((*it)->visible))
					it++;
				
				if (it != end1 ){
				
					min = (*it)->timestamp;
					max = (*it)->timestamp;
					
					it++;
				
					for (; it != end1; it++){
				
						if (!(*it)->visible)
							continue;
							
						shared_ptr<Timestamp> t = (*it)->timestamp;
				
						if ( (*t) < (*min))
							min = t;
						
						if ( (*max) < (*t) )
							max = t;
						
					}
				}
				
				const double ret = ( max->since( (*min) ).getMilliseconds() )/(end1-results.begin() );
				
				if (ret < 0.1)
					return 0.099;
			
				return ret;
			}
			
			double ret = 0;
			double size = 0;
			
			for (auto it = results.begin(); it != results.end(); it++){
			
				if (!(*it)->visible)
					continue;
			
				ret += (*it)->get(s2);
				size++;
			}
			
			return ret / size;
		}
		
		if (UtilString::startsWith(s, "max ")){
		
			const string s2 = UtilString::afterSplit(s, "max ");
			
			double ret = UtilMath::largestNegativeDouble();
			
			for (auto it = results.begin(); it != results.end(); it++){
			
				if (!(*it)->visible)
					continue;
			
				ret = UtilMath::maxVal<double>(ret, (*it)->get(s2));
			}
			
			return ret;
		}
		
		if (UtilString::startsWith(s, "min ")){
		
			const string s2 = UtilString::afterSplit(s, "min ");
			
			double ret = UtilMath::largestDouble();
			
			for (auto it = results.begin(); it != results.end(); it++){
			
				if (!(*it)->visible)
					continue;
			
				ret = UtilMath::minVal<double>(ret, (*it)->get(s2));
			}
			
			return ret;
		}
		
		return -1;
	}
	
	inline void add(shared_ptr<SummaryResult> result){
		
		results.push_back(std::move(result));
	}

};


class VLDBPlot : public Plot{

	stringstream axis;
	
	bool isLogLog = false;
	bool yq = false;

		
		

public:
	inline VLDBPlot(const string title, const string x, const string y) : Plot(title, x, y){
	
		add(Plots::EQUIDEPTH, Plots::EQUIDEPTH_STYLE);
		add(Plots::DYADIC, Plots::DYADIC_STYLE);
		add(Plots::DIGITHIST, Plots::DIGITHIST_STYLE);
		add(Plots::SLICEHIST, Plots::SLICEHIST_STYLE);
		
		add(Plots::SLICEHIST+"2", "red, dotted, thick");
		add(Plots::SLICEHIST+"3", "red, dashed, thick");
		add(Plots::SLICEHIST+"4", Plots::SLICEHIST_STYLE);
		add(Plots::SLICEHIST+"5", Plots::SLICEHIST_STYLE);
		
		add(Plots::INDSLICEHIST, Plots::INDSLICEHIST_STYLE);
		add(Plots::RANKHIST, Plots::RANKHIST_STYLE);
		
		add(Plots::RANDOMSAMPLING, Plots::RANDOMSAMPLING_STYLE);
		
		add("sampling (cf = 99\\%)", "dotted");
		
		add("verifiedsampling", "violet, mark=x");
	}
	
	inline VLDBPlot& loglog(){
	
		axis << ",xmode=log, ymode=log";
		
		isLogLog = true;
		
		return (*this);
	}
	
	
	inline VLDBPlot& outerLegend(){
	
		axis << ",legend pos=outer north east";
		
		return (*this);
	}
	
	inline VLDBPlot& xAxisPercent(){
	
		axis << ",xtick={0.0000001,0.00001,0.0001,0.001,0.01,0.1,1.0,10.0,100.0,1000,10000}";
		axis << ",xticklabels={$\\frac{1}{10^6}\\%$,$\\frac{1}{10^5}\\%$,$\\frac{1}{10^4}\\%$,$\\frac{1}{10^3}\\%$,$0.01\\%$,$0.1\\%$,$1\\%$,$10\\%$,$100\\%$,$1000\\%$,$10000\\%$}";
		
		return (*this);
	}
	
	inline VLDBPlot& yAxisPercent(){
	
		axis << ",ytick={0.0000001,0.00001,0.0001,0.001,0.01,0.1,1.0,10.0,100.0,1000,10000}";
		axis << ",yticklabels={$\\frac{1}{10^6}\\%$,$\\frac{1}{10^5}\\%$,$\\frac{1}{10^4}\\%$,$\\frac{1}{10^3}\\%$,$0.01\\%$,$0.1\\%$,$1\\%$,$10\\%$,$100\\%$,$1000\\%$,$10000\\%$}";
		
		return (*this);
	}
	
	inline VLDBPlot& yAxisHours(){
	
		axis << ",ytick={0.0002777778,0.01666667,0.5,1,5,10,24,48,72}";
		axis << ",yticklabels={$1 sec$,$1 min$,$0.5h$,$1h$,$5h$,$10h$,$24h$,$48h$,$3 days$},";
		
		return (*this);
	}
	
	inline VLDBPlot& yAxisQ(){
	
		axis << ",ytick={1.0,1.25, 1.5, 1.75, 2, 3,4,5,10.0,100.0,1000,10000}";
		axis << ",yticklabels={$1$,$1.25$,$1.5$,$1.75$,$2$, $3$,$4$,$5$,$\\ge 10$,$100$,$1000$,$10000$}";
		yq = true;
		
		return (*this);
	}
	
	inline VLDBPlot& xAxisSize(){
	
		axis << ",xtick={0.001,0.01,0.1,1,10,100,1000,10000,100000, 1000000,10000000,100000000,1000000000,1000000000},xticklabels={$KB$,{$_{10}$},{$_{100}$},$MB$,{$_{10}$},{$_{100}$},$GB$,{$_{10}$},{$_{100}$}, $TB$,{$_{10}$},{$_{100}$}, $EB$}";
		return (*this);
	}
	
	inline VLDBPlot& yAxisSize(){

		axis << ",ytick={0.001,0.01,0.1,1,10,100,1000,10000,100000, 1000000,10000000,100000000,1000000000,1000000000},yticklabels={$KB$,{$_{10}$},{$_{100}$},$MB$,{$_{10}$},{$_{100}$},$GB$,{$_{10}$},{$_{100}$}, $TB$,{$_{10}$},{$_{100}$}, $EB$}";
		//axis << ",ytick={1,1000,1000000,1000000000},yticklabels={$MB$,$GB$,$TB$,$EB$}";
		return (*this);
	}
	
	inline VLDBPlot& xMax(long double j){
	
		axis << ",xmax=" << j;
		return (*this);
	}
	
	inline VLDBPlot& scale(long double x, long double y){
	
		axis << ",xscale=" << x;
		axis << ",yscale=" << y;
		return (*this);
	}
	
	inline const string toString(){
	
		if (allPoints.size() == 0)
			return getPgfString(axis.str());
	
		if (isLogLog){
		
			axis << ",xmin=" << pow(10.0, floor(log10(this->xmin)));
			axis << ",xmax=" << pow(10.0, ceil(log10(this->xmax)));
			
			if (!yq)
				axis << ",ymin=" << pow(10.0, floor(log10(this->ymin)));
			
			if (yq){
			
			
				if (allPoints.size() > 0){
					
					axis << ",ymin=1";
				
					double yymax = 0;
				
					for (auto it = allPoints.begin(); it != allPoints.end(); it++){
				
						if ( (*it)[1] <= 10)
							yymax = UtilMath::maxVal(yymax, (*it)[1] );
					}
					
					double newymax = 10;
					
					vector<double> ticks = {1.0,1.25, 1.5, 1.75, 2, 3,4,5,10.0,100.0,1000,10000};
					
					for (int j = ticks.size()-1; j >= 0; j--)
						if ( yymax <= ticks[j])
							newymax = ticks[j];
					
					axis << ",ymax=" << newymax;
				}
			}
			
			if (!yq)	
				axis << ",ymax=" << pow(10.0, ceil(log10(this->ymax)));
			
		} else {
		
		
			axis << ",xmin=" << pow(10.0, floor(log10(this->xmin)));
			axis << ",xmax=" << pow(10.0, ceil(log10(this->xmax)));
			
			axis << ",ymin=" << pow(10.0, floor(log10(this->ymin)));
			axis << ",ymax=" << pow(10.0, ceil(log10(this->ymax)));
		}
		
	
		return getPgfString(axis.str() );
	}
	
	inline void saveStandaloneFile(const string path){
	
		stringstream s;
		
		s << "\\documentclass[class=minimal,border=0pt]{standalone}" << endl;
		
		s << "\\usepackage{tikz}" << endl;
		s << "\\usepackage{pgfplots}" << endl;
		s << "\\begin{document}" << endl;
		
		s << getPgfString(axis.str() );
		
		s << "\\end{document}" << endl;
		
		const string folder = UtilString::getUnixFolder(path);
		const string filename = UtilString::getUnixFileName(path);
		const string ending = UtilString::getFileEnding(path);
		
		UtilString::writeToFile( folder+filename+".tex", s.str() );
		
		if (UtilString::contains(ending, "pdf")){
				
			const string syscall = "pdflatex -output-directory \""+folder+"\" \""+folder+filename+".tex\"";
			
			std::system(syscall.c_str() );
		}
		
	}
	
	inline void saveFile(const string path){
	
		const string folder = UtilString::getUnixFolder(path);
		const string filename = UtilString::getUnixFileName(path);
		const string ending = UtilString::getFileEnding(path);
		
		UtilString::writeToFile( folder+filename+".tex", getPgfString(axis.str()) );
		
		if (UtilString::contains(ending, "pdf")){
				
			const string syscall = "pdflatex -output-directory \""+folder+"\" \""+folder+filename+".tex\"";
			
			std::system(syscall.c_str() );
		}
	}
	
};

class Experiments{

public:	
	vector<shared_ptr<SummaryResults>> results;
	
	
	vector<string> summaries;
	
	inline Experiments(){
	
		
	}
	
	const bool REMOVE_RANKHIST = false;
	
	
	inline void add(const Experiments& e){
	
		const string idAttr = "mb";
	
		for (auto it = e.results.begin(); it != e.results.end(); it++){
		
			bool merged = false;
		
			for (auto it2 = results.begin(); it2 != results.end(); it2++){
			
				if ( UtilString::equals( (*it)->summaryName, (*it2)->summaryName ) && ((*it)->get(idAttr) == (*it2)->get(idAttr) ) ){
				
					(*it2)->merge( (*(*it)) );
					merged = true;
				}
			}
			
			if (!merged)
				results.push_back( (*it) );
		}
	}
	
	/*
	inline bool isLinDominated(const string name, const string xAttr, const string yAttr, double x, double y){
	
	
		vector<Point
	
		for (auto it1 = results.begin(); it1 != results.end(); it1++){
				
			if (!UtilString::contains( (*it1)->summaryName, name))
				continue;
				
			
			
			
			if ( (*it)->get("mb") < 0.1)
					continue;
				
				if ( ((*it)->get("mb") <= (*it2)->get("mb")) && ((*it)->get("eps") <= (*it2)->get("eps")) ){
				
					if ((*it2)->get("mb") < k2Threshold)
						k2Threshold = (*it2)->get("mb");
			}
		}
	}*/
	
	
	inline void filterSelectivity(double selmin, double selmax){
	
	
		for (auto it = results.begin(); it != results.end(); it++){
			
			//const double sel = (*it)->get("sel");
			//(*it)->visible = (selmin <= sel) && (sel <= selmax);
			
			(*it)->filterSelectivity(selmin, selmax);	
		}
	}
	
	inline void selectOnly(const string contains, const string hasToContain){
	
	
		for (auto it = results.begin(); it != results.end(); it++){
			
			if ( UtilString::contains( (*it)->summaryName, contains) )
				(*it)->visible = UtilString::contains( (*it)->summaryName, hasToContain);
		}
	}
	
	inline void selectRandom(long m){
	
	
		for (auto it = results.begin(); it != results.end(); it++){
			
			(*it)->selectRandom(m);	
		}
	}
	
	
	inline shared_ptr<SummaryResults> getSummary(const string name, const string attr, double attrmax) {
	
		shared_ptr<SummaryResults> ret;
		bool first = true;
	
		for (auto it = results.begin(); it != results.end(); it++){
			
			const string s = getName( *(*it) );
			
			if (UtilString::contains(s, name) && (*it)->get(attr) <= attrmax && (first || ((*it)->get(attr)) >= ret->get(attr))){
				ret = (*it);
				first = false;
			}
		}
	
		return ret;
	}
	
	inline shared_ptr<SummaryResults> getSummary2(const string name, const string attr, double attrmin) {
	
		shared_ptr<SummaryResults> ret;
		bool first = true;
	
		for (auto it = results.begin(); it != results.end(); it++){
			
			const string s = getName( *(*it) );
			
			if (UtilString::contains(s, name) && (*it)->get(attr) >= attrmin && (first || ((*it)->get(attr)) <= ret->get(attr))){
				ret = (*it);
				first = false;
			}
		}
	
		return ret;
	}
	
	inline shared_ptr<SummaryResults> getSummaryClosestInSize(const string name, const string attr, double attrmax) {
	
		shared_ptr<SummaryResults> ret;
		bool first = true;
		
		double dist = attrmax;
	
		for (auto it = results.begin(); it != results.end(); it++){
			
			if (!(*it)->visible)
				continue;
			
			
			const string s = getName( *(*it) );
			
			if (UtilString::contains(s, name) ){
				
				const double v = (*it)->get(attr);
				const double dist1 = v > attrmax ? v/attrmax : attrmax/v;
				
				if (first || (dist1 < dist) ){
				
					ret = (*it);
					dist = dist1;
				}
				
				first = false;	
			}
		}
	
		return ret;
	}
	
	inline void filterSize(double mbmax){
		
		unordered_map<string, double> map;
	
		for (auto it = results.begin(); it != results.end(); it++){
			
			const string s = getName( *(*it) );
			
			const double mb = (*it)->get("mb");
			
			if ( mb <= 2*mbmax){
				
				if (map.count(s) == 0 ||mb > map[s])
					map[s] = mb;
			}
		}
		
		for (auto it = map.begin(); it != map.end(); it++){
		
			cout << it->first << " " << it->second << endl;
		}
		
		for (auto it = results.begin(); it != results.end(); it++){
		
			(*it)->visible = (*it)->get("mb") == map[getName( *(*it)) ];
		}
	}
	
	inline void filter(bool b){
	
		for (auto it = results.begin(); it != results.end(); it++)
			(*it)->filter(b);
		
	}
	
	inline void selectMinY(const string a, const string b, const string yAttr){
	
		vector<Point> as;
		vector<Point> bs;
		
		const string xAttr = "mb";
		
		for (auto it = results.begin(); it != results.end(); it++){
			
			if (!(*it)->visible)
				continue;
		
			if (UtilString::contains((*it)->summaryName, a))
				as.push_back (Point(true, (*it)->get(xAttr), (*it)->get(yAttr) ));
				
			if (UtilString::contains((*it)->summaryName, b))
				bs.push_back (Point(true, (*it)->get(xAttr), (*it)->get(yAttr) ));
		}
		
		if (as.size() > 0){
			PointComparator pc(0, true);
		
			std::sort( as.begin(), as.end(), pc);
			std::sort( bs.begin(), bs.end(), pc);
		
			for (auto it = results.begin(); it != results.end(); it++){
		
				if (!(*it)->visible)
					continue;
		
				bool isA = UtilString::contains((*it)->summaryName, a);
				bool isB = UtilString::contains((*it)->summaryName, b);
		
				if ( (!isA) && (!isB) )
					continue;
				
				vector<Point>& vec = isA ? bs : as;
				
				const double x0 = (*it)->get(xAttr);
				const double y0 = (*it)->get(yAttr);
		
				for (long j = 0; j < vec.size(); j++){
				
					const double x1 = vec[j][0];
					const double y1 = vec[j][1];
					
					if ( (j == vec.size()-1) && (x1 <= x0) && (y1 <= y0) ){
						(*it)->visible = false;	
						break;
					}
				
					for (long k = j+1; k < vec.size(); k++){
					
						const double x2 = vec[k][0];
						const double y2 = vec[k][1];
						
						assert (x1 <= x2);
					
						if (y2 > y1)
							continue;
						
						if ( (x1 <= x0) && (x0 <= x2)){
					
							const double y3 = (log(x0)-log(x1))/( log(x2)-log(x1) )*( log(y2)-log(y1) )+log(y1);
						
							if (y3 < log(y0) ){
								(*it)->visible = false;
								break;
							}
						}
					
					}
					
					if (!(*it)->visible)
						break;
					
				}
			}
		}
	}
	
	inline const string getPlotName(const string name) const{
	
		constvector<string> summarynames;
		constvector<string> summarytitles;
		
		summarynames.push_back("rankhist");
		summarytitles.push_back("\\RH");
		
		summarynames.push_back("slicehist");
		summarytitles.push_back("\\SH");
			
		summarynames.push_back("digithist");
		summarytitles.push_back("\\DH");
		
		summarynames.push_back("equidepth");
		summarytitles.push_back("\\ED");
		
		summarynames.push_back("dyadichist");
		summarytitles.push_back("\\DY");
		
		summarynames.push_back("randomsampling");
		summarytitles.push_back("\\RS");
		
		summarynames.push_back("sampling");
		summarytitles.push_back("\\SMP");
		
		for (int j = 0; j < summarynames.size(); j++){
			
			if (UtilString::contains(name, summarynames[j] )){
				
				return summarytitles[j];
			}
		}
		
		return "_";
	}
	
	
	
	inline const string getShortName(const string name) const{
	
		constvector<string> summarynames;
		constvector<string> summarytitles;
		
		summarynames.push_back("rankhist");
		summarytitles.push_back("RH");
		
		summarynames.push_back("slicehist");
		summarytitles.push_back("SH");
			
		summarynames.push_back("digithist");
		summarytitles.push_back("DH");
		
		summarynames.push_back("equidepth");
		summarytitles.push_back("ED");
		
		summarynames.push_back("dyadichist");
		summarytitles.push_back("DY");
		
		summarynames.push_back("randomsampling");
		summarytitles.push_back("RS");
		
		summarynames.push_back("sampling");
		summarytitles.push_back("SMP");
		
		for (int j = 0; j < summarynames.size(); j++){
			
			if (UtilString::contains(name, summarynames[j] )){
				
				return summarytitles[j];
			}
		}
		
		return "_";
	}
	
	
	inline const string getName(const SummaryResults& r) const{
	
		constvector<string> summarynames;
		constvector<string> summarytitles;
		
		summarynames.push_back("rank");
		//summarytitles.push_back(REMOVE_RANKHIST ? "slicehist" : "rankhist");
		summarytitles.push_back("slicehist");
		
		summarynames.push_back("slice");
		summarytitles.push_back("slicehist");
			
		summarynames.push_back("digit");
		summarytitles.push_back("digithist");
		
		summarynames.push_back("dyadic");
		summarytitles.push_back("dyadichist");
		
		summarynames.push_back("equidepth");
		summarytitles.push_back("equidepth");
		
		
		summarynames.push_back("vsampling");
		summarytitles.push_back("sampling");
		
		summarynames.push_back("randomsampling");
		summarytitles.push_back("sampling");
		
		summarynames.push_back("verified");
		summarytitles.push_back("skip");
		
		const string s0 = r.summaryName;
		
		for (int j = 0; j < summarynames.size(); j++){
			
			if (UtilString::contains(s0, summarynames[j] )){
				
				return summarytitles[j];
			}
		}
		
		return "skip";
	}
	
	
	
	inline const string getBestTable(){
	
		constvector<string> names;
		
		names.push_back("slicehist");
		names.push_back("dyadichist");
		names.push_back("sampling");
		names.push_back("equidepth");
		names.push_back("digithist");
		
		vector<double> sizes = {1,10,100, 1000};
		
		constvector<string> qs;
		
		qs.push_back("p75 qerr");
		qs.push_back("p75 qwid");
		
		vector<double> sels = {0.01, 0.1, 1.0, 5.0, 10.0};
		
		stringstream ss;
		
		ss << "%";
		
		for (int q = 0 ; q < qs.size(); q++){
		
			for (int s = 0; s < sels.size(); s++){
				
				for (int b = 1; b <= 1; b++){
				
					ss << "("<< qs[q] <<","<< sels[s] <<","<< b	<<")";
				}
			}
		}
		
		ss << endl;
		
		for (int j = 0 ; j < sizes.size(); j++){
		
			ss << "$" << sizes[j] << "$ MB";
			
			for (int q = 0 ; q < qs.size(); q++){
		
				for (int s = 0; s < sels.size(); s++){
				
					for (int b = 1; b <= 1; b++){
			
						ss << "\t & \t";
					
						double best = -1;
						int bestInd = -1;
						
						double secondBest = -1;
						int secondBestInd = -1;
					
						for (int n = 0; n < names.size(); n++){
					
							bool isSH = (n == 0);
							
							shared_ptr<SummaryResults> res = getSummary(names[n], "mb", isSH ? sizes[j]*1.1 : sizes[j]*1.5);
							
							if (res){
							
								// 0...999 1000...1999 2000...2999
							
								for (long m = 0; m < res->results.size(); m++)
									res->results.at(m)->visible = false;
							
								for (long m = 0; m < 1000; m++)
									res->results.at( (b == 1 ? 0: 5000)+s*1000 )->visible = true;
							
								const double v = res->get(qs[q]);
							
								if (bestInd == -1 || v < best){
						
									bestInd = n;
									best = v;
								
								} else if (secondBestInd == -1 || v < secondBest){
							
									secondBest = v;
									secondBestInd = n;
								}
							
							}		
						}
						
						if (bestInd == 0)
							ss << "\\first{} ";
							
						
						if (bestInd == -1)
							ss << "NONE";
						else
							ss << getShortName(names[bestInd]);
					
						if (secondBestInd != -1)
							ss << " \\better{" << secondBest/best << "}";
					}
				}
			}
			
			ss << "\\\\" << endl;
		}
		
		return ss.str();
		
	}
	
	inline const string getTable() {
	
		constvector<string> names;
		
		names.push_back("slicehist");
		names.push_back("dyadichist");
		names.push_back("sampling");
		names.push_back("equidepth");
		names.push_back("digithist");
		
		vector<double> sizes = {100}; // {1, 10, 100, 1000};
		
		constvector<string> sizenames;
		
		//sizenames.push_back("1MB");
		//sizenames.push_back("10MB");
		sizenames.push_back("100MB");
		//sizenames.push_back("1GB");
		
		constvector<string> met;
		
		met.push_back("mb");
		
		met.push_back("eps");
		met.push_back("construction hours");
		met.push_back("avg query time");
		
		constvector<string> qmet;
		
		//qmet.push_back("avg qerr");
		//qmet.push_back("avg qwid");
		
		qmet.push_back("p75 qerr");
		qmet.push_back("p75 qwid");
		
		
		constvector<string> metnames;
		
		metnames.push_back("summary size");
		
		metnames.push_back("$\\varepsilon$-error guarantee");
		
		metnames.push_back("construction time");
		
		metnames.push_back("avg query time");
		
		constvector<string> qmetnames;
		
		qmetnames.push_back("75th percentile q-error");
		qmetnames.push_back("75th percentile q-bounds");
		
		//if (REMOVE_RANKHIST)
		//	selectMinY("slice", "rank", "eps");
		
		
		vector<double> sels = {0.01, 0.1, 1.0, 5.0, 10.0};
		
		stringstream s;
		
		const bool BODY = false;
		if (BODY){
			s << "\\documentclass{article}" << endl;
			s << "\\begin{document}" << endl;
		}
		
		const bool B = true; // orientation of table
		
		s << "{\\tiny" << endl;
		s << "\\begin{tabular}{l";
		
		for (int k = 0; k < sizes.size(); k++)
			for (int j = 0; j < names.size(); j++)
				s << "l";
		
		s << "}" << endl;
		
		if (B){
		
		
		
			//for (int j = 0; j < sizes.size(); j++)
			//	s << "& \\multicolumn{"<< names.size() << "}{l}{" << sizenames[j] << "}";
			
			s << "\\\\" << endl;
			
			for (int k = 0; k < sizes.size(); k++)
				for (int j = 0; j < names.size(); j++)
					s << "&" << getPlotName( names[j] );
			
			s << "\\\\" << endl;
			
			for (int j = 0; j < met.size(); j++){
		
				s << metnames[j];
			
				for (int l = 0; l < sizes.size(); l++){
				
					const double minSize = getSummary(names[0], "mb", sizes[l])->get("mb");
					
						for (int k = 0; k < names.size(); k++){
							
							shared_ptr<SummaryResults> res = getSummary2(names[k], "mb", minSize);
							
			
						s << "\t&\t";
						
						
			
						if (res){
						
							s << "$";
							
							const double v = res->get(met[j]);
						
							if (UtilString::contains(met[j], "mb")){
							
								s << UtilString::bytesToString( v*1024*1024);
							
						
							} else if (UtilString::contains(met[j], "eps")){
							
								s << UtilString::numToStr( v, 3) << "\\%";
							
							} else if (UtilString::contains(met[j], "query")){
							
								if (v < 0.1){
								
									s << "< 0.1 ms";
									
								} else{
								
									s << UtilString::msToShortString( v );
								}
							
								
								
							} else if (UtilString::contains(met[j], "construction")){
							
								s << UtilString::msToShortString( v*3600*1000 );
							
							} else {
							
								if (v == UtilErr::MAX_QERR)
									s << "\\ge";
								
								if (v > 1 && v < 1.1){
								
									if (v < 1.00001){
									
										s << "< 1.00001";
									} else {
									
										s << UtilString::numToStr( v, 5);
									}
								
								} else {
								
									s << UtilString::numToStr( v, 3);
								}
							}
							
							s << "$";
						}
					}
				}
				
				s << "\\\\" << endl;
			}
			
			
			for (int biased = 1; biased <= 1; biased++)
			for (int sel = 0; sel <= 2; sel++){
			
				for (int j = 0; j < qmet.size(); j++){
		
					s << getTitle(qmet[j]) << " sel " << sels[sel] << " biased " << biased;
					
					s << qmetnames[j];
			
					for (int l = 0; l < sizes.size(); l++){
					
						const double minSize = getSummary(names[0], "mb", sizes[l])->get("mb");
					
						for (int k = 0; k < names.size(); k++){
							
							shared_ptr<SummaryResults> res = getSummary2(names[k], "mb", minSize);
							
							s << "&";
			
							if (res){
							
								
							
								for (long m = 0; m < res->results.size(); m++)
									res->results.at(m)->visible = false;
									
								
								for (long m = 0; m < 1000; m++)
									res->results.at( (biased == 1 ? 0: 5000)+sel*1000 )->visible = true;
									
								
								res->filterSelectivity(sels[sel]*0.5*0.01, sels[sel]*1.5*0.01);
								
								
								const double v = res->get(qmet[j]);
								
								if (v > 1 && v < 1.1){
								
									if (v < 1.00001){
									
										s << "< 1.00001";
									} else {
									
										s << UtilString::numToStr( v, 5);
									}
								
								} else {
								
									s << UtilString::numToStr( v, 1);
								}
								
								//s << UtilString::numToStr( res->get(qmet[j]), 1);
								
								res->filter(false);
							}
						}
					}
			
					s << "\\\\" << endl;
				}
			}

		
		} else {
		
			for (int j = 0; j < names.size(); j++)
				s << "& \\multicolumn{"<< sizes.size() << "}{l}{" << getShortName( names[j] ) << "}";
			
			s << "\\\\" << endl;
		
		
			for (int j = 0; j < names.size(); j++)
				for (int k = 0; k < sizes.size(); k++)
					s << "&" << sizenames[k];
		
			s << "\\\\" << endl;
		
			for (int j = 0; j < met.size(); j++){
		
				s << getTitle(met[j]);
			
				for (int k = 0; k < names.size(); k++){
				
					for (int l = 0; l < sizes.size(); l++){
			
			
						shared_ptr<SummaryResults> res = getSummary(names[k], "mb", sizes[l]*3);
						
						s << "&";						
						
						if (res)
							s << UtilString::numToStr( res->get(met[j]), 1);
					}
				}
			
				s << "\\\\" << endl;
			}
		}
		s << "\\end{tabular}" << endl;
		
		s << "}" << endl;
		
		if (BODY){
			
			
			s << "\\end{document}" << endl;
		}
		
		return s.str();
	
	}
	
	
	inline const string getTitle(const string name) const{
	
		constvector<string> names;
		constvector<string> titles;
		
		names.push_back("eps");
		titles.push_back("{$\\varepsilon$-error guarantee}");
		
		names.push_back("mb");
		titles.push_back("{summary size [MB]}");
		
		names.push_back("selgroup");
		titles.push_back("{selectivity [$\\%$]}");
		
		names.push_back("volgroup");
		titles.push_back("{data volume [$\\%$]}");
		
		names.push_back("densgroup");
		titles.push_back("{data density}");

		names.push_back("avg qerr");
		titles.push_back("{avg q-error}");
		
		names.push_back("avg qwid");
		titles.push_back("{avg upper bound for q-error}");
		
		names.push_back("p95 qerr");
		titles.push_back("{95th percentile q-error}");
		
		names.push_back("p95 qwid");
		titles.push_back("{95th percentile upper bound for q-error}");
		
		names.push_back("p75 qerr");
		titles.push_back("{75th percentile q-error}");
		
		names.push_back("p75 qwid");
		titles.push_back("{75th percentile q-bounds}");
		
		names.push_back("max qerr");
		titles.push_back("{max q-error}");
		
		names.push_back("max qwid");
		titles.push_back("{max upper bound for q-error}");
	
		for (int j = 0; j < names.size(); j++){
		
			if (UtilString::equals(name, names[j]))
				return titles[j];
		}		
		
		return name;
	}
	
	inline shared_ptr<VLDBPlot> newPlot(const string title, const string xAttr, const string yAttr) const{
	
		
		const string x = getTitle(xAttr);
		const string y = getTitle(yAttr);
		
		
		shared_ptr<VLDBPlot> ret( new VLDBPlot( title, x, y ));
		
		if (!UtilString::contains(yAttr, "viol"))
			ret->loglog();
		
		
		bool ERR =  UtilString::contains(xAttr, "err") || UtilString::contains(xAttr, "wid");
		bool XQ = UtilString::contains(xAttr, "q");
		bool YQ = UtilString::contains(yAttr, "q");
		
		if (UtilString::equals(xAttr, "eps") || (ERR & (!XQ) ) ||UtilString::contains(xAttr, "sel") )
			ret->xAxisPercent();
		
		if (UtilString::equals(yAttr, "eps") || (ERR & (!YQ) ) )
			ret->yAxisPercent();
			
		if (UtilString::contains(yAttr, "hours"))
			ret->yAxisHours();
			
		
		if (YQ)
			ret->yAxisQ();
		
		if (UtilString::equals(xAttr, "mb"))
			ret->xAxisSize();
			
		if (UtilString::equals(yAttr, "mb"))
			ret->yAxisSize();
			
		
		
		//if (REMOVE_RANKHIST)
		//	selectMinY("slice", "rank", "eps");
		
		
		bool aggregate = true;
		int aggregationMode = -123;
		
		if (UtilString::contains(yAttr, "p95") ){
			aggregationMode = 95;
			aggregate = true;
		}
		
		if (UtilString::contains(yAttr, "p75") ){
			aggregationMode = 75;
			aggregate = true;
		}
		
		if (UtilString::contains(yAttr, "avg") ){
			aggregationMode = 0;
			aggregate = true;
		}
		
		if (UtilString::contains(yAttr, "max") ){
			aggregationMode = 1;
			aggregate = true;
		}
		
		if (UtilString::contains(yAttr, "min") ){
			aggregationMode = -1;
			aggregate = true;
		}
		
		//if (UtilString::contains(yAttr, "query"))
		//	aggregate = false;
		
		
		
		for (auto it = results.begin(); it != results.end(); it++){
			
			if ( !(*it)->visible)
				continue;
			
			const string s0 = getName( *(*it) );
			
			if (UtilString::contains(s0, "skip"))
				continue;
			
			
			std::regex r0("\\s+.*");
					
			const string s = std::regex_replace(s0, r0, "");
		  
		  	
			if (aggregate){
			
				vector<double> xs = (*it)->getVector(xAttr);
				vector<double> ys = (*it)->getVector(yAttr);
				
				assert (ys.size() > 0);
				
				for (int k = 0; k < ys.size(); k++){
					
					Point p(2);
					
					p[0] = xs[k];
					p[1] = ys[k];
					ret->add(s, p);
				}

			} else {

				Point p(2);
				p[0] = (*it)->get(xAttr);
				p[1] = (*it)->get(yAttr);
				ret->add(s, p);
			}
		}
		
		if (aggregate)
			ret->aggregate(aggregationMode);
		else
			ret->sort();
		
		return ret;
	}
	
	inline shared_ptr<unordered_map<string, vector<Point> >> getPlot( const string xAttr, const string yAttr) const{
	
		shared_ptr<unordered_map<string, vector<Point> >> ret( new unordered_map<string, vector<Point> >());
	
		unordered_map<string, vector<Point>>& map = (*ret);	
		
		for (auto it = results.begin(); it != results.end(); it++){
		
			const string s0 = (*it)->summaryName;
			
			std::regex r0("\\s+.*");
					
			const string s = std::regex_replace(s0, r0, "");
		
			if (map.count(s) == 0){
			
				vector<Point> v;
				map[s] = v;
			}
			
			Point p(2);
			p[0] = (*it)->get(xAttr);
			p[1] = (*it)->get(yAttr);
			
			//cout << p << endl;
			
			map[s].push_back(p);
		}
	
		PointComparator pc(0, true);
		
		for (auto it = map.begin(); it != map.end(); it++)
			std::sort( (it->second).begin(), (it->second).end(), pc);
			
		return ret;
	}
	
	inline double get(const string summaryName, const string s) const{
	
		for (auto it = results.begin(); it != results.end(); it++){
		
			if ( UtilString::startsWith((*it)->summaryName, summaryName) )
				return (*it)->get(s);
		}
		
		return -1;
	}
	
	inline double getLast(const string s) const{
		
		double last = -1;
	
		for (auto it = results.begin(); it != results.end(); it++)
			last = (*it)->get(s);
		
		return last;
	}
	
	inline vector<Point> get(const string summaryName, const string xAttr, const string yAttr) const{
	
		vector<Point> ret;
		
		for (auto it = results.begin(); it != results.end(); it++){
		
			if ( UtilString::startsWith((*it)->summaryName, summaryName) ){
			
				Point p(2);
				p[0] = (*it)->get(xAttr);
				p[1] = (*it)->get(yAttr);
				
				ret.push_back(p);
			}
		}
		
		PointComparator pc(0, true);
		
		std::sort(ret.begin(), ret.end(), pc);
		
		return ret;
	}
	
	inline void fromJson(const json& j){
	
		for (int i = 0; i < j.size(); i++){
		
			const string s = j[i]["summary"];
			const string dn = j[i]["dataname"];
			
			SummaryObserver obs;
			
			shared_ptr<SummaryResults> r(new SummaryResults(obs, s,dn) );
			
			r->fromJson(j[i]);
			
			results.push_back( std::move(r) );
		}
	}
	
	inline const string toJson() const{
	
		vector<string> values;
		
		for (auto it = results.begin(); it != results.end(); it++)
			values.push_back( (*it)->toJson() );
	
		return UtilJson::toJson(values);
	}
};





class Task{

public:
	
	json j;
	
	inline Task(const string s){
	
		parse(fromParamsToJson(s));
		execute();
	}
	
	inline Task(){
	
	}
	
	inline const string fromParamsToJson(const string& s){

		const string sep = ";";

		std::regex r0("["+sep+"]");
		
		const string s2 = std::regex_replace(s, r0, sep+sep);
	
		std::regex r1("([("+sep+"])\\s*([^)"+sep+"\"']+)\\s*([)"+sep+"])");
	
		const string s3 = std::regex_replace(s2, r1, "$1\"$2\"$3");
	
		std::regex r2("["+sep+"]["+sep+"]");
	
		const string s4 = std::regex_replace(s3, r2, ",");
	
		std::regex r3("[(]([^)]+)[)]");
	
		const string s5 = std::regex_replace(s4, r3, "[$1]");
	
		std::regex r4("[-][-]([^\\s-]+)");
	
		const string s6 = std::regex_replace(s5, r4, ",\"$1\":");
	
		if (s6.length() == 0)
			return "{}";
	
		return "{"+s6.substr(1)+"}";
	}
	
	inline void parse(const string& s){
	
		j = UtilJson::parse(s);
		
		constvector<string> vs = variants("{data,ranges,summary}");
		
		for (int i = 0; i < vs.size(); i++){
		
			const string s = vs[i];
			
			if (j.count(s) > 0){
			
				vector<string> v = j[s];
				j[s] = variants2(v);
			}
		}
	}
	
	shared_ptr<SummaryResult> evaluate(const SummaryResults& ret, const DataSet& data, const DataSummary& summary, const Box& b, const string& src){	
			
		shared_ptr<SummaryResult> r(new SummaryResult() );
		
		r->datadims = ret.datadims;
		r->datasize = ret.datasize;
		
		r->source = src;
		
		r->query = b;
		r->timestamp.reset(new Timestamp() );
		
		//try {
		
			r->lowerbound = summary.boxCount(b.min, b.max, QueryMode::LB);
			r->upperbound = summary.boxCount(b.min, b.max, QueryMode::UB);
			r->estimate = summary.boxCount(b.min, b.max, QueryMode::EST);
			
		//} catch (...){
		
		//}
		
		if (r->lowerbound > b.count){
		
			
		}
		
		if (r->upperbound < b.count){
		
			
		}
		
		return r;
	}
	
	
		
	inline void execute(){
	
		const string DATA = "data";
		
		const string SAMPLE = "sample";
		
		const string DIMSEL = "selectdims";
		
		const string RANGES = "ranges";
		
		
		const string PLOTS = "plots";
		
		const string SAVEPLOTS = "saveplots";
		
		const string SAVERESULTS = "saveresults";
		const string LOADRESULTS = "loadresults";
		
		const string FILTER = "filter";
		
		const string DEBUG = "debug";
		
		const string SUMMARY = "summary";
		
		const string SAVE_DATA = "savedata";
		const string SAVE_RANGES = "saveranges";
		const string SAVE_SUMMARY = "savesummaries";
		
		const string VIS = "vis";
		
		const string VIZ = "viz";
		
		const string SAVE_VIS = "savevis";
		
		const string SAVE_VIZ = "saveviz";
		
		const string FILE_PREFIX = "file ";
		const string CSV_ENDING = ".csv";
		
		const string IMPORT = "import";
		
		
		if (j.count("summaries") > 0)
			j[SUMMARY] = j["summaries"];
			
		if (j.count("savesummary") > 0)
			j[SAVE_SUMMARY] = j["savesummary"];
		
		if (j.count("savequeries") > 0)
			j[SAVE_RANGES] = j["savequeries"];
		
		if (j.count("queries") > 0)
			j[RANGES] = j["queries"];
			
		if (j.count("workload") > 0)
			j[RANGES] = j["workload"];
		
		
		
		
		constvector<string> vs;
		
		vs.push_back(DATA);
		vs.push_back(RANGES);
		vs.push_back(SUMMARY);
		vs.push_back(PLOTS);
		vs.push_back(LOADRESULTS);
		vs.push_back(FILTER);
		
		for (int vsi = 0; vsi < vs.size(); vsi++){
			
			const string ss = vs[vsi];
			
			if (j.count(ss) > 0){
			
				vector<string> v = j[ss];
				vector<string> v2 = variants2(v);
				j[ss] = v2;
			}
		}
		
		cout << "json interpretation " << j << endl;
		
		shared_ptr<Experiments> results(new Experiments() );
		
		
		
		if (j.count(IMPORT) > 0){
		
			const string selstr = j[DIMSEL][0];
			
			const IntSeq seq(selstr);
			
			const vector<int>& selectedDims = seq.seq;
			
			DataReader dr;
			
			for (int j1 = 0; j1 < j[IMPORT].size(); j1++){
			
				const string s = j[IMPORT][j1];
				
				constvector<string> v = glob(s);
				
				for (int i = 0; i < v.size(); i++){
				
					const string ss = v[i];
				
					cout << "add import source " << ss << endl;
					dr.addSource(ss, selectedDims);
				}
			}
		
			if (j.count(SAVE_DATA) > 0){
		
				const string jsonfile = j[SAVE_DATA][0];
				
				const string datfile = UtilString::replaceAll(jsonfile, ".json", ".dat");
		
				DatWriter dat(datfile, selectedDims);
				
				long size = 0;
		
				for (vector<double> p; dr.nextPoint(p); ){
				
					if (dat.writeDirect(p))
						size++;
					
					if (size % 1000000 == 0)
						cout << "written " << size << " points so far ..." << endl;
					
				}
			
				dr.close();
				
				dat.close();
			
				cout << "written data file " << datfile << " with " << size << " points" << endl;
				
				DataDescription desc(selectedDims.size());
		
				desc.setDataFile(datfile);
			
				desc.setSize(size);
				
				Point domainMin(dat.min.size() );
				Point domainMax(dat.max.size() );
		
				for (int i = 0; i < dat.min.size(); i++)
					domainMin[i] = dat.min[i];
		
				for (int i = 0; i < dat.max.size(); i++)
					domainMax[i] = dat.max[i];
		
				Box box(domainMin, domainMax);
		
				desc.setBoundingBox(box);
		
				UtilString::writeToFile( jsonfile, desc.toJson() );
			
			}
		
			return;
		}
		
		
		
		
		if (j.count(DATA) > 0){
		
			if (j.count(DIMSEL) == 0){
		
				vector<string> s = {"all"};
		
				j[DIMSEL] = s;
			}
		
			for (int j1 = 0; j1 < j[DATA].size(); j1++){
			
				const string dataString = j[DATA][j1];
		
			
				const bool CSV = UtilString::startsWith(dataString, FILE_PREFIX) && UtilString::endsWith(dataString, CSV_ENDING);
			
				int dims = 0;
				
				if (CSV){
			
					const string import = dataString.substr(FILE_PREFIX.length() );
					CSVReader csv(import);
					dims = csv.getDimNum();
				
				} else {
				
					shared_ptr<DataSet> datatemp = newDataSet(dataString);
					dims = datatemp->getDims();
				}
			
				for (int j5 = 0; j5 < j[DIMSEL].size(); j5++){
				
					string dataString2 = dataString;
					
					DimSelection dimSel;
				
					dimSel.setDataDims("0:"+to_string(dims-1) );
				
					const string s = j[DIMSEL][j5];
				
					cout << "--selectdims (" << s << ")" <<  endl;
					
					if (UtilString::equals(s, "all")){
				
						dimSel.selectDims("0:"+to_string(dims-1));
				
					} else {
				
						dimSel.selectDims(s);
					}
				
					if (CSV){
				
						cout << "import csv " << dataString2 << endl;
			
						const string import = dataString2.substr(FILE_PREFIX.length() );
				
						CSVReader csv(import);
				
						const string dataName = UtilString::getFileName(import)+"_"+dimSel.selectedDims.getFileString();
					
					 	
					
						const string folder = j.count(SAVE_DATA) > 0 ? j[SAVE_DATA][0] : "";
					
						if (UtilString::endsWith(folder, ".json")){
						
							const string path =  UtilString::replaceAll(folder, ".json", "");
						
							cout << "write to " << path << endl;
					
							csv.writeDataFiles(path, dimSel.ind);
					
							dataString2 = "file "+folder;
						
						} else {
						
							const string path = folder+dataName;
					
							cout << "write to " << path << endl;
					
							csv.writeDataFiles(path, dimSel.ind);
					
							dataString2 = "file "+path+".json";
				
						}
					
					
						const string select = dimSel.selectedDims.toString();
					
						dimSel.setDataDims( select );
						dimSel.selectDims( select );
					}
				
				
				
					cout << "--"<< DATA << "(" << dataString2 << ")" <<  endl;
					
					shared_ptr<DataSet> data = newDataSet(dataString2);
					
					//data->selectDims(dimSel.selectedDims.toString() );
					
					if (j.count(SAMPLE) > 0){
						
						const string s = j[SAMPLE][0];
						const Params params(s);
						
						long SAMPLE_POINTS = params.get("num");
						
						data = data->getRandomSubset(SAMPLE_POINTS, "123");
					}
					
					if (j.count(VIZ) > 0){
					
						const string s = j[VIZ][0];
						
						const Params params(s);
						
						long XRES = params.get("xres");
						long YRES = params.get("yres");
						
						array<float,3> black = {0,0,0};
						
						long VIZ_POINTS = params.get("num");
						
						constvector<string> folder;
						
						if (j.count(SAVE_VIZ) > 0){
						
							const string f = j[SAVE_VIZ][0];
							folder.push_back(f);
							
						} else {
						
							folder.push_back("");
						}
						
						const string filename = "_"+data->getFileName()+".png"; //dataName+".png";
						
						{
							
							shared_ptr<DataSet> sample = data->getRandomSubset(VIZ_POINTS);
							createImage( (*sample), folder[0]+"random_"+to_string(VIZ_POINTS)+filename, XRES, YRES, true  );
						}
	
						{
	
							shared_ptr<DataSet> sample = getSample( (*data), VIZ_POINTS);
							createImage( (*sample), folder[0]+"deterministic_"+to_string(VIZ_POINTS)+filename, XRES, YRES, true  );
						}
					}
					
					if (j.count(VIS) > 0){
					
						/*X11
						const string s = j[VIS][0];
						
						const Params params(s);
						
						long XRES = params.get("xres");
						long YRES = params.get("yres");
						
						DataImage img(XRES*4,YRES*4);
						array<float,3> black = {0,0,0};
							
						cout << "--" << VIS << "(" << s << ")" <<  endl;
						
						long VIS_POINTS = params.get("num");
						
						shared_ptr<DataSet> sample = VIS_POINTS == 0 ? data : data->getRandomSubset(VIS_POINTS, "abc");
						
						const long size = sample->getSize();
						
						double w = 1.0/size;
						
						Point p(2);
						
						Timer t("visdata");
						
						t.start(60, 10000, size);
						
						for (auto it2 = sample->begin(); it2 != sample->end(); it2++){
								
							p = (*it2);
							img.drawPoint(p, black, w);
							
							t.tick();
						}
						
						t.end();
						
						img.resize(XRES,YRES);
						
						if (j.count(SAVE_VIS) > 0){
							
							const string folder = j[SAVE_VIS][0];
							
							const string dataName = UtilString::getFileName(dataString2)+"_"+dimSel.selectedDims.getFileString();
							
							const string path = folder+dataName+".png";
							img.save(path);
						}*/
						
					}
				
					if (j.count(SAVE_DATA) > 0){
				
						const string folder = j[SAVE_DATA][0];
				
						cout << "--" << SAVE_DATA << "(" << folder << ")" <<  endl;
				
						if (CSV){
					
							// ALREADY SAVED DATA
					
						} else {
				
							if (UtilString::endsWith(folder, ".json")){
						
								const string jsonfile = folder;
								const string datafile =  UtilString::replaceAll(jsonfile, ".json", ".dat");
								
				
								UtilString::writeToFile( jsonfile, data->getDescription().toJson(datafile) );
								data->writeToFile(datafile);
							} else {
						
								const string datafile = folder+data->getFileName()+".dat";
								const string jsonfile = folder+data->getFileName()+".json";
				
								UtilString::writeToFile( jsonfile, data->getDescription().toJson(datafile) );
								data->writeToFile(datafile);
							}
				
						
						
						}
					}
				
			
					
					vector<shared_ptr<RangeQueries>> workloads;
					
			
					if (j.count(RANGES) > 0){
			
						for (int j2 = 0; j2 < j[RANGES].size(); j2++){
			
							const string rangeString = j[RANGES][j2];
			
							cout << "--" << RANGES << "(" << rangeString << ")" << endl;
							
							shared_ptr<RangeQueries> workload = newWorkload( *data, rangeString );
							
							if (j.count(SAVE_RANGES) > 0){
							
								const string path = j[SAVE_RANGES][0];
								
								cout << "--" << SAVE_RANGES << "(" << path << ")" << endl;
							
								if (UtilString::endsWith(path, ".json") ){
							
									workload->write(path);
								
								} else {
											
									workload->write( path+workload->getName()+".json");
								}
							}
							
							workloads.push_back( std::move(workload) );
						}
					}
					
					if (j.count(SUMMARY) > 0){
					
						unordered_map<double, shared_ptr<DataSummary>> verifyMapK1;
						unordered_map<double, shared_ptr<DataSummary>> verifyMapK2;
						
						for (int j2 = 0; j2 < j[SUMMARY].size(); j2++){
			
							const string summaryString = j[SUMMARY][j2];
			
							Params params(summaryString);
			
							const double vkb = params.get("vkb");
			
							if (vkb > 0 && verifyMapK1.count(vkb) == 0){
			
								Timestamp t1;
			
								verifyMapK1[vkb] = newDataSummary(*data, "slicehist -k 1 -kb "+to_string(vkb) ); 
								
								Timestamp t2;
								
								if (verifyMapK1.count(vkb) > 0)
									verifyMapK1[vkb]->observer["constructionhours"] = t2.since(t1).getHours();
								
								TheoreticalSliceHist tsh1(data->getDims() );
								TheoreticalSliceHist tsh2(data->getDims() );
								
								tsh1.setK(1);
								tsh2.setK(2);
								
								tsh1.setBytes(vkb*1000);
								tsh2.setBytes(vkb*1000);
								
								if (tsh2.getEps() < tsh1.getEps() ){
								
									verifyMapK2[vkb] = newDataSummary(*data, "slicehist -k 2 -kb "+to_string(vkb) );
									
									Timestamp t3;	
									verifyMapK2[vkb]->observer["constructionhours"] = t3.since(t2).getHours();
								}
								
							}
						}
						
						for (int j2 = 0; j2 < j[SUMMARY].size(); j2++){
			
							const string summaryString = j[SUMMARY][j2];
					
							cout << "--" << SUMMARY << "(" << summaryString << ")" << endl;
				
							const Params p(summaryString);
				
							shared_ptr<Timestamp> t1(new Timestamp() );
				
							shared_ptr<DataSummary> summary;
							
							/*try{ 
							*/
								summary = newDataSummary(*data, summaryString);
							/*} catch (...){
							
								cout << "!!! exception during summary creation " << endl;
							}*/
							
							if (RandomSampleSummary* rs = dynamic_cast<RandomSampleSummary*>(summary.get() ) ){
							
								Params params(summaryString);
								
								const double vkb = params.get("vkb");
								
								
								if (vkb > 0){
									assert (verifyMapK1.count(vkb) > 0);
									assert (verifyMapK2.count(vkb) > 0);
								
									if (verifyMapK1.count(vkb) > 0)
									if (SliceHistSummary* sh1 = dynamic_cast<SliceHistSummary*>(verifyMapK1[vkb].get() ) )
										rs->verify( (*sh1) );
								
									if (verifyMapK2.count(vkb) > 0)		
									if (SliceHistSummary* sh2 = dynamic_cast<SliceHistSummary*>(verifyMapK2[vkb].get() ) )
										rs->verify( (*sh2) );
								}
							}
							
							const string summaryName = summary ? summary->toString() : p.toString();
					
							shared_ptr<Timestamp> t2(new Timestamp() );
				
							if (summary){
							
								cout << "created " << summaryName << " [ " << UtilString::bytesToString(summary->getBytes()) << "] in " << UtilString::msToShortString( t2->since((*t1)).getMilliseconds() ) << endl;
							
								if (j.count(VIS) > 0){
								
									/*
					
									const string s = j[VIS][0];
									const string path = s+summaryName+".png";
					
									cout << "--" << VIS << "(" << s << ")" << endl;
									
									const Params params(s);
									
									const long VIS_POINTS = params.get("num");
									
									long XRES = params.get("xres");
									long YRES = params.get("yres");
									
									
									DataImage img(XRES*4,YRES*4);
									array<float,3> black = {0,0,0};
									
									RandomGenerator rg("abc",0,1);
									
									const long size = VIS_POINTS;
									
									double w = 1.0/size;
						
									Point p(2);
									
									Timer t("vissummary");
						
									t.start(60, 10000, size);
									
									for (long j = 0; j < size; j++){
										
										summary->drawRandomPoint(rg, p);
										img.drawPoint(p, black, w);
										
										t.tick();
									}
									
									t.end();
									
									img.resize(XRES,YRES);
									
									if (j.count(SAVE_VIS) > 0){
										
										const string folder = j[SAVE_VIS][0];
										const string path = folder+summaryName+".png";
										img.save(path);
									}
									*/
								}
								
								if (j.count(DEBUG) > 0){
								
									//cout << "DEBUG SUMMARY " << summaryString << endl;
									//summary->debug();
									//cout << "/DEBUG SUMMARY " << summaryString << endl;
								}
								
								
								if (j.count(RANGES) > 0){
					
									shared_ptr<SummaryResults> ret(new SummaryResults(summary->observer, summaryName, data->getName() ) );
				
									ret->constructionStart = t1;
									ret->constructionEnd = t2;
						
									ret->datasize = data->getSize();
									ret->datadims = data->getDims();
						
									ret->bytes = summary->getBytes();
									
									if (j.count(DEBUG) > 0){
										
										
										cout << "DEBUG evaluate summary " << summaryString << endl;
										
										vector<Point> debugPoints;
										
										long jjj = 0;
										
										for (auto itt = workloads.begin(); itt != workloads.end(); itt++){
										
											for (auto it = (*itt)->ranges.begin(); it != (*itt)->ranges.end(); it++){
										
												jjj++;
										
												const Box& b = (*it);
												
												const double LB = summary->boxCount(b.min, b.max, QueryMode::LB);
												const double UB = summary->boxCount(b.min, b.max, QueryMode::UB);
												double T = b.count;	
												
												if ( (LB > T) || (T > UB)){
												
													long sanitycount = 0;
													
													for (auto it2 = data->begin(); it2 != data->end(); it2++){
													
														if (b.contains( (*it2) ))
															sanitycount++;
													}
													
													if (sanitycount != T){
													
														cout << "CORRECTED QUERY COUNT " << T << " != " << sanitycount << " TO ACTUAL COUNT " << endl;
													
													}
													
													T = sanitycount;
													(*it).count = T;
													
												}
												
												if ( (LB > T) || (T > UB)){
												
													cout << "query" << (jjj) << " LB " << LB;
												
													if (LB <= T)
														cout << " <= ";
													else
														cout << " !!!>!!! ";
				
													cout << T;
			
													if (T <= UB)
														cout << " <= ";
													else
														cout << " !!!>!!! ";
													
													cout << UB << " UB" << endl;
													
												}
												
												if (LB > T){
												
													//cout << "QUERY " << b << endl;
													//cout << "LB " << LB << " !!!>!!! " << T << endl;
													
													
													if (debugPoints.size() == 0){
														for (auto it = data->begin(); it != data->end(); it++)
															debugPoints.push_back( (*it) );
													}
													
													summary->activateDebug(&debugPoints);
												
													summary->boxCount(b.min, b.max, QueryMode::LB);
													
													summary->activateDebug(NULL);
												}
												
												if (UB < T){
												
													//cout << "QUERY " << b << endl;
													//cout << "UB " << UB << " !!!<!!! " << T << endl;
													
													if (debugPoints.size() == 0){
														for (auto it = data->begin(); it != data->end(); it++)
															debugPoints.push_back( (*it) );
													}
													
													summary->activateDebug(&debugPoints);
												
													summary->boxCount(b.min, b.max, QueryMode::UB);
													summary->activateDebug(NULL);
												}
												
											}
										}
										
										cout << "/DEBUG evaluate summary " << summaryString << endl;
									}
				
									cout << "evaluate summary " << summaryString << endl;
					
									for (auto itt = workloads.begin(); itt != workloads.end(); itt++){
									
										auto it2 = (*itt)->sources.begin();
										
										
										
										for (auto it = (*itt)->ranges.begin(); it != (*itt)->ranges.end(); it++){
											
											string source;
											
											if (it2 != (*itt)->sources.end()){
											
												source = (*it2);
												it2++;
											} else {
											
												source = "unknown";
											}
											
											ret->add( evaluate( (*ret), (*data), (*summary), (*it), source ) );
											
											
											
											
										}
									}
									
									
									cout << "finished summary evaluation" << endl;
									
									results->results.push_back( std::move(ret) );
									
									if (j.count(PLOTS) > 0){
		
										const string cmd = "";//j.dump();

										for (int i = 0; i < j[PLOTS].size(); i++){
											
											const string plotString = j[PLOTS][i];
		
											//cout << "relation (" << plotString << ")" << endl;
		
											std::regex r0("\\s*[-][>]\\s*");
		
											const string jsonString = "[\""+std::regex_replace(plotString, r0, "\",\"")+"\"]";
		
											cout << jsonString << endl;
		
											json jj = UtilJson::parse(jsonString);
		
											const string xAttr = jj[0];
											const string yAttr = jj[1];
											
											const double X = results->getLast(xAttr);
											const double Y = results->getLast(yAttr);
											
											
											cout << X << " -> " << Y << "( " << xAttr << " -> " << yAttr << ")" << endl;
										}
									}
									
									try {
									
										if (j.count(SAVERESULTS) > 0){
										
											const string resultPath = j[SAVERESULTS][0];
											UtilString::writeToFile( resultPath+".incomplete", results->toJson() );
										}
									
									} catch (...){
									
										cout << "exception during save results" << endl;
									}
									
									try {
									
										if (j.count(PLOTS) > 0){
									
											stringstream plotstream;
		
											const string cmd = "";
				
											for (int i = 0; i < j[PLOTS].size(); i++){
			
												const string plotString = j[PLOTS][i];
		
												std::regex r0("\\s*[-][>]\\s*");
		
												const string jsonString = "[\""+std::regex_replace(plotString, r0, "\",\"")+"\"]";
		
												json jj = UtilJson::parse(jsonString);
		
												const string xAttr = jj[0];
												const string yAttr = jj[1];
											
												shared_ptr<VLDBPlot> plot = results->newPlot(cmd, xAttr, yAttr);
												
												plotstream << plot->toString() << endl << endl;
											}
										
											if (j.count(SAVEPLOTS) > 0){
											
												const string path = j[SAVEPLOTS][0];
												const string folder = UtilString::getUnixFolder(path);
												const string filename = UtilString::getUnixFileName(path);
												const string ending = UtilString::getFileEnding(path);
											
												UtilString::writeToFile(folder+filename+".tex" , getPgfLatexString(plotstream.str()) );
												
											}
										}
										
									} catch (...){
									
										cout << "exception during plots" << endl;
									
									}
									
									
									
								}
								
								
								
							}
						}
					}
				}
			}
		}
		
		if (j.count(LOADRESULTS) > 0){
		
			
		
			for (int i = 0; i < j[LOADRESULTS].size(); i++){
			
				cout << "--" << LOADRESULTS << "(" << j[LOADRESULTS][i] << ")" << endl;
			
				const string s = j[LOADRESULTS][i];
				
				Experiments res;
				
				ifstream in(s);
				
				const bool OPEN = in.is_open();
				
				in.close();
				
				if (!OPEN)
					continue;
				
				try{ 
							
					json j;
				
					UtilJson::readFile(j, s);
				
					res.fromJson(j);
				
					results->add(res);
					
				} catch (...){
							
					cout << "!!! exception during load results " << endl;
				}
				
				
			}
			
			
		}
		
			
		if (j.count(SAVERESULTS) > 0){
		
			cout << "--" << SAVERESULTS << "(" << j[SAVERESULTS][0] << ")" << endl;
			
			UtilString::writeToFile( j[SAVERESULTS][0], results->toJson() );
		}
		
		if (j.count(FILTER) > 0){
		
			const string filter = j[FILTER][0];
		
			Params p(filter);
			
			if (p.get("selectivity") > 0){
			
				const double selmin = p.get("selectivity")*0.5;
				const double selmax = p.get("selectivity")*1.5;
			
				results->filterSelectivity(selmin, selmax);
			}
		
			
		
			if (p.get("mb") > 0){
			
				results->filterSize(p.get("mb"));
			}
			
			if (p.get("sample") > 0){
			
				results->selectRandom(p.get("sample"));
			}
			
			if (p.get("selectverified") > 0){
			
				results->selectOnly("sampl", "vkb");
				
			}
		}
		
		stringstream plotstream;
		
		
		
		if (j.count(PLOTS) > 0){
		
			const string cmd = "";//j.dump();
			
			cout << "create plots ... " << endl;

			for (int i = 0; i < j[PLOTS].size(); i++){
			
				const string plotString = j[PLOTS][i];
		
				cout << "--" << PLOTS << "(" << plotString << ")" << endl;
		
				std::regex r0("\\s*[-][>]\\s*");
		
				const string jsonString = "[\""+std::regex_replace(plotString, r0, "\",\"")+"\"]";
				
				cout << jsonString << endl;
		
				json jj = UtilJson::parse(jsonString);
		
				const string xAttr = jj[0];
				const string yAttr = jj[1];
		
		
		
				shared_ptr<VLDBPlot> plot = results->newPlot(cmd, xAttr, yAttr);
			
				cout << plot->getLatexString() << endl;
				
				if (j.count(SAVEPLOTS) > 0){
				
				
					const string path = j[SAVEPLOTS][0];
					const string folder = UtilString::getUnixFolder(path);
					const string filename = UtilString::getUnixFileName(path);
					const string ending = UtilString::getFileEnding(path);
					
					const string x = UtilString::removeAll(xAttr, ' ');
					const string y = UtilString::removeAll(yAttr, ' ');
					
					const string out = folder+filename+"_"+x+"_"+y+".tex";
					
					plot->saveStandaloneFile(out);
					
					if (UtilString::contains(ending, "pdf")){
					
 						const string syscall = "pdflatex -output-directory \""+folder+"\" \""+out+"\"";
						std::system(syscall.c_str() );
					}
											
					/*
												UtilString::writeToFile(folder+file+".tex" , getPgfLatexString(plotstream.str()) );
				
					stringstream f;
					
					f << "\\documentclass[class=minimal,border=0pt]{standalone}" << endl;
					f << "\\usepackage{tikz}" << endl;
					f << "\\usepackage{pgfplots}" << endl;
					f << "\\begin{document}" << endl;
					f << plot->toString() << endl;								
					f << "\\end{document}" << endl;
					
					const string savePlots = j[SAVEPLOTS][0];
					const string filename = UtilString::replaceAll(savePlots, ".tex", "_"+xAttr+"_"+yAttr+".tex");
													
					UtilString::writeToFile(filename, f.str() );*/
													
				}
				
				plot->outerLegend();
				
				plotstream << plot->toString() << endl << endl;
			
				/*
				auto map = results->getPlot(xAttr, yAttr);
		
				for (auto it1 = map->begin(); it1 != map->end(); it1++){
			
					cout << "plot for " << it1->first << endl;
			
					for (auto itt = it1->second.begin(); itt != it1->second.end(); itt++)
						cout << "p " << (*itt) << endl;
				}*/
			}
		}
		
		if (j.count(SAVEPLOTS) > 0){
		
			cout << "--" << SAVEPLOTS << "(" << j[SAVEPLOTS][0] << ")" << endl;
		
			const string path = j[SAVEPLOTS][0];
			const string folder = UtilString::getUnixFolder(path);
			const string filename = UtilString::getUnixFileName(path);
			const string ending = UtilString::getFileEnding(path);
			
			if (UtilString::contains(path, "table")){
			
				UtilString::writeToFile( folder+filename+"_table.tex", results->getTable() );
			
				UtilString::writeToFile( folder+filename+"_besttable.tex", results->getBestTable() );
			
			} else {
			
				UtilString::writeToFile( folder+filename+".tex", getPgfLatexString(plotstream.str()) );
			
				if (UtilString::contains(ending, "pdf")){
					
					const string syscall = "pdflatex -output-directory \""+folder+"\" \""+folder+filename+".tex\"" ;
					std::system(syscall.c_str() );
				}
			
			}
		}
		
	}
};


#endif