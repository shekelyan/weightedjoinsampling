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

#ifndef HEADER_DATASET_H
#define HEADER_DATASET_H

#include <approxdata/data/data.h>




// meta information about dataset to help read it etc
class DataDescription{

	private:
		vector<double> domainAdd;
		vector<double> domainMul;
		vector<double> domainDiv;
		
		vector<string> name = {"noname"};
		vector<string> dataFile = {"nofile"};
		
		vector<long> domainDistinct;
		vector<string> types;
		
		long size = 0;
		
		int dims = 0;
		
		Box boundingBox;
		

	public:
		
		inline DataDescription(int dims, string type="real"): dims(dims){
			
			setName("noname");
		}
		
		inline long getSize() const{
		
			return size;
		}
		
		inline void setSize(long size){
		
			this->size = size;
		}
		
		inline int getDims() const{
		
			return dims;
		}
		
		inline void setName(const string name){
		
			this->name.clear();
			this->name.push_back(name);
		}
		
		inline void setDataFile(const string file){
		
			this->dataFile.clear();
			this->dataFile.push_back(file);
		}
		
		inline DataDescription(const string file){
			
			
			if (UtilString::endsWith(file, ".json")){
				readJSON(file);
				return;
			}
			
			
			cout << file << endl;
			
			assert (false);
			/*
			setName("data");
			Point domainMin;
			Point domainMax;
			
			if (UtilString::endsWith(file, "ini")){
			
				IniReader r(file);
				dims = r.getInt("description", "dims");
			
				assert (dims > 0);
			
				size = r.getInt("description", "size");
				setDataFile( r.getString("description", "file") );
				
				
				for (int i = 0; i < dims; i++){
			
					domainMin.push_back(r.getDouble("boundaries", "min"+to_string(i)));
					domainMax.push_back(r.getDouble("boundaries", "max"+to_string(i)));
				
					domainDistinct.push_back(r.getInt("boundaries", "distinct"+to_string(i)));
				}
			}
			
			
			Box box;
			
			box.setMin(domainMin);
			box.setMax(domainMax);
			
			setBoundingBox(box);*/
		}
		
		
		
		inline void readJSON(const json& j){

			dims = j["data"][0]["dimensionality"];
			size = j["data"][0]["size"];
	
			setName(j["name"]);
	
			setDataFile( j["data"][0]["file"] );
	
			Point domainMin(dims);
			Point domainMax(dims);
	
			for (int i = 0; i < dims; i++){
	
				domainMin[i] = j["data"][0]["dimensions"][i]["boundaries"][0];
				domainMax[i] = j["data"][0]["dimensions"][i]["boundaries"][1];
		
				domainDistinct.push_back(numeric_limits<long>::max() );
			}
			
			setBoundingBox(Box(domainMin, domainMax, -1));
		}
		
		inline void readJSON(const string& file){
			
			json j;
			
			UtilJson::readFile(j, file);
			
			readJSON(j);
		}
		
		inline const string toJson(const string datafile="") const{
		
			vector<string> attr1;
			vector<string> vals1;
			
			{
	
				attr1.push_back("name");
				vals1.push_back(UtilJson::toJson(name[0]) );
					
	
				vector<string> attr2;
				vector<string> vals2;
				
				{
		
					attr2.push_back("size");
					vals2.push_back(UtilJson::toJson(size) );
		
					attr2.push_back("dimensionality");
					vals2.push_back(UtilJson::toJson(dims) );
		
					
					attr2.push_back("file");
					vals2.push_back(datafile.length() > 0 ? UtilJson::toJson(datafile) : UtilJson::toJson(dataFile[0]) );
		
					vector<string> dims;
					
					{
						for (int i = 0; i < boundingBox.min.size(); i++){
		
							vector<string> attr3;
							vector<string> vals3;
							
							{
								vector<string> boundaries;
			
								{
									boundaries.push_back(UtilJson::toJson( boundingBox.min[i] ) );
									boundaries.push_back(UtilJson::toJson( boundingBox.max[i] ) );
								}
							
								attr3.push_back("boundaries");
								vals3.push_back( UtilJson::toJson(boundaries) );
							
							}
							
							dims.push_back( UtilJson::toJson(attr3, vals3));
						}
					}
					
					attr2.push_back("dimensions");
					vals2.push_back(UtilJson::toJson(dims));
				}
				
				
				vector<string> datas;
				{
					datas.push_back( UtilJson::toJson(attr2, vals2));
				}
				
				attr1.push_back("data");
				vals1.push_back(UtilJson::toJson(datas) );
			}
			
			
		
			return UtilJson::toJson(attr1, vals1);
		}
		
		inline void writeJSON(ofstream& out) const{
			
			out << toJson() << endl;
		}
		
		inline const string& getDataFile() const{
		
			return dataFile[0];
		}
		
		inline void setBoundingBox(const Box& b){
			
			boundingBox = b;
			
			domainAdd.clear();
			domainMul.clear();
		
			for (int i = 0; i < dims; i++){
			
				domainAdd.push_back(-boundingBox.min[i]);
				domainMul.push_back( 1.0/(boundingBox.max[i]-boundingBox.min[i]) );
			}	
		}
		
		inline const Box& getBoundingBox() const{
		
			return boundingBox;
		}
		
		inline const string& getName() const{
		
			return name[0];
		}
		
		inline const string& getFileName() const{
		
			return name[0];
		}
		
		inline double normalise(double d, int i) const{
			
			d += domainAdd[i];
			d *= domainMul[i];
			
			if (d < 0)
				return 0;
				
			if (d > 1)
				return 1;
			
			return d;
		}
		
		inline void normalise(Point& d) const{
		
			for (int i = 0; i < dims; i++){
			
				d[i] += domainAdd[i];
				d[i] *= domainMul[i];
				
				if (d[i] < 0)
					d[i] = 0;
				
				if (d[i] > 1)
					d[i] = 1;
			}
			
			
		}
		
		inline double denormalise(double d, int i) const{
			
			d /= domainMul[i];
			d -= domainAdd[i];
			
			return d;
		}
		
		inline void normalise(Point& d, const vector<int>& dimVec) const{
			
			int k = 0;
			for (auto it = dimVec.begin(); it != dimVec.end(); it++){
			
				const int i = (*it);
				d[k] += domainAdd[i];
				d[k] *= domainMul[i];
				
				if (d[k] < 0)
					d[k] = 0;
				
				if (d[k] > 1)
					d[k] = 1;
				
				k++;
			}
		}
		
		inline void print() const{
		
			cout << "dataFile " << getDataFile() << endl;
			cout << "dims " << getDims() << endl;	
			cout << "size " << getSize() << endl;	
			
		}
		
		inline void processTuple(Point& t){
		
			boundingBox.enclose(t);	
			size++;
		}
		
		inline string toString() const{
		
			std::stringstream ss;
	
			ss << "data \t file" << getDataFile() << " dims " << getDims() << " size " << getSize() << endl;
			
			ss << boundingBox << endl;
			
			return ss.str();
		}
	
		friend inline std::ostream& operator<<(std::ostream &os, const DataDescription& m) { 
			
			return os << m.toString();
		}
};


class DatWriter{
public:
	ofstream dat;
	std::vector<int> selectedDims;
	
	std::vector<double> p;
		
	std::vector<double> min;
	std::vector<double> max;
	
	char buffer[16184];
	
	inline DatWriter(const string file, std::vector<int> dimSel){
	
		dat.open(file, std::ios::binary|std::ios::out);
		
		selectedDims = dimSel;
		
		for (int j = 0; j < selectedDims.size(); j++){
		
			p.push_back(0);
			
			min.push_back(UtilMath::largestDouble() );	
			max.push_back(UtilMath::largestNegativeDouble() );
		}
	}
	
	
	inline DatWriter(const string file, int dims){
	
		dat.open(file, std::ios::binary|std::ios::out);
		
		for (int j = 0; j < dims; j++){
		
			p.push_back(0);
			
			min.push_back(UtilMath::largestDouble() );	
			max.push_back(UtilMath::largestNegativeDouble() );
		}
	}
	
	inline bool writeDirect(const vector<double>& p){
	
		for (int i = 0; i < p.size(); i++){
			
			double v = p[i];
			
			if (std::isnan(v))
				return false;
			
			if (v < min[i])
				min[i] = v;
				
			if (v > max[i])
				max[i] = v;
				
			dat.write( reinterpret_cast<char*>( &v ), sizeof(v) );
		}
		
		return true;
	}
	
	inline bool write(const vector<double>& p){
	
		for (int j = 0; j < selectedDims.size(); j++){
					
			const int i = selectedDims[j];
			
			const double v = p[i];
			
			if (std::isnan(v))
				return false;
			
			if (v < min[i])
				min[i] = v;
				
			if (v > max[i])
				max[i] = v;
				
			dat << v;
		}
		
		return true;
	}
	
	inline void close(){
	
		dat.flush();
		dat.close();
	}
		
	

	
};


class IntSeq{

public:

	vector<int> seq;

	inline IntSeq(){
	
	
	}

	inline IntSeq(const IntSeq& s){
	
		seq = s.seq;
	}

	inline IntSeq(const vector<int>& s){
	
		seq = s;
	}
	
	inline IntSeq(IntSeq& s){
	
		seq = s.seq;
	}

	inline IntSeq(const string s){
		
		stringstream commaSplit(s);
		string commaPart;
		
		while ( getline(commaSplit, commaPart, ',') ){
		
			stringstream rangeSplit(commaPart);
			string rangePart;
		
			int a = -1;
			int b = -1;
			
			if (getline(rangeSplit, rangePart, ':')){
			
				istringstream numStream(rangePart);
				numStream >> a;
				b = a;
				
			}
			
			if (getline(rangeSplit, rangePart, ':')){
			
				istringstream numStream(rangePart);
				numStream >> b;
			}
			
			assert (a != -1 && b != -1);
			
			if (a <= b){
			
				for (int i = a; i <= b; i++)
					seq.push_back(i);
				
			} else {
			
				for (long i = a; i >= b; i--)
					seq.push_back(i);
			}
		}
		
	}
	
	inline int get(int i) const{
	
		return seq.at(i);
	}
	
	inline bool operator==(const IntSeq& s) const{
	
		if (seq.size() != s.seq.size() )
			return false;
	
		vector<int> a = seq;
		vector<int> b = s.seq;
		
		sort(a.begin(), a.end());
		sort(b.begin(), b.end());
		
		for (int j = 0; j < a.size(); j++){
		
			if (a.at(j) != b.at(j) )
				return false;
		}
		
		return true;
	}
	
	inline bool contains(int i) const{
	
		for (auto it = seq.begin(); it != seq.end(); it++)
			if ( (*it) == i)
				return true;
				
		return false;
	}
	
	
	inline int size() const{
	
		return seq.size();
	}
	
	inline string toString() const{
		
		std::stringstream ss;
	
		for (int i = 0; i < seq.size(); i++){
			if (i > 0)
				ss << ",";
			
			ss << to_string(seq[i]);
		}
		
		return ss.str();
	}
	
	inline string getFileString() const {
	
		std::stringstream ss;
	
		for (int i = 0; i < seq.size(); i++){
			if (i > 0)
				ss << "_";
			
			ss << to_string(seq[i]);
		}
		
		return ss.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const IntSeq& m) { 
    	
    	return os << m.toString();
	}

};

class DimSelection{


	public:
	
	
		IntSeq dataDims;
		IntSeq selectedDims;
		
		vector<int> ind;
		
		int dims = 0;
		
		bool eq = false;
		
		// selected dimensions (first dimension is "1") given as comma-separated string, e.g., "1,4,5"
		
		
		inline void setDataDims(const string s){
		
			dataDims = IntSeq(s);
			
			selectDims(s);
			
			assert (dataDims == selectedDims);
		}
		
		inline int getDataDims() const{
		
			return dataDims.size();
		}
		
		inline int getDims() const{
		
			return selectedDims.size();
		}
		
		
		inline void selectDims(const IntSeq& seq){
		
			selectedDims = seq;
			
			dims = selectedDims.size();
			
			ind.clear();
			
			int dim = 0;
			
			for (int j = 0; j < dataDims.size(); j++){
			
				if (selectedDims.contains( dataDims.get(j) ))
					ind.push_back(dim);
				
				dim++;
			}
			
			eq = (selectedDims == dataDims);
			
			assert (ind.size() == selectedDims.size() );
		}
		
		inline void selectDims(const string& s){
		
			IntSeq seq(s);
		
			selectDims(seq);
		}
		
		inline void selectDims(const vector<int>& s){
		
			IntSeq seq(s);
		
			selectDims(seq);
		}
		
		
		inline void print() const{
		
			cout << "dims: ";
			UtilString::printList(ind);
		}
		
		inline string getFileString() const{
		
			std::stringstream ss;
		
			for (int i = 0; i < dims; i++){
				if (i > 0)
					ss << "_";
				
				ss << to_string(ind[i]);
			}
			
			return "dims"+ss.str();
		}
		
		inline string getString() const{
		
			std::stringstream ss;
		
			for (int i = 0; i < dims; i++){
				if (i > 0)
					ss << ",";
				
				ss << to_string(ind[i]+1);
			}
			
			return ss.str();
		}
		
		inline void translate(Point& highdim, Point& lowdim) const{
			
			if (eq){
			
				lowdim = highdim;
				
			} else {
			
				for (int i = 0; i < dims; i++)
					lowdim[i] = highdim[ind[i]];
					
				lowdim.id = highdim.id;
			}
			
			
			
		}
};


class CSVReader{
public:
	ifstream f;
	
    char sep;
    
    int dims = 0;
    bool fileopen = false;

	const string file;
	

	
    inline void open(){
    	
    	if (!fileopen){
    		f.open(file, std::ios::in | std::ios::binary);
    		fileopen = true;
    	}
    	
    	assert (f.is_open() );
    }
    
    inline void close(){
    	
    	assert (fileopen);
    	fileopen = false;
    	f.close();
    }

    inline CSVReader(const string file) : file(file){
    
    	// AUTOMATICALLY FIND SEPARATING CHARACTER
    
    	sep = 'X';
    	
    	std::bitset<256> isNumChar;
    	
    	for (int i = 0; i < 256; i++){
			
			const char c = i;
			
			isNumChar[i] = false;
			
			if (UtilMath::isIn(c, 'e','+','-', '.'))
				isNumChar[i] = true;
		
			if (UtilMath::isIn(c, '1','2','3','4','5'))
				isNumChar[i] = true;
				
			if (UtilMath::isIn(c, '6','7','8','9','0'))
				isNumChar[i] = true;
		}
    	
    	open();
    	{
			vector<int> hist;
			
			for (int i = 0; i < 256; i++)
				hist.push_back(0);
		
			if (f.is_open()) {
				
				long size = 100;
				
				for (std::string line; std::getline(f, line); ) {
					
       				size--;
					
					if (size < 100)
					for (int i = 0; i < line.length(); i++){
					
						const int c = line[i];
						
						if (!isNumChar[c])
							hist[c]++;
					}
					
					if (size == 0)
						break;
    			}
    			
				long sepCount = 0;
				
				for (int i = 0; i < hist.size(); i++){
				
					if (hist.at(i) > sepCount){
						
						sepCount = hist.at(i);
						sep = (char) i;
					}
				}
			}
		}
		close();
		
		cout << "sep:" << sep << endl;
    	
    	open();
    	{	
			long size = 1;
				
				
			std::string line;
			std::getline(f, line);
    		
    		dims = 1;
    			
    		for (int i = 0; i < line.length(); i++)
    			if (line[i] == sep)
    				dims++;
    	}
    	close();
    	
    	cout << "dims:" << dims << endl;
    	
    }
    
    inline int getDimNum() const {
    
    	return dims;
    }
    
    inline void writeDataFiles(const string out, const vector<int>& selectedDims){
	
		int dims = getDimNum();
		
		vector<int> isSelected;
		
		for (int i = 0; i < dims; i++)
			isSelected.push_back(0);
		
		
		cout << "dims ";
		for (int j = 0; j < selectedDims.size(); j++){
			
			const int i = selectedDims[j];
			cout << i << ",";
			
			isSelected[i] = 1;
			
		}
		
		cout << endl;
		
		DatWriter dat(out+".dat", selectedDims);
	
		long size = 0;
		
		std::string::size_type sz;
		
		vector<double> p;
		
		open();
    	{	
					
    		for (std::string line; std::getline(f, line); ) {
    	
    			p.clear();
    	
    			stringstream symbs(line);
    	
    			bool error = false;
    	
    			int i = 0;
    	
    			for (std::string symb; std::getline(symbs, symb, sep); i++) {
    			
    				if (isSelected[i] != 1)
    					continue;
    				
    				try {
    					double d = std::stod(symb,&sz);
    					
    					if (d < -10e+12)
    						error = true;
    					
    					if (d > 10e+12)
    						error = true;
    					
    					p.push_back(d);
    					
    				} catch (std::invalid_argument& e){
    					
    					p.push_back(NAN);
    					
    					error = true;
    				}
    			}
    	
    			/*
    			for (std::string symb; std::getline(symbs, symb, sep); ) {
    				
    				try {
    					double d = std::stod(symb,&sz);
    					
    					if (d < -10e+12)
    						error = true;
    					
    					if (d > 10e+12)
    						error = true;
    					
    					p.push_back(d);
    					
    				} catch (std::invalid_argument& e){
    					
    					p.push_back(NAN);
    					
    					error = true;
    				}
    			}*/
    			
    			if (error)
    				cout << line << endl;
    			
    			if (dat.writeDirect(p))
					size++;
					
				if (size % 1000000 == 0)
					cout << "written " << size << " points so far ..." << endl;
    		}
    		
    	}
    	close();
		
		dat.close();
		
		
		cout << "written data file " << out << ".dat" << " with " << size << " points" << endl;
		
		DataDescription desc(selectedDims.size());
		
		desc.setDataFile(out+".dat");
		
		desc.setSize(size);
		
		Point domainMin(dat.min.size() );
		Point domainMax(dat.max.size() );
		
		for (int i = 0; i < dat.min.size(); i++)
			domainMin[i] = dat.min[i];
		
		for (int i = 0; i < dat.max.size(); i++)
			domainMax[i] = dat.max[i];
		
		Box box(domainMin, domainMax);
		
		desc.setBoundingBox(box);
		
		UtilString::writeToFile( out+".json", desc.toJson() );
		
		cout << "written metadata file " << out << ".json" << endl;
	}
};
/*
class CSVReader{
public:
	ifstream f;
		
	const int buflen = 8096;
	
	char* buf = NULL;
    
	char* readptr = NULL;
	
	char* fullbufend = NULL;
	
    char* bufend = NULL;
    char* numptr = NULL;
    char* cur = NULL;
    
    char sep;
    
    int dims = 0;
    
    std::bitset<256> isNumChar;
    
    std::bitset<256> isWhitespaceChar;
    
   	bool numOpen = false;
   	
   	int numlen = 0;
   	
   	bool eof = false;

	bool fileopen = false;

	const string file;

	
    inline void open(){
    	
    	if (!fileopen){
    		f.open(file, std::ios::in | std::ios::binary);
    		fileopen = true;
    	}
    	
    	assert (f.is_open() );
    }
    
    inline void close(){
    	
    	assert (fileopen);
    	fileopen = false;
    	f.close();
    }

    inline CSVReader(const string file) : file(file){
    
    	sep = 'X';
    	
    	for (int i = 0; i < 256; i++){
			
			const char c = i;
			
			isNumChar[i] = false;
			
			if (UtilMath::isIn(c, 'e','+','-', '.'))
				isNumChar[i] = true;
		
			if (UtilMath::isIn(c, '1','2','3','4','5'))
				isNumChar[i] = true;
				
			if (UtilMath::isIn(c, '6','7','8','9','0'))
				isNumChar[i] = true;
		}
    	
    	open();
    	{
			vector<int> hist;
			
			for (int i = 0; i < 256; i++)
				hist.push_back(0);
		
			if (f.is_open()) {
				
				long size = 100;
				
				for (std::string line; std::getline(f, line); ) {
					
       				size--;
					
					if (size < 100)
					for (int i = 0; i < line.length(); i++){
					
						const int c = line[i];
						
						if (!isNumChar[c])
							hist[c]++;
					}
					
					if (size == 0)
						break;
    			}
    			
				long sepCount = 0;
				
				for (int i = 0; i < hist.size(); i++){
				
					if (hist.at(i) > sepCount){
						
						sepCount = hist.at(i);
						sep = (char) i;
					}
				}
			}
		}
		close();
		
		cout << "sep:" << sep << endl;
    	
    	
    	open();
    	{	
			long size = 1;
				
    		for (std::string line; std::getline(f, line); ) {
    	
    			size--;
    			
    			dims = 1;
    			
    			for (int i = 0; i < line.length(); i++)
    				if (line[i] == sep)
    					dims++;
    			
    			if (size == 0)
					break;
    		}
    	}
    	close();
    	
    	cout << "dims:" << dims << endl;
    	
    	buf = new char[buflen*2];
    	
    	fullbufend = buf+buflen*2;
    	
    	readptr = buf;
    	
    	bufend = buf;
    	cur = bufend;
    }
    
    inline ~CSVReader(){
    
    	delete[] buf;
    }
    
    
    
    inline bool readTuple(vector<double>& v){
    	
    	v.clear();
    
    	if (cur == bufend && eof)
    		return false;
    		
    	assert (!numOpen);
    	
    	stringstream numstr;
    	double d;
    	
    	while(true){
    		
    		if (cur == bufend){
    			
    			readptr = buf;
    			
    			if (readptr+buflen >= fullbufend)
    				return false;
    			
    			f.read(readptr, buflen);
    			
				bufend = readptr+f.gcount();
				
				if (bufend == readptr){
					eof = true;
					
					if (numOpen){
    				
    					if (numstr << d)
    						v.push_back(d);
    					else
    						v.push_back( NAN  );
    					
    					numOpen = false;
    					numstr.clear();
						numstr.str(std::string());
    					
    				} else {
    				
    					v.push_back( NAN  );
    				}
					
					return true;
				}
				
				cur = buf;
    		}
    		
    		const char c = (*cur);
    		
    		cur++;
    		
    		const int ii = c;
    		
    		if (c == sep || c == '\n'){
    		
    			if (numOpen){
    				
					if (numstr << d)
    					v.push_back(d);
    				else
    					v.push_back( NAN  );
    				
					numOpen = false;
					numstr.clear();
					numstr.str(std::string());
    				
    			} else {
    			
    				v.push_back( NAN  );
    			}
    			
    			if (c == '\n'){
    			
    				if (v.size() != dims){
    				
    					numOpen = false;
    					numstr.clear();
						numstr.str(std::string());
    					v.clear();
    					
    					continue;
    				}
    				
    				return true;
    			}
    		
    		} else if (isNumChar[ii]){
    				
    			numOpen = true;
    			numstr << c;
    		}
    		
    	}
    	
    	assert (false);
    }
    
    inline int getDimNum() const {
    
    	return dims;
    }
    
    inline void writeDataFiles(const string out, const vector<int>& selectedDims){
	
		int dims = getDimNum();
		
		DatWriter dat(out+".dat", selectedDims);
	
		long size = 0;
		
		std::stringstream stream;
		
		open();
		
		for (vector<double> p; readTuple(p); ){
		
			if (dat.write(p))
				size++;
			
			if (size % 1000000 == 0)
				cout << "written " << size << " points so far ..." << endl;
		}
		
		//cout << "dims " << dims << " =?= " << dat.min.size() << " min.size()" << endl; 
		
		close();
		
		
		dat.close();
		
		
		cout << "written data file " << out << ".dat" << " with " << size << " points" << endl;
		
		DataDescription desc(selectedDims.size());
		
		desc.setDataFile(out+".dat");
		
		desc.setSize(size);
		
		Point domainMin(dat.min.size() );
		Point domainMax(dat.max.size() );
		
		for (int i = 0; i < dat.min.size(); i++)
			domainMin[i] = dat.min[i];
		
		for (int i = 0; i < dat.max.size(); i++)
			domainMax[i] = dat.max[i];
		
		Box box(domainMin, domainMax);
		
		desc.setBoundingBox(box);
		
		UtilString::writeToFile( out+".json", desc.toJson() );
		
		cout << "written metadata file " << out << ".json" << endl;
	}
    
};
*/


/*
class CSVReader{
public:
	
	const string file;
	
	bool skipNan = false;
	
	std::bitset<256> isNumChar;
	
	inline CSVReader(const string f): file(f){
		
		for (int i = 0; i < 256; i++){
			
			const char c = i;
			
			isNumChar[i] = false;
			
			if (UtilMath::isIn(c, 'e','+','-', '.'))
				isNumChar[i] = true;
		
			if (UtilMath::isIn(c, '1','2','3','4','5'))
				isNumChar[i] = true;
				
			if (UtilMath::isIn(c, '6','7','8','9','0'))
				isNumChar[i] = true;
		}
	}
	
	inline void readLine(const string& line, vector<double>& v, char* buffer) const{
		
		bool open = false;
		
		bool nanOpen = true;
		
		char* ptr = buffer;
		
		for (int i = 0; i < line.length(); i++){
			
			const char c = line[i];
			
			if (isNumChar[(int)c] ){
				
				open = true;
				
				nanOpen = false;
				
				(*ptr) = c;
				ptr++;

			} else {
				
				if (open){
					
					v.push_back(std::strtod(buffer, &ptr) );
					nanOpen = false;
					
				} else if (!skipNan){
					
					if (c == ','){
					
						if (nanOpen){
							
							v.push_back(NAN);
							nanOpen = false;
						}
						
					} else {
						
						nanOpen = true;
					}
				}
				
				ptr = buffer;
				open = false;
			}
		}
		
		if (nanOpen)
			v.push_back(NAN);
		
		if (open)
			v.push_back(std::strtod(buffer, &ptr) );
		
		
	}
	
	inline int getDimNum() const {
		
		std::ifstream in(file);
	
		char buffer[16184];
		char csvbuffer[16184];		
		
		in.rdbuf()->pubsetbuf(buffer, 16184);
		
		int dims = -1;
		long size = 0;
		
		vector<double> v;
		
		std::stringstream stream;
		
		if (in.is_open()) {
				
			std::string line;
			
			while (std::getline(in, line)){
				
				readLine(line, v, &(csvbuffer[0]));
				
				const int dimCount = v.size();
				
				v.clear();
				
				cout << dimCount << "::" << line << endl;
					
				if (size == 0){
				
					size++;
					continue;
				}
				
				
				
				
				if (dimCount > 0){
					
					if (dims == -1)
						dims = dimCount;
					
					size++;
					
					if (size > 10)
						return dims;
					
					assert (dimCount == dims);
				}
			}
		}
		
		in.close();
		
		assert (dims != -1);
		
		return dims;
	}
	
	inline void writeDataFiles(const string out, const vector<int>& selectedDims) const{
	
		std::vector<double> p;
		
		int dims = getDimNum();
		
		std::ifstream in(file);
		
		char buffer[16184];
		
		char csvbuffer[16184];		
		
		in.rdbuf()->pubsetbuf(buffer, 16184);
		
		DatWriter dat(out+".dat", selectedDims);
	
		long size = 0;
		
		std::stringstream stream;
		
		if (in.is_open()) {
				
			std::string line;
		
			while (std::getline(in, line)){
				
				readLine(line, p, &(csvbuffer[0]));
				
				if (p.size() > 0){
					
					if (p.size() != dims){
						
						p.clear();
						continue;
						cout << "p.size() " << p.size() << " != " << dims << " dims" << endl;
					}
					
					if (dat.write(p))
						size++;
					
					if (size % 1000000 == 0)
						cout << "written " << size << " points so far ..." << endl;						
					
					assert (p.size() == dims);
				}
				
				p.clear();
				
			}
		}
		
		//cout << "dims " << dims << " =?= " << dat.min.size() << " min.size()" << endl; 
		
		dat.close();
		
		in.close();
		
		
		cout << "written data file " << out << ".dat" << " with " << size << " points" << endl;
		
		DataDescription desc(selectedDims.size());
		
		desc.setDataFile(out+".dat");
		
		desc.setSize(size);
		
		Point domainMin(dat.min.size() );
		Point domainMax(dat.max.size() );
		
		for (int i = 0; i < dat.min.size(); i++)
			domainMin[i] = dat.min[i];
		
		for (int i = 0; i < dat.max.size(); i++)
			domainMax[i] = dat.max[i];
		
		Box box(domainMin, domainMax);
		
		desc.setBoundingBox(box);
		
		UtilString::writeToFile( out+".json", desc.toJson() );
		
		cout << "written metadata file " << out << ".json" << endl;
		
	}
};*/


class DataSet{

public:
	class iterator{
	
	private:
	
		std::unique_ptr<PointStream> rawDataStream;
		std::unique_ptr<Timer> t;
		
		const DataSet* ds;
		
		Point rawPoint;
		Point current;
		long pos = 0;
		const long datasize;
		
		inline void read(){
		
			assert(!sortingIterator);
			
			rawDataStream->read(rawPoint);
			ds->dimSel.translate(rawPoint, current);
			ds->getDescription().normalise(current, ds->dimSel.ind);
		}
	
		// sorting related
		
		const bool verbose;
		const bool sortingIterator;
		
		int currentBlock = -1;
		long blockPos = 0;
		
		vector<Point> block;
		vector<Point> borders;
		
		vector<int> lex;
	
		const int DIMS;
		
		bool blockSorted = false;
		
		const PointComparator comp;
	
	public:
		
		
		inline iterator() : sortingIterator(false), comp(0), DIMS(0), datasize(0), block(0), borders(0), lex(0), verbose(false){
		
			
		}
		
		inline iterator(const DataSet* ds, bool verb) : sortingIterator(false), comp(0), verbose(verb), block(0), borders(0), lex(0), ds(ds), pos(0), rawPoint(ds->getRawDims() ), DIMS(ds->getDims() ),current(ds->getDims() ), datasize(ds->getSize() ){
		
			if (verbose)
				t.reset(new Timer(ds->getDescription().getName() ));
		
			if (verbose)
				t->start(10, 1000, datasize);
			
			rawDataStream.reset(ds->getStream() );
			
			read();
			
			pos = 0;
		}
		
		inline iterator(const DataSet* ds, const PointComparator& comp, long blockSize) : verbose(ds->getVerbose() ), comp(comp), sortingIterator(true), ds(ds), DIMS(ds->getDims() ), datasize(ds->getSize() ){
		
			block.reserve(blockSize);
			
			if (verbose){
				t.reset(new Timer(ds->getDescription().getName()+"sorted" ));
				t->start(10, 100, datasize);
			}
			
			
			shared_ptr<vector<Point>> s = std::move( ds->getRandomSubsetPoints(100000, "abc") );
			
			vector<Point>& sample = (*s);
			
			const long sampleBlockSize = round(0.9*blockSize*1.0/datasize*sample.size() );
			
			sort(sample.begin(), sample.end(), comp);
			
			{
				Point first(DIMS);
				first.setAll(std::numeric_limits<double>::lowest() );
				borders.push_back(first);
			}
			
			
			for (long j = sampleBlockSize; j < (sample.size()-sampleBlockSize); j += sampleBlockSize){
			
				assert(sample.at(j).size() == DIMS);
			
				borders.push_back(sample.at(j));
				
			}
			
			{
				Point last(DIMS);
				last.setAll(std::numeric_limits<double>::max() );
				borders.push_back(last);
			}
			
			(*this) += 0;
		}
		
		inline Point& operator*(){
		
        	return sortingIterator ? block[blockPos] : current;
    	}
    	
    	inline void close(){
    	
    		if (pos < datasize)
    			if (verbose)
    				t->end();
    	}
    	
    	inline Point& getPoint(){
		
			return sortingIterator ? block[blockPos] : current;
    	}
    	
    	inline Point& getRawPoint(){
		
			assert (!sortingIterator);
			
        	return rawPoint;
    	}
    	
    	inline const long operator()() const{
    		
    		return 0;
    	}
    	
    	inline long getPos() const{
    	
    		return pos;
    	}
    	
		inline bool operator !=(const long b) const{
		
			if (pos < datasize){
			
				return true;
				
			} else {
		
				if (verbose)
					t->end();
				
				return false;
			}
			
		}
		
		inline void loadNextBlock(){
	
			currentBlock++;
			
			block.clear();
			
			if (currentBlock >= (borders.size()-1) ){
			
				blockSorted = true;
				return;
			}
			
			PointStream* ps = ds->getStream();
			
			const Point& min = borders.at(currentBlock);
			const Point& max = borders.at(currentBlock+1);
			
			const int i = comp.getMainDim();
			
			const double mini = min[i];
			const double maxi = max[i];
			
			{
				
				for (auto it = ds->beginUnsorted(false); it != ds->end(); it++){
				
					const Point& p = (*it);
					const double pi = p[i];
			
					if (pi > maxi)
						continue;
			
					if (pi < mini)
						continue;
					
					if (comp.lessOrEqual(min, p) && comp.less(p, max) )
						block.push_back(p);
				}
				
			}
			
			//cout << "mini " << mini << " maxi " << maxi << endl;
			//cout << "block.size() " << block.size() << endl;
			
			blockSorted = false;
		}
		
		inline iterator& operator +=(const long n){
	
			if (sortingIterator){
			
				blockPos += n;
				pos += n;
			
				if (pos < datasize){
			
					while (blockPos >= block.size()){
			
						blockPos -= block.size();
						loadNextBlock();
					}
				
					if (!blockSorted){
			
						//UtilPoint::sort(lex, block);
						
						sort(block.begin(), block.end(), comp);
						
						//assert (UtilPoint::isSorted(lex, block) );
						blockSorted = true;
					}
				}
				
				if (verbose)
					t->tick();
				
				return (*this);
			}
	
			if (n <= 0)
				return (*this);
	
			pos += n;
			
			if (pos >= datasize)
				return (*this);
			
			if (verbose)
				t->tick();
			
			if (n > 1)
				rawDataStream->skip(n-1);
			
			read();
			
			return (*this);
		}
		
		inline void print(){
			
			const Point p = getPoint();
		
			cout << p << endl;
		}
		
		inline iterator& operator++(){
			
			return ( (*this) += 1 );
		}
		
		inline iterator& operator++(int n){
			
			return ( (*this) += 1 );
		}
		
		
	};


protected:
	
	inline DataSet(const string s){// : t(s){
	
		
	}
	
	inline virtual PointStream* getStream() const = 0;

	bool sortingIterator = false;
		
	//vector<int> sortingLex = {};
	
	PointComparator comp = 0;
	
	long sortingIteratorBuffer = 0;
	
	bool verbose = false;

public:

	virtual ~DataSet(){};

	DimSelection dimSel;
	
	iterator end;
	
	inline virtual void enableSorting(const PointComparator& comp1, long sb= 1000000){
	
		comp = comp1;
		sortingIterator = true;
		sortingIteratorBuffer = sb;
	}
	
	inline void disableSorting(){
	
		sortingIterator = false;
	}
	
	virtual inline const string getName() const{
	
		return "dataset";
	}
	
	virtual inline const string getFileName() const{
	
		return "dataset";
	}
	
	iterator begin() const{
	
		const DataSet* ds = this;
		
		if (sortingIterator)
			return iterator(ds, comp, sortingIteratorBuffer);
		
		return iterator(ds, getVerbose());
	}
	
	shared_ptr<iterator> getIterator() const{
	
		const DataSet* ds = this;
		
		shared_ptr<iterator> ret;
		
		if (sortingIterator)
			ret.reset(new iterator(ds, comp, sortingIteratorBuffer));
		else
			ret.reset(new iterator(ds, getVerbose() ) );
		
		return ret;
	}
	
	shared_ptr<iterator> getEnd() const{
	
		shared_ptr<iterator> ret(new iterator() );
		return ret;
	}
	
	
	iterator beginUnsorted(bool verb=true) const{
	
		const DataSet* ds = this;
		
		return iterator(ds, getVerbose() && verb);
	}
	
	inline virtual const DataDescription& getDescription() const = 0;
	
	inline void setVerbose(bool verbose){
	
		this->verbose = verbose;
	}
	
	inline bool getVerbose() const{
		
		return verbose;
	}

	inline virtual long getSize() const = 0;
	
	inline void add(Point& p){
	
		throw "modification not supported!";
	}
	
	inline virtual long boxCount(const Point& qmin, const Point& qmax) const{
	
		long ret = 0;
	
		for (auto it = beginUnsorted(); it != end(); it++){
		
			if (UtilPoint::contains(qmin, qmax, (*it) ))
				ret++;
		}
		
		return ret;
	}


	inline void selectDims(const string& s){
	
		dimSel.selectDims(s);
		
	}
	
	inline void setDataDims(const string& s){
	
		dimSel.setDataDims(s);
	}
	
	inline void selectAllDims(){
		
		selectDims("1:"+to_string(getDescription().getDims() ));
	}
	
	inline void print() const {
	
		getDescription().print();
		dimSel.print();
		
		
		int i = 0;
		
		for (auto it = beginUnsorted(false); it != end(); it++){
		
			it.print();
			
			if (i == 4)
				it += getSize()-10;
				
			i++;
		}
	}
	
	inline void writeToFile(const string& file) const {
	
		const int DIMS = getDims();
		
		
		ofstream fdata(file, std::ios::binary );
		
		assert (fdata.is_open() );
		
		
	
		for (auto it = begin(); it != end(); it++){
		
			const Point& p = (*it);
		
			for (int i = 0; i < DIMS; i++){
			
				double d = p[i];
				fdata.write( reinterpret_cast<char*>( &d ), sizeof(d) );
			}
		}
		
		fdata.close();
		
		cout << "WRITTEN DATA TO FILE: " << file << endl;
		
	}
	
	inline int getDims() const{
	
		return dimSel.getDims();
	}
	
	inline int getRawDims() const{
	
		return dimSel.getDataDims();
	}
	
	inline virtual shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=1000000) const = 0;
	inline virtual shared_ptr<DataSet> getRandomSubset(long size, const string seed="abc123") const = 0;
	inline virtual shared_ptr<DataSet> getRandomPoints(long size, const string seed="abc123") const = 0;
	inline virtual shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const = 0;
	
	inline shared_ptr<DataSet> getEmptyCopy(long size=1) const{
		
		shared_ptr<vector<Point>> v(new vector<Point>());
		
		if (size >= 0)
			v->reserve(size);
		
		return getCopy(v);
	}
	
	/*
	inline virtual shared_ptr< vector<Point> > getPointsInRange(const Box& b, long size = getSize() ) const{
		
		shared_ptr<vector<Point>> vec(new vector<Point>() );
		
		for (auto it = beginUnsorted(false); it != end(); it++){
			
			if (b.contains( (*it) ))
				vec->push_back( (*it) );
		}
		
		return vec;
	}*/
	
	
	inline virtual shared_ptr< vector<Point> > getRandomSubsetPoints(long size, const string seed="123") const{
		
		
		shared_ptr<vector<Point>> vec(new vector<Point>() );
		
		assert (size <= getSize() );
		
		if (size >= getSize() ){
		
			for (auto it = beginUnsorted(false); it != end(); it++)
				vec->push_back( (*it) );
			
			std::random_shuffle(vec->begin(), vec->end() );
			
			return vec;
		}
		
		shared_ptr<vector<long>> inds = UtilCollection::getRandomIndices( getSize(), size, seed);
		
		vec->reserve(size);
		
		auto it1 = beginUnsorted(false);
		long last = 0;
		
		long k = 0;
		
		for (auto it = inds->begin(); it != inds->end(); it++){
		
			const long ind = (*it);
		
			it1 += ind-last;
			
			Point p = (*it1);
			
			vec->push_back(p);
			
			last = ind;
		}
		
		it1.close();
		
		std::random_shuffle(vec->begin(), vec->end() );
		
		
		return vec;
	}
	
	inline virtual shared_ptr< vector<Point> > getRandomUnifPoints(long size, string seed="123") const{
		
		shared_ptr<vector<Point>> vec(new vector<Point>() );
		
		vec->reserve(size);
		
		UniformPoints unif(getRawDims(), seed );
		
		for (long j = 0; j < size; j++){
		
			Point p;
			
			unif.randomPoint(getDescription().getBoundingBox(), p);
			
			vec->push_back(p);
		}
		
		return vec;
	}
	
	
	
	
	
	const virtual Point operator[](std::size_t idx) const{
		
		auto it = begin();
		
		it += idx;
		
		if (it != end() )
			return (*it);
		
		throw string("out of bounds!");
	}
	
	inline double getCount(double selPercent) const{
	
		return selPercent/100.0*getSize();
	}
	
};


////////////////////////////////////////////////////////////////////////////////////////////////

class VectorDataSet : public DataSet{

private:
		
	DataDescription desc;

protected:
	
	inline PointStream* getStream() const override{
	
		return new VectorPointStream(vec);
	}

public:
	
	shared_ptr<vector<Point>> vec;
	
	inline VectorDataSet(int dims) : desc(dims), vec(new vector<Point>() ), DataSet("vector"){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
		
		Point pmin(dims);
		Point pmax(dims);
		
		for (int i = 0; i < dims; i++){
		
			pmin[i] = 0;
			pmax[i] = 1;
		}
		
		Box b(pmin, pmax);
		
		desc.setSize(0);
		desc.setBoundingBox(b);
	}
	
	inline VectorDataSet(const Box& b, shared_ptr<vector<Point>> vec, int dims) : desc(dims), vec(vec), DataSet("vector"){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
		
		desc.setSize(vec->size());
		desc.setBoundingBox(b);
	}
	
	inline void add(Point& p){
	
		vec->push_back(p);
		desc.setSize(vec->size());
	}
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return vec->size();
	}
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
	
		Box b(pmin, pmax);	
		
		shared_ptr<DataSet> ret( new VectorDataSet(b, v , getDims() ));
		
		ret->setDataDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<vector<Point>> vecCopy(new vector<Point>() );
		
		vecCopy->reserve(vec->size() );
		
		for (auto it = vec->begin(); it != vec->end(); it++)
			vecCopy->push_back( (*it) );
		
		shared_ptr<DataSet> ret( new VectorDataSet(desc.getBoundingBox(), vecCopy, getDims() ));
		ret->setDataDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp);
		
		return ret;
	}
	
	inline void enableSorting(const PointComparator& comp, long memory=0) override{
	
		sort( vec->begin(), vec->end(), comp );
		sortingIterator = false;
		sortingIteratorBuffer = memory;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed) );
	}	
};

////////////////////////////////////////////////////////////////////////////////////////////////

class FileDataSet : public DataSet{

private:
		
	DataDescription desc;
	const string descFile;
		
protected:
	
	inline PointStream* getStream() const override{
	
		if (UtilString::contains(desc.getDataFile(), ".csv") )
			return new CSVPointStream(desc.getDataFile());
			
		if (UtilString::contains(desc.getDataFile(), ".pdat") ){	
		
			return new BinaryPointStream(desc.getDataFile(), desc.getDims(), 16);
		}
			
		if (UtilString::contains(desc.getDataFile(), ".dat") ){
		
			return new BinaryPointStream(desc.getDataFile(), desc.getDims());
		}
		
		return NULL;
	}

public:

	inline const string getName() const override{
	
		return "file "+descFile;
	}
	
	inline const string getFileName() const override{
	
		return UtilString::getFileName(desc.getDataFile())+"_"+dimSel.getFileString();
	};

	inline FileDataSet(const string descFile) : descFile(descFile), desc(descFile), DataSet(descFile){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<DataSet> ret( new FileDataSet(descFile));
		
		ret->setDataDims(dimSel.dataDims.toString() );
		ret->selectDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp, memory);
		
		return ret;
	}
	
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return desc.getSize();
	}
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
	
		Box b(pmin, pmax);	
		shared_ptr<DataSet> ret( new VectorDataSet(b, v, getDims() ));
		
		ret->selectDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed ) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed ) );
	}
};


//

class TestDataSet : public DataSet{

private:
		
	DataDescription desc;
	const string descFile;
	
	const int DIMS;
	const long points;
	const int rep;
	
protected:
	
	inline PointStream* getStream() const override{
		
		return new TestPointStream(DIMS, points/rep);
	}

public:
	
	inline TestDataSet(int dims, long pointNum, int rep=1) : DIMS(dims), rep(rep), points(pointNum*rep), descFile(""), desc(dims), DataSet("test"){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
		
		desc.setSize(points);
		
		Box b1(DIMS);
		desc.setBoundingBox(b1);
		
		for (auto it = this->begin(); it != this->end(); it++)
			b1.enclose( (*it) );
		
		desc.setBoundingBox(b1);
	}
	
	
	inline long boxCount(const Point& qmin, const Point& qmax) const override{
	
		
		double ret = 1;
		
		for (int i = 0; i < DIMS; i++){
			
			const double frac = (pow(qmax[i], i+1)-pow(qmin[i], i+1) );
			ret *= frac;
		}
		
		return ret*points;
	}
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return desc.getSize();
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<DataSet> ret( new TestDataSet(DIMS, points));
		
		ret->setDataDims(dimSel.dataDims.toString() );
		ret->selectDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp, memory);
		
		return ret;
	}
	
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
		
		Box b(pmin, pmax);	
		shared_ptr<DataSet> ret( new VectorDataSet(b, v, getDims() ));
		
		ret->selectDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed ) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed ) );
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////

class ZipfDataSet : public DataSet{

private:
		
	DataDescription desc;
	const string descFile;
	
	const string seed;
	const int DIMS;
	const long points;
	const long clusters;
	
protected:
	
	inline PointStream* getStream() const override{
		
		return new ZipfPointStream(seed, DIMS, clusters, points);
	}

public:

	inline ZipfDataSet(int dims, long pointNum, long clusterNum, const string s="abc") : seed(s), DIMS(dims), clusters(clusterNum), points(pointNum), descFile(""), desc(dims), DataSet("zipf"){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
		
		desc.setSize(pointNum);
		
		Box b1(DIMS);
		desc.setBoundingBox(b1);
		
		for (auto it = this->begin(); it != this->end(); it++)
			b1.enclose( (*it) );
		
		desc.setBoundingBox(b1);
	}
	
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return desc.getSize();
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<DataSet> ret( new ZipfDataSet(DIMS, points, clusters, seed));
		
		ret->setDataDims(dimSel.dataDims.toString() );
		ret->selectDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp, memory);
		
		return ret;
	}
	
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
		
		Box b(pmin, pmax);	
		shared_ptr<DataSet> ret( new VectorDataSet(b, v, getDims() ));
		
		ret->selectDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed ) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed ) );
	}
	
	inline const string getName() const override{
	
		return "zipf -dims "+to_string(DIMS)+" -num "+to_string(points)+" -clusters "+to_string(clusters);
	}
};


class HyperplaneDataSet : public DataSet{

private:
		
	DataDescription desc;
	const string descFile;
	
	const string seed;
	const int DIMS;
	
	long hyperplanes;
	long points;

protected:

	inline PointStream* getStream() const override{
		
		return new HyperplanePointStream(seed, DIMS, hyperplanes, points);
	}

public:
	
	inline HyperplaneDataSet(int dims, long pointNum, long hyperNum, const string s="abc") : seed(s), DIMS(dims), hyperplanes(hyperNum), points(pointNum), descFile(""), desc(dims), DataSet("hyperplanes"){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
		
		desc.setSize(pointNum);
		
		Point p1(DIMS, 0);
		Point p2(DIMS, 1);
		
		Box b1(p1, p2);
		desc.setBoundingBox(b1);
		
	}
	
	
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return desc.getSize();
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<DataSet> ret( new HyperplaneDataSet(DIMS, points, hyperplanes, seed));
		
		ret->setDataDims(dimSel.dataDims.toString() );
		ret->selectDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp, memory);
		
		return ret;
	}
	
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
		
		Box b(pmin, pmax);	
		shared_ptr<DataSet> ret( new VectorDataSet(b, v, getDims() ));
		
		ret->selectDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed ) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed ) );
	}
	
	
	inline const string getFileName() const override{
	
		return "hyperplanes_"+to_string(DIMS)+"d_num"+UtilString::integerToShortString(points)+"_hyper"+UtilString::integerToShortString(hyperplanes)+"_"+dimSel.getFileString();
	}
	
	inline const string getName() const override{
	
		return "hyperplanes -dims "+to_string(DIMS)+" -num "+UtilString::integerToShortString(points)+" -planes "+UtilString::integerToShortString(hyperplanes);
	}
};


class ChainPointStream  : public PointStream{
	
	vector<shared_ptr<DataSet>> datasets;
	
	
	shared_ptr<DataSet::iterator> end;
	shared_ptr<DataSet::iterator> it;
	
	long countdown;
	
	long id = 0;

public:
	inline ChainPointStream(const vector<shared_ptr<DataSet>>& ds){
	
		for (auto it = ds.begin(); it != ds.end(); it++)
			datasets.push_back( (*it) );
		
		nextDataSet();
	}
	
	inline bool read(Point& p) override{
	
		if (it){
		
		} else {
		
			return false;
		}
		
		p = (*(*it));
		
		p.id = id;
		
		id++;
		
		countdown--;
		
		if (countdown == 0)
			nextDataSet();
		else
			(*it)++;
		
		return true;
	}
	
	inline void nextDataSet(){
	
		it.reset();
		end.reset();
		
		if (datasets.size() > 0){
	
			it = std::move( datasets.back()->getIterator() );
			end = std::move( datasets.back()->getEnd() );
		
			countdown = datasets.back()->getSize()-1;
			
			datasets.pop_back();
		}
		
		
		while (countdown == 0 && datasets.size() > 0){
		
			nextDataSet();
		}
		
	}
	
	inline bool skip(long off) override{
	
		while (off > countdown){
		
			nextDataSet();
			off -= countdown;
		}
		
		(*it) += off;
		
		return true;
	}
	
	inline void close() override{
			
		(*it) += countdown;
		(*it) != (*end)();	
	}
	

};

class ChainDataSet : public DataSet{

private:
		
	DataDescription desc;
	
	vector<shared_ptr<DataSet>> datasets;

protected:

	inline PointStream* getStream() const override{
		
		return new ChainPointStream(datasets);
	}

public:
	
	inline ChainDataSet(const vector<shared_ptr<DataSet>>& ds) : desc( (*ds[0]).getDims() ), DataSet("chain"){
	
		for (auto it = ds.begin(); it != ds.end(); it++){
		
			datasets.push_back( (*it) );
		}
	
		setDataDims("1:"+to_string(desc.getDims() ) );
		
		long size = 0;
		
		for (long j = 0; j < ds.size(); j++)
			size += (*ds[j]).getSize();
		
		
		desc.setSize(size);
		
		Point p1(desc.getDims(), 0);
		Point p2(desc.getDims(), 1);
		
		Box b1(p1, p2);
		desc.setBoundingBox(b1);
	}
	
	
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return desc.getSize();
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<DataSet> ret( new ChainDataSet(datasets));
		
		ret->setDataDims(dimSel.dataDims.toString() );
		ret->selectDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp, memory);
		
		return ret;
	}
	
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
		
		Box b(pmin, pmax);	
		shared_ptr<DataSet> ret( new VectorDataSet(b, v, getDims() ));
		
		ret->selectDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed ) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed ) );
	}
	
	
	inline const string getFileName() const override{
	
		return "chain";
	}
	
	inline const string getName() const override{
	
		return "chain";
	}
};

class SphereDataSet : public DataSet{

private:
		
	DataDescription desc;
	const string descFile;
	
	const string seed;
	const int DIMS;
	
	long points;
	int mode;

protected:

	inline PointStream* getStream() const override{
		
		return new SpherePointStream(seed, DIMS, mode);
	}

public:
	
	inline SphereDataSet(int dims, long pointNum, int mode, const string s="abc") : mode(mode), seed(s), DIMS(dims), points(pointNum), descFile(""), desc(dims), DataSet("sphere"){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
		
		desc.setSize(pointNum);
		
		Point p1(DIMS, 0);
		Point p2(DIMS, 1);
		
		Box b1(p1, p2);
		desc.setBoundingBox(b1);
		
	}
	
	
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return desc.getSize();
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<DataSet> ret( new SphereDataSet(DIMS, points, mode, seed));
		
		ret->setDataDims(dimSel.dataDims.toString() );
		ret->selectDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp, memory);
		
		return ret;
	}
	
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
		
		Box b(pmin, pmax);	
		shared_ptr<DataSet> ret( new VectorDataSet(b, v, getDims() ));
		
		ret->selectDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed ) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed ) );
	}
	
	
	inline const string getFileName() const override{
	
		return "sphere_"+to_string(DIMS)+"d_num"+UtilString::integerToShortString(points)+"_"+dimSel.getFileString();
	}
	
	inline const string getName() const override{
	
		return "sphere -dims "+to_string(DIMS)+" -num "+UtilString::integerToShortString(points);
	}
};

inline shared_ptr<DataSet> newDataSet(const string s){

	if ( UtilString::contains(s, "#") ){
	
		shared_ptr<DataSet> ret;
	
		stringstream f(s);
	
		for (std::string ss; std::getline(f, ss,'#'); ){
		
			if (ret)
				return ret->getRandomSubset( UtilString::strToInt(ss) );
			
			ret = newDataSet(ss);
		}
		
		assert (false);
	}
	

	shared_ptr<DataSet> ret;

	if (UtilString::contains(s, "+")){
	
		stringstream f(s);
	
		vector<shared_ptr<DataSet>> datas;
		
		for (std::string ss; std::getline(f, ss,'+'); )
			datas.push_back(newDataSet(UtilString::trim(ss) ));
		
		ret.reset(new ChainDataSet(datas));
		
		return ret;
	}

	cout << "newDataSet " << s << endl;

	

	const string FILE_PREFIX = "file ";

	

	if (UtilString::startsWith(s, FILE_PREFIX)){
	
		ret.reset(new FileDataSet(s.substr(FILE_PREFIX.length()) ));
		
	} else {

	

		Params p(translateUnits(s));
		
		if (UtilString::startsWith(p.name, "uniform"))
			ret.reset(new SphereDataSet(p.get("dims"), p.get("num"), 0));
		
		
		if (UtilString::startsWith(p.name, "sphere"))
			ret.reset(new SphereDataSet(p.get("dims"), p.get("num"), 1));
		
		if (UtilString::startsWith(p.name, "spiral"))
			ret.reset(new SphereDataSet(p.get("dims"), p.get("num"), 2 ));
		
		
		if (UtilString::startsWith(p.name, "hyperplanes"))
			ret.reset(new HyperplaneDataSet(p.get("dims"), p.get("num"), p.get("planes")  ));
	}
	
	if (ret){
	
		return ret;
	}
		
	assert (false);
}


#endif
