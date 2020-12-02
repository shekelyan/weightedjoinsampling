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

#ifndef SHEKELYAN_DATAIMPORT_H
#define SHEKELYAN_DATAIMPORT_H

#include <approxdata/testing/testing.h>

#include <glob.h>
#include <gzip/compress.hpp>
#include <gzip/decompress.hpp>
#include <gzip/utils.hpp>

constvector<std::string> glob(const std::string& pattern) {
    using namespace std;

    // glob struct resides on the stack
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // do the glob operation
    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if(return_value != 0) {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    // collect all the filenames into a std::list<std::string>
    constvector<string> filenames;
    for(size_t i = 0; i < glob_result.gl_pathc; ++i) {
        filenames.push_back(string(glob_result.gl_pathv[i]));
    }

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}

class LineReader{

	
	
protected:
	ifstream in;
	inline LineReader(){}
	
public:
	
	inline LineReader(const string& filename) : in(filename){
		
		
	}
	
	inline bool readLine(string& line){
	
		if (std::getline(in, line))
			return true;
		else
			return false;
	}
	
	inline void close(){
	
		in.close();
	}
};

class GzipLineReader{
	
	shared_ptr<std::stringstream> stream;

public:
	inline GzipLineReader(const string& filename){
	
		ifstream in(filename, std::ios_base::in | std::ios_base::binary);
		
		assert(in.is_open() );
		
		string zip((std::istreambuf_iterator<char>(in.rdbuf())),
								   std::istreambuf_iterator<char>());
		in.close();
		
		std::size_t limit = 500 * 1024 * 1024; 
		// file should be about 500 mb uncompressed
		gzip::Decompressor decomp(limit);
		
		string output;
		
		decomp.decompress(output, zip.data(), zip.size());	
		
		stream.reset(new std::stringstream(output) );
	}

	inline bool readLine(string& line){
	
		if (std::getline( (*stream), line))
			return true;
		else
			return false;
	}
	
	inline void close(){
	
	}
	
	
};



class DataReader{


	vector<string> sources;
	
	shared_ptr<LineReader> lineReader;
	shared_ptr<GzipLineReader> gzipLineReader;
	
	vector<vector<int>> dimSels;
	
	std::bitset<256> sel;
	
	int type = -1;
	
	const int CSV = 1;
	const int JSON = 8;
	const int CSVGZ = 2;
	const int PLY = 4;
	
	char sep = ',';
	
	vector<char> seps;
	
	bool firstSource = true;
	
	string line;
	string symb;
	
	size_t sz;

public:	
	inline void addSource(const string& filename, const vector<int>& dimSel, char sep = ','){
	
		sources.push_back(filename);
		dimSels.push_back(dimSel);
		seps.push_back(sep);
	}
	
	inline bool nextSource(){
	
		if (!firstSource){
		
			close();
			
			sources.pop_back();
			dimSels.pop_back();
			seps.pop_back();
		}
		
		
		firstSource = false;
		
		if (sources.size() == 0)
			return false;
		
		
		const vector<int>& v = dimSels.back();
		const string& s = sources.back();
		sep = seps.back();
		
		sel.reset();
		
		for (int i = 0; i < v.size(); i++)
			sel.set( v[i] );
		
		if (UtilString::endsWith(s, ".csv"))
			type = CSV;
			
		if (UtilString::endsWith(s, ".ply"))
			type = PLY;
		
		if (UtilString::endsWith(s, ".csv.gz"))
			type = CSVGZ;
			
		assert (type != -1);
		
		if ( type & (CSV | PLY)  )
			lineReader.reset(new LineReader(s) );
		
		if (type & CSVGZ)
			gzipLineReader.reset(new GzipLineReader(s) );
		
		return true;
	}
	
	inline bool nextPoint(vector<double>& p){
	
		if (type & (CSV | PLY | CSVGZ) ){
			
			const bool hasNextLine = firstSource ? false : (type & (CSV | PLY) ) ? lineReader->readLine(line) : (type & CSVGZ) ? gzipLineReader->readLine(line) : false;
			
			if (!hasNextLine){
			
				const bool hasNextSource = nextSource();
				
				if (!hasNextSource)
					return false;
			}
			
			std::stringstream ss(line);
			
			p.clear();
			
			bool error = false;
			
			for (int i = 0; std::getline(ss, symb, sep); i++){
			
				if (!sel[i])
					continue;
				
				bool errorSymb = false;
				
				try {
					double d = std::stod(symb,&sz);
			
					if (d < -10e+12)
						errorSymb = true;
			
					if (d > 10e+12)
						errorSymb = true;
			
					p.push_back(d);
				
					//cout << symb << ", ";
			
				} catch (std::invalid_argument& e){
					
					p.push_back(NAN);
			
					errorSymb = true;
				
					//cout << symb << ", ";
				}
				
				if (errorSymb){
					cout << "error symb: " << symb << endl;
					error = true;
				}
				
			}
			
			if (error)
				cout << "error line: " << line << endl;
			
			return true;
		}
	
		return false;
	}

	inline void close(){
	
		if (lineReader){
			lineReader->close();
			lineReader.reset();
		}
			
		if (gzipLineReader){
			gzipLineReader->close();
			gzipLineReader.reset();
		}
	}

};


#endif