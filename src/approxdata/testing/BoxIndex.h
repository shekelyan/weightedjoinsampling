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

#ifndef SHEKELYAN_BOXINDEX_H
#define SHEKELYAN_BOXINDEX_H

#include <approxdata/testing/testing.h>

class ValueIndex{


private:


	// 
	
	inline long ind(long ind, int rel) const{
	
		return (ind << 1)+1+rel; 
	}

public:
	vector<double> v;
	
	
	inline void addVal(double x){
	
		v.push_back(x);
	}
	
	inline void finalise(){
	
		sort( v.begin(), v.end() );
		v.erase( unique( v.begin(), v.end() ), v.end() );
	}
	
	inline long largestIndex() const{
	
		return ind(v.size()-1, +1);
	}
	
	
	inline long getIndex(double val) const{
	
		if (v.size() == 0)
			return -1;
	
		if (val < v[0])
			return ind(0, -1);
		
		if (val == v[0])
			return ind(0, 0);
		
		if (val == v[v.size()-1])
			return ind(v.size()-1, 0);
		
		if (val > v[v.size()-1])
			return ind(v.size()-1, +1);
		
		auto low = std::lower_bound (v.begin(), v.end(), val);  // first >= 
		
		const long lb = low-v.begin();
		
		if ( val == (*low))
			return ind(lb, 0);
		
		if ( val < (*low))
			return ind(lb, -1);
			
		if ( val > (*low))
			return ind(lb, +1);
			
		return -1;
	}


};

class BoxIndex{

protected:

	// for each dimension the query boundary values
	vector< ValueIndex > discLut;
	
	// for each dimension and for each query boundary bucket set of box-indices containing it
	vector< vector< vector<int> >> intersectedLut;
	
	inline void prepare(){
	
		for (int i = 0; i < DIMS; i++){
		
			vector<double> v;
			
			
			ValueIndex vi;
	
			for (auto it = boxes.begin(); it != boxes.end(); it++){
	
				vi.addVal(it->min[i]);
				vi.addVal(it->max[i]);
			}
			
			vi.finalise();
			
			discLut.push_back(vi);
			
			vector<vector<int>> lutDim;
			
			for (long j = 0; j <= vi.largestIndex(); j++){
			
				vector<int> intersectedIndices;
				
				for (long k = 0; k < boxes.size(); k++){
			
					const long minIndex = getIndex( boxes[k].min[i], i );
					const long maxIndex = getIndex( boxes[k].max[i], i );
					
					if (minIndex <= j && maxIndex >= j)
						intersectedIndices.push_back(k);
				}
				
				lutDim.push_back(intersectedIndices);
			}
			
			intersectedLut.push_back(lutDim);	
		}
	}

	
	inline long getIndex(double val, int dim) const{
	
		return discLut[dim].getIndex(val);
	}
	
public:
	
	const int DIMS;
	
	vector<Box>& boxes;
	
	const bool fast;
	
	inline BoxIndex(int dims, vector<Box>& boxesPtr, int mode) : boxes(boxesPtr), fast(mode != 0), DIMS(dims){
		
		if (fast)
			prepare();
			
		reset();
	}
	
	inline void reset(){
		
		for (auto it = boxes.begin(); it != boxes.end(); it++)
			it->count = 0;
	}
	
	inline void print(){
		
	}
	
	inline void count(Point& p){
	
		if (!fast){
			
			for (auto it = boxes.begin(); it != boxes.end(); it++){
		
				if ( it->contains(p) )
					it->count++;
			}
			
			return;
		}
		
		vector<int>* b = NULL;
		
		// find most selective dimension
		for (int i = 0; i < DIMS; i++){
		
			vector<int>* c = &(intersectedLut[i][getIndex(p[i], i)]);
			
			if (b == NULL || (c->size() < b->size()) )
				b = c;
		}
		
		for (auto it = b->begin(); it != b->end(); it++){
			
			const int j = (*it);
			
			if ( boxes[j].contains(p) )
				boxes[j].count++;
		}
	}
	

};







#endif