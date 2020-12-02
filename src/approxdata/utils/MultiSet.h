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

#ifndef SHEKELYAN_MULTISET_H
#define SHEKELYAN_MULTISET_H

#include <approxdata/utils/utils.h>

// Class for multisets over integers


 // E = element type (any type supporting < and == operators), C = count type (any type that can be set to 0 and added)
template<typename E, typename C>
class MultiSet{


 
protected:

	vector<E> unsortedElems;
	vector<C> unsortedCounts;
	
	vector<E> sortedElems;
	vector<C> sortedCounts;
	
	vector<E> tempElems;
	vector<C> tempCounts;
	
	bool unsortedIsSorted = true;
	
	
	template<typename A,typename B>
	inline void clear(vector<A>& v1, vector<B>& v2){
	
		v1.clear();
		v2.clear();
	}
	
	template<typename A,typename B>
	inline void shrink(vector<A>& v1, vector<B>& v2){
	
		if (memorySaving){
			v1.shrink_to_fit();
			v2.shrink_to_fit();
		}
	}
	
public:

	virtual ~MultiSet(){};
	
	vector<E> elems;
	vector<C> counts;
	
	bool memorySaving = false;
	
	long batchSize = 1000000;
	
	C total = 0;
	
	inline MultiSet(long size=100, long batchSize=-1, bool memSave = true){
	
		this->memorySaving = memSave;
		
		if (batchSize == -1)
			batchSize = size/4;
		
		unsortedElems.reserve(batchSize);
		
		total = 0;
		
		if (!memorySaving){
			
			sortedElems.reserve(batchSize);
			sortedCounts.reserve(batchSize);
			
			elems.reserve(size);
			counts.reserve(size);
			
			tempElems.reserve(size);
			tempCounts.reserve(size);
		}
	}
	
	inline long getDistinctNum() const{
	
		return elems.size();
	}
	
	inline C getTotalNum() const{
	
		return total;
	}
	
	inline void print() const{
	
		for (long k = 0; k < elems.size(); k++){
		
			cout << elems.at(k) << " " << counts.at(k) << endl;
		}
	}
	
	inline void mergeDuplicates(
	const vector<E>& srcElems, const vector<C>& srcCounts, vector<E>& dstElems, vector<C>& dstCounts){
		
		clear(dstElems, dstCounts);
		
		if (srcElems.size() == 0)
			return;
		
		dstElems.reserve(srcElems.size() );
		dstCounts.reserve(srcElems.size() );
		
		const bool COUNTS = srcCounts.size() > 0;
		
		E elem = srcElems[0];
		C count = COUNTS ? srcCounts[0] : 1;
		
		assert (count > 0);
		
		for (long j = 1; j < srcElems.size(); j++){
	
			const E e = srcElems[j];
			const C c = COUNTS ? srcCounts[j] : 1;
			
			assert (c > 0);
			
			if (e == elem){
				
				count += c;
			
			} else {
			
				dstElems.push_back(elem);
				dstCounts.push_back(count);
				
				elem = e;
				count = c;
			}
		}
		
		dstElems.push_back(elem);
		dstCounts.push_back(count);
		
	}
	
	inline void processUnsorted(){
	
		if (unsortedElems.size() == 0)
			return;
		
		if (!unsortedIsSorted){
			
			assert (unsortedCounts.size() == 0);
			std::sort(unsortedElems.begin(), unsortedElems.end() );
		}
		
		mergeDuplicates(unsortedElems, unsortedCounts, sortedElems, sortedCounts);
		
		clear(unsortedElems, unsortedCounts);
		
		unsortedIsSorted = true;
	}
	
	
	
	
	inline void processSorted(){
	
		if (sortedElems.size() == 0)
			return;
		
		auto elemsIt = elems.begin();
		auto countsIt = counts.begin();
		auto elemsEnd = elems.end();
		
		auto sortedElemsIt = sortedElems.begin();
		auto sortedCountsIt = sortedCounts.begin();
		auto sortedElemsEnd = sortedElems.end();
		
		clear(tempElems, tempCounts);
		
		tempElems.reserve( elems.size()+sortedElems.size() );
		tempCounts.reserve( counts.size()+sortedCounts.size() );
		
		while (sortedElemsIt != sortedElemsEnd){
			
			const E oldElem = (*sortedElemsIt);
			const C oldCount = (*sortedCountsIt);
				
			if (elemsIt == elemsEnd){
			
				tempElems.push_back(oldElem);
				tempCounts.push_back(oldCount);
				
				sortedElemsIt++;
				sortedCountsIt++;
				
			} else {
			
				const E newElem = (*elemsIt);
				const C newCount = (*countsIt);
			
				if ( oldElem < newElem){
			
					// take from old
				
					tempElems.push_back(oldElem);
					tempCounts.push_back(oldCount);
				
					sortedElemsIt++;
					sortedCountsIt++;
				
				
				} else if (oldElem == newElem){
			
					// take from both
			
					tempElems.push_back(newElem);
					tempCounts.push_back(oldCount+newCount);
				
					elemsIt++;
					sortedElemsIt++;
					sortedCountsIt++;
					countsIt++;
				
				} else {
			
					// take from new
				
					const E newElem = (*elemsIt);
					const C newCount = *countsIt;
			
					tempElems.push_back(newElem);
					tempCounts.push_back(newCount);
				
					elemsIt++;
					countsIt++;
				}
			}
		}
		
		while (elemsIt != elemsEnd){
		
			const E newElem = (*elemsIt);
			const C newCount = (*countsIt);
			
			tempElems.push_back(newElem);
			tempCounts.push_back(newCount);
				
			elemsIt++;
			countsIt++;
		}
		
		clear(elems, counts);
		clear(sortedElems, sortedCounts);
		
		elems = tempElems;
		counts = tempCounts;
		
		clear(tempElems, tempCounts);
		
		shrink(elems, counts);
		shrink(tempElems, tempCounts);
	}
	
	inline void addSorted(E e){
		
		total++;
		
		unsortedElems.push_back(e);
		
		if (unsortedElems.size() >= batchSize){
		
			processUnsorted();
			processSorted();
		}
	}
	
	inline void addSorted(E e, C c=1){
		
		total += c;
		
		if ( (c != 1) || (unsortedCounts.size() > 0) ){
			
			for (long j = unsortedCounts.size(); j < unsortedElems.size(); j++)
				unsortedCounts.push_back(1);
			
			unsortedCounts.push_back(c);
		}
		
		unsortedElems.push_back(e);
		
		if (unsortedElems.size() >= batchSize){
		
			processUnsorted();
			processSorted();
		}
	}
	
	
	inline void finalize(){
	
		processUnsorted();
		processSorted();
	}
	
	inline void add(E e, C c=1){
		
		if (unsortedIsSorted)
			unsortedIsSorted = false;
		
		addSorted(e, c);
	}

};





#endif