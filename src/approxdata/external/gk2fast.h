/*
 * Implementation of the GK algorithm
 * Copyright (c) 2013 Lu Wang <coolwanglu@gmail.com>
 */

/*
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/*
shekelyan: changes made to support ordered points
*/

//#ifndef GK_H__ // shekelyan: REMOVED
//#define GK_H__ // shekelyan: REMOVED

class FastEntry{

	double* ptr;
	unsigned int g,delta;
};


class FastEntryComparator : public std::binary_function<const FastEntry*, const FastEntry*, bool>{

	int dims;
	int mainDim;

	inline bool operator()(const FastEntry*& a, const FastEntry*& b) const{
		
		return less(a,b);
	}	

};
bool cmp(const FastEntry* a, const FastEntry* b) {
    return ...;
}

class FastEntries : public std::binary_function<const FastEntry*, const FastEntry*, bool>{

	vector<FastEntry*> available;
	
	multiset<FastEntry*> set;

	vector<FastEntry> entries;
	vector<double> doubles;
	
	inline FastEntries(long n, int dims) : doubles(n*dims), entries(n){
		
		double* dbeg = doubles.data();
		double* dend = dbeg+doubles.size();
		
		long j = 0;
		
		for (double* d = dbeg; d != dend; d += dims){
			entries[j].ptr = d;
			
			available.push_back(&entries[j]);
			
			j++;
		}
		
		auto cmp = [](int a, int b) { return ... };
		std::set<int, decltype(cmp)> s(cmp);
	}
	
	inline void insertBefore(vector<double>& vec, unsigned int g, unsigned int delta){
	
		FastEntry* f = available.back();
		available.pop_back();
		
		for (int i = 0; i < vec.size(); i++)
			f->ptr[i] = vec[i];
		
		f->g = g;
		f->delta = delta;
		
		set.insert(f);
	}
	
	inline void remove(FastEntry* g){
	
		available.push_back(g);
		
		set.erase(g);
	}
};

/*
class FastEntries{

	vector<double> data;

	vector<FastEntry> vec;
	
	long size = 0;
	
	inline void remove(long k){
		
		vec[k].ptr = vec[k-1].ptr;
		vec[k].removed = true;
	}
	
	inline void set(long a, long b){
	
		vec[a].ptr = vec[b].ptr;
		vec[a].g = vec[b].g;
		vec[a].delta = vec[b].delta;
	}
	
	inline void insertBefore(long k, double* ptr, unsigned int g, unsigned int delta){
	
		long nextEmpty = k;
		
		while (!vec[nextEmpty].removed)
			nextEmpty++;
		
		for (long j = nextEmpty; j > k; j--)
			set(j, j-1);
		
		vec[k].ptr = ptr;
		vec[k].g = g;
		vec[k].delta = delta;
	}
	
	inline void compress(){
	
		
		long toRemove = 0;
		
		for (long j = 0; j < vec.size(); j++)
			if (vec[j].removed)
				toRemove++;
		
		
	
		long j = 0;
		
		size = 0;
		
		
		while (j < vec.size() ){
		
			if (vec[j].removed){
				j++;
				
				
			}
				
			set(size, j);
			size++;
		}
		
		
	}
	
	
	
};*/


#ifndef GK2_H__ //shekelyan: ADDED
#define GK2_H__	//shekelyan: ADDED

#include <map>
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <unordered_map>
#include <queue>
#include <cassert>

// template<class IntType> // shekelyan: REMOVED
template<class IntType, class IntComp> // shekelyan: ADDED
//class GK // shekelyan: REMOVED
class GK2 // shekelyan: ADDED
{
public:
    long long last_n;
    double EPS;
    IntType min_v; // #ADD
    IntType max_v; // #MODIFY

    class entry{
    public:
        entry () { }
        entry (unsigned int _g, unsigned int _d) : g(_g), delta(_d) { }
        unsigned int g,delta;
    };

	// typedef std::multimap<IntType, entry> map_t; // shekelyan: REMOVED
    typedef std::multimap<IntType, entry, IntComp> map_t; // shekelyan: ADDED
    map_t entries_map;
    unsigned int max_gd;

    class SummaryEntry
    {
    public:
        IntType v;
        // long long g, delta; // shekelyan: REMOVED
        double g; // shekelyan: ADDED
        long long delta; // shekelyan: ADDED
    };
    std::vector<SummaryEntry> summary;

	// GK(double eps, uint64_t max_v = std::numeric_limits<IntType>::max()) // shekelyan: REMOVED
    GK2(double eps, const IntType& min_v, const IntType& max_v, const IntComp& comp) // shekelyan: ADDED
        :last_n(0)
        ,EPS(eps) 
        ,min_v(min_v)// shekelyan: ADDED
        ,max_v(max_v)
        ,entries_map(comp)// shekelyan: ADDED
    {
        entry e(1,0);
        entries_map.insert(std::make_pair(max_v, e));
    }

    void finalize () {
        summary.clear();
        summary.reserve(entries_map.size());
        SummaryEntry se;
        se.g = 0;
        max_gd = 0;
        for(auto & p : entries_map)
        {
            const auto & e = p.second;
            unsigned int gd = e.g + e.delta;
            if(gd > max_gd) max_gd = gd;

            se.v = p.first;
            se.g += e.g;
            se.delta = e.delta;
            summary.push_back(se);
        }
        max_gd /= 2;
        
        entries_map.clear(); // shekelyan: ADDED
        entries_map.shrink_to_fit(); // shekelyan: ADDED
    }

	// IntType query_for_value(double rank) const // shekelyan: REMOVED
    const IntType query_for_value(double rank) const // shekelyan: ADDED
    {
        SummaryEntry se;
        se.g = rank*last_n+max_gd;
        se.delta = 0;
        auto iter = upper_bound(summary.begin(), summary.end(),
                se,
                [](const SummaryEntry & e1, const SummaryEntry & e2)->bool{
                    return (e1.g + e1.delta) < (e2.g + e2.delta);
                });

        if(iter == summary.begin())
        {
        	// return 0; // shekelyan: REMOVED
            return min_v; // shekelyan: ADDED
        }

        return (--iter)->v;
    } 
    // record when an item might be removed
    class ThresholdEntry
    {
    public:
        unsigned int threshold; // the item could be removed only when the global threshold is no less than this value
        typename map_t::iterator map_iter;

        bool operator > (const ThresholdEntry & te) const {
            return threshold > te.threshold;
        }
    };

    std::priority_queue<ThresholdEntry, std::vector<ThresholdEntry>, std::greater<ThresholdEntry>> compress_heap;

	// void feed(IntType v) { // shekelyan: REMOVED
    void feed(const IntType& v) { // shekelyan: ADDED
        ++last_n;
		
		// if((uint64_t)v == max_v) return; // shekelyan: REMOVED
        if(v == max_v) return; // shekelyan: ADDED
		
        const auto iter = entries_map.upper_bound(v);
        //iter cannot be end()
		
        entry & ecur = iter->second;
        const unsigned int threshold = (unsigned int)std::floor(EPS * last_n * 2);
        unsigned int tmp = ecur.g + ecur.delta;
        if(tmp < threshold)
        {
            //no need to insert
            ++(ecur.g);

            // now the heap entries for recur and previous one (if there is) are affected
            // the threshold should have been increased by 1
            // we will fix this when trying to remove one
        }
        else
        {
            auto iter2 = entries_map.insert(iter, std::make_pair(v, entry(1, tmp-1)));

            // note that the heap entry for the previous one (if there is) is not affected!
            compress_heap.push({tmp + 1, iter2});

            // try to remove one
            while(true)
            {
                auto top_entry = compress_heap.top();
                if(top_entry.threshold > threshold)
                {
                    // all threshold in the heap are no less than the real threshold of the entry
                    // so we cannnot remove any tuple right now
                    break;
                }
				
                compress_heap.pop();
                // check if it is true
                auto map_iter2 = top_entry.map_iter; 
                auto map_iter1 = map_iter2++; 
                assert(map_iter2 != entries_map.end());

                auto & e1 = map_iter1->second;
                auto & e2 = map_iter2->second;

                auto real_threshold = e1.g + e2.g + e2.delta;

                if(real_threshold <= threshold) {
                    e2.g += e1.g;
                    entries_map.erase(map_iter1);
                    break;
                }
                else
                {
                    top_entry.threshold = real_threshold;
                    compress_heap.push(top_entry);
                }
            }
        }
    }
};

#endif // GK_H__