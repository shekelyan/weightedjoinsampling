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

#ifndef SHEKELYAN_ZCURVECOORDS_H
#define SHEKELYAN_ZCURVECOORDS_H

#include <approxdata/utils/utils.h>

// higher significance => lower level of detail
// msb stores coordinate of lod==0
// lsb stores coordinate of lod==lodMax

class ZCurveCoords{

	public:
		int lodMax = 0;
		
		const int DIMS;
		
		vector<long> dimMasks;
		vector<int> dimMaskBits;
		
		vector<long> dimCellNum;
		vector<long> lookupInterleave;
		
		inline ZCurveCoords(int dims, int lodMax=62) : DIMS(dims){
		
			init(lodMax);
		}
		
		inline int getLod() const{
		
			return lodMax;
		}
		
		inline int getLod(int i) const{
		
			return dimMaskBits[i];
		}
		
		inline void init(int newLod){
			
			lodMax = newLod;
			
			vector<long> splitsLeft;
			
			for (int i = 0; i < DIMS; i++)
				splitsLeft.push_back( std::numeric_limits<long>::max() );
			
			{
				long dimSkipped = 0;
			
				dimMasks.clear();
			
				for (int i = 0; i < DIMS; i++)
					dimMasks.push_back(0);
			
				for (int lod = 1; lod <= lodMax; lod++){
			
					if (lod > 1)
						for (int k = 0; k < DIMS; k++)
							dimMasks[k] <<= 1;
				
					long j = -1;
				
					while (true){
				
						j = (lod+dimSkipped) % DIMS;
					
						if (splitsLeft[j] > 0){
					
							splitsLeft[j]--;
							break;
						
						} else {
					
							dimSkipped++;
						}
					}
				
					dimMasks[j] |= 1L;
				}
			}
			
			{
				dimMaskBits.clear();
			
				for (int i = 0; i < DIMS; i++)
					dimMaskBits.push_back(UtilBits::countBits(dimMasks[i]));
			}
			
			{
				dimCellNum.clear();
			
				for (int i = 0; i < DIMS; i++)
					dimCellNum.push_back(1L << dimMaskBits[i]);
			}
			
			{
				lookupInterleave.clear();
				
				long mask = 0;
				
				for (int k = 0; k < 8; k++)
					mask |= 1L << (k*DIMS);
				
				for (long j = 0; j < 256; j++){
					lookupInterleave.push_back( UtilBits::spreadBitsOverMask(j, mask) );
				}
			}
			
		}
		
		inline long getDimCoordWithMask(long x, long mask) const{
		
			return UtilBits::collectBitsFromMask(x, mask);
		}
		
		
		inline long getCellsPerDim(int i) const{
		
			return dimCellNum[i];
		}
		
		inline bool contains(long imin, long imax, long x) const{
		
			if (imin == -1)
				return false;
				
			if (imax == -1)
				return false;
			
			if (x == -1)
				return false;
			
			for (int i = 0; i < DIMS; i++){
			
				const long m = dimMasks[i];
				
				const long xm = x & m;
				
				if ( xm < (imin & m) )
					return false;
				
				if ( xm > (imax & m) )
					return false;
			}
			
			return true;
		}
		
		
		inline long getDimCoord(long x, int dim) const{
		
			return getDimCoordWithMask(x, dimMasks[dim] );
		}
		
		
		inline long sameDimCoord(long x, long y, int dim) const{
		
			const long m = dimMasks[dim];
		
			return (x & m) == (y & m);
		}
		
		
		inline long decAllCoords(long x) const{
		
			long ret = x;
			
			for (int i = 0; i < DIMS; i++){
				
				const long m = dimMasks[i];
				
				if ( (ret & m) == 0)
					return -1;
				
				ret = UtilBits::decrementInMask(ret, m);
			}
		
			return ret;
		}
		
		inline long incAllCoords(long x) const{
		
			long ret = x;
			
			for (int i = 0; i < DIMS; i++){
				
				const long m = dimMasks[i];
				
				if ( (ret & m) == m)
					return -1;
				
				ret = UtilBits::incrementInMask(ret, m);
			}
			
			return ret;
		}
		
		inline long setDimCoord(long ret, long bin, int dim) const {
		
			if (ret == -1)
				return -1;
		
			if (bin < 0)
				return -1;
			
			if (DIMS == 1){
				
				return ret |bin;
				
			} else if (DIMS == 2){ // 2d coordinates <= 2^64 
			
				//assert(bin <= (1L << 63) );
				
				return ret |(
				(lookupInterleave[bin & 255L] ) |
				(lookupInterleave[(bin >> 8) & 255L] << (1*(8+(DIMS-1)*8)  ) ) |
				(lookupInterleave[(bin >> 16) & 255L] << (2*(8+(DIMS-1)*8) ) ) |
				(lookupInterleave[(bin >> 24) & 255L] << (3*(8+(DIMS-1)*8) ) ) 
				) << ((lodMax-dim) % DIMS);
				
			} else if (DIMS == 3){
			
				//assert(bin < (1L << 32) );
			
				return ret |(
				(lookupInterleave[bin & 255L] ) |
				(lookupInterleave[(bin >> 8) & 255L] << (1*(8+(DIMS-1)*8)  ) ) |
				(lookupInterleave[(bin >> 16) & 255L] << (2*(8+(DIMS-1)*8) ) )
				) << ((lodMax-dim) % DIMS);
				
			} else if (DIMS < 8){ // 4d coordinates <= 2^16
				
				//assert(bin < (1L << 16) );
				
				return ret |(
				(lookupInterleave[bin & 255L] ) |
				(lookupInterleave[(bin >> 8) & 255L] << (1*(8+(DIMS-1)*8)  ) )
				) << ((lodMax-dim) % DIMS);
				
			}  else { //if (DIMS >= 8){ // 8d coordinates <= 2^8
			
				//assert(bin < (1L << 8) );
				
				return ret |(lookupInterleave[bin & 255L]) << ((lodMax-dim) % DIMS);	
			}
			
			assert (false);
			return -1;
		}
		
		
		inline long getIndex(const vector<long>& p) const{
		
			long ret = 0;
			
			for (int i = 0; i < DIMS; i++)
				ret = setDimCoord(ret, p[i], i);
			
			return ret;
			
		}
		
};

#endif