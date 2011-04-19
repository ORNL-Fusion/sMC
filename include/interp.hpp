#ifndef INTERP_HPP_
#define INTERP_HPP_

#include "constants.hpp"
#include "array2D.hpp"
#include <cmath>

class CinterpSpans {

    public:

        REAL rfront, rback, zfront, zback;
        unsigned int rsize, zsize;

        CinterpSpans() : rfront(0), rback(0), rsize(0), zfront(0), zback(0), zsize(0) {};
        CinterpSpans(
            const REAL _rfront, const REAL _rback, const unsigned int _rsize, 
            const REAL _zfront, const REAL _zback, const unsigned int _zsize ) {
                rfront = _rfront;
                rback = _rback;
                zfront = _zfront;
                zback = _zback;
                rsize = _rsize;
                zsize = _zsize;
            }  
};

class CinterpIndex {

	public: 
		REAL i, j;
		int i1, i2, j1, j2;
        unsigned int stat;

#ifdef __CUDACC__        
        __host__ __device__
#endif
            CinterpIndex () { i=0.0; j=0.0; i1=0; i2=0; j1=0; j2=0; stat=0; };
};

namespace {

#ifdef __CUDACC__
__host__ __device__
#endif
CinterpIndex get_index ( const REAL rIn, const REAL zIn, const CinterpSpans &spans ) {

    CinterpIndex index;

	index.j = (rIn - spans.rfront) / ( spans.rback - spans.rfront ) * (spans.rsize-1);
	index.i = (zIn - spans.zfront) / ( spans.zback - spans.zfront ) * (spans.zsize-1);

	index.i1 = floor(index.i);
	index.i2 = ceil(index.i);
	index.j1 = floor(index.j);
	index.j2 = ceil(index.j);

    // Check if particle is off grid	
    if( index.i1<0 || index.i2>=(spans.zsize-1) || index.j1<0 || index.j2>=(spans.rsize-1) ) {
        index.stat += 1;
    }

	return index;
}

// bi-linear interpolation
// see wikipedia ;)

#ifdef __CUDACC__
__host__ __device__
#endif
REAL bilinear_interp 
    ( const CinterpIndex &index , const array2D<REAL,BCHECK> &data ) {

	REAL f11 = data(index.i1,index.j1);
	REAL f21 = data(index.i2,index.j1);
	REAL f12 = data(index.i1,index.j2);
	REAL f22 = data(index.i2,index.j2);

	// (x2-x1)(y2-y1) == 1 since i'm using indices

	REAL dataOut = f11 * (index.i2-index.i)*(index.j2-index.j)
			+ f21 * (index.i-index.i1)*(index.j2-index.j)
			+ f12 * (index.i2-index.i)*(index.j-index.j1)
			+ f22 * (index.i-index.i1)*(index.j-index.j1); 

    return dataOut;
}

} // end namespace
#endif
