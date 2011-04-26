#ifndef INTERP_HPP_
#define INTERP_HPP_

#include "constants.hpp"
#include "array2D.hpp"
#include <cmath>

//#ifdef __CUDACC__
//#include "C/src/simplePrintf/cuPrintf.cu"
//#endif

class CinterpSpans {

    public:

        REAL mfront, mback, nfront, nback;
        unsigned int NCOL, NROW;

        CinterpSpans() : mfront(0), mback(0), NROW(0), nfront(0), nback(0), NCOL(0) {};
        CinterpSpans(
            const REAL _mfront, const REAL _mback, const unsigned int _NROW, 
            const REAL _nfront, const REAL _nback, const unsigned int _NCOL ) {
                mfront = _mfront;
                mback = _mback;
                nfront = _nfront;
                nback = _nback;
                NCOL = _NCOL;
                NROW = _NROW;
            }  
};

class CinterpIndex {

	public: 
		REAL m, n;
		int m1, m2, n1, n2;
        unsigned int stat;

#ifdef __CUDACC__        
        __host__ __device__
#endif
            CinterpIndex () { m=0.0; n=0.0; m1=0; m2=0; n1=0; n2=0; stat=0; };
};

namespace {

#ifdef __CUDACC__
__host__ __device__
#endif
CinterpIndex get_index ( const REAL mIn, const REAL nIn, const CinterpSpans &spans ) {

    CinterpIndex index;

	index.m = (mIn - spans.mfront) / ( spans.mback - spans.mfront ) * (spans.NROW-1);
	index.n = (nIn - spans.nfront) / ( spans.nback - spans.nfront ) * (spans.NCOL-1);

	// Fudge if particle identically @ grid point
	if(index.m-floor(index.m)==0.0)
			index.m+=0.0001;
	if(index.n-floor(index.n)==0.0)
			index.n+=0.0001;

	index.m1 = floor(index.m);
	index.m2 = ceil(index.m);
	index.n1 = floor(index.n);
	index.n2 = ceil(index.n);

    // Check if particle is off grid	
    if( index.m1<0 || index.m2>=(spans.NROW-1) || index.n1<0 || index.n2>=(spans.NCOL-1) ) {
        index.stat+=1;
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

	REAL f11 = data(index.m1,index.n1);
	REAL f21 = data(index.m2,index.n1);
	REAL f12 = data(index.m1,index.n2);
	REAL f22 = data(index.m2,index.n2);

//#ifdef __CUDA_ARCH__
//	cuPrintf("f11: %f, f21: %f, f12: %f, f22: %f\n", f11,f21,f12,f22);
//#endif

	// (x2-x1)(y2-y1) == 1 since i'm using indices

	REAL dataOut = f11 * (index.m2-index.m)*(index.n2-index.n)
			+ f21 * (index.m-index.m1)*(index.n2-index.n)
			+ f12 * (index.m2-index.m)*(index.n-index.n1)
			+ f22 * (index.m-index.m1)*(index.n-index.n1); 

    return dataOut;
}

} // end namespace
#endif
