#ifndef INTERP_HPP_
#define INTERP_HPP_

#include "constants.hpp"
#include "array2D.hpp"
#include <cmath>

#ifdef __CUDA_ARCH__
#define PRINT cuPrintf 
#else
#define PRINT printf
#endif

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
#if DEBUGLEVEL >= 4
		PRINT("\t\t\t\t%s line: %i\n",__FILE__,__LINE__);
		PRINT("\t\t\t\tmIn: %f\n",mIn);
		PRINT("\t\t\t\tnIn: %f\n",nIn);
		PRINT("\t\t\t\tindex.m: %f\n",index.m);
		PRINT("\t\t\t\tindex.n: %f\n",index.n);
		PRINT("\t\t\t\tindex.m1: %i\n",index.m1);
		PRINT("\t\t\t\tindex.m2: %i\n",index.m2);
		PRINT("\t\t\t\tindex.n1: %i\n",index.n1);
		PRINT("\t\t\t\tindex.n2: %i\n",index.n2);
#endif
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

	REAL x = index.n - index.n1;
	REAL y = index.m2 - index.m;

	REAL f00 = data(index.m2,index.n1);
	REAL f10 = data(index.m2,index.n2);
	REAL f01 = data(index.m1,index.n1);
	REAL f11 = data(index.m1,index.n2);

	REAL dataOut = f00*(1-x)*(1-y)+f10*x*(1-y)+f01*(1-x)*y+f11*x*y;

    return dataOut;
}

} // end namespace
#endif
