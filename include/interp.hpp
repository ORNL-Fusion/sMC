#ifndef INTERP_HPP_
#define INTERP_HPP_

#include "constants.hpp"
#include "array2D.hpp"

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
        
        __host__ __device__
            CinterpIndex () { i=0.0; j=0.0; i1=0; i2=0; j1=0; j2=0; stat=0; };
};

__host__ __device__
CinterpIndex get_index ( const REAL rIn, const REAL zIn, const CinterpSpans &spans );

__host__ __device__
REAL bilinear_interp 
    ( const CinterpIndex &index , const array2D<REAL,BCHECK> &data );


#endif
