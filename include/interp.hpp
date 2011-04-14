#ifndef INTERP_HPP_
#define INTERP_HPP_

#include "constants.hpp"
//#include "boost/multi_array.hpp"
#include "array2D.hpp"

class interpSpans {

    public:

        REAL rfront, rback, zfront, zback;
        unsigned int rsize, zsize;

        interpSpans() : rfront(0), rback(0), rsize(0), zfront(0), zback(0), zsize(0) {};
        interpSpans(
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

class interpIndex {

	public: 
		REAL i, j;
		int i1, i2, j1, j2;
        unsigned int stat;
        
        __host__ __device__
            interpIndex () { i=0.0; j=0.0; i1=0; i2=0; j1=0; j2=0; stat=0; };
};

__host__ __device__
interpIndex get_index ( const REAL rIn, const REAL zIn, const interpSpans &spans );

__host__ __device__
REAL bilinear_interp 
    ( const interpIndex &index , const array2D<REAL,BCHECK> &data );


#endif
