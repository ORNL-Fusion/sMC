#include "constants.hpp"
//#include "boost/multi_array.hpp"
#include "interp.hpp"
#include "array2D.hpp"

__host__ __device__
interpIndex get_index ( const REAL rIn, const REAL zIn, 
	  const REAL rfront, const REAL rback, const unsigned int rsize,
   	  const REAL zfront, const REAL zback, const unsigned int zsize ) {

    interpIndex index;

	index.j = (rIn - rfront) / ( rback - rfront ) * (rsize-1);
	index.i = (zIn - zfront) / ( zback - zfront ) * (zsize-1);

	index.i1 = floor(index.i);
	index.i2 = ceil(index.i);
	index.j1 = floor(index.j);
	index.j2 = ceil(index.j);

    // Check if particle is off grid	
    if( index.i1<0 || index.i2>=(zsize-1) || index.j1<0 || index.j2>=(rsize-1) ) {
        index.stat += 1;
    }

	return index;
}

// bi-linear interpolation
// see wikipedia ;)
__host__ __device__
REAL bilinear_interp 
    ( const interpIndex &index , const array2D<REAL,BCHECK> &data ) {

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

