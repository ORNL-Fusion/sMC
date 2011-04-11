#ifndef INTERP_HPP_
#define INTERP_HPP_

#include "constants.hpp"
#include "boost/multi_array.hpp"

class interpIndex {

	public: 
		REAL i, j;
		int i1, i2, j1, j2;
        unsigned int stat;
        
        __host__ __device__
            interpIndex () { i=0.0; j=0.0; i1=0; i2=0; j1=0; j2=0; stat=0; };
};

__host__ __device__
interpIndex get_index ( const REAL rIn, const REAL zIn, 
	  const REAL rfront, const REAL rback, const unsigned int rsize,
   	  const REAL zfront, const REAL zback, const unsigned int zsize );

__host__ __device__
REAL bilinear_interp 
    ( const interpIndex &index , const REAL **data );


#endif
