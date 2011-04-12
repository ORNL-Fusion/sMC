#ifndef CUDA_WRAP_H_
#define CUDA_WRAP_H_

#include "particle.hpp"
#include "constants.hpp"
#include "eqdsk.hpp"
#include <vector>
//#include "boost/multi_array.hpp"
#include "array2d.hpp"

struct cu_ptr_pitch {
    REAL *ptr;
    size_t pitch;
};

struct cu_ptrs {

    REAL *r, *z;
    cu_ptr_pitch 
        bmag, b_r, b_p, b_z, 
        bCurv_r, bCurv_p, bCurv_z,
        bGrad_r, bGrad_p, bGrad_z,
        bDotGradB;
};

class cu_interpIndex {

	public: 
		REAL i, j;
		int i1, i2, j1, j2;
		unsigned int stat;
};


int 
copy_particles_to_device (std::vector<Cgc_particle> &H_particles);

cu_ptr_pitch 
copy_2D_to_device ( array2D<REAL,BCHECK> &data2D, const unsigned int nRow, const unsigned int nCol );

REAL* 
copy_1D_to_device ( std::vector<REAL> &h_data1D, const unsigned int n );

int 
cu_test_cuda ( const cu_ptrs &d_ptrs, const int nRow, const int nCol );

#endif

