#ifndef CUDA_WRAP_H_
#define CUDA_WRAP_H_

#include "particle.hpp"
#include "constants.hpp"
#include "eqdsk.hpp"

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
		unsigned int i1, i2, j1, j2;
};


int copy_particles_to_device (std::vector<Cgc_particle> &H_particles);
cu_ptr_pitch copy_2D_to_device ( eqdsk::arr2D_ &data2D, const unsigned int nRow, const unsigned int nCol );
REAL* copy_1D_to_device ( eqdsk::arr1D_ &h_data1D, const unsigned int n );
int cu_test_cuda ( cu_ptrs &d_ptrs, int nRow, int nCol );

#endif

