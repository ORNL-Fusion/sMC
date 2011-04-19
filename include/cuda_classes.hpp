#ifndef CUDA_CLASSES_
#define CUDA_CLASSES_

#include "constants.hpp"

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

#endif
