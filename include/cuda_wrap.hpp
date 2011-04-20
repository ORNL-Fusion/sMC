#ifndef CUDA_WRAP_H_
#define CUDA_WRAP_H_

#include "particle.hpp"
#include "constants.hpp"
#include "eqdsk.hpp"
#include <vector>
#include "array2d.hpp"
#include "interp.hpp"
#include "cuda_classes.hpp"

Ccu_ptr_pitch copy_2D_to_device ( const array2D<REAL,BCHECK> &data2D, 
    const unsigned int nRow, const unsigned int nCol );

REAL* copy_1D_to_device ( const std::vector<REAL> &h_data1D, const unsigned int n );

extern "C" 
int cu_test_cuda ( std::vector<Cgc_particle> &H_particles, 
    const int nRow, const int nCol, const CinterpSpans &spans, const Ceqdsk &eqdsk );

#endif

