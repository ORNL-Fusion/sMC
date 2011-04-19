#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>
#include "eqdsk.hpp"
#include "constants.hpp"
#include "cuda_classes.hpp"
#include "array2D.hpp"
#include "gc_eqom.hpp"
#include "rkf.hpp"
#include "interp.hpp"
#include <vector>
#include "particle.hpp"
#include "cuPrintf.cu"

__global__ void cu_move_particles (Cgc_particle *const particles, 
                    const cu_ptrs dataPtrs, const CinterpSpans spans) {

    Ctextures textures;
    
    textures.bmag.ptr = dataPtrs.bmag.ptr;
    textures.bmag.pitch = dataPtrs.bmag.pitch;

    textures.bDotGradB.ptr   = dataPtrs.bDotGradB.ptr;
    textures.bDotGradB.pitch = dataPtrs.bDotGradB.pitch;

    textures.bCurv_r.ptr     = dataPtrs.bCurv_r.ptr;
    textures.bCurv_r.pitch   = dataPtrs.bCurv_r.pitch;
    textures.bCurv_p.ptr     = dataPtrs.bCurv_p.ptr;
    textures.bCurv_p.pitch   = dataPtrs.bCurv_p.pitch;
    textures.bCurv_z.ptr     = dataPtrs.bCurv_z.ptr;
    textures.bCurv_z.pitch   = dataPtrs.bCurv_z.pitch;

    textures.bGrad_r.ptr     = dataPtrs.bGrad_r.ptr;
    textures.bGrad_r.pitch   = dataPtrs.bGrad_r.pitch;
    textures.bGrad_p.ptr     = dataPtrs.bGrad_p.ptr;
    textures.bGrad_p.pitch   = dataPtrs.bGrad_p.pitch;
    textures.bGrad_z.ptr     = dataPtrs.bGrad_z.ptr;
    textures.bGrad_z.pitch   = dataPtrs.bGrad_z.pitch;

    textures.b_r.ptr     = dataPtrs.b_r.ptr;
    textures.b_r.pitch   = dataPtrs.b_r.pitch;
    textures.b_p.ptr     = dataPtrs.b_p.ptr;
    textures.b_p.pitch   = dataPtrs.b_p.pitch;
    textures.b_z.ptr     = dataPtrs.b_z.ptr;
    textures.b_z.pitch   = dataPtrs.b_z.pitch;

    int stat;
	for(unsigned int p=0;p<10;p++) {

		if(!particles[p].status) {

			//cuPrintf("Particle No. %i\n",p) ;
			//cuPrintf("eV: %f\n",particles[p].energy_eV);

            stat = move_particle ( particles[p], textures, spans, p );
			
		} // end if(!particle[p].status)
	} // end for(p)

}

__global__ void check2Dcpy ( const cu_ptr_pitch data_, 
				const unsigned int nRow, const unsigned int nCol ) {

    array2D<REAL,BCHECK> data2D;
    data2D.ptr = data_.ptr;
    data2D.pitch = data_.pitch;

	for (int r=0;r<nRow;++r) {
			REAL *row = (REAL*)((char*)data_.ptr + r*data_.pitch);
			for (int c=0;c<nCol;++c) {
					REAL element = row[c];
                    cuPrintf("%i %i %f\n", r, c, element);
                    cuPrintf("%i %i %f\n", r, c, data2D(r,c));
			}
	}
}

__global__ void check1Dcpy ( REAL *data1D, const unsigned int n ) {

	for(int i=0;i<n;++i) {
			//cuPrintf("%i %f\n", i, data1D[i]);
	}
}

cu_ptr_pitch copy_2D_to_device
( const array2D<REAL,BCHECK> &data2D, const unsigned int M, const unsigned int N ) {

    cu_ptr_pitch out;
	size_t size = N * sizeof(REAL);

	cudaMallocPitch ( (void**)&out.ptr, &out.pitch, size, M );
	cudaMemcpy2D ( out.ptr, out.pitch, &data2D(0,0), 
					size, size, M, cudaMemcpyHostToDevice );

	return out;
}

REAL* copy_1D_to_device 
( const std::vector<REAL> &h_data1D, const unsigned int n ) {

	REAL *d_data1D;
	size_t size = n * sizeof(REAL);
	cudaMalloc ( (void**)&d_data1D, size );
	cudaMemcpy ( d_data1D, &h_data1D[0], size, cudaMemcpyHostToDevice);

	return d_data1D;
}

extern "C" 
int cu_test_cuda 
( std::vector<Cgc_particle> &H_particles, 
    const int nRow, const int nCol, const CinterpSpans &spans, const Ceqdsk &eqdsk ) {

    cudaPrintfInit();

    // Copy particles to device
    thrust::device_vector<Cgc_particle> D_particles ( H_particles.begin(), H_particles.end() );
    Cgc_particle * D_particles_ptr = thrust::raw_pointer_cast(&D_particles[0]);

    // Copy background data to device
    cu_ptrs D_ptrs;

	D_ptrs.r = copy_1D_to_device (eqdsk.r,eqdsk.nCol);
	D_ptrs.z = copy_1D_to_device (eqdsk.z,eqdsk.nRow);

	D_ptrs.bmag = copy_2D_to_device (eqdsk.bmag,eqdsk.nRow,eqdsk.nCol);
    
	D_ptrs.b_r = copy_2D_to_device (eqdsk.br,eqdsk.nRow,eqdsk.nCol);
	D_ptrs.b_p = copy_2D_to_device (eqdsk.bp,eqdsk.nRow,eqdsk.nCol);
	D_ptrs.b_z = copy_2D_to_device (eqdsk.bz,eqdsk.nRow,eqdsk.nCol);

	D_ptrs.bCurv_r = copy_2D_to_device (eqdsk.bCurvature_r,eqdsk.nRow,eqdsk.nCol);
	D_ptrs.bCurv_p = copy_2D_to_device (eqdsk.bCurvature_p,eqdsk.nRow,eqdsk.nCol);
	D_ptrs.bCurv_z = copy_2D_to_device (eqdsk.bCurvature_z,eqdsk.nRow,eqdsk.nCol);

	D_ptrs.bGrad_r = copy_2D_to_device (eqdsk.bGradient_r,eqdsk.nRow,eqdsk.nCol);
	D_ptrs.bGrad_p = copy_2D_to_device (eqdsk.bGradient_p,eqdsk.nRow,eqdsk.nCol);
	D_ptrs.bGrad_z = copy_2D_to_device (eqdsk.bGradient_z,eqdsk.nRow,eqdsk.nCol);

    std::cout << "Testing 1D memcopy ..." << std::endl;
	check1Dcpy<<<1,1>>>( D_ptrs.z, nRow );
    std::cout << "Testing 2D memcopy ..." << std::endl;
	check2Dcpy<<<1,1>>>( D_ptrs.bmag, nRow, nCol );

    std::cout << "Moving particles ..." << std::endl;
	cu_move_particles<<<1,1>>>( D_particles_ptr, D_ptrs, spans );

    cudaPrintfDisplay (stdout, true);
    cudaPrintfEnd();

    return 0;
}

