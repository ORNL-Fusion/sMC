#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>
#include "rk_gc_particle.hpp"
#include "eqdsk.hpp"
#include "constants.hpp"
#include "cuda_wrap.h"
#include "cuPrintf.cu"

using namespace std;
/*
__device__ cu_interpIndex cu_get_index ( const REAL rIn, const REAL zIn, 
	  const REAL rfront, const REAL rback, const unsigned int rsize,
   	  const REAL zfront, const REAL zback, const unsigned int zsize ) {

	cu_interpIndex index;

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

__device__ REAL cu_bilinear_interp ( cu_interpIndex index, REAL *data ) {

		REAL dataOut;

		return dataOut;
}
*/

__device__ Crk cu_vGC () {

		Crk vGC;

		return vGC;
}

__global__ void cu_push () {
}

__global__ void check2Dcpy ( REAL *data2D, 
				size_t pitch, const unsigned int nRow, const unsigned int nCol ) {

	for (int r=0;r<nRow;++r) {
			REAL *row = (REAL*)((char*)data2D + r*pitch);
			for (int c=0;c<nCol;++c) {
					REAL element = row[c];
                    cuPrintf("%i %i %f\n", r, c, element);
			}
	}
}

__global__ void check1Dcpy ( REAL *data1D, const unsigned int n ) {

	for(int i=0;i<n;++i) {
			cuPrintf("%i %f\n", i, data1D[i]);
	}
}

int copy_particles_to_device (vector<Cgc_particle> &H_particles) {

    cout << "First CUDA call :)" << endl;

    cout << "\tCopying particle list from HOST -> DEVICE ... ";
    thrust::device_vector<Cgc_particle> D_particles ( H_particles.begin(), H_particles.end() );
    cout << "DONE" << endl;

    //cout << "\tCopying particle list DEVICE -> HOST ... ";
    //thrust::host_vector<Cgc_particle> H_tmp ( D_particles.begin(), D_particles.end() );
    //cout << "DONE" << endl;

    //float **D_test2D;
    //cudaMalloc( (void**) &D_test2D, eqdsk.nRow * eqdsk.nCol * sizeof(float *));

	return 0;
}

cu_ptr_pitch copy_2D_to_device 
( boost::multi_array<REAL,2> &data2D, const unsigned int M, const unsigned int N ) {

    cu_ptr_pitch out;
	size_t size = N * sizeof(REAL);

	cudaMallocPitch ( (void**)&out.ptr, &out.pitch, size, M );
	cudaMemcpy2D ( out.ptr, out.pitch, &data2D[0][0], 
					size, size, M, cudaMemcpyHostToDevice );

	return out;
}

REAL* copy_1D_to_device 
( std::vector<REAL> &h_data1D, const unsigned int n ) {

	REAL *d_data1D;
	size_t size = n * sizeof(REAL);
	cudaMalloc ( (void**)&d_data1D, size );
	cudaMemcpy ( d_data1D, &h_data1D[0], size, cudaMemcpyHostToDevice);

	return d_data1D;
}

int cu_test_cuda ( const cu_ptrs &d_ptrs, const int nRow, const int nCol ) {

    cudaPrintfInit();

    cout << "Launching testKernel ..." << endl;

	check1Dcpy<<<1,1>>>( d_ptrs.z, nRow );
	check2Dcpy<<<1,1>>>( d_ptrs.bmag.ptr, d_ptrs.bmag.pitch, nRow, nCol );

    cudaPrintfDisplay (stdout, true);
    cudaPrintfEnd();

    return 0;
}

