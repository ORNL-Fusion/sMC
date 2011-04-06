#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>
#include "particle.hpp"
#include "eqdsk.hpp"
#include "constants.hpp"
#include "cuda_wrap.h"
#include <cuda.h>
#include "cuPrintf.cu"

using namespace std;

/*
__device__ int cu_get_index 
	( const REAL rIn, const REAL zIn, cu_interpIndex &index ) const {

	index.j = (rIn - r.front()) / ( r.back() - r.front() ) * (r.size()-1.0);
	index.i = (zIn - z.front()) / ( z.back() - z.front() ) * (z.size()-1.0);

	index.i1 = floor(index.i);
	index.i2 = ceil(index.i);
	index.j1 = floor(index.j);
	index.j2 = ceil(index.j);

    // Check if particle is off grid	
    if( index.i1<0 || index.i2>=(z.size()-1) ||
        index.j1<0 || index.j2>=(r.size()-1) ) {
        return 1;
    }

	return 0;
}
*/

__global__ void testKernel ( REAL *D_data2D, 
				size_t pitch, const unsigned int nRow, const unsigned int nCol ) {
	int tid;
    tid = blockIdx.x * blockDim.x + threadIdx.x;
    cuPrintf("%d\n", tid);	

	for (int r=0;r<nRow;++r) {
			REAL *row = (REAL*)((char*)D_data2D + r*pitch);
			for (int c=0;c<nCol;++c) {
					REAL element = row[c];
                    cuPrintf("%i %i %d\n", r, c, element);
			}
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

cu_ptr_pitch copy_2D_to_device ( eqdsk::arr2D_ &data2D, const unsigned int nRow, const unsigned int nCol ) {

    cu_ptr_pitch out;

	cudaMallocPitch ( (void**)&out.ptr, &out.pitch, nCol * sizeof(REAL), nRow );
	cudaMemcpy2D ( out.ptr, out.pitch, &data2D[0][0], nCol * sizeof(REAL), nCol, nRow, cudaMemcpyHostToDevice );
	//testKernel<<<100,512>>>(D_data2D,pitch,nRow,nCol);

	return out;
}

REAL* copy_1D_to_device ( eqdsk::arr1D_ &h_data1D, const unsigned int n ) {

	REAL *d_data1D;
	size_t size = n * sizeof(REAL);
	cudaMalloc ( &d_data1D, size );
	cudaMemcpy ( d_data1D, &h_data1D[0], size, cudaMemcpyHostToDevice);

	return d_data1D;
}

int cu_test_cuda ( cu_ptrs &d_ptrs, int nRow, int nCol ) {

    cudaPrintfInit();

    cout << "Launching testKernel ..." << endl;

    testKernel<<<100,512>>>(d_ptrs.bmag.ptr, d_ptrs.bmag.pitch, nRow, nCol);

    cudaPrintfDisplay (stdout, true);
    cudaPrintfEnd();

    return 0;
}

