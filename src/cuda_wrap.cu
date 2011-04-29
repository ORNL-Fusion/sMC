#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>
#include "constants.hpp"
#include "cuda_classes.hpp"
#include "array2D.hpp"
#include "gc_eqom.hpp"
#include "rkf.hpp"
#include "interp.hpp"
#include <vector>
#include "particle.hpp"
#include "C/src/simplePrintf/cuPrintf.cu"
#include <ctime>

#define NTHREADS 128 
#define NTHREADBLOCKS 256 

#define SETUP_TEXTURE(texRef,array2D,NROW,NCOL) \
	(texRef).normalized = 0; \
  	(texRef).filterMode = cudaFilterModeLinear; \
  	(texRef).addressMode[0] = cudaAddressModeClamp; \
  	(texRef).addressMode[1] = cudaAddressModeClamp; \
	cudaBindTexture2D(0,&(texRef),(array2D).ptr,\
					&channelDesc,(NCOL),(NROW),(array2D).pitchBytes)

__global__ void cu_move_particles (Cgc_particle *const particles, 
                    Ctextures *const textures, const CinterpSpans spans, const unsigned int nP) {

	unsigned int my_idx = blockDim.x * blockIdx.x + threadIdx.x;

    if(textures[0].isValid()) {

        int stat;
	    for(unsigned int p=my_idx;p<nP;p+=NTHREADS*NTHREADBLOCKS) {
			
	    	if(!particles[p].status) {

                stat = move_particle ( particles[p], textures[0], spans, p );
	    		
	    	} // end if(!particle[p].status)
	    } // end for(p)

    } 
    else {

        cuPrintf("ERROR: dataPtrs not valid");

    } // end if(dataPtrs.isValid()) == false
}

__global__ void checkTextures ( Ctextures *const textures ) {

	unsigned int my_idx = blockDim.x * blockIdx.x + threadIdx.x;
	cuPrintf("%i %i %f\n", my_idx, 40, textures[0].bmag(my_idx,40));
	float x,y;
	x = 40;
	y = my_idx;
	cuPrintf("%i %i %f\n", my_idx, 40, tex2D(texRef_bmag,x,y) );

}

__global__ void check2Dcpy ( Ctextures *const textures, 
				const unsigned int nRow, const unsigned int nCol ) {

	for (int r=0;r<10;++r) {
			REAL *row = (REAL*)((char*)textures[0].bmag.ptr 
							+ r*textures[0].bmag.pitchBytes);
			for (int c=40;c<41;++c) {
					REAL element = row[c];
                    cuPrintf("%i %i %f\n", r, c, element);
                    cuPrintf("%i %i %f\n", r, c, textures[0].bmag(r,c));
					float x,y;
					x = c;
					y = r;
                    cuPrintf("%i %i %f\n", r, c, tex2D(texRef_bmag,x,y));
			}
	}
}

__global__ void check1Dcpy ( REAL *data1D, const unsigned int n ) {

	for(int i=0;i<n;++i) {
			//cuPrintf("%i %f\n", i, data1D[i]);
	}
}

void copy_2D_to_device
( const array2D<REAL,BCHECK> &h_data2D, array2D<REAL,BCHECK> &d_data2D ) {

	size_t width = h_data2D.N * sizeof(REAL);
	size_t height = h_data2D.M;
	size_t spitchBytes = h_data2D.N * sizeof(REAL);

	cudaMallocPitch(&d_data2D.ptr, &d_data2D.pitchBytes, width, height);
	cudaMemcpy2D(d_data2D.ptr, d_data2D.pitchBytes, &h_data2D(0,0), 
					spitchBytes, width, height, cudaMemcpyHostToDevice);
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
    const int nRow, const int nCol, const CinterpSpans &spans, const Ctextures &h_textures ) {

    // Copy particles to device

    //std::cout << "\tCopying particle list to device ... ";
    //thrust::device_vector<Cgc_particle> D_particles ( H_particles.begin(), H_particles.end() );
    //Cgc_particle * D_particles_ptr = thrust::raw_pointer_cast(&D_particles[0]);
    //std::cout << "DONE" << std::endl;

	Cgc_particle *D_particles;
	size_t sizeParticles = H_particles.size() * sizeof(Cgc_particle);
	cudaMalloc ( (void**)&D_particles, sizeParticles );
	cudaMemcpy ( D_particles, &H_particles[0], sizeParticles, cudaMemcpyHostToDevice);

    // Copy background data to device

    std::cout << "\tCopying textures to device ... ";

    Ctextures h_d_textures;

	copy_2D_to_device (h_textures.bmag,h_d_textures.bmag);

	copy_2D_to_device (h_textures.b_r,h_d_textures.b_r);
	copy_2D_to_device (h_textures.b_p,h_d_textures.b_p);
	copy_2D_to_device (h_textures.b_z,h_d_textures.b_z);

	copy_2D_to_device (h_textures.bCurv_r,h_d_textures.bCurv_r);
	copy_2D_to_device (h_textures.bCurv_p,h_d_textures.bCurv_p);
	copy_2D_to_device (h_textures.bCurv_z,h_d_textures.bCurv_z);

	copy_2D_to_device (h_textures.bGrad_r,h_d_textures.bGrad_r);
	copy_2D_to_device (h_textures.bGrad_p,h_d_textures.bGrad_p);
	copy_2D_to_device (h_textures.bGrad_z,h_d_textures.bGrad_z);

	copy_2D_to_device (h_textures.bDotGradB,h_d_textures.bDotGradB);

    std::cout << "DONE" << std::endl;

	if(h_d_textures.isValid()) {
		std::cout << "\tTextures are valid :)" << std::endl;
		std::cout << "\tWith pitchBytes of: " << h_d_textures.bmag.pitchBytes << std::endl;
	}
	else 
		std::cout << "\tTextures are NOT valid :)" << std::endl;

	Ctextures *d_textures;
	size_t size = sizeof(Ctextures);
	cudaMalloc ( (void**)&d_textures, size );
	cudaMemcpy ( d_textures, &h_d_textures, size, cudaMemcpyHostToDevice);

	// Setup hardware textures
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	SETUP_TEXTURE(texRef_bmag,h_d_textures.bmag,nRow,nCol);
	SETUP_TEXTURE(texRef_bDotGradB,h_d_textures.bDotGradB,nRow,nCol);

	SETUP_TEXTURE(texRef_bCurv_r,h_d_textures.bCurv_r,nRow,nCol);
	SETUP_TEXTURE(texRef_bCurv_p,h_d_textures.bCurv_p,nRow,nCol);
	SETUP_TEXTURE(texRef_bCurv_z,h_d_textures.bCurv_z,nRow,nCol);

	SETUP_TEXTURE(texRef_bGrad_r,h_d_textures.bGrad_r,nRow,nCol);
	SETUP_TEXTURE(texRef_bGrad_p,h_d_textures.bGrad_p,nRow,nCol);
	SETUP_TEXTURE(texRef_bGrad_z,h_d_textures.bGrad_z,nRow,nCol);

	SETUP_TEXTURE(texRef_b_r,h_d_textures.b_r,nRow,nCol);
	SETUP_TEXTURE(texRef_b_p,h_d_textures.b_p,nRow,nCol);
	SETUP_TEXTURE(texRef_b_z,h_d_textures.b_z,nRow,nCol);

    cudaPrintfInit();

    //std::cout << "Testing 1D memcopy ..." << std::endl;
	//check1Dcpy<<<1,1>>>( D_ptrs.z, nRow );
    //std::cout << "Testing 2D memcopy ..." << std::endl;
	//check2Dcpy<<<1,1>>>( d_textures, nRow, nCol );
    //std::cout << "Testing textures ..." << std::endl;
	//checkTextures<<<1,10>>>( d_textures );

    std::cout << "Moving particles ..." << std::endl;
	time_t startTime, endTime;
	startTime = time ( NULL );

	unsigned int nP;
	nP = H_particles.size();
	//nP = 10;
	cu_move_particles<<<NTHREADBLOCKS,NTHREADS>>>
			( D_particles, d_textures, spans, nP);

    cudaPrintfDisplay (stdout, true);
    cudaPrintfEnd();

	endTime = time ( NULL );
	std::cout << "CUDA Run took: " << difftime ( endTime, startTime ) << std::endl;

    return 0;
}

