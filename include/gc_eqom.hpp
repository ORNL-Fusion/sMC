#ifndef GC_EQOM_HPP_
#define GC_EQOM_HPP_

#include "rk_gc_particle.hpp"
#include "eqdsk.hpp"
#include "constants.hpp"
#include "interp.hpp"
#include "array2D.hpp"
#include <iostream>

#ifdef __CUDA_ARCH__
#include "C/src/simplePrintf/cuPrintf.cu"
#endif

class Ctextures {
    public:
        array2D<REAL,BCHECK> bmag,
            bCurv_r, bCurv_p, bCurv_z,
		    bDotGradB, 
		    bGrad_r, bGrad_p, bGrad_z,
            b_r, b_p, b_z;

#ifdef __CUDACC__
        __host__ __device__
#endif
        bool isValid() const {
       
            bool err = true; 

            if(!bmag.ptr || bmag.pitchBytes<=1) { 
#ifndef __CUDA_ARCH__
				std::cout << "bmag NOT valid" << std::endl;
				std::cout << bmag.ptr << "  " << bmag.pitchBytes << std::endl;
#else
				cuPrintf("bmag NOT valid\n");
#endif
                err = false;
			}

            if(!b_r.ptr || b_r.pitchBytes<=1) {
#ifndef __CUDA_ARCH__
				std::cout << "b_r NOT valid" << std::endl;
#else
				cuPrintf("b_r NOT valid\n");
#endif
                err = false;
			}
            if(!b_p.ptr || b_p.pitchBytes<=1) 
                err = false;
            if(!b_z.ptr || b_z.pitchBytes<=1) 
                err = false;

            if(!bCurv_r.ptr || bCurv_r.pitchBytes<=1) 
                err = false;
            if(!bCurv_p.ptr || bCurv_p.pitchBytes<=1) 
                err = false;
            if(!bCurv_z.ptr || bCurv_z.pitchBytes<=1) 
                err = false;

            if(!bGrad_r.ptr || bGrad_r.pitchBytes<=1) 
                err = false;
            if(!bGrad_p.ptr || bGrad_p.pitchBytes<=1) 
                err = false;
            if(!bGrad_z.ptr || bGrad_z.pitchBytes<=1) 
                err = false;

            if(!bDotGradB.ptr || bDotGradB.pitchBytes<=1) 
                err = false;

            return err;

        }
};

// Calculate vGC given position, mu and vPar.

#ifdef __CUDACC__
__host__ __device__
#endif
Crk vGC ( const REAL dt, const Crk &p, const REAL mu, 
            const Ctextures &textures, const CinterpSpans &spans, int &err ) {

	Crk vGC;

    if(!err) {

	    // get background data(s) at particle location
        
	    CinterpIndex index;
	    index = get_index (p.z,p.r,spans);
        
	    if(index.stat>=1) {
            err++;
#ifndef __CUDA_ARCH__
			std::cout << "\tget_index() result off grid:" << std::endl;
			printf("p.r: %f, p.z: %f\n",p.r,p.z);
			printf("spans.mfront: %f, spans.mback: %f, spans.nfront: %f, spans.nback %f, spans.NROW: %i, spans.NCOL: %i\n", 
							spans.mfront, spans.mback, spans.nfront, spans.nback, 
							spans.NROW, spans.NCOL);
#else
			cuPrintf("get_index() result off grid: %i %i %i %i %f %f\n", 
							index.m1, index.m2, index.n1, index.n2, index.m, index.n);
			cuPrintf("p.r: %f, p.z: %f\n",p.r,p.z);
			cuPrintf("spans.mfront: %f, spans.mback: %f, spans.nfront: %f, spans.nback %f, spans.NROW: %i, spans.NCOL: %i\n", 
							spans.mfront, spans.mback, spans.nfront, spans.nback, 
							spans.NROW, spans.NCOL);
#endif
            return vGC;
        }
		else {
#ifdef __CUDA_ARCH__
			cuPrintf("get_index() GOOD: %i %i %i %i %f %f\n", 
							index.m1, index.m2, index.n1, index.n2, index.m, index.n);
#else
			printf("get_index() GOOD: %i %i %i %i %f %f\n", 
							index.m1, index.m2, index.n1, index.n2, index.m, index.n);
#endif
		}

        REAL bmag = bilinear_interp ( index, textures.bmag );
#ifdef __CUDA_ARCH__
			cuPrintf("bmag: %f\n", bmag);
#else
			printf("bmag: %f\n", bmag);
#endif
	
        REAL b_r = bilinear_interp ( index, textures.b_r );
        REAL b_p = bilinear_interp ( index, textures.b_p );
        REAL b_z = bilinear_interp ( index, textures.b_z );

        REAL bCurv_r = bilinear_interp ( index, textures.bCurv_r );
        REAL bCurv_p = bilinear_interp ( index, textures.bCurv_p );
        REAL bCurv_z = bilinear_interp ( index, textures.bCurv_z );

        REAL bGrad_r = bilinear_interp ( index, textures.bGrad_r );
        REAL bGrad_p = bilinear_interp ( index, textures.bGrad_p );
        REAL bGrad_z = bilinear_interp ( index, textures.bGrad_z );

        REAL bDotGradB = bilinear_interp ( index, textures.bDotGradB );

	    REAL unitb_r = b_r / bmag;
	    REAL unitb_p = b_p / bmag;
	    REAL unitb_z = b_z / bmag;

	    // vPer
	    REAL vPer = sqrtf ( 2.0 * mu * bmag / _mi );
	    // dvPar_dt
	    REAL dvPar_dt = -mu / _mi * bDotGradB;
	    // Here vGC is a dvGC and so vGC.vPar is really a dvPar.
	    // I'm just using the Crk class as containers for x/dt and v/dt quantities.
	    vGC.vPar = dvPar_dt; 
	    REAL vPar = p.vPar;

	    // vGC
	    vGC.r = vPar * unitb_r + pow(vPer,2) * bGrad_r + pow(vPar,2) * bCurv_r;
	    vGC.p = vPar * unitb_p + pow(vPer,2) * bGrad_p + pow(vPar,2) * bCurv_p;
	    vGC.z = vPar * unitb_z + pow(vPer,2) * bGrad_z + pow(vPar,2) * bCurv_z;

//#ifdef __CUDA_ARCH__
//
//	cuPrintf("end of vGC vPar: %f\n",vGC.vPar);
//	cuPrintf("end of vPer: %f\n",vPer);
//	cuPrintf("end of vPar: %f\n",vPar);
//  
//	cuPrintf("end of vGC bDotGradB: %e\n",bDotGradB);
//
//	cuPrintf("end of vGC r: %f\n",vGC.r);
//	cuPrintf("end of vGC p: %f\n",vGC.p);
//	cuPrintf("end of vGC z: %f\n",vGC.z);
//
//	cuPrintf("end of vGC unitb_r: %f\n",unitb_r);
//	cuPrintf("end of vGC unitb_p: %f\n",unitb_p);
//	cuPrintf("end of vGC unitb_z: %f\n",unitb_z);
//
//	cuPrintf("end of vGC bGrad_r: %e\n",bGrad_r);
//	cuPrintf("end of vGC bGrad_p: %e\n",bGrad_p);
//	cuPrintf("end of vGC bGrad_z: %e\n",bGrad_z);
//
//	cuPrintf("end of vGC bCurv_r: %e\n",bCurv_r);
//	cuPrintf("end of vGC bCurv_p: %e\n",bCurv_p);
//	cuPrintf("end of vGC bCurv_z: %e\n",bCurv_z);
//#else
//	printf("end of vGC vPar: %f\n",vGC.vPar);
//	printf("end of vPer: %f\n",vPer);
//	printf("end of vPar: %f\n",vPar);
//  
//	printf("end of vGC bDotGradB: %e\n",bDotGradB);
//
//	printf("end of vGC r: %f\n",vGC.r);
//	printf("end of vGC p: %f\n",vGC.p);
//	printf("end of vGC z: %f\n",vGC.z);
//
//	printf("end of vGC unitb_r: %f\n",unitb_r);
//	printf("end of vGC unitb_p: %f\n",unitb_p);
//	printf("end of vGC unitb_z: %f\n",unitb_z);
//
//	printf("end of vGC bGrad_r: %e\n",bGrad_r);
//	printf("end of vGC bGrad_p: %e\n",bGrad_p);
//	printf("end of vGC bGrad_z: %e\n",bGrad_z);
//
//	printf("end of vGC bCurv_r: %e\n",bCurv_r);
//	printf("end of vGC bCurv_p: %e\n",bCurv_p);
//	printf("end of vGC bCurv_z: %e\n",bCurv_z);
//#endif
	
    } // End if(!err) 

	return vGC;
}

#endif
