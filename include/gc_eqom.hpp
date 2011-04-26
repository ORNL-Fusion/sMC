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

            if(!bmag.ptr || bmag.pitch<=1) { 
#ifndef __CUDA_ARCH__
				std::cout << "bmag NOT valid" << std::endl;
				std::cout << bmag.ptr << "  " << bmag.pitch << std::endl;
#else
				cuPrintf("bmag NOT valid\n");
#endif
                err = false;
			}

            if(!b_r.ptr || b_r.pitch<=1) {
#ifndef __CUDA_ARCH__
				std::cout << "b_r NOT valid" << std::endl;
#else
				cuPrintf("b_r NOT valid\n");
#endif
                err = false;
			}
            if(!b_p.ptr || b_p.pitch<=1) 
                err = false;
            if(!b_z.ptr || b_z.pitch<=1) 
                err = false;

            if(!bCurv_r.ptr || bCurv_r.pitch<=1) 
                err = false;
            if(!bCurv_p.ptr || bCurv_p.pitch<=1) 
                err = false;
            if(!bCurv_z.ptr || bCurv_z.pitch<=1) 
                err = false;

            if(!bGrad_r.ptr || bGrad_r.pitch<=1) 
                err = false;
            if(!bGrad_p.ptr || bGrad_p.pitch<=1) 
                err = false;
            if(!bGrad_z.ptr || bGrad_z.pitch<=1) 
                err = false;

            if(!bDotGradB.ptr || bDotGradB.pitch<=1) 
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
#else
			cuPrintf("get_index() result off grid: %i %i %i %i %f %f\n", 
							index.m1, index.m2, index.n1, index.n2, index.m, index.n);
			cuPrintf("\t p.r: %f, p.z: %f\n",p.r,p.z);
#endif
            return vGC;
        }

        REAL bmag = bilinear_interp ( index, textures.bmag );

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
	    // I'm just using the Crk class as containers for x/dx and v/dv quantities.
	    vGC.vPar = dvPar_dt; 
	    REAL vPar = p.vPar;

	    // vGC
	    vGC.r = vPar * unitb_r + pow(vPer,2) * bGrad_r + pow(vPar,2) * bCurv_r;
	    vGC.p = vPar * unitb_p + pow(vPer,2) * bGrad_p + pow(vPar,2) * bCurv_p;
	    vGC.z = vPar * unitb_z + pow(vPer,2) * bGrad_z + pow(vPar,2) * bCurv_z;

#ifdef __CUDA_ARCH__
//
//	cuPrintf("end of vGC vPar: %f\n",vGC.vPar);
//	cuPrintf("end of vGC vPer: %f\n",vPer);
//	cuPrintf("end of vGC vPar: %f\n",vPar);
//	
//	cuPrintf("end of vGC bmag: %f\n",bmag);
//	cuPrintf("vGC spans: %f %f %i %f %f %i\n",
//					spans.mfront, spans.mback, spans.NROW, 
//					spans.nfront, spans.nback, spans.NCOL);
//	cuPrintf("vGC index: %i %i %f %i %i %f %i\n",
//					index.m1, index.m2, index.m, 
//					index.n1, index.n2, index.n, index.stat);

//	cuPrintf("end of vGC bDotGradB: %f\n",bDotGradB);
//
//	cuPrintf("end of vGC r: %f\n",vGC.r);
//	cuPrintf("end of vGC p: %f\n",vGC.p);
//	cuPrintf("end of vGC z: %f\n",vGC.z);
//
//	cuPrintf("end of vGC unitb_r: %f\n",unitb_r);
//	cuPrintf("end of vGC unitb_p: %f\n",unitb_p);
//	cuPrintf("end of vGC unitb_z: %f\n",unitb_z);
//
//	cuPrintf("end of vGC bGrad_r: %f\n",bGrad_r);
//	cuPrintf("end of vGC bGrad_p: %f\n",bGrad_p);
//	cuPrintf("end of vGC bGrad_z: %f\n",bGrad_z);
//
//	cuPrintf("end of vGC bCurv_r: %f\n",bCurv_r);
//	cuPrintf("end of vGC bCurv_p: %f\n",bCurv_p);
//	cuPrintf("end of vGC bCurv_z: %f\n",bCurv_z);
//
#endif
	
    } // End if(!err) 

	return vGC;
}

#endif
