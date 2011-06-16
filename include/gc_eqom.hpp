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

#ifdef __CUDACC__
texture<float,cudaTextureType2D,cudaReadModeElementType>
	texRef_bmag, texRef_bDotGradB,
	texRef_bCurv_r, texRef_bCurv_p, texRef_bCurv_z,
	texRef_bGrad_r, texRef_bGrad_p, texRef_bGrad_z,
	texRef_b_r, texRef_b_p, texRef_b_z;
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

// Calculate vPer given position and mu

#ifdef __CUDACC__
__host__ __device__
#endif
REAL get_vPer ( const Crk &p, const REAL mu, const Ctextures &textures, const CinterpSpans &spans ) {

	    CinterpIndex index;
	    index = get_index (p.z,p.r,spans);
        REAL bmag = bilinear_interp ( index, textures.bmag );
	    return sqrtf ( 2.0 * mu * bmag / _mi );
}

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
///#ifdef __PROFILING__
///		index.stat = 0;
///#endif
        
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
//		else {
//#ifdef __CUDA_ARCH__
//			cuPrintf("get_index() GOOD: %i %i %i %i %f %f\n", 
//							index.m1, index.m2, index.n1, index.n2, index.m, index.n);
//#else
//			printf("get_index() GOOD: %i %i %i %i %f %f\n", 
//							index.m1, index.m2, index.n1, index.n2, index.m, index.n);
//#endif
//		}

#ifndef __CUDA_ARCH__	
        REAL bmag = bilinear_interp ( index, textures.bmag );
        REAL bDotGradB = bilinear_interp ( index, textures.bDotGradB );

        REAL b_r = bilinear_interp ( index, textures.b_r );
        REAL b_p = bilinear_interp ( index, textures.b_p );
        REAL b_z = bilinear_interp ( index, textures.b_z );

        REAL bCurv_r = bilinear_interp ( index, textures.bCurv_r );
        REAL bCurv_p = bilinear_interp ( index, textures.bCurv_p );
        REAL bCurv_z = bilinear_interp ( index, textures.bCurv_z );

        REAL bGrad_r = bilinear_interp ( index, textures.bGrad_r );
        REAL bGrad_p = bilinear_interp ( index, textures.bGrad_p );
        REAL bGrad_z = bilinear_interp ( index, textures.bGrad_z );

//#ifndef __CUDA_ARCH__	

#else
		float bmag = tex2D(texRef_bmag,index.n+0.5f,index.m+0.5f);
		float bDotGradB = tex2D(texRef_bDotGradB,index.n+0.5f,index.m+0.5f);

		float b_r = tex2D(texRef_b_r,index.n+0.5f,index.m+0.5f);
		float b_p = tex2D(texRef_b_p,index.n+0.5f,index.m+0.5f);
		float b_z = tex2D(texRef_b_z,index.n+0.5f,index.m+0.5f);

		float bCurv_r = tex2D(texRef_bCurv_r,index.n+0.5f,index.m+0.5f);
		float bCurv_p = tex2D(texRef_bCurv_p,index.n+0.5f,index.m+0.5f);
		float bCurv_z = tex2D(texRef_bCurv_z,index.n+0.5f,index.m+0.5f);

		float bGrad_r = tex2D(texRef_bGrad_r,index.n+0.5f,index.m+0.5f);
		float bGrad_p = tex2D(texRef_bGrad_p,index.n+0.5f,index.m+0.5f);
		float bGrad_z = tex2D(texRef_bGrad_z,index.n+0.5f,index.m+0.5f);

#endif

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

	    // vGC
	    vGC.r = p.vPar * unitb_r + pow(vPer,2) * bGrad_r + pow(p.vPar,2) * bCurv_r;
	    vGC.p = p.vPar * unitb_p + pow(vPer,2) * bGrad_p + pow(p.vPar,2) * bCurv_p;
	    vGC.z = p.vPar * unitb_z + pow(vPer,2) * bGrad_z + pow(p.vPar,2) * bCurv_z;

#if DEBUGLEVEL >= 4
		printf("v_r: %f, v_p: %f, v_z: %f\n", vGC.r, vGC.p, vGC.z);
		printf("vPer: %f, vPar: %f\n", vPer, p.vPar);
#endif


    } // End if(!err) 

	return vGC;
}

#endif
