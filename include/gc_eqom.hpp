#ifndef GC_EQOM_HPP_
#define GC_EQOM_HPP_

#include "rk_gc_particle.hpp"
#include "eqdsk.hpp"
#include "constants.hpp"
#include "interp.hpp"
#include "array2D.hpp"
#include <iostream>

#ifdef __CUDA_ARCH__
#include "cuPrintf.cu"
#endif

class Ctextures {
    public:
        array2D<REAL,BCHECK> bmag,
            bCurv_r, bCurv_p, bCurv_z,
		    bDotGradB, 
		    bGrad_r, bGrad_p, bGrad_z,
            b_r, b_p, b_z;
};

// Calculate vGC given position, mu and vPar.
#ifdef __CUDACC__
__host__ __device__
#endif
Crk vGC ( const REAL dt, const Crk &p0, const REAL mu, 
            const Ctextures &textures, const CinterpSpans &spans, int &err ) {

	Crk vGC;

    if(!err) {

	    // get background data(s) at particle location
        
	    CinterpIndex index;
	    index = get_index ( p0.r, p0.z,	spans );
        
	    if(index.stat>=1) {
            err++;
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
	    REAL vPar = p0.vPar;

	    // vGC
	    vGC.r = vPar * unitb_r + pow(vPer,2) * bGrad_r + pow(vPar,2) * bCurv_r;
	    vGC.p = vPar * unitb_p + pow(vPer,2) * bGrad_p + pow(vPar,2) * bCurv_p;
	    vGC.z = vPar * unitb_z + pow(vPer,2) * bGrad_z + pow(vPar,2) * bCurv_z;

    } // End if(!err) 

	return vGC;
}

#endif
