#include "eqdsk.hpp"
#include "particle.hpp"
#include "constants.hpp"
#include <iostream>

using namespace constants;
using namespace std;

int average_vGC 
	(const C_rkGCparticle &p1, 
	 const C_rkGCparticle &p2, 
	 const C_rkGCparticle &p3, 
	 const C_rkGCparticle &p4, 
	 C_rkGCparticle &pf, const REAL dt) {

	pf.v_r = ( p1.v_r + 2.0 * p2.v_r + 2.0 * p3.v_r + p4.v_r ) / 6.0;
	pf.v_p = ( p1.v_p + 2.0 * p2.v_p + 2.0 * p3.v_p + p4.v_p ) / 6.0;
	pf.v_z = ( p1.v_z + 2.0 * p2.v_z + 2.0 * p3.v_z + p4.v_z ) / 6.0;

	pf.r	= p1.r + dt * pf.v_r;
	pf.p	= p1.p + dt * pf.v_p;
	pf.z	= p1.z + dt * pf.v_z;

	pf.vPar = ( p1.vPar + 2.0 * p2.vPar + 2.0 * p3.vPar + p4.vPar ) / 6.0;
	
	pf.mu = p1.mu;

	return 0;	
}

// Fill in position, mu and vPar for the next rk4 step.
int euler ( const C_rkGCparticle &p1, C_rkGCparticle &p2, const REAL dt ) {

	p2.mu = p1.mu;

	p2.r = p1.r + p1.v_r * dt; 
	p2.p = p1.p + p1.v_p * dt; 
	p2.z = p1.z + p1.v_z * dt; 

	p2.vPar = p1.vPar + p1.dvPar_dt * dt;

	return 0;
}

// Calculate vGC given position, mu and vPar.
Crk vGC ( const REAL dt, const Crk &p0, const REAL mu, const REAL vPar0, Ceqdsk &eqdsk, int err ) {

	// get background data(s) at particle location
	
	Ceqdsk::interpIndex index;
	err = eqdsk.get_index ( p0.r, p0.z, index );

	//cout << "\t"<< __FILE__ << endl;
	//cout << "\tREAL i: "<<index.i << endl;
	//cout << "\tREAL j: "<<index.j << endl;

    REAL bmag = eqdsk.bilinear_interp ( index, eqdsk.bmag );

    REAL b_r = eqdsk.bilinear_interp ( index, eqdsk.br );
    REAL b_p = eqdsk.bilinear_interp ( index, eqdsk.bp );
    REAL b_z = eqdsk.bilinear_interp ( index, eqdsk.bz );

    REAL bCurv_r = eqdsk.bilinear_interp ( index, eqdsk.bCurvature_r );
    REAL bCurv_p = eqdsk.bilinear_interp ( index, eqdsk.bCurvature_p );
    REAL bCurv_z = eqdsk.bilinear_interp ( index, eqdsk.bCurvature_z );

    REAL bGrad_r = eqdsk.bilinear_interp ( index, eqdsk.bGradient_r );
    REAL bGrad_p = eqdsk.bilinear_interp ( index, eqdsk.bGradient_p );
    REAL bGrad_z = eqdsk.bilinear_interp ( index, eqdsk.bGradient_z );

    REAL bDotGradB = eqdsk.bilinear_interp ( index, eqdsk.bDotGradB );

	REAL unitb_r = b_r / bmag;
	REAL unitb_p = b_p / bmag;
	REAL unitb_z = b_z / bmag;

	// vPer
	REAL vPer = sqrt ( 2.0 * mu * bmag / _mi );
	// dvPar_dt
	REAL dvPar_dt = -mu / _mi * bDotGradB;
	REAL vPar = vPar0 + dvPar_dt * dt;

	// vGC
	Crk vGC;
	vGC.r = vPar * unitb_r + pow(vPer,2) * bGrad_r + pow(vPar,2) * bCurv_r;
	vGC.p = vPar * unitb_p + pow(vPer,2) * bGrad_p + pow(vPar,2) * bCurv_p;
	vGC.z = vPar * unitb_z + pow(vPer,2) * bGrad_z + pow(vPar,2) * bCurv_z;

	return vGC;
}


