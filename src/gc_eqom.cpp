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
	 C_rkGCparticle &pf, const float dt) {

	pf.v_r = ( p1.v_r + 2.0 * p2.v_r + 2.0 * p3.v_r + p4.v_r ) / 6.0;
	pf.v_p = ( p1.v_p + 2.0 * p2.v_p + 2.0 * p3.v_p + p4.v_p ) / 6.0;
	pf.v_z = ( p1.v_z + 2.0 * p2.v_z + 2.0 * p3.v_z + p4.v_z ) / 6.0;

	pf.r	= p1.r + dt * pf.v_r;
	pf.p	= p1.p + dt * pf.v_p;
	pf.z	= p1.z + dt * pf.v_z;

	return 0;	
}

// Fill in position, mu and vPar for the next rk4 step.
int euler ( const C_rkGCparticle &p1, C_rkGCparticle &p2, const float dt ) {

	p2.mu = p1.mu;

	p2.r = p1.r + p1.v_r * dt / 2.0; 
	p2.p = p1.p + p1.v_p * dt / 2.0; 
	p2.z = p1.z + p1.v_z * dt / 2.0; 

	p2.vPar = p1.dvPar_dt * dt / 2.0;

	return 0;
}

// Calculate vGC given position, mu and vPar.
int vGC ( C_rkGCparticle &p0, Ceqdsk &eqdsk ) {

	// get background data(s) at particle location

	Ceqdsk::interpIndex index;

	int stat = eqdsk.get_index ( p0.r, p0.z, index );

	if(stat) { return 1; }

	cout << "\tfloat i: "<<index.i << endl;
	cout << "\tfloat j: "<<index.j << endl;

    p0.bmag = eqdsk.bilinear_interp ( index, eqdsk.bmag );

    p0.b_r = eqdsk.bilinear_interp ( index, eqdsk.br );
    p0.b_p = eqdsk.bilinear_interp ( index, eqdsk.bp );
    p0.b_z = eqdsk.bilinear_interp ( index, eqdsk.bz );

    p0.bCurv_r = eqdsk.bilinear_interp ( index, eqdsk.bCurvature_r );
    p0.bCurv_p = eqdsk.bilinear_interp ( index, eqdsk.bCurvature_p );
    p0.bCurv_z = eqdsk.bilinear_interp ( index, eqdsk.bCurvature_z );

    p0.bGrad_r = eqdsk.bilinear_interp ( index, eqdsk.bGradient_r );
    p0.bGrad_p = eqdsk.bilinear_interp ( index, eqdsk.bGradient_p );
    p0.bGrad_z = eqdsk.bilinear_interp ( index, eqdsk.bGradient_z );

    p0.bDotGradB = eqdsk.bilinear_interp ( index, eqdsk.bDotGradB );

	p0.unitb_r = p0.b_r / p0.bmag;
	p0.unitb_p = p0.b_p / p0.bmag;
	p0.unitb_z = p0.b_z / p0.bmag;

	// vPer
	p0.vPer = sqrt ( 2.0 * p0.mu * p0.bmag / _mi );

	// vGC
	p0.v_r = p0.vPar * p0.unitb_r + pow(p0.vPer,2) * p0.bGrad_r + pow(p0.vPar,2) * p0.bCurv_r;
	p0.v_p = p0.vPar * p0.unitb_p + pow(p0.vPer,2) * p0.bGrad_p + pow(p0.vPar,2) * p0.bCurv_p;
	p0.v_z = p0.vPar * p0.unitb_z + pow(p0.vPer,2) * p0.bGrad_z + pow(p0.vPar,2) * p0.bCurv_z;

	// dvPar_dt
	p0.dvPar_dt = -p0.mu / _mi * p0.bDotGradB;

	return 0;
}


