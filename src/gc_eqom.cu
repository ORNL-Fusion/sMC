#include "eqdsk.hpp"
#include "rk_gc_particle.hpp"
#include "constants.hpp"
#include <iostream>
#include "interp.hpp"

using namespace constants;
using namespace std;


// Calculate vGC given position, mu and vPar.
Crk vGC ( const REAL dt, const Crk &p0, const REAL mu, Ceqdsk &eqdsk, int &err ) {

	Crk vGC;

	// get background data(s) at particle location
	
	interpIndex index;
	index = get_index ( p0.r, p0.z,	
                eqdsk.r.front(), eqdsk.r.back(), eqdsk.r.size(),
                eqdsk.z.front(), eqdsk.z.back(), eqdsk.z.size() );

	if(index.stat) {

	    cout << "\t" << __FILE__  << "\tREAL i: "<<index.i << endl;
	    cout << "\t" << __FILE__  << "\tREAL j: "<<index.j << endl;
	    cout << "\t" << __FILE__  << "\tREAL r: "<< p0.r << endl;
	    cout << "\t" << __FILE__  << "\tREAL z: "<< p0.z << endl;

        return vGC;
    }

    REAL bmag = bilinear_interp ( index, eqdsk.bmag );

    REAL b_r = bilinear_interp ( index, eqdsk.br );
    REAL b_p = bilinear_interp ( index, eqdsk.bp );
    REAL b_z = bilinear_interp ( index, eqdsk.bz );

    REAL bCurv_r = bilinear_interp ( index, eqdsk.bCurvature_r );
    REAL bCurv_p = bilinear_interp ( index, eqdsk.bCurvature_p );
    REAL bCurv_z = bilinear_interp ( index, eqdsk.bCurvature_z );

    REAL bGrad_r = bilinear_interp ( index, eqdsk.bGradient_r );
    REAL bGrad_p = bilinear_interp ( index, eqdsk.bGradient_p );
    REAL bGrad_z = bilinear_interp ( index, eqdsk.bGradient_z );

    REAL bDotGradB = bilinear_interp ( index, eqdsk.bDotGradB );

	REAL unitb_r = b_r / bmag;
	REAL unitb_p = b_p / bmag;
	REAL unitb_z = b_z / bmag;

	//cout << index.i<<"\t"<<index.j<<"\t"<<eqdsk.br[index.i][index.j]<<"\t"<<eqdsk.bp[index.i][index.j]<<endl;
	//cout << b_r <<"\t"<<b_p<<"\t"<<b_z<<endl;
	//cout << unitb_r <<"\t"<<unitb_p<<"\t"<<unitb_z<<endl;
	//cout << bCurv_r <<"\t"<<bCurv_p<<"\t"<<bCurv_z<<endl;
	//cout << bGrad_r <<"\t"<<bGrad_p<<"\t"<<bGrad_z<<endl;

	// vPer
	REAL vPer = sqrt ( 2.0 * mu * bmag / _mi );
	//cout << "\tvPer: "<<vPer<<endl;
	//cout << "\t"<<p0.r<<" "<<p0.p<<" "<<p0.z<<endl;
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

	return vGC;
}


