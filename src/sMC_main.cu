#include <iostream>
#include "eqdsk.hpp"
#include "particle.hpp"
#include "read_pl.hpp"
#include "constants.hpp"
#include <vector>
#include <ctime>
#include "cuda_wrap.h"
#include "interp.hpp"
#include "rkf.hpp"

using namespace std;

int main ()
{

	cout << "sMC ... c++ & CUDA version :)" << endl;

    int _Z = 1;
    int amu = 1;
	unsigned int nR = 64, nZ = 128;

	string fName = "data/eqdsk";
	//string fName = "data/g122080.03100";

	Ceqdsk eqdsk;	
	int stat;
	stat = eqdsk.read_file ( fName );
	stat = eqdsk.calc_b ( nZ, nR );
	stat = eqdsk.write_ncfile ( "output/bField.nc" );
	stat = eqdsk.bForceTerms ( _Z, amu );

    // Create the interpSpans container
    CinterpSpans spans (eqdsk.r.front(), eqdsk.r.back(), eqdsk.r.size(), 
            eqdsk.z.front(), eqdsk.z.back(), eqdsk.z.size() );

	vector<Cgc_particle> particles;	

	string plName = "data/pList.dav.nc.000";
	stat = read_pl ( plName, particles );

	// Get bMag @ particle locations and calculate mu.
	CinterpIndex index;
	REAL bmag_p;
	for(unsigned int p=0;p<particles.size();p++){
		index = get_index ( particles[p].r, particles[p].z, spans );
		bmag_p = bilinear_interp ( index, eqdsk.bmag );
		particles[p].mu = ( amu * _mi ) * pow(particles[p].vPer,2) / ( 2.0 * bmag_p );
		particles[p].energy_eV = 0.5 * ( amu * _mi ) * 
				( pow(particles[p].vPer,2) + pow(particles[p].vPar,2) ) / _e;
	} // end for(unsigned int p=0;p<particles.size();p++)	


	time_t startTime, endTime;
	startTime = time ( NULL );

	for(unsigned int p=0;p<10;p++) {

		if(!particles[p].status) {

			cout << "Particle No. " << p << endl;
			cout << "eV: " << particles[p].energy_eV << endl;

            stat = move_particle ( particles[p], eqdsk, spans, p );
			
		} // end if(!particle[p].status)
	} // end for(p)

	endTime = time ( NULL );

	cout << "Run took: " << difftime ( endTime, startTime ) << endl;

    cout << "*** CUDA ***" << endl;

	stat = copy_particles_to_device ( particles );

    cu_ptrs d_ptrs;

	d_ptrs.r = copy_1D_to_device (eqdsk.r,eqdsk.nCol);
	d_ptrs.z = copy_1D_to_device (eqdsk.z,eqdsk.nRow);

	d_ptrs.bmag = copy_2D_to_device (eqdsk.bmag,eqdsk.nRow,eqdsk.nCol);
    
	d_ptrs.b_r = copy_2D_to_device (eqdsk.br,eqdsk.nRow,eqdsk.nCol);
	d_ptrs.b_p = copy_2D_to_device (eqdsk.bp,eqdsk.nRow,eqdsk.nCol);
	d_ptrs.b_z = copy_2D_to_device (eqdsk.bz,eqdsk.nRow,eqdsk.nCol);

	d_ptrs.bCurv_r = copy_2D_to_device (eqdsk.bCurvature_r,eqdsk.nRow,eqdsk.nCol);
	d_ptrs.bCurv_p = copy_2D_to_device (eqdsk.bCurvature_p,eqdsk.nRow,eqdsk.nCol);
	d_ptrs.bCurv_z = copy_2D_to_device (eqdsk.bCurvature_z,eqdsk.nRow,eqdsk.nCol);

	d_ptrs.bGrad_r = copy_2D_to_device (eqdsk.bGradient_r,eqdsk.nRow,eqdsk.nCol);
	d_ptrs.bGrad_p = copy_2D_to_device (eqdsk.bGradient_p,eqdsk.nRow,eqdsk.nCol);
	d_ptrs.bGrad_z = copy_2D_to_device (eqdsk.bGradient_z,eqdsk.nRow,eqdsk.nCol);

    stat = cu_test_cuda ( d_ptrs, eqdsk.nRow, eqdsk.nCol );

	cout << "End of program :)" << endl;

	return 0;
}
