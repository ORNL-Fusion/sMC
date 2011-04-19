#include <iostream>
#include "eqdsk.hpp"
#include "particle.hpp"
#include "read_pl.hpp"
#include "constants.hpp"
#include <vector>
#include <ctime>
#include "interp.hpp"
#include "rkf.hpp"

#include "cuda_wrap.hpp"

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

    // Create the textures container
    Ctextures textures;
    textures.bmag = eqdsk.bmag;
    textures.bDotGradB = eqdsk.bDotGradB;

    textures.b_r = eqdsk.br;
    textures.b_p = eqdsk.bp;
    textures.b_z = eqdsk.bz;

    textures.bCurv_r = eqdsk.bCurvature_r;
    textures.bCurv_p = eqdsk.bCurvature_p;
    textures.bCurv_z = eqdsk.bCurvature_z;

    textures.bGrad_r = eqdsk.bGradient_r;
    textures.bGrad_p = eqdsk.bGradient_p;
    textures.bGrad_z = eqdsk.bGradient_z;

    // Create the particle list
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

            stat = move_particle ( particles[p], textures, spans, p );
			
		} // end if(!particle[p].status)
	} // end for(p)

	endTime = time ( NULL );

	cout << "Run took: " << difftime ( endTime, startTime ) << endl;

    cout << "*** CUDA ***" << endl;

    stat = cu_test_cuda ( particles, eqdsk.nRow, eqdsk.nCol, spans, eqdsk );

	cout << "End of program :)" << endl;

	return 0;
}
