#include <iostream>
#include "eqdsk.hpp"
#include "particle.hpp"
#include "gc_eqom.hpp"
#include "read_pl.hpp"
#include "constants.hpp"

using namespace std;
using namespace constants;

int main ()
{

	cout << "sMC ... c++ & CUDA version :)" << endl;


	string fName = "data/g129x129_1051206002.01120.cmod";
	//string fName = "data/g122080.03100";

	Ceqdsk eqdsk;	
	int stat;
	stat = eqdsk.read_file ( fName );
	stat = eqdsk.write_ncfile ( "output/bField.nc" );
	stat = eqdsk.bForceTerms ();

	vector<C_GCparticle> particles;	

	string plName = "data/pList.dav.nc.000";
	stat = read_pl ( plName, particles );

	// Get bMag @ particle locations and calculate mu.
	Ceqdsk::interpIndex index;
	float bmag_p;
	for(unsigned int p=0;p<particles.size();p++){
		stat = eqdsk.get_index(particles[p].r,particles[p].z,index);
		bmag_p = eqdsk.bilinear_interp ( index, eqdsk.bmag );
		particles[p].mu = _mi * (particles[p].vPer,2) / ( 2.0 * bmag_p );
		cout<<particles[p].mu<<endl;
	}	

	float dt = 1e-8;

	C_rkGCparticle rkp1, rkp2, rkp3, rkp4, rkpf;

	for(int p=0;p<10;p++) {

		// Initialize the first step integrator particle
		rkp1.r = particles[p].r; 
		rkp1.p = particles[p].r; 
		rkp1.z = particles[p].r; 

		rkp1.vPer = particles[p].vPer;
		rkp1.vPar = particles[p].vPar;
		
		rkp1.mu = particles[p].mu;

		for(int t=0;t<1;t++) {

			// Given a position, mu and vPar calculate vGC
			stat = vGC ( rkp1, eqdsk );

			stat = euler ( rkp1, rkp2, dt / 2.0);
			stat = vGC ( rkp2, eqdsk );

			stat = euler ( rkp2, rkp3, dt / 2.0);
			stat = vGC ( rkp3, eqdsk );

			stat = euler ( rkp3, rkp4, dt);
			stat = vGC ( rkp4, eqdsk );

			stat = average_vGC ( rkp1, rkp2, rkp3, rkp4, rkpf, dt );

		}
	}

	cout << "End of program :)" << endl;

	return 0;
}
