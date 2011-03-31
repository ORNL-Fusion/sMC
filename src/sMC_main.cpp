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
	REAL bmag_p;
	for(unsigned int p=0;p<particles.size();p++){
		stat = eqdsk.get_index(particles[p].r,particles[p].z,index);
		bmag_p = eqdsk.bilinear_interp ( index, eqdsk.bmag );
		particles[p].mu = _mi * (particles[p].vPer,2) / ( 2.0 * bmag_p );
	}	

	REAL dt = 1e-8;
	int err = 0;

	//C_rkGCparticle rkp1, rkp2, rkp3, rkp4, rkpf;
	Crk K1, K2, K3, K4, K5, K6, w;

	for(int p=0;p<1;p++) {

		if(!particles[p].status) {

			cout << "Particle No. " << p << endl;

			// Initial position w

			w.r = particles[p].r;
			w.p = particles[p].p;
			w.z = particles[p].z;

			for(int t=0;t<1;t++) {

				// Given a position, mu and vPar calculate vGC
				K1 = dt * vGC ( 0.0, 
						w, particles[p].mu, particles[p].vPar, eqdsk, err );

				K2 = dt * vGC ( 1.0/4.0 * dt, 
						w + K1 * 1.0/4.0, 
                        particles[p].mu, particles[p].vPar, eqdsk, err );

				K3 = dt * vGC ( 3.0/8.0 * dt, 
				        w + K1 * 3.0/32.0 + K2 * 9.0/32.0,
                        particles[p].mu, particles[p].vPar, eqdsk, err );

				K4 = dt * vGC ( 12.0/13.0 * dt, 
				        w + K1 * 1932.0/2179.0 + K2 * (-7200.0/2197.0) + K3 * 7296.0/2197.0,
						particles[p].mu, particles[p].vPar, eqdsk, err );

				K5 = dt * vGC ( dt, 
				        w + K1 * 439.0/216.0 + K2 * (-8.0) + K3 * 3680.0/513.9 + K4 * (-845.0/4104.0),
						particles[p].mu, particles[p].vPar, eqdsk, err );

				K6 = dt * vGC ( 0.5 * dt, 
				        w + K1 * (-8.0/27.0) + K2 * 2.0 + K3 * (-3544.0/2565.0) + K4 * 1859.0/4104.0 + K5 * (-11.0/40.0),
						particles[p].mu, particles[p].vPar, eqdsk, err );

				if(particles[p].status) {
					cout << "\twall" << endl;
					break;
				}
			}
		}
	}

	cout << "End of program :)" << endl;

	return 0;
}
