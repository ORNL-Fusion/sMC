#include <iostream>
#include "eqdsk.hpp"
#include "particle.hpp"
#include "gc_eqom.hpp"
#include "read_pl.hpp"
#include "constants.hpp"
#include <vector>
#include <netcdfcpp.h>

using namespace std;
using namespace constants;

int main ()
{

	cout << "sMC ... c++ & CUDA version :)" << endl;

    int _Z = 1;
    int amu = 1;

	string fName = "data/eqdsk";
	//string fName = "data/g122080.03100";

	Ceqdsk eqdsk;	
	int stat;
	stat = eqdsk.read_file ( fName );
	stat = eqdsk.write_ncfile ( "output/bField.nc" );
	stat = eqdsk.bForceTerms ( _Z, amu );

	vector<C_GCparticle> particles;	

	string plName = "data/pList.dav.nc.000";
	stat = read_pl ( plName, particles );

	// Get bMag @ particle locations and calculate mu.
	Ceqdsk::interpIndex index;
	REAL bmag_p;
	for(unsigned int p=0;p<particles.size();p++){
		stat = eqdsk.get_index(particles[p].r,particles[p].z,index);
		bmag_p = eqdsk.bilinear_interp ( index, eqdsk.bmag );
		particles[p].mu = ( amu * _mi ) * (particles[p].vPer,2) / ( 2.0 * bmag_p );
	}	


	// Runge-Kutta-Fehlberg integrator
	// pg. 254 Burden and Faires

	unsigned int ii = 0;	
	int err = 0;
	REAL t = 0.0;
	REAL runTime = 1e-4;
	REAL dt = 1e-8;
	REAL dtMax = 1e-6;
	REAL dtMin = 1e-10;
	REAL TOL = 1e-4;
	unsigned int FLAG;

	Crk K1, K2, K3, K4, K5, K6, w, R;
	REAL dvPar1=0, dvPar2=0, dvPar3=0, dvPar4=0, dvPar5=0, dvPar6=0;

	for(int p=0;p<100;p++) {

		if(!particles[p].status) {

			cout << "Particle No. " << p << endl;

			// Initial position w

			w.r = particles[p].r;
			w.p = particles[p].p;
			w.z = particles[p].z;
			w.vPar = particles[p].vPar;

			vector<REAL> rOut(1,particles[p].r);	
			vector<REAL> pOut(1,particles[p].p);	
			vector<REAL> zOut(1,particles[p].z);	

			FLAG = 1;
            ii = 0;
            err = 0;
            t = 0.0;

			while(FLAG==1 && ii < 10000) {

				cout << endl;
				cout << "\tStep: " << ii << endl;
				cout << "\tdt: " << dt << endl;
				cout << "\tt: " << t << endl;
				cout << "\tmu_mi: " << particles[p].mu / _mi << endl;
                cout << "\terr: " << err << endl;

				w.print();

				// Given a position, mu and vPar calculate vGC
				K1 = dt * vGC ( 0.0, 
						w, 
						particles[p].mu, eqdsk, err );

				K2 = dt * vGC ( 1.0/4.0 * dt, 
						w + K1 * 1.0/4.0, 
                        particles[p].mu, eqdsk, err );

				K3 = dt * vGC ( 3.0/8.0 * dt, 
				        w + K1 * 3.0/32.0 + K2 * 9.0/32.0,
                        particles[p].mu, eqdsk, err );

				K4 = dt * vGC ( 12.0/13.0 * dt, 
				        w + K1 * 1932.0/2197.0 + K2 * (-7200.0/2197.0) + K3 * 7296.0/2197.0,
						particles[p].mu, eqdsk, err );

				K5 = dt * vGC ( dt, 
				        w + K1 * 439.0/216.0 + K2 * (-8.0) + K3 * 3680.0/513.0 + K4 * (-845.0/4104.0),
						particles[p].mu, eqdsk, err );

				K6 = dt * vGC ( 0.5 * dt, 
				        w + K1 * (-8.0/27.0) + K2 * 2.0 + K3 * (-3544.0/2565.0) + K4 * 1859.0/4104.0 + K5 * (-11.0/40.0),
						particles[p].mu, eqdsk, err );

				R = Kabs ( 1.0/360.0*K1 - 128.0/4275.0*K3 - 2197.0/75240.0*K4 + 1.0/50.0*K5 + 2.0/55.0*K6 ) / dt;

				REAL R_ = Kmax ( R );

				REAL delta = 0.84 * pow ( TOL / R_, 1.0/4.0 );

				cout << "\tdelta: " << delta << endl;

				if(R_<=TOL) {
						// Approximation accepted
						t += dt;
						w += 25.0/216.0*K1 + 1408.0/2565.0*K3 + 2197.0/4104.0*K4 - 1.0/5.0*K5;
				}

				// Adjust time step
				if(delta<=0.1) {
					dt = 0.1 * dt;
				}
				else if(delta>=4.0) {
					dt = 4.0 * dt;
				}
				else {
					dt = delta * dt;
				}

				// Catch max dt
				if(dt>dtMax) {
					cout << "\tdtMax reached: " << dt <<" "<< dtMax << endl;
					dt = dtMax;
				}

				// End of desired time
				if(t>=runTime || err) {
					FLAG = 0;
				}
				else if(t+dt>runTime) {
					// Make sure integration ends == runTime
					dt = runTime - t;
				}
				else if(dt<dtMin) {
					cout << "\tdtMin reached: " << dt <<" "<< dtMin << endl;
					dt = dtMin;
					FLAG = 0;
				}

				//if(particles[p].status) {
				//	cout << "\twall" << endl;
				//	break;
				//}
			
				ii++;

				rOut.push_back(w.r);
				pOut.push_back(w.p);
				zOut.push_back(w.z);

			} // end while(FLAG=1)

			string orb_fName = "output/orbit.nc";
			NcFile dataFile ( &orb_fName[0], NcFile::Replace );
			if (!dataFile.is_valid()) {
				cout << "ERROR: Could not open nc file for writing." << endl;
				return 1;
			}	

			unsigned int nSteps = rOut.size();
			NcDim *nDim = dataFile.add_dim("n",nSteps);

			NcVar *nc_rOrb = dataFile.add_var ("rOrb", ncFloat, nDim );
			NcVar *nc_pOrb = dataFile.add_var ("pOrb", ncFloat, nDim );
			NcVar *nc_zOrb = dataFile.add_var ("zOrb", ncFloat, nDim );

			nc_rOrb->put(&rOut[0],nSteps);
			nc_pOrb->put(&pOut[0],nSteps);
			nc_zOrb->put(&zOut[0],nSteps);

		} // end if(!particle[p].status)
	} // end for(p)


	cout << "End of program :)" << endl;

	return 0;
}
