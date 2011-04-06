#include <iostream>
#include "eqdsk.hpp"
#include "particle.hpp"
#include "gc_eqom.hpp"
#include "read_pl.hpp"
#include "constants.hpp"
#include <vector>
#include <netcdfcpp.h>
#include <sstream>
#include <ctime>
#include "cuda_wrap.h"

#define __SAVE_ORBITS__

using namespace std;
using namespace constants;

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

	vector<Cgc_particle> particles;	

	string plName = "data/pList.dav.nc.000";
	stat = read_pl ( plName, particles );

	// Get bMag @ particle locations and calculate mu.
	Ceqdsk::interpIndex index;
	REAL bmag_p;
	for(unsigned int p=0;p<particles.size();p++){
		stat = eqdsk.get_index(particles[p].r,particles[p].z,index);
		bmag_p = eqdsk.bilinear_interp ( index, eqdsk.bmag );
		particles[p].mu = ( amu * _mi ) * pow(particles[p].vPer,2) / ( 2.0 * bmag_p );
		particles[p].energy_eV = 0.5 * ( amu * _mi ) * 
				( pow(particles[p].vPer,2) + pow(particles[p].vPar,2) ) / _e;
	}	


	// Runge-Kutta-Fehlberg integrator
	// pg. 254 Burden and Faires

	unsigned int ii = 0;	
	int err = 0;
	REAL t = 0.0;
	REAL runTime = 1e-3;
	REAL dt;
	REAL dtMax = 1e-4;
	REAL dtMin = 1e-9;
	REAL EPS = 1e-3;
	unsigned int FLAG;


	time_t startTime, endTime;
	startTime = time ( NULL );

	for(int p=0;p<10;p++) {

		if(!particles[p].status) {

			cout << "Particle No. " << p << endl;
			cout << "eV: " << particles[p].energy_eV << endl;

			Crk K1, K2, K3, K4, K5, K6, w, R;

			// Initial position w

			w.r = particles[p].r;
			w.p = particles[p].p;
			w.z = particles[p].z;
			w.vPar = particles[p].vPar;

#ifdef __SAVE_ORBITS__
			vector<REAL> rOut(1,particles[p].r);	
			vector<REAL> pOut(1,particles[p].p);	
			vector<REAL> zOut(1,particles[p].z);	
#endif

			// Reset counters for new paricle 
			FLAG = 1;
            ii = 0;
            err = 0;
            t = 0.0;
			dt = dtMin;

			while(FLAG==1) {

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

				if(err) {
					particles[p].status = err;
					break;
				}

				R = Kabs ( 1.0/360.0*K1 - 128.0/4275.0*K3 - 2197.0/75240.0*K4 + 1.0/50.0*K5 + 2.0/55.0*K6 ) / dt;

				REAL R_ = Kmax ( R );

				// Error control is scaled with energy
				REAL TOL = EPS * particles[p].energy_eV;
				REAL delta = 0.84 * pow ( TOL / R_, 0.25 );

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
				if(t>=runTime) {
					FLAG = 0;

					cout << "\teV: "<<particles[p].energy_eV<<endl;
				}
				else if(t+dt>runTime) {
					// Make sure integration ends == runTime
					dt = runTime - t;
				}
				else if(dt<dtMin) {
					cout << "\tdtMin reached: " << dt <<" "<< dtMin << endl;
					cout << "\tdelta: "<<delta<<endl;
					cout << "\tR_: "<<R_<<endl;
					cout << "\tTOL: "<<TOL<<endl;
					cout << "\tvPer: "<<particles[p].vPer<<endl;
					cout << "\teV: "<<particles[p].energy_eV<<endl;
					dt = dtMin;
					FLAG = 0;
				}

				ii++;

#ifdef __SAVE_ORBITS__
				rOut.push_back(w.r);
				pOut.push_back(w.p);
				zOut.push_back(w.z);
#endif
			} // end while(FLAG=1)

			particles[p].r = w.r; 
			particles[p].p = w.p; 
			particles[p].z = w.z; 
			particles[p].vPar = w.vPar;

			cout << "\t nSteps: " << ii << endl;

#ifdef __SAVE_ORBITS__
			stringstream orb_fName;
			orb_fName << setfill('0');
		   	orb_fName << "output/" << setw(4) << p << "orbit.nc";
			NcFile dataFile ( &orb_fName.str()[0], NcFile::Replace );
			if (!dataFile.is_valid()) {
				cout << "ERROR: Could not open nc file for writing." << endl;
				return 1;
			}	

			unsigned int nSteps = rOut.size();
			NcDim *nDim = dataFile.add_dim("n",nSteps);
			NcDim *sDim = dataFile.add_dim("s",1);

			NcVar *nc_rOrb = dataFile.add_var ("rOrb", ncFloat, nDim );
			NcVar *nc_pOrb = dataFile.add_var ("pOrb", ncFloat, nDim );
			NcVar *nc_zOrb = dataFile.add_var ("zOrb", ncFloat, nDim );

			NcVar *nc_vPar = dataFile.add_var ("vPar", ncFloat, sDim );

			NcVar *nc_stat = dataFile.add_var ("stat", ncInt, sDim );

			nc_rOrb->put(&rOut[0],nSteps);
			nc_pOrb->put(&pOut[0],nSteps);
			nc_zOrb->put(&zOut[0],nSteps);
			nc_stat->put(&particles[p].status,1);
			nc_vPar->put(&particles[p].vPar,1);
#endif
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
