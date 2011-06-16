#ifndef RKF_HPP_
#define RKF_HPP_

#include <vector>
#include <iostream>
#include <sstream>
#include "rk_gc_particle.hpp"
#include "particle.hpp"
#include "gc_eqom.hpp"
#include "interp.hpp"
#include <netcdfcpp.h>
#include "constants.hpp"

#ifdef __CUDA_ARCH__
#include "C/src/simplePrintf/cuPrintf.cu"
#endif

#define __SAVE_ORBITS__

// Runge-Kutta-Fehlberg integrator
// pg. 254 Burden and Faires

namespace {

#ifdef __CUDACC__
__host__ __device__
#endif
int move_particle ( Cgc_particle &p, const Ctextures &textures, 
        const CinterpSpans &spans, const unsigned int pp ) {

    // Reset counters for new paricle 
 	unsigned int ii = 0, jj = 0;	
	int err = 0;
	REAL t = 0.0;
	REAL runTime = 1e-3;
	REAL dt;
	REAL dtMax = 1e-4;
	REAL dtMin = 1e-9;
	REAL EPS = 1e-5;
	unsigned int FLAG = 1;

	Crk K1, K2, K3, K4, K5, K6, R, v_ms;

	dt = dtMin;

	// Initial position w
	Crk w(p.r,p.p,p.z,p.vPar);

#ifndef __CUDA_ARCH__
#ifdef __SAVE_ORBITS__
	std::vector<REAL> rOut(1,p.r);	
	std::vector<REAL> pOut(1,p.p);	
	std::vector<REAL> zOut(1,p.z);	
#endif
#endif

	while(FLAG==1) {

		// Given a position, mu and vPar calculate vGC
		// K1,K2, etc contain dx (position) and dvPar.
		
		K1 = dt * vGC ( 0.0, 
				w, 
				p.mu, textures, spans, err );

		K2 = dt * vGC ( 1.0/4.0 * dt, 
				w + K1 * 1.0/4.0, 
                p.mu, textures, spans, err );
		
		K3 = dt * vGC ( 3.0/8.0 * dt, 
		        w + K1 * 3.0/32.0 + K2 * 9.0/32.0,
                p.mu, textures, spans, err );

		K4 = dt * vGC ( 12.0/13.0 * dt, 
		        w + K1 * 1932.0/2197.0 + K2 * (-7200.0/2197.0) + K3 * 7296.0/2197.0,
				p.mu, textures, spans, err );

		K5 = dt * vGC ( dt, 
		        w + K1 * 439.0/216.0 + K2 * (-8.0) + K3 * 3680.0/513.0 + K4 * (-845.0/4104.0),
				p.mu, textures, spans, err );

		K6 = dt * vGC ( 0.5 * dt, 
		        w + K1 * (-8.0/27.0) + K2 * 2.0 + K3 * (-3544.0/2565.0) + K4 * 1859.0/4104.0 + K5 * (-11.0/40.0),
				p.mu, textures, spans, err );

		if(err) {
			p.status = err;
#ifndef __CUDA_ARCH__
			std::cout << "\tParticle lost." << std::endl;
#else
			cuPrintf("\tParticle lost, err = %i\n", err);
#endif
            return err;
		}

		R = Kabs ( 1.0/360.0*K1 - 128.0/4275.0*K3 - 2197.0/75240.0*K4 + 1.0/50.0*K5 + 2.0/55.0*K6 ) / dt;

		REAL R_ = Kmax ( R );

		// Error control is scaled with energy
		REAL TOL = EPS * p.energy_eV;
		REAL delta = 0.84 * powf ( TOL / R_, 0.25 );

		if(R_<=TOL) {
				// Approximation accepted
				t += dt;
				w += 25.0/216.0*K1 + 1408.0/2565.0*K3 + 2197.0/4104.0*K4 - 1.0/5.0*K5;
				jj++;

				// Apply pitch angle scatter
				// 
				// Monte-Carlo evaluation of transport coefficients
				// Boozer and Kuo-Petravic
				// Phys. Fluids, 24, 851 (1981)

				v_ms = vGC ( 0.0, w, p.mu, textures, spans, err );

				REAL vPer = get_vPer ( w, p.mu, textures, spans );

				double vMag_ms = sqrt( pow(w.vPar,2) + pow(vPer,2) );

				double pitch_0 = w.vPar / vMag_ms;

				float coulomb_log = 23.0; 
				double E_background_eV = 2.0e3;
				double T_background_eV = 2.0/3.0 * E_background_eV;
				double n_background_cm3 = 1e13;

				// vTh is here done in SI units since the x variable is dimensionless
				double vTh_ms = sqrt ( 2 * T_background_eV * _e / (p.amu*_mi) );

				// cgs units here following Boozer
				double nu_B = coulomb_log / 10.0 / 3e6 * sqrt(2.0/p.amu) * 
						n_background_cm3 / pow(T_background_eV,1.5);

				// x is the dimensionless ratio of v to vTh
				double x =  vMag_ms / vTh_ms;

				double phi_x = erf ( x );	
				// Here we use the analytic erf derivative (Phi prime)
				double phi_x_dx = 2 / sqrt (_pi) * exp ( -pow(x,2) );
				double psi_x = ( phi_x - x * phi_x_dx ) / ( 2 * pow(x,2) );

				double nu_D = 3 * sqrt (_pi/2) * nu_B * (phi_x-psi_x)/pow(x,3);

				// random number generation

				srand(time(0));
				int pm = rand() % 2; // 0 or 1, i.e., odd or even
				if(!pm) pm = -1;

				// apply pitch kick
				
				double pitch_1 = pitch_0 * (1-nu_D*dt) * pm * sqrt ( (1-pow(pitch_0,2)) * nu_D*dt );
#if DEBUGLEVEL >= 3
#ifndef __CUDA_ARCH__
			printf  ("\tvMag = %f\n", sqrt(pow(p.vPer,2)+pow(p.vPar,2)));
			printf  ("\tvPar = %f\n", w.vPar);
			printf  ("\tvPer = %f\n", vPer);
			printf  ("\tv = %f\n", vMag_ms);
			printf  ("\tvTh = %f\n", vTh_ms);
			printf  ("\tv/vTh = %f\n", x);
			printf  ("\tphi_x = %f\n", phi_x);
			printf  ("\tpsi_x = %f\n", psi_x);
			printf  ("\tphi_x_dx = %e\n", phi_x_dx);
			printf  ("\tnu_D/nu_B = %e\n", nu_D/nu_B);
			printf  ("\tnu_D = %e\n", nu_D);
			printf  ("\tnu_B = %e\n", nu_B);
			printf  ("\tnu_D*dt = %e\n", nu_D*dt);
			printf  ("\tdt: %e, dtMin: %e\n", dt, 0.1/nu_D);
#else
			cuPrintf("\tv/vTh = %f\n", x);
			cuPrintf("\tnu_D/nu_B = %e\n", nu_D/nu_B);
			cuPrintf("\tnu_D*dt = %e\n", nu_D*dt);
			cuPrintf  ("\tdt: %e, dtMin: %e\n", dt, 0.1/nu_D);
#endif
			FLAG = 0; 
#endif
		}

//#ifndef __PROFILING__
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
//#endif
		// Catch max dt
		if(dt>dtMax) {
#ifndef __CUDA_ARCH__
			std::cout << "\tdtMax reached: " << dt <<" "<< dtMax << std::endl;
#else
			cuPrintf("dtMax reached: %f, %f\n",dt,dtMax);
#endif
		dt = dtMax;
		}

		// End of desired time
		if(t>=runTime) {
			FLAG = 0;
		}
		else if(t+dt>runTime) {
			// Make sure integration ends == runTime
			dt = runTime - t;
		}
		else if(dt<dtMin || dt!=dt) {
#ifndef __CUDA_ARCH__
			printf("dtMin reached: %e, %e\n",dt,dtMin);
			printf("delta: %f\n",delta);
			printf("R_: %e\n",R_);
			printf("TOL: %f\n",TOL);
			printf("vPer: %f\n",p.vPer);
			printf("eV: %f\n",p.energy_eV);
#else
			cuPrintf("dtMin reached: %e, %e\n",dt,dtMin);
			cuPrintf("delta: %f\n",delta);
			cuPrintf("R_: %e\n",R_);
			cuPrintf("TOL: %f\n",TOL);
			cuPrintf("vPer: %f\n",p.vPer);
			cuPrintf("eV: %f\n",p.energy_eV);
#endif
			dt = dtMin;
			FLAG = 0;
		}

		ii++;

#ifndef __CUDA_ARCH__
#ifdef __SAVE_ORBITS__
		rOut.push_back(w.r);
		pOut.push_back(w.p);
		zOut.push_back(w.z);
#endif
#endif

	} // end while(FLAG=1)

	p.r = w.r; 
	p.p = w.p; 
	p.z = w.z; 
	p.vPar = w.vPar;

#ifndef __CUDA_ARCH__
	std::cout << "\t nSteps: " << ii <<" with "<< jj << " accepted."<< std::endl;
#else
    cuPrintf("nSteps: %i with %i accepted.\n", ii,jj);
#endif

#ifndef __CUDA_ARCH__
#ifdef __SAVE_ORBITS__
	std::stringstream orb_fName;
	orb_fName << std::setfill('0');
   	orb_fName << "output/" << std::setw(4) << pp << "orbit.nc";
	NcFile dataFile ( &orb_fName.str()[0], NcFile::Replace );
	if (!dataFile.is_valid()) {
		std::cout << "ERROR: Could not open nc file for writing." << std::endl;
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
	nc_stat->put(&p.status,1);
	nc_vPar->put(&p.vPar,1);
#endif
#endif

    return err;

} // end void move_particle()

} // end namespace

#endif
