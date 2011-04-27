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
 	unsigned int ii = 0;	
	int err = 0;
	REAL t = 0.0;
	REAL runTime = 1e-3;
	REAL dt;
	REAL dtMax = 1e-4;
	REAL dtMin = 1e-9;
	REAL EPS = 1e-3;
	unsigned int FLAG = 1;

	Crk K1, K2, K3, K4, K5, K6, R;

	dt = dtMin;

	// Initial position w
	Crk w(p.r,p.p,p.z,p.vPar);

//#ifdef __CUDA_ARCH__
//	cuPrintf("move_particle eV: %f\n",p.energy_eV);
//	cuPrintf("move_particle r: %f\n",p.r);
//	cuPrintf("move_particle p: %f\n",p.p);
//	cuPrintf("move_particle z: %f\n",p.z);
//	cuPrintf("move_particle vPar: %f\n",p.vPar);
//
//	cuPrintf("move_particle r: %f\n",w.r);
//	cuPrintf("move_particle p: %f\n",w.p);
//	cuPrintf("move_particle z: %f\n",w.z);
//	cuPrintf("move_particle vPar: %f\n",w.vPar);
//
//	cuPrintf("move_particle spans: %f %f %i %f %f %i\n",
//					spans.rfront, spans.rback, spans.rsize, 
//					spans.zfront, spans.zback, spans.zsize);
//#endif

#ifndef __CUDA_ARCH__
#ifdef __SAVE_ORBITS__
	std::vector<REAL> rOut(1,p.r);	
	std::vector<REAL> pOut(1,p.p);	
	std::vector<REAL> zOut(1,p.z);	
#endif
#endif

	while(FLAG==1) {

		// Given a position, mu and vPar calculate vGC
		K1 = dt * vGC ( 0.0, 
				w, 
				p.mu, textures, spans, err );

		if(!K1.isValid()) {
			err++;
#ifndef __CUDA_ARCH__
			std::cout << "\tK1 not valid" << std::endl;
#else
			cuPrintf("\tK1 not valid\n");
#endif
			return err;
		}	

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
			std::cout << "\tdtMin reached: " << dt <<" "<< dtMin << std::endl;
			std::cout << "\tdelta: "<<delta<<std::endl;
			std::cout << "\tR_: "<<R_<<std::endl;
			std::cout << "\tTOL: "<<TOL<<std::endl;
			std::cout << "\tvPer: "<<p.vPer<<std::endl;
			std::cout << "\teV: "<<p.energy_eV<<std::endl;
#else
			cuPrintf("dtMin reached: %f, %f\n",dt,dtMin);
			cuPrintf("delta: %f\n",delta);
			cuPrintf("R_: %f\n",R_);
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
	std::cout << "\t nSteps: " << ii << std::endl;
#else
    cuPrintf("nSteps: %i\n", ii);
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
