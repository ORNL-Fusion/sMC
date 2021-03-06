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

#ifdef __CUDA_ARCH__
#define PRINT cuPrintf 
#else
#define PRINT printf
#endif

// Runge-Kutta-Fehlberg integrator
// pg. 254 Burden and Faires

namespace {

#ifdef __CUDACC__
__host__ __device__
#endif
float newPitch ( float pitch_0, float nu_D_dt ) {
	int pm = rand() % 2; // 0 or 1, i.e., odd or even
	if(!pm) pm = -1;
	double pitch_1 = pitch_0 * (1-nu_D_dt) 
			+ pm * sqrt ( (1-pow(pitch_0,2)) * nu_D_dt );
#if DEBUGLEVEL >= 3
	if(pitch_1 > 1 || pitch_1 < -1) {
			PRINT("\t\t\t*** Error in pitch angle update!!!\n");
			PRINT("\t\t\t*** pitch_0: %f, pitch_1: %f\n", pitch_0, pitch_1);
			exit(0);
	}
#endif
	return  pitch_1;
}

#ifdef __CUDACC__
__host__ __device__
#endif
float newEnergy ( float energy_0_eV, float vTh_ms, int amu, float T_background_eV, float nu_B, float nu_E, float nu_E_dt ) {

	int pm = rand() % 2; // 0 or 1, i.e., odd or even
	if(!pm) pm = -1;

	// Numerical evaluation of the dnu_E_dE term
	// See comments in the nu_D & nu_E section
	// I think that a once-per-orbit evaluation of this
	// term should be OK, and at worst if we include gradients
	// in T_background_eV / vTh_ms it would become an average.
	
	double dE_eV = energy_0_eV * 0.05;
	double v1 = sqrt ( 2 * (energy_0_eV - dE_eV) * _e / (amu*_mi) );
	double v2 = sqrt ( 2 * (energy_0_eV + dE_eV) * _e / (amu*_mi) );
	double x1 =  v1 / vTh_ms;
	double x2 =  v2 / vTh_ms;

	double phi_x1 = erf ( x1 );	
	double phi_x1_dx = 2 / sqrt (_pi) * exp ( -pow(x1,2) );
	double psi_x1 = ( phi_x1 - x1 * phi_x1_dx ) / ( 2 * pow(x1,2) );

	double phi_x2 = erf ( x2 );	
	double phi_x2_dx = 2 / sqrt (_pi) * exp ( -pow(x2,2) );
	double psi_x2 = ( phi_x2 - x2 * phi_x2_dx ) / ( 2 * pow(x2,2) );

	double nu_E1 = 3 * sqrt (_pi/2) * nu_B * psi_x1/x1;
	double nu_E2 = 3 * sqrt (_pi/2) * nu_B * psi_x2/x2;

	double dnu_E_dE = (nu_E2 - nu_E1) / ( 2 * dE_eV );

	double energy_1_eV = energy_0_eV - (2*nu_E_dt) * 
			(energy_0_eV - (3.0/2.0 + energy_0_eV / nu_E * dnu_E_dE)*T_background_eV)
			+ pm*2*sqrt(T_background_eV*energy_0_eV*(nu_E_dt));

#if DEBUGLEVEL >= 3
	if(energy_1_eV <= 0) {
			PRINT("\t\t\t*** Negative energy update!!!\n");
			exit(0);
	}
#endif

	return energy_1_eV;
}

// Calculate vPer given position and mu

#ifdef __CUDACC__
__host__ __device__
#endif
REAL get_vPer_ms ( const Crk &w, const Cgc_particle &p, const Ctextures &textures, const CinterpSpans &spans ) {

	    CinterpIndex index;
	    index = get_index (w.z,w.r,spans);
        REAL bmag = bilinear_interp ( index, textures.bmag );
	    return sqrtf ( 2.0 * p.mu * bmag / (p.amu*_mi) );
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
void update_vPar_and_mu ( float pitch_1, float energy_1_eV, const float vPer_ms_0, Crk &w, Cgc_particle &p ) {

#if DEBUGLEVEL >= 5
	double vPar_0 = w.vPar;
	double mu_0 = p.mu;
#endif

	// Update vPar
	double vMag_ms_1 = sqrt ( 2 * energy_1_eV * _e / (p.amu*_mi) );
	w.vPar = vMag_ms_1 * pitch_1;

	// Update mu
	double vPerSq_ms_1 = pow(vMag_ms_1,2)-pow(w.vPar,2);
	// Use the old vPer and mu to calculate B_T
	float B_T = (p.amu*_mi*pow(vPer_ms_0,2))/(p.mu * 2);

	p.mu = (p.amu*_mi) * vPerSq_ms_1 / ( 2 * B_T );
	p.energy_eV = energy_1_eV;
	//p.vPar = w.vPar;
	//p.vPer = sqrt ( vPerSq_ms_1 );

#if DEBUGLEVEL >= 5
	PRINT("\t\t\tvPar_0: %e, vPar_1: %e\n", vPar_0, w.vPar);
	PRINT("\t\t\tmu_0: %e, mu_1: %e\n", mu_0, p.mu);
#endif
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
void update_p ( const Crk &w, Cgc_particle &p ) {

	p.r = w.r; 
	p.p = w.p; 
	p.z = w.z; 
	p.vPar = w.vPar;

	// Remember that p.energy_eV is updated in the integration
	// as it is needed to scale the accuracy.
	double vMag_ms = sqrt ( 2 * p.energy_eV * _e / (p.amu*_mi) );
	double vPerSq_ms = pow(vMag_ms,2)-pow(p.vPar,2);
	p.vPer = sqrt ( vPerSq_ms );

}

#ifdef __CUDACC__
__host__ __device__
#endif
int move_particle ( Cgc_particle &p, const Ctextures &textures, 
        const CinterpSpans &spans, const unsigned int pp, const int seed ) {

    // Reset counters for new paricle 
 	unsigned int ii = 0, jj = 0;	
	int err = 0;
	REAL t = 0.0;
	REAL runTime = 2e-1;
	REAL dt;
	REAL dtMax = 1e-4;
	REAL dtMin = 1e-10;
	REAL EPS = 1e-3;
	unsigned int FLAG = 1;

#if GOOSE >= 1
	// Single poloidal orbit vars
	float poloidalDistanceOld = 0;
	float poloidalDistanceNew = 0;
	float poloidalDistanceClosed = 1e-3;
	Crk poloidalStart(p.r,p.p,p.z,p.vPar);
	int gooseFac = 1;
	bool goose = true;
	int poloidalPoints = 0;
	int nPoloidalOrbits = 0;
	float t_poloidalStart = 0;
#endif

	// Diffusion vars
	float nu_dt_max = 0.1;
	float nu_B = 0;
	float nu_D = 0, nu_D_dt_orbitTotal=0;
	float nu_E = 0, nu_E_dt_orbitTotal=0;
	float dt_nu = -1;
	float t_orbitTotal = 0;
	double vMag_ms = 0, vTh_ms;
	float coulomb_log = 23.0; 
	double E_background_eV = 2.0e3;
	double T_background_eV = 2.0/3.0 * E_background_eV;
	double n_background_cm3 = 5e13;

	// random number generation
	srand(seed);
	int pm = 0;

	Crk K1, K2, K3, K4, K5, K6, R, v_ms;

	dt = dtMin;

	// Initial position w
	Crk w(p.r,p.p,p.z,p.vPar);

#ifndef __CUDA_ARCH__
#ifdef __SAVE_ORBITS__
	std::vector<REAL> rOut(1,p.r);	
	std::vector<REAL> pOut(1,p.p);	
	std::vector<REAL> zOut(1,p.z);	
	std::vector<REAL> E_eVOut(1,p.energy_eV);	
	std::vector<REAL> pitchOut(1,p.vPar/sqrt(pow(p.vPer,2)+pow(p.vPar,2)));	
	std::vector<REAL> dtOut(1,dtMin);	
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
			FLAG = 0;
			PRINT("\tParticle lost, err = %i\n", err);
			break;
		}

		R = Kabs ( 1.0/360.0*K1 - 128.0/4275.0*K3 - 2197.0/75240.0*K4 + 1.0/50.0*K5 + 2.0/55.0*K6 ) / dt;

		REAL R_ = Kmax ( R );

		// Error control is scaled with vMag
		REAL TOL = EPS * sqrt(p.energy_eV);
		REAL delta = 0.84 * powf ( TOL / R_, 0.25 );

		if(R_<=TOL) {
			// Approximation accepted
			t += dt;

			w += 25.0/216.0*K1 + 1408.0/2565.0*K3 + 2197.0/4104.0*K4 - 1.0/5.0*K5;
			jj++;
#if GOOSE >= 1
			t_orbitTotal += dt;
			poloidalPoints++;
			poloidalDistanceOld = poloidalDistanceNew;
			poloidalDistanceNew = sqrt ( pow(poloidalStart.r-w.r,2)+pow(poloidalStart.z-w.z,2) );
#if DEBUGLEVEL >= 4
			PRINT("\t\tPoloidalDistance: %f\n", poloidalDistanceNew );
#endif
#endif
			// Calculate pitch and energy diffusion coefficients 
			// 
			// Monte-Carlo evaluation of transport coefficients
			// Boozer and Kuo-Petravic
			// Phys. Fluids, 24, 851 (1981)

			v_ms = vGC ( 0.0, w, p.mu, textures, spans, err );

			REAL vPer_ms_0 = get_vPer_ms ( w, p, textures, spans );

			vMag_ms = sqrt( pow(w.vPar,2) + pow(vPer_ms_0,2) );

			// vTh is here done in SI units since the x variable is dimensionless
			vTh_ms = sqrt ( 2 * T_background_eV * _e / (p.amu*_mi) );

			// cgs units here following Boozer
			nu_B = coulomb_log / 10.0 / 3e6 * sqrt(2.0/p.amu) * 
					n_background_cm3 / pow(T_background_eV,1.5);

			// x is the dimensionless ratio of v to vTh
			double x =  vMag_ms / vTh_ms;

			double phi_x = erf ( x );	
			// Here we use the analytic erf derivative (Phi prime)
			double phi_x_dx = 2 / sqrt (_pi) * exp ( -pow(x,2) );
			double psi_x = ( phi_x - x * phi_x_dx ) / ( 2 * pow(x,2) );

			nu_D = 3 * sqrt (_pi/2) * nu_B * (phi_x-psi_x)/pow(x,3);
			nu_E = 3 * sqrt (_pi/2) * nu_B * psi_x/x;

			// Check if the particle energy has dropped so low as to
			// require reducing dt to keep nu * dt << 1
			
			if(nu_D*dt>nu_dt_max || nu_E*dt>nu_dt_max) {
				// Unaccept the approximation and reduce dt
				t -= dt;
				w -= 25.0/216.0*K1 + 1408.0/2565.0*K3 + 2197.0/4104.0*K4 - 1.0/5.0*K5;
				jj--;
				dt_nu = nu_dt_max / nu_D * 0.9;
				if(nu_dt_max / nu_E * 0.9 < dt_nu) 
						dt_nu = nu_dt_max / nu_E * 0.9;
#if GOOSE >= 1
				goose = false;
				t_orbitTotal -= dt;
				poloidalPoints--;
#endif
#if DEBUGLEVEL >= 4 
				PRINT("\t\tPariticle so slow that nu*dt !<< 1\n");
				PRINT("\t\tdt: %e, dt_nu: %e\n",dt, dt_nu );
				PRINT("\t\tnu_D*dt: %f, nu_E*dt: %f\n", nu_D*dt, nu_E*dt );
				PRINT("\t\tv/vTh: %f\n", x);
#endif
			} 

#if GOOSE >= 1
			nu_D_dt_orbitTotal += nu_D*dt;
			nu_E_dt_orbitTotal += nu_E*dt;

			if(!goose && dt_nu<0) {
#else
			if(dt_nu<0) {
#endif
				double pitch_0 = w.vPar / vMag_ms;
#if PITCH_SCATTERING >= 1
				double pitch_1 = newPitch ( pitch_0, nu_D*dt );
#else
				double pitch_1 = pitch_0;
#endif

				double energy_0_eV = 0.5 * p.amu*_mi * pow(vMag_ms,2) / _e;
#if ENERGY_SCATTERING >= 1
				double energy_1_eV = 
						newEnergy ( energy_0_eV, vTh_ms, p.amu, T_background_eV, nu_B, nu_E, 
										nu_E*dt );
#else
				double energy_1_eV = energy_0_eV;
#endif
				update_vPar_and_mu ( pitch_1, energy_1_eV, vPer_ms_0, w, p );

#ifndef __CUDA_ARCH__
#ifdef __SAVE_ORBITS__
				rOut.push_back(w.r);
				pOut.push_back(w.p);
				zOut.push_back(w.z);
				E_eVOut.push_back(energy_1_eV);
				pitchOut.push_back(pitch_1);
				dtOut.push_back(dt);
#endif
#endif

			}

#if DEBUGLEVEL >= 4
			PRINT("\tvMag = %f\n", sqrt(pow(p.vPer,2)+pow(p.vPar,2)));
			PRINT("\tvPar = %f\n", w.vPar);
			PRINT("\tvPer = %f\n", vPer_ms_0);
			PRINT("\tv = %f\n", vMag_ms);
			PRINT("\tvTh = %f\n", vTh_ms);
			PRINT("\tv/vTh = %f\n", x);
			PRINT("\tphi_x = %f\n", phi_x);
			PRINT("\tpsi_x = %f\n", psi_x);
			PRINT("\tphi_x_dx = %e\n", phi_x_dx);
			PRINT("\tnu_D/nu_B = %e\n", nu_D/nu_B);
			PRINT("\tnu_D = %e\n", nu_D);
			PRINT("\tnu_B = %e\n", nu_B);
			PRINT("\tnu_D*dt = %e\n", nu_D*dt);
#if GOOSE >= 1
			PRINT("\tdt: %e, dtMin: %e\n", dt, nu_dt_max/nu_D);
#endif
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
			PRINT("dtMax reached: %f, %f\n",dt,dtMax);
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
			PRINT("\t\tdtMin reached: %e, %e\n",dt,dtMin);
			PRINT("\t\tdelta: %f\n",delta);
			PRINT("\t\tR_: %e\n",R_);
			PRINT("\t\tTOL: %f\n",TOL);
			PRINT("\t\tvPer: %f\n",p.vPer);
			PRINT("\t\teV: %f\n",p.energy_eV);
			dt = dtMin;
			FLAG = 0;
		}

#if GOOSE >= 1
		// Adjust dt such that the poloidal steps are within poloidalDistanceClosed
		if(fabs(poloidalDistanceNew-poloidalDistanceOld)>=poloidalDistanceClosed && goose) {
			dt = 0.95 * dt * poloidalDistanceClosed/fabs(poloidalDistanceNew-poloidalDistanceOld);
		}

		// Accelerate with goose-ing
		if(R_<=TOL && goose) {
		if(poloidalDistanceNew-poloidalDistanceOld<0 // approaching
						&& w.vPar*poloidalStart.vPar>0 // pick the right orbit side for trapped particles
						&& poloidalDistanceNew<=dt*sqrt(pow(v_ms.r,2)+pow(v_ms.z,2)) ) { // if we done one more dt step it will be too far 
			if(poloidalDistanceNew<poloidalDistanceClosed) { // close enough to call it a complete orbit

				// Calculate the goose-ing factor
				gooseFac = nu_dt_max/(nu_D_dt_orbitTotal);
				if(nu_dt_max/(nu_E_dt_orbitTotal)<gooseFac)
						gooseFac = nu_dt_max/(nu_E_dt_orbitTotal);

				// Catch for orbits that cannot be accelerated
				if(gooseFac<1) {
#if DEBUGLEVEL >=2 
					PRINT("\t\tOrbit cannot be accelerated, gooseFac: %i\n", gooseFac);
					PRINT("\t\tnu_D_dt_orbitTotal: %f\n", nu_D_dt_orbitTotal);
					PRINT("\t\tnu_E_dt_orbitTotal: %f\n", nu_E_dt_orbitTotal);
#endif
					goose=false;

					// Reset state to that of the start of this poloidal orbit
					// since this orbit cannot be accelerated and we need to 
					// apply the diffusion as we go along the orbit. This is
					// usually the impact of a particle thermalizing.
					
					w = poloidalStart;
					t = t_poloidalStart;
					dt = dtMin;
				} 
				else {

					// Ensure the goose-ing does not push the particle past "runTime"
					if(t+(gooseFac-1)*t_orbitTotal>runTime) {
						int gooseFacAdj = (runTime-t)/t_orbitTotal;
#if DEBUGLEVEL >= 3
						PRINT("\t\n");
						PRINT("\tAdjusting goose-ing factor from %i to %i\n",gooseFac,gooseFacAdj);
						PRINT("\tt_orbitTotal: %e\n",t_orbitTotal);
						PRINT("\tdt: %e\n",dt);
						PRINT("\tpoloidalDistanceOld: %f\n",poloidalDistanceOld);
						PRINT("\tpoloidalDistanceNew: %f\n",poloidalDistanceNew);
#endif
						if(t_orbitTotal<dt) exit(0);
						gooseFac = gooseFacAdj;
					}

					if(gooseFac<1) {
#if DEBUGLEVEL >= 3
						PRINT("\tAdjusting goose-ing factor from %i to %i\n",gooseFac,1);
						PRINT("\tdt: %e\n",dt);
#endif
						gooseFac = 1;
					}
					t += (gooseFac-1)*t_orbitTotal;
					t_poloidalStart = t;
					nPoloidalOrbits++;

					double pitch_0 = w.vPar / vMag_ms;
#if PITCH_SCATTERING >= 1
					double pitch_1 = newPitch ( pitch_0, nu_D_dt_orbitTotal*gooseFac );
#else
					double pitch_1 = pitch_0;
#endif

					double energy_0_eV = 0.5 * p.amu*_mi * pow(vMag_ms,2) / _e;
#if ENERGY_SCATTERING >= 1
					double energy_1_eV = 
							newEnergy ( energy_0_eV, vTh_ms, p.amu, T_background_eV, nu_B, nu_E, 
											nu_E_dt_orbitTotal*gooseFac );
#else
					double energy_1_eV = energy_0_eV;
#endif

#if DEBUGLEVEL >= 3
					PRINT("\t\n");
					PRINT("\tt: %f\n",t);
					PRINT("\tdt: %e\n",dt);
					PRINT("\tt_orbitTotal: %e\n",t_orbitTotal);
					PRINT("\tnu_D_dt: %e\n",nu_D_dt_orbitTotal);
					PRINT("\tnu_E_dt: %e\n",nu_E_dt_orbitTotal);
					PRINT("\tgooseFac: %i\n", gooseFac);
					PRINT("\tenergy0[eV]: %e\n", energy_0_eV);
					PRINT("\tenergy1[eV]: %e\n", energy_1_eV);
					PRINT("\tpitch0: %f\n", pitch_0);
					PRINT("\tpitch1: %f\n", pitch_1);
					PRINT("\tpoloidalDistanceOld: %f\n",poloidalDistanceOld);
					PRINT("\tpoloidalDistanceNew: %f\n",poloidalDistanceNew);
					PRINT("\tpoloidalPoints: %i\n",poloidalPoints);
					PRINT("\tnu_D_dt_orbitTotal*gooseFac: %f\n", nu_D_dt_orbitTotal*gooseFac);
					PRINT("\tnu_E_dt_orbitTotal*gooseFac: %f\n", nu_E_dt_orbitTotal*gooseFac);
					PRINT("\tpitch change: %f\%\n", fabs((pitch_1-pitch_0)/2.0)*100);
					PRINT("\tenergy change: %f\%\n", fabs((energy_0_eV-energy_1_eV)/energy_0_eV)*100);
					PRINT("\tv/vTh: %f\%\n",sqrt(energy_1_eV)/sqrt(3.0/2.0*T_background_eV)*100);
#endif
					REAL vPer_ms_0 = get_vPer_ms ( w, p, textures, spans );
					update_vPar_and_mu ( pitch_1, energy_1_eV, vPer_ms_0, w, p );

#ifndef __CUDA_ARCH__
#ifdef __SAVE_ORBITS__
					rOut.push_back(w.r);
					pOut.push_back(w.p);
					zOut.push_back(w.z);
					E_eVOut.push_back(energy_1_eV);
					pitchOut.push_back(pitch_1);
					dtOut.push_back(dt);
#endif
#endif
					t_orbitTotal = 0;
					nu_D_dt_orbitTotal = 0;
					nu_E_dt_orbitTotal = 0;
					dt = dtMin;

					// Update the position with which respect to we test for poloidal orbit completion

					poloidalStart = w;
					poloidalDistanceOld = poloidalDistanceNew;
					poloidalDistanceNew = 0;
					poloidalPoints = 0;
				} // End of else (gooseFac<1) 
			}
			else {
				dt=0.95*poloidalDistanceNew/sqrt(pow(v_ms.r,2)+pow(v_ms.z,2));
			} // End of else (poloidalDistanceNew<poloidalDistanceClosed)

		}
		} // End of (R_<=TOL && goose)
#endif

		// Update dt for the case where nu * dt is not << 1 even
		// without goosing, i.e., uber cool particles E/T << 1
#if GOOSE >= 1
		if(!goose && dt_nu>0) {
#else
		if(dt_nu>0) {
#endif
			dt = dt_nu;
			dt_nu = -1;
		}

		ii++;

	} // end while(FLAG=1)

	update_p ( w, p );

#if DEBUGLEVEL >= 1
	PRINT("\n");
    PRINT("nSteps: %i with %i accepted.\n", ii,jj);
#endif

#if GOOSE >= 1
#if DEBUGLEVEL >= 2
	PRINT("\tnPoloidalOrbits: %i\n", nPoloidalOrbits);
#endif
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
	NcVar *nc_E_eV = dataFile.add_var ("E_eV", ncFloat, nDim );
	NcVar *nc_pitch = dataFile.add_var ("pitch", ncFloat, nDim );
	NcVar *nc_dt = dataFile.add_var ("dt", ncFloat, nDim );

	NcVar *nc_vPar = dataFile.add_var ("vPar", ncFloat, sDim );
	NcVar *nc_nu_B = dataFile.add_var ("nu_B", ncFloat, sDim );
	NcVar *nc_stat = dataFile.add_var ("stat", ncInt, sDim );

	nc_rOrb->put(&rOut[0],nSteps);
	nc_pOrb->put(&pOut[0],nSteps);
	nc_zOrb->put(&zOut[0],nSteps);
	nc_E_eV->put(&E_eVOut[0],nSteps);
	nc_pitch->put(&pitchOut[0],nSteps);
	nc_dt->put(&dtOut[0],nSteps);

	nc_stat->put(&p.status,1);
	nc_vPar->put(&p.vPar,1);
	nc_nu_B->put(&nu_B,1);

#endif
#endif

    return err;

} // end void move_particle()

} // end namespace

#endif
