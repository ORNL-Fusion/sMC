#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include <iostream>
#include "constants.hpp"
#include <iomanip>

class Cgc_particle {

	public:
		REAL r, p, z, vPer, vPar, weight, mu, energy_eV;
		int status, amu, Z;
};

#endif
