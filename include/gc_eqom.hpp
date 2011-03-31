#ifndef GC_EQOM_HPP_
#define GC_EQOM_HPP_

#include "particle.hpp"
#include "eqdsk.hpp"
#include "constants.hpp"

Crk vGC ( const REAL dt, const Crk &p0, const REAL mu, Ceqdsk &eqdsk, int &err );
int euler ( const C_rkGCparticle &p1, C_rkGCparticle &p2, const REAL dt );
int average_vGC 
	(const C_rkGCparticle &p1, 
	 const C_rkGCparticle &p2, 
	 const C_rkGCparticle &p3, 
	 const C_rkGCparticle &p4, 
	 C_rkGCparticle &pf, const REAL dt);

#endif
