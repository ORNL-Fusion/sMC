#ifndef GC_EQOM_HPP_
#define GC_EQOM_HPP_

#include "particle.hpp"
#include "eqdsk.hpp"
#include "constants.hpp"

CK& vGC ( REAL dt, CK &p0, REAL mu, REAL vPar0, Ceqdsk &eqdsk, int err );
//int vGC ( C_rkGCparticle &p, Ceqdsk &eqdsk);
int euler ( const C_rkGCparticle &p1, C_rkGCparticle &p2, const REAL dt );
int average_vGC 
	(const C_rkGCparticle &p1, 
	 const C_rkGCparticle &p2, 
	 const C_rkGCparticle &p3, 
	 const C_rkGCparticle &p4, 
	 C_rkGCparticle &pf, const REAL dt);

#endif
