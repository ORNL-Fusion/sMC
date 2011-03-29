#ifndef GC_EQOM_HPP_
#define GC_EQOM_HPP_

#include "particle.hpp"
#include "eqdsk.hpp"

int vGC ( C_rkGCparticle &p, Ceqdsk &eqdsk);
int euler ( const C_rkGCparticle &p1, C_rkGCparticle &p2, const float dt );
int average_vGC 
	(const C_rkGCparticle &p1, 
	 const C_rkGCparticle &p2, 
	 const C_rkGCparticle &p3, 
	 const C_rkGCparticle &p4, 
	 C_rkGCparticle &pf, const float dt);

#endif
