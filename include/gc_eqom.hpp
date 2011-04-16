#ifndef GC_EQOM_HPP_
#define GC_EQOM_HPP_

#include "rk_gc_particle.hpp"
#include "eqdsk.hpp"
#include "constants.hpp"
#include "interp.hpp"

Crk vGC ( const REAL dt, const Crk &p0, const REAL mu, 
            const Ceqdsk &eqdsk, const CinterpSpans &span, int &err );

#endif
