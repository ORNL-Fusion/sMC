#include "eqdsk.hpp"
#include "particle.hpp"
#include <iostream>

Cparticle eq_o_motion ( Cparticle p0, Ceqdsk &eqdsk ) {

	Cparticle p1;

	// get b at particle location
	
	float bHere;		

	p0.row = (p0.r - eqdsk.r.front()) / ( eqdsk.r.back() - eqdsk.r.front() )
		 		* eqdsk.r.size();
	p0.col = (p0.z - eqdsk.z.front()) / ( eqdsk.z.back() - eqdsk.z.front() )
		 		* eqdsk.z.size();

	// bi-linear interpolation
	// see wikipedia ;)
		
	float fq11 = eqdsk.bmag[floor(p0.row)][floor(p0.col)];
	float fq21 = eqdsk.bmag[ceil(p0.row)][floor(p0.col)];
	float fq12 = eqdsk.bmag[floor(p0.row)][ceil(p0.col)];
	float fq22 = eqdsk.bmag[ceil(p0.row)][ceil(p0.col)];

	// (x2-x1)(y2-y1) == 1 since i'm using indices

	float x1 = floor(p0.row);
	float x2 = ceil(p0.row);
	float y1 = floor(p0.col);
	float y2 = ceil(p0.col);
	float x = p0.row;
	float y = p0.col;

	bHere = fq11 * (x2-x)*(y2-y)
			+ fq21 * (x-x1)*(y2-y)
			+ fq12 * (x2-x)*(y-y1)
			+ fq22 * (x-x1)*(y-y1); 

	std::cout << "p0.row: "<<p0.row<<std::endl;
	std::cout << "p0.col: "<<p0.col<<std::endl;
	std::cout << "bVals: "<<fq11<<" "<<fq21<<" "<<fq12<<" "<<fq22<<std::endl;
	std::cout << "bHere: "<<bHere<<std::endl;

	return p1;
}

