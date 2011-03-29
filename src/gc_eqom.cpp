#include "eqdsk.hpp"
#include "particle.hpp"
#include <iostream>


// bi-linear interpolation
// see wikipedia ;)

float bilinear_interp_data 
    ( float i, float j, int i1, int i2, int j1, int j2, eqdsk::arr2D_ &data ) {

	float f11 = data[i1][j1];
	float f21 = data[i2][j1];
	float f12 = data[i1][j2];
	float f22 = data[i2][j2];

	// (x2-x1)(y2-y1) == 1 since i'm using indices

	float dataOut = f11 * (i2-i)*(j2-j)
			+ f21 * (i-i1)*(j2-j)
			+ f12 * (i2-i)*(j-j1)
			+ f22 * (i-i1)*(j-j1); 

    return dataOut;
}

int eq_o_motion ( Cparticle &p0, Cparticle &p1, Ceqdsk &eqdsk ) {

	// get background data(s) at particle location
	
	p0.row = (p0.r - eqdsk.r.front()) / ( eqdsk.r.back() - eqdsk.r.front() )
		 		* eqdsk.r.size();
	p0.col = (p0.z - eqdsk.z.front()) / ( eqdsk.z.back() - eqdsk.z.front() )
		 		* eqdsk.z.size();

	int x1 = floor(p0.row);
	int x2 = ceil(p0.row);
	int y1 = floor(p0.col);
	int y2 = ceil(p0.col);

    // Check if particle is off grid	
    if( x1<0 || x2>=eqdsk.nRow_ ||
        y1<0 || y2>=eqdsk.nCol_ ) {

        std::cout << "\tERROR: position outside eqdsk grid." << std::endl;
        std::cout << "\tx1: " << x1 << std::endl;
        std::cout << "\tx2: " << x2 << std::endl;
        std::cout << "\ty1: " << y1 << std::endl;
        std::cout << "\ty2: " << y2 << std::endl;

        return 1;
    }

	std::cout << "\tp0.row: "<<p0.row<<std::endl;
	std::cout << "\tp0.col: "<<p0.col<<std::endl;

    float bmag_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.bmag );

    float br_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.br );
    float bp_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.bp );
    float bz_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.bz );

    float bCurvature_r_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.bCurvature_r );
    float bCurvature_p_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.bCurvature_p );
    float bCurvature_z_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.bCurvature_z );

    float bGradient_r_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.bGradient_r );
    float bGradient_p_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.bGradient_p );
    float bGradient_z_p0 = bilinear_interp_data ( p0.row, p0.col, x1, x2, y1, y2, eqdsk.bGradient_z );

	return 0;
}


