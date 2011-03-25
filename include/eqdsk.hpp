#ifndef EQDSK_HPP_
#define EQDSK_HPP_

#include <vector>
#include "data2D.hpp"

class Ceqdsk {

	public:

        // variables
		char *header;
		int idum, nCol, nRow, nbbbs, limitr;
		float rdim,zdim,rcentr,rleft,zmid;
        float rmaxis,zmaxis,simag,sibry,bcentr;
        float current;
		std::vector<float> fpol, pres, ffprim, pprime,
			r, z, fluxGrid, qpsi, rbbbs, zbbbs, rlim, zlim;
		std::vector<std::vector<float> > psizr, br, bz, bp, bmag, fpolzr;
        float dr, dz;
        bool ascending_flux;

        // functions 
		int read_file ( std::string );
		int write_ncfile ( std::string );
		int bCurvature ();

	private:

		float xdum;

};

#endif


