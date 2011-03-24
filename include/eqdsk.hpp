#ifndef EQDSK_HPP_
#define EQDSK_HPP_

#include "Data.hpp"
#include <vector>

using namespace std;

class Ceqdsk {

	public:

        // variables
		char *header;
		int idum, nCol, nRow, nbbbs, limitr;
		float rdim,zdim,rcentr,rleft,zmid;
        float rmaxis,zmaxis,simag,sibry,bcentr;
        float current;
		vector<float> fpol, pres, ffprim, pprime,
			r, z, fluxGrid, qpsi, rbbbs, zbbbs, rlim, zlim;
		vector<vector<float> > psizr, br, bz, bp, bmag, fpolzr;
        float dr, dz;
        bool ascending_flux;

        // functions 
		int read_file ( string );
		int write_ncfile ( string );

	private:

		float xdum;

};

#endif


