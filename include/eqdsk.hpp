#ifndef EQDSK_HPP_
#define EQDSK_HPP_

#include <vector>
#include "boost/multi_array.hpp"

class Ceqdsk {
		
	private:

		float xdum;
		int nCol_, nRow_;

		typedef boost::multi_array<float,2> arr2D_;
		typedef std::vector<float> arr1D_;

	public:

        // variables
		char *header;
		int idum,  nbbbs, limitr;
		float rdim,zdim,rcentr,rleft,zmid;
        float rmaxis,zmaxis,simag,sibry,bcentr;
        float current;

		arr2D_ psizr, br, bz, bp, bmag, fpolzr;
		arr1D_ fpol, pres, ffprim, pprime,
			r, z, fluxGrid, qpsi, rbbbs, zbbbs, rlim, zlim;

        float dr, dz;
        bool ascending_flux;

		arr2D_ bCurvature_r, bCurvature_p, bCurvature_z,
			   gradB_r, gradB_z, bDotGradB, 
			   bGradient_r, bGradient_p, bGradient_z;

        // functions 
		int read_file ( std::string );
		int write_ncfile ( std::string );
		int bForceTerms ();

};

#endif


