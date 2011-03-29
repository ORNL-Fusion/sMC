#ifndef EQDSK_HPP_
#define EQDSK_HPP_

#include <vector>
#include "boost/multi_array.hpp"

namespace eqdsk {
    typedef boost::multi_array<float,2> arr2D_;
    typedef std::vector<float> arr1D_;
}

using namespace eqdsk;

class Ceqdsk {
		
	private:

		float xdum;

	public:

	    int nCol_, nRow_;

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

		class interpIndex {

			public: 
				float i, j;
				int i1, i2, j1, j2;
		};

        // functions 
		int read_file ( std::string );
		int write_ncfile ( std::string );
		int bForceTerms ();
		int get_index 
			( const float r, const float z, interpIndex &index );
		float bilinear_interp 
    		( const interpIndex &index , const eqdsk::arr2D_ &data );
};

#endif


