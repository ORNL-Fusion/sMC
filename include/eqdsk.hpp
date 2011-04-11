#ifndef EQDSK_HPP_
#define EQDSK_HPP_

#include <vector>
#include "boost/multi_array.hpp"

#include "constants.hpp"

//namespace eqdsk {
//    typedef boost::multi_array<REAL,2> boost::multi_array<REAL,2>;
//    typedef std::vector<REAL> std::vector<REAL>;
//}

//using namespace eqdsk;

class Ceqdsk {
		
	private:

		REAL xdum;
	    int nCol_, nRow_;
		boost::multi_array<REAL,2> psizr_;
		std::vector<REAL> r_, z_, fpol_,fluxGrid_;
		REAL dr_, dz_;

	public:

        // variables

		int nRow, nCol;

		char *header;
		int idum,  nbbbs, limitr;
		REAL rdim,zdim,rcentr,rleft,zmid;
        REAL rmaxis,zmaxis,simag,sibry,bcentr;
        REAL current;

		boost::multi_array<REAL,2> psizr, br, bz, bp, bmag, fpolzr;
		std::vector<REAL> pres, ffprim, pprime,
			r, z, qpsi, rbbbs, zbbbs, rlim, zlim;

        REAL dr, dz;
        bool ascending_flux;

		boost::multi_array<REAL,2> bCurvature_r, bCurvature_p, bCurvature_z,
			   gradB_r, gradB_z, bDotGradB, 
			   bGradient_r, bGradient_p, bGradient_z;

		// Default constructor
		Ceqdsk () {}; 

		// Copy constructor
		Ceqdsk ( const Ceqdsk &eqdsk ) {*this = eqdsk;}	

        // functions 
		int read_file ( std::string );
		int write_ncfile ( std::string );
		int bForceTerms ( const int _Z, const int amu);
		int calc_b ( const unsigned int nrow, const unsigned int ncol );
};

#endif


