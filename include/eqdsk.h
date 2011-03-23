#ifndef EQDSK_H
#define EQDSK_H

using namespace std;

class Ceqdsk {

	public:

        // variables
		char *header;
		int idum, nw, nh, nbbbs, limitr;
		float rdim,zdim,rcentr,rleft,zmid;
        float rmaxis,zmaxis,simag,sibry,bcentr;
        float current;
		float *fpol, *pres, *ffprim, *pprime;
		float **psizr, *r, *z, *fluxGrid;
        float *qpsi, *rbbbs, *zbbbs, *rlim, *zlim;
        float dr, dz, **br, **bz, **bp, **bmag;
        bool ascending_flux;

        // functions 
		int read_file ( string );

	private:

		float xdum;

};

#endif


