#ifndef EQDSK_H
#define EQDSK_H

#include <vector>

using namespace std;

class Ceqdsk {

	public:

		char *header;
		int idum, nw, nh;
		float rdim,zdim,rcentr,rleft,zmid;
        float rmaxis,zmaxis,simag,sibry,bcentr;
        float current;

		vector<float> fpol, pres, ffprim, pprime;
		vector<vector<float>> psizr;
        float *qpsi, *r, *z, *fluxGrid, *fpol_, *fluxGrid_;
 
		int read_file ( string );

	private:

		float xdum;

};

#endif


