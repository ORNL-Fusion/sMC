#include <iostream>
#include <netcdfcpp.h>
#include "particle.hpp"
#include <vector>

using namespace std;

int read_pl ( const string fName, vector<C_GCparticle> &particles )  {

	cout << "Reading particle list file " << fName << endl;

	NcFile dataFile ( fName.c_str(), NcFile::ReadOnly );

	if(!dataFile.is_valid()) {
			cout << "Couldn't open file.\n";
			return 1;
	}

	NcDim *nc_nP;
	nc_nP = dataFile.get_dim("nP");

	int nP = nc_nP->size();
	particles.resize(nP);

	cout << "\tnP: " << nP << endl;

	NcVar *nc_r, *nc_z, *nc_vPer, *nc_vPar, *nc_weight, *nc_status;

	nc_r = dataFile.get_var("R");
	nc_z = dataFile.get_var("z");
	nc_vPer = dataFile.get_var("vPer");
	nc_vPar = dataFile.get_var("vPar");
	nc_weight = dataFile.get_var("weight");
	nc_status = dataFile.get_var("status");

	float *tmp_r = new float[nP];
	float *tmp_z = new float[nP];
	float *tmp_vPer = new float[nP];
	float *tmp_vPar = new float[nP];
	float *tmp_weight = new float[nP];
	int *tmp_status = new int[nP];

	nc_r->get(&tmp_r[0],nP);
	nc_z->get(&tmp_z[0],nP);
	nc_vPer->get(&tmp_vPer[0],nP);
	nc_vPar->get(&tmp_vPar[0],nP);
	nc_weight->get(&tmp_weight[0],nP);
	nc_status->get(&tmp_status[0],nP);

	for(int p=0;p<nP;p++) {

			particles[p].r = tmp_r[p];
			particles[p].z = tmp_z[p];
			particles[p].vPer = tmp_vPer[p];
			particles[p].vPar = tmp_vPar[p];
			particles[p].weight = tmp_weight[p];
			particles[p].status = tmp_status[p];
	}

	cout << "DONE" << endl;

	return 0;
}
