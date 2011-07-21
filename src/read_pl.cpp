#include <iostream>
#include <netcdfcpp.h>
#include "particle.hpp"
#include <vector>

using namespace std;

int read_pl ( const string fName, vector<Cgc_particle> &particles )  {

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

	NcVar *nc_r, *nc_p, *nc_z, *nc_vPer, 
		  *nc_vPar, *nc_weight, *nc_status, 
		  *nc_energy_eV, *nc_mu, *nc_amu, *nc_Z;

	nc_r = dataFile.get_var("r");
	nc_p = dataFile.get_var("p");
	nc_z = dataFile.get_var("z");
	nc_vPer = dataFile.get_var("vPer");
	nc_vPar = dataFile.get_var("vPar");
	nc_weight = dataFile.get_var("weight");
	nc_status = dataFile.get_var("status");
	nc_energy_eV = dataFile.get_var("E_eV");
	nc_mu = dataFile.get_var("mu");
	nc_amu = dataFile.get_var("amu");
	nc_Z = dataFile.get_var("Z");

	float *tmp_r = new float[nP];
	float *tmp_p = new float[nP];
	float *tmp_z = new float[nP];
	float *tmp_vPer = new float[nP];
	float *tmp_vPar = new float[nP];
	float *tmp_weight = new float[nP];
	int *tmp_status = new int[nP];
	float *tmp_energy_eV = new float[nP];
	float *tmp_mu = new float[nP];
	int *tmp_amu = new int[nP];
	int *tmp_Z = new int[nP];

	nc_r->get(&tmp_r[0],nP);
	nc_p->get(&tmp_p[0],nP);
	nc_z->get(&tmp_z[0],nP);
	nc_vPer->get(&tmp_vPer[0],nP);
	nc_vPar->get(&tmp_vPar[0],nP);
	nc_weight->get(&tmp_weight[0],nP);
	nc_status->get(&tmp_status[0],nP);
	nc_energy_eV->get(&tmp_energy_eV[0],nP);
	nc_mu->get(&tmp_mu[0],nP);
	nc_amu->get(&tmp_amu[0],nP);
	nc_Z->get(&tmp_Z[0],nP);

	for(int p=0;p<nP;p++) {

			particles[p].r = tmp_r[p];
			particles[p].p = tmp_p[p];
			particles[p].z = tmp_z[p];
			particles[p].vPer = tmp_vPer[p];
			particles[p].vPar = tmp_vPar[p];
			particles[p].weight = tmp_weight[p];
			particles[p].status = tmp_status[p];
			particles[p].energy_eV = tmp_energy_eV[p];
			particles[p].mu = tmp_mu[p];
			particles[p].amu = tmp_amu[p];
			particles[p].Z = tmp_Z[p];

	}

	cout << "DONE" << endl;

	return 0;
}

int write_pl ( const string fName, vector<Cgc_particle> &particles )  {

	cout << "Writing particle list file " << fName << endl;

	NcFile dataFile ( fName.c_str(), NcFile::Replace );

	if(!dataFile.is_valid()) {
			cout << "Couldn't open file.\n";
			return 1;
	}

	NcDim *nc_nP = dataFile.add_dim("nP", particles.size() );

	int nP = particles.size();

	cout << "\tnP: " << nP << endl;

	NcVar *nc_r, *nc_p, *nc_z, 
		  *nc_vPer, *nc_vPar, *nc_weight, 
		  *nc_status, *nc_energy_eV, *nc_amu, *nc_Z, *nc_mu;

	nc_r = dataFile.add_var("r", ncFloat, nc_nP);
	nc_p = dataFile.add_var("p", ncFloat, nc_nP);
	nc_z = dataFile.add_var("z", ncFloat, nc_nP);
	nc_vPer = dataFile.add_var("vPer", ncFloat, nc_nP);
	nc_vPar = dataFile.add_var("vPar", ncFloat, nc_nP);
	nc_weight = dataFile.add_var("weight", ncFloat, nc_nP);
	nc_status = dataFile.add_var("status", ncInt, nc_nP);
	nc_energy_eV = dataFile.add_var("E_eV", ncInt, nc_nP);
	nc_amu = dataFile.add_var("amu", ncInt, nc_nP);
	nc_Z = dataFile.add_var("Z", ncInt, nc_nP);
	nc_mu = dataFile.add_var("mu", ncInt, nc_nP);

	float *tmp_r = new float[nP];
	float *tmp_p = new float[nP];
	float *tmp_z = new float[nP];
	float *tmp_vPer = new float[nP];
	float *tmp_vPar = new float[nP];
	float *tmp_weight = new float[nP];
	int *tmp_status = new int[nP];
	float *tmp_energy_eV = new float[nP];
	float *tmp_mu = new float[nP];
	int *tmp_amu = new int[nP];
	int *tmp_Z = new int[nP];

	for(int p=0;p<nP;p++) {

			tmp_r[p] = particles[p].r;
			tmp_p[p] = particles[p].p;
			tmp_z[p] = particles[p].z;
			tmp_vPer[p] = particles[p].vPer;
			tmp_vPar[p] = particles[p].vPar;
			tmp_weight[p] = particles[p].weight;
			tmp_status[p] = particles[p].status;
			tmp_energy_eV[p] = particles[p].energy_eV;
			tmp_mu[p] = particles[p].mu;
			tmp_amu[p] = particles[p].amu;
			tmp_Z[p] = particles[p].Z;

	}

	nc_r->put(&tmp_r[0],nP);
	nc_p->put(&tmp_p[0],nP);
	nc_z->put(&tmp_z[0],nP);
	nc_vPer->put(&tmp_vPer[0],nP);
	nc_vPar->put(&tmp_vPar[0],nP);
	nc_weight->put(&tmp_weight[0],nP);
	nc_status->put(&tmp_status[0],nP);
	nc_energy_eV->put(&tmp_energy_eV[0],nP);
	nc_mu->put(&tmp_mu[0],nP);
	nc_amu->put(&tmp_amu[0],nP);
	nc_Z->put(&tmp_Z[0],nP);

	cout << "DONE" << endl;

	return 0;
}
