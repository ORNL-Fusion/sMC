#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "particle.hpp"
#include <vector>
#include <iostream>

using namespace std;

int copy_data ( vector<Cgc_particle> &H_particles ) {


    //cout << "First CUDA call :)" << endl;

    //cout << "\tCopying particle list from HOST -> DEVICE ... ";
    //thrust::device_vector<Cgc_particle> D_particles ( H_particles.begin(), H_particles.end() );
    //cout << "DONE" << endl;

    //cout << "\tCopying particle list DEVICE -> HOST ... ";
    //thrust::host_vector<Cgc_particle> H_tmp ( D_particles.begin(), D_particles.end() );
    //cout << "DONE" << endl;

    //float *D_test2D;

    //cudaMalloc( (void**) &D_test2D, eqdsk.nRow * eqdsk.nCol * sizeof(float *));


    return 0;
}

