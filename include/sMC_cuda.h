#ifndef SMC_CUDA_H_
#define SMC_CUDA_H_

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "particle.hpp"
#include <vector>

int copy_data ( std::vector<Cgc_particle> &H_particles );

#endif
