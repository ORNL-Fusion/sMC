#include "particle.hpp"
#include "constants.hpp"
#include "eqdsk.hpp"

int copy_particles_to_device (std::vector<Cgc_particle> &H_particles);
REAL* copy_2D_to_device ( eqdsk::arr2D_ &data2D, const unsigned int nRow, const unsigned int nCol );
REAL* copy_1D_to_device ( eqdsk::arr1D_ &h_data1D, const unsigned int n );


