#ifndef READ_PL_HPP_
#define READ_PL_HPP_

#include "particle.hpp"
#include <vector>

int read_pl ( const std::string fName, std::vector<Cgc_particle> &particles );
int write_pl ( const std::string fName, std::vector<Cgc_particle> &particles );

#endif
