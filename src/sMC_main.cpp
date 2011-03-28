#include <iostream>
#include "eqdsk.hpp"
#include "particle.hpp"
#include "gc_eqom.hpp"

using namespace std;


int main ()
{

	cout << "sMC ... c++ & CUDA version :)" << endl;


	string fName = "data/g129x129_1051206002.01120.cmod";
	Ceqdsk eqdsk;	
	int stat;
	stat = eqdsk.read_file ( fName );
	stat = eqdsk.write_ncfile ( "output/bField.nc" );
	stat = eqdsk.bForceTerms ();

	Cparticle p0, p1;

	p0.r = 0.6;
	p0.p = 0.0;
	p0.z = 0.0;

	p1 = eq_o_motion ( p0, eqdsk );


	cout << "End of program :)" << endl;

	return 0;
}
