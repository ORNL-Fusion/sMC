#include <iostream>
#include "eqdsk.hpp"

using namespace std;

int main ()
{

	cout << "sMC ... c++ & CUDA version :)" << endl;


	string fName = "data/g129x129_1051206002.01120.cmod";
	Ceqdsk eqdsk;	
	int stat = eqdsk.read_file ( fName );
	stat = eqdsk.write_ncfile ( "output/bField.nc" );

	cout << "End of program :)" << endl;

	return 0;
}
