#include <iostream>
#include "eqdsk.h"

using namespace std;

int main ()
{

	cout << "sMC ... c++ & CUDA version :)" << endl;

	cout << "Read g-eqdsk data file ..." << endl;

	string fName = "data/g129x129_1051206002.01120.cmod";
	Ceqdsk eqdsk;	
	int stat = eqdsk.read_file ( fName );

	cout << "DONE" << endl;

	return 0;
}
