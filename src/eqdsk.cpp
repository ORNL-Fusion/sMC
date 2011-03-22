#include <iostream>
#include <fstream>
#include "eqdsk.h"

using namespace std;

int Ceqdsk::read_file ( string fName ) {

	ifstream inFile ( fName.c_str() );

	if ( inFile.good() ) {

		int headerLength = 48;
		header = new char [headerLength];

		inFile.read ( header, headerLength );

		inFile >> idum >> nw >> nh;	
		inFile >> rdim >> zdim >> rcentr >> rleft >> zmid;
        inFile >> rmaxis >> zmaxis >> simag >> sibry >> bcentr;
        inFile >> current >> simag >> xdum >> rmaxis >> xdum; 
        inFile >> zmaxis >> xdum >> sibry >> xdum >> xdum; 

		cout << "header: " << header << endl;
		cout << "idum: " << idum << endl;
		cout << "nw: " << nw << endl;
		cout << "nh: " << nh << endl;

		cout << "rdim: " << rdim << endl;
		cout << "zdim: " << zdim << endl;
		cout << "rcentr: " << rcentr << endl;
		cout << "rleft: " << rleft << endl;

		cout << "rmaxis: " << rmaxis << endl;
		cout << "zmaxis: " << zmaxis << endl;
		cout << "simag: " << simag << endl;
		cout << "sibry: " << sibry << endl;
		cout << "bcentr: " << bcentr << endl;

		for (int i=0;i<nw;i++)
		{
			float number;
			inFile >> number;
			fpol.push_back(number);
		}

		for (int i=0;i<nw;i++)
		{
			float number;
			inFile >> number;
			pres.push_back(number);
		}

		for (int i=0;i<nw;i++)
		{
			float number;
			inFile >> number;
			ffprim.push_back(number);
		}

		for (int i=0;i<nw;i++)
		{
			float number;
			inFile >> number;
			pprime.push_back(number);
		}

		psizr.reserve(nh);
		for (int j=0;j<nh;j++)
		{
			vector<float> temp;
			for (int i=0;i<nh;i++)
			{
				float number;
				inFile >> number;
				temp.push_back(number);
			}
			psizr.push_back(temp);
		}


		/*for (int j=0;j<nh;j++)
		{
			for (int i=0;i<nw;i++)
			{
				inFile >> psizr[i][j];
			}
		}*/	 
		/*
		qpsi ( nw ), &
            r ( nw ), z ( nh ), fluxGrid ( nw ), fpol_(nw), &
            fluxGrid_(nw) )
			*/
 
	}
	else {

		cout << "ERROR: file '" << fName << "' does not exist?" << endl;
		return 1;

	}

	return 0;

}
