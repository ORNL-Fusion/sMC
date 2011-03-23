#include <iostream>
#include <fstream>
#include "eqdsk.h"

using namespace std;

int Ceqdsk::read_file ( string fName ) {

	ifstream inFile ( fName.c_str() );

	if ( inFile.good() ) {

        // read data from file

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

        fpol = new float[nw];
        pres = new float[nw];
        ffprim = new float[nw];
        pprime = new float[nw];

		for (int i=0;i<nw;i++)
			inFile >> fpol[i];

		for (int i=0;i<nw;i++)
			inFile >> pres[i];

		for (int i=0;i<nw;i++)
			inFile >> ffprim[i];

		for (int i=0;i<nw;i++)
			inFile >> pprime[i];

        psizr = new float *[nw];
        for (int i=0;i<nw;i++)
            psizr[i] = new float[nh];

		for (int j=0;j<nh;j++)
		{
			for (int i=0;i<nw;i++)
			{
				inFile >> psizr[i][j];
			}
		}	 

        qpsi = new float[nw];

		for (int i=0;i<nw;i++)
			inFile >> qpsi[i];

        inFile >> nbbbs >> limitr;

        rbbbs = new float[nbbbs];
        zbbbs = new float[nbbbs];

        for (int i=0;i<nbbbs;i++)
        {
            inFile >> rbbbs[i] >> zbbbs[i];
        }

        rlim = new float[limitr];
        zlim = new float[limitr];

        for (int i=0;i<limitr;i++)
        {
            inFile >> rlim[i] >> zlim[i];
        }

        inFile.close ();


        // Calculate other quantities from read data

        dr   = rdim / ( nw - 1 );
        dz   = zdim / ( nh - 1 );

        ascending_flux = true;
        if ( sibry > simag )
            ascending_flux = false;

        float fStep = ( sibry - simag ) / ( nw - 1 );

        r = new float[nw];
        z = new float[nh];

        for (int i=0;i<nw;i++)
            r[i] = i * dr + rleft;

        for (int j=0;j<nh;j++)
            z[j] = j * dz + zmid - zdim / 2.0;

        fluxGrid = new float[nw];
        for (int i=0;i<nw;i++)
            fluxGrid[i] = i * fStep + simag;

        br = new float *[nw];
        for (int i=0;i<nw;i++)
            br[i] = new float[nh];

        bz = new float *[nw];
        for (int i=0;i<nw;i++)
            bz[i] = new float[nh];

        bp = new float *[nw];
        for (int i=0;i<nw;i++)
            bp[i] = new float[nh];

        bmag = new float *[nw];
        for (int i=0;i<nw;i++)
            bmag[i] = new float[nh];

	}
	else {

		cout << "ERROR: file '" << fName << "' does not exist?" << endl;
		return 1;

	}

	return 0;

}
