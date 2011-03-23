#include <iostream>
#include <fstream>
#include <vector>
#include "eqdsk.hpp"
#include "deriv.hpp"

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

		fpol.resize(nw);
		pres.resize(nw);
		ffprim.resize(nw);
		pprime.resize(nw);

		for (int i=0;i<nw;i++)
			inFile >> fpol[i];

		for (int i=0;i<nw;i++)
			inFile >> pres[i];

		for (int i=0;i<nw;i++)
			inFile >> ffprim[i];

		for (int i=0;i<nw;i++)
			inFile >> pprime[i];

        psizr.resize(nw);
        for (int i=0;i<nw;i++)
            psizr[i].resize(nh);

		for (int j=0;j<nh;j++)
		{
			for (int i=0;i<nw;i++)
			{
				inFile >> psizr[i][j];
			}
		}	 

        qpsi.resize(nw);

		for (int i=0;i<nw;i++)
			inFile >> qpsi[i];

        inFile >> nbbbs >> limitr;

        rbbbs.resize(nbbbs);
        zbbbs.resize(nbbbs);

        for (int i=0;i<nbbbs;i++)
        {
            inFile >> rbbbs[i] >> zbbbs[i];
        }

        rlim.resize(limitr);
        zlim.resize(limitr);

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

        r.resize(nw);
        z.resize(nh);

        for (int i=0;i<nw;i++)
            r[i] = i * dr + rleft;

        for (int j=0;j<nh;j++)
            z[j] = j * dz + zmid - zdim / 2.0;

        fluxGrid.resize(nw);
        for (int i=0;i<nw;i++)
            fluxGrid[i] = i * fStep + simag;

        br.resize(nw);
        for (int i=0;i<nw;i++)
            br[i].resize(nh);

        bz.resize(nw);
        for (int i=0;i<nw;i++)
            bz[i].resize(nh);

        bp.resize(nw);
        for (int i=0;i<nw;i++)
            bp[i].resize(nh);

        bmag.resize(nw);
        for (int i=0;i<nw;i++)
            bmag[i].resize(nh);

		for (int i=0;i<nw;i++)
		{
			vector<float> tmpData (nh);
			vector<float> tmpRes (nw);
			for (int j=0;j<nh;j++) tmpData[j] = psizr[i][j];

        	tmpRes = deriv ( tmpData, dz );
		}

	}
	else {

		cout << "ERROR: file '" << fName << "' does not exist?" << endl;
		return 1;

	}

	return 0;

}
