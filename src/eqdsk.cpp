#include <iostream>
#include <fstream>
#include <vector>
#include "eqdsk.hpp"
#include "deriv.hpp"
#include "interpolation.h"
#include <cmath>
#include <netcdfcpp.h>

using namespace std;

int Ceqdsk::read_file ( string fName ) {

	cout << "Reading g-eqdsk file " << fName << endl;

	ifstream inFile ( fName.c_str() );

	if ( inFile.good() ) {

        // read data from file

		int headerLength = 48;
		header = new char [headerLength];

		inFile.read ( header, headerLength );

		inFile >> idum >> nCol >> nRow;	
		inFile >> rdim >> zdim >> rcentr >> rleft >> zmid;
        inFile >> rmaxis >> zmaxis >> simag >> sibry >> bcentr;
        inFile >> current >> simag >> xdum >> rmaxis >> xdum; 
        inFile >> zmaxis >> xdum >> sibry >> xdum >> xdum; 

		cout << "\t header: " << header << endl;
		cout << "\t idum: " << idum << endl;
		cout << "\t nCol: " << nCol << endl;
		cout << "\t nRow: " << nRow << endl;

		cout << "\t rdim: " << rdim << endl;
		cout << "\t zdim: " << zdim << endl;
		cout << "\t rcentr: " << rcentr << endl;
		cout << "\t rleft: " << rleft << endl;

		cout << "\t rmaxis: " << rmaxis << endl;
		cout << "\t zmaxis: " << zmaxis << endl;
		cout << "\t simag: " << simag << endl;
		cout << "\t sibry: " << sibry << endl;
		cout << "\t bcentr: " << bcentr << endl;

		fpol.resize(nCol);
		pres.resize(nCol);
		ffprim.resize(nCol);
		pprime.resize(nCol);

		for (int j=0;j<nCol;j++)
			inFile >> fpol[j];

		for (int j=0;j<nCol;j++)
			inFile >> pres[j];

		for (int j=0;j<nCol;j++)
			inFile >> ffprim[j];

		for (int j=0;j<nCol;j++)
			inFile >> pprime[j];

        psizr.resize(nCol);
        for (int j=0;j<nCol;j++)
            psizr[j].resize(nRow);

		for (int j=0;j<nCol;j++)
		{
			for (int i=0;i<nRow;i++)
			{
				inFile >> psizr[i][j];
			}
		}	 

        qpsi.resize(nCol);

		for (int j=0;j<nCol;j++)
			inFile >> qpsi[j];

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

        dr   = rdim / ( nCol - 1 );
        dz   = zdim / ( nRow - 1 );

        ascending_flux = true;
        if ( sibry > simag )
            ascending_flux = false;

        float fStep = ( sibry - simag ) / ( nCol - 1 );

        r.resize(nCol);
        z.resize(nRow);

        for (int j=0;j<nCol;j++)
            r[j] = j * dr + rleft;

        for (int i=0;i<nRow;i++)
            z[i] = i * dz + zmid - zdim / 2.0;

        fluxGrid.resize(nCol);
        for (int j=0;j<nCol;j++)
            fluxGrid[j] = j * fStep + simag;

        br.resize(nCol);
        for (int j=0;j<nCol;j++)
            br[j].resize(nRow);

        bz.resize(nCol);
        for (int j=0;j<nCol;j++)
            bz[j].resize(nRow);

        bp.resize(nCol);
        for (int j=0;j<nCol;j++)
            bp[j].resize(nRow);

        bmag.resize(nCol);
        for (int j=0;j<nCol;j++)
            bmag[j].resize(nRow);

		// br = -dpsi/dz * 1/r
		for (int j=0;j<nCol;j++)
		{
			vector<float> tmpData (nRow);
			vector<float> tmpRes (nCol);
			for (int i=0;i<nRow;i++) 
				tmpData[i] = psizr[i][j];

        	tmpRes = deriv ( tmpData, dz );

			for (int i=0;i<nRow;i++)
				br[i][j] = -tmpRes[i] / r[j];

		}

		cout << "\t dr: " << dr << endl;
		cout << "\t dz: " << dz << endl;

		// bz = dpsi/dr * 1/r
		for (int i=0;i<nRow;i++)
		{
			vector<float> tmpData (nCol);
			vector<float> tmpRes (nRow);
			for (int j=0;j<nCol;j++) 
				tmpData[j] = psizr[i][j];

        	tmpRes = deriv ( tmpData, dr );
			for (int j=0;j<nCol;j++)
				bz[i][j] = tmpRes[j] / r[j];

		}


		// Interpolate fpol from fluxGrid to r,z 2D space
		// Using ALGLIB

		// Initialize the AGLIB data arrays
		alglib::real_1d_array AG_fluxGrid;
		alglib::real_1d_array AG_fpol;
		alglib::spline1dinterpolant AG_s;

		// AGLIB is double only, so copy float vectors to double
		std::vector<double> fluxGrid_dbl(fluxGrid.begin(),fluxGrid.end());
		std::vector<double> fpol_dbl(fpol.begin(),fpol.end());

		// Population the ALGLIB arrays
		AG_fluxGrid.setcontent(nCol,&fluxGrid_dbl[0]);
		AG_fpol.setcontent(nCol,&fpol_dbl[0]);

		// Build the spline
		alglib::spline1dbuildcubic ( AG_fluxGrid, AG_fpol, AG_s );

		// Calculate fpol on the 2D r,z mesh
		fpolzr.resize(nCol);
		for (int j=0;j<nCol;j++)
			fpolzr[j].resize(nRow);

		for (int j=0;j<nCol;j++) {
			for(int i=0;i<nRow;i++) {

				fpolzr[i][j] = alglib::spline1dcalc(AG_s,psizr[i][j]);
				bp[i][j] = fpolzr[i][j] / r[i];

			}
		}

		// Magnitude of b
		for (int j=0;j<nCol;j++) {
			for (int i=0;i<nRow;i++) {
				bmag[i][j] = 
					sqrt ( pow(br[i][j],2)+pow(bp[i][j],2)+pow(bz[i][j],2) );
				//cout <<i<<" "<<j<<" "<<bmag[i][j] <<" "<<psizr[i][j]<<endl;

			}
		}

	}
	else {

		cout << "ERROR: file '" << fName << "' does not exist?" << endl;
		return 1;

	}

	cout << "DONE." << endl;

	return 0;

}

int Ceqdsk::write_ncfile ( const string fName ) {

	cout << "Writing file " << fName << endl;

	NcFile dataFile ( &fName[0], NcFile::Replace );
	if (!dataFile.is_valid()) {
		cout << "ERROR: Could not open nc file for writing." << endl;
		return 1;
	}	

	NcDim *rDim = dataFile.add_dim("nr",nCol);
	NcDim *zDim = dataFile.add_dim("nz",nRow);

	NcVar *nc_r = dataFile.add_var ("r", ncFloat, rDim );
	NcVar *nc_z = dataFile.add_var ("z", ncFloat, zDim );

	nc_r->put(&r[0],nCol);
	nc_z->put(&z[0],nRow);

	NcVar *nc_br = dataFile.add_var ("br", ncFloat, zDim, rDim );
	NcVar *nc_bp = dataFile.add_var ("bp", ncFloat, zDim, rDim );
	NcVar *nc_bz = dataFile.add_var ("bz", ncFloat, zDim, rDim );
	NcVar *nc_bmag = dataFile.add_var ("bmag", ncFloat, zDim, rDim );
	NcVar *nc_psizr = dataFile.add_var ("psizr", ncFloat, zDim, rDim );

	nc_br->put(&br[0][0],nRow,nCol);
	nc_bp->put(&bp[0][0],nRow,nCol);
	nc_bz->put(&bz[0][0],nRow,nCol);
	nc_bmag->put(&bmag[0][0],nRow,nCol);
	nc_psizr->put(&psizr[0][0],nRow,nCol);

	NcVar *nc_fpol = dataFile.add_var ("fpol", ncFloat, rDim );
	NcVar *nc_fluxGrid = dataFile.add_var ("fluxGrid", ncFloat, rDim );

	nc_fpol->put(&fpol[0],nCol);
	nc_fluxGrid->put(&fluxGrid[0],nCol);

	cout << "DONE." << endl;

	return 0;
}
