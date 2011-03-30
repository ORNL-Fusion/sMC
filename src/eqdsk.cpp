#include <iostream>
#include <fstream>
#include <vector>
#include "eqdsk.hpp"
#include "deriv.hpp"
#include "interpolation.h"
#include <cmath>
#include <netcdfcpp.h>
#include "boost/multi_array.hpp"
#include "constants.hpp"

using namespace std;
using namespace constants;

int Ceqdsk::get_index 
	( const REAL rIn, const REAL zIn, Ceqdsk::interpIndex &index ) {

	index.i = (rIn - r.front()) / ( r.back() - r.front() ) * r.size();
	index.j = (zIn - z.front()) / ( z.back() - z.front() ) * z.size();

	index.i1 = floor(index.i);
	index.i2 = ceil(index.i);
	index.j1 = floor(index.j);
	index.j2 = ceil(index.j);

    // Check if particle is off grid	
    if( index.i1<0 || index.i2>=(nRow_-1) ||
        index.j1<0 || index.j2>=(nCol_-1) ) {
        return 1;
    }

	return 0;
}

// bi-linear interpolation
// see wikipedia ;)
REAL Ceqdsk::bilinear_interp 
    ( const Ceqdsk::interpIndex &index , const eqdsk::arr2D_ &data ) {

	REAL f11 = data[index.i1][index.j1];
	REAL f21 = data[index.i2][index.j1];
	REAL f12 = data[index.i1][index.j2];
	REAL f22 = data[index.i2][index.j2];

	// (x2-x1)(y2-y1) == 1 since i'm using indices

	REAL dataOut = f11 * (index.i2-index.i)*(index.j2-index.j)
			+ f21 * (index.i-index.i1)*(index.j2-index.j)
			+ f12 * (index.i2-index.i)*(index.j-index.j1)
			+ f22 * (index.i-index.i1)*(index.j-index.j1); 

    return dataOut;
}

int Ceqdsk::read_file ( string fName ) {

	cout << "Reading g-eqdsk file " << fName << endl;

	ifstream inFile ( fName.c_str() );

	if ( inFile.good() ) {

        // read data from file

		int headerLength = 48;
		header = new char [headerLength];

		inFile.read ( header, headerLength );

		inFile >> idum >> nCol_ >> nRow_;	
		inFile >> rdim >> zdim >> rcentr >> rleft >> zmid;
        inFile >> rmaxis >> zmaxis >> simag >> sibry >> bcentr;
        inFile >> current >> simag >> xdum >> rmaxis >> xdum; 
        inFile >> zmaxis >> xdum >> sibry >> xdum >> xdum; 

		cout << "\t header: " << header << endl;
		cout << "\t idum: " << idum << endl;
		cout << "\t nCol_: " << nCol_ << endl;
		cout << "\t nRow_: " << nRow_ << endl;

		cout << "\t rdim: " << rdim << endl;
		cout << "\t zdim: " << zdim << endl;
		cout << "\t rcentr: " << rcentr << endl;
		cout << "\t rleft: " << rleft << endl;

		cout << "\t rmaxis: " << rmaxis << endl;
		cout << "\t zmaxis: " << zmaxis << endl;
		cout << "\t simag: " << simag << endl;
		cout << "\t sibry: " << sibry << endl;
		cout << "\t bcentr: " << bcentr << endl;

		fpol.resize(nCol_);
		pres.resize(nCol_);
		ffprim.resize(nCol_);
		pprime.resize(nCol_);

		for (int j=0;j<nCol_;j++)
			inFile >> fpol[j];

		for (int j=0;j<nCol_;j++)
			inFile >> pres[j];

		for (int j=0;j<nCol_;j++)
			inFile >> ffprim[j];

		for (int j=0;j<nCol_;j++)
			inFile >> pprime[j];

		psizr.resize(boost::extents[nRow_][nCol_]);

		for (int i=0;i<nRow_;i++)
		{
			for (int j=0;j<nCol_;j++)
			{
				inFile >> psizr[i][j];
			}
		}	 

        qpsi.resize(nCol_);

		for (int j=0;j<nCol_;j++)
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

        dr   = rdim / ( nCol_ - 1 );
        dz   = zdim / ( nRow_ - 1 );

        ascending_flux = true;
        if ( sibry > simag )
            ascending_flux = false;

        REAL fStep = ( sibry - simag ) / ( nCol_ - 1 );

        r.resize(nCol_);
        z.resize(nRow_);

        for (int j=0;j<nCol_;j++)
            r[j] = j * dr + rleft;

        for (int i=0;i<nRow_;i++)
            z[i] = i * dz + zmid - zdim / 2.0;

        fluxGrid.resize(nCol_);
        for (int j=0;j<nCol_;j++)
            fluxGrid[j] = j * fStep + simag;

        br.resize(boost::extents[nRow_][nCol_]);
        bz.resize(boost::extents[nRow_][nCol_]);
        bp.resize(boost::extents[nRow_][nCol_]);
        bmag.resize(boost::extents[nRow_][nCol_]);

		// br = -dpsi/dz * 1/r
		for (int j=0;j<nCol_;j++)
		{
			vector<REAL> tmpData (nRow_);
			vector<REAL> tmpRes (nCol_);
			for (int i=0;i<nRow_;i++) 
				tmpData[i] = psizr[i][j];

        	tmpRes = deriv ( tmpData, dz, 4 );

			for (int i=0;i<nRow_;i++)
				br[i][j] = -tmpRes[i] / r[j];

		}

		cout << "\t dr: " << dr << endl;
		cout << "\t dz: " << dz << endl;

		// bz = dpsi/dr * 1/r
		for (int i=0;i<nRow_;i++)
		{
			vector<REAL> tmpData (nCol_);
			vector<REAL> tmpRes (nRow_);
			for (int j=0;j<nCol_;j++) 
				tmpData[j] = psizr[i][j];

        	tmpRes = deriv ( tmpData, dr, 4 );
			for (int j=0;j<nCol_;j++)
				bz[i][j] = tmpRes[j] / r[j];

		}


		// Interpolate fpol from fluxGrid to r,z 2D space
		// Using ALGLIB

		// Initialize the AGLIB data arrays
		alglib::real_1d_array AG_fluxGrid;
		alglib::real_1d_array AG_fpol;
		alglib::spline1dinterpolant AG_s;

		// AGLIB is double only, so copy REAL vectors to double
		std::vector<double> fluxGrid_dbl(fluxGrid.begin(),fluxGrid.end());
		std::vector<double> fpol_dbl(fpol.begin(),fpol.end());

		// Population the ALGLIB arrays
		AG_fluxGrid.setcontent(nCol_,&fluxGrid_dbl[0]);
		AG_fpol.setcontent(nCol_,&fpol_dbl[0]);

		// Build the spline
		alglib::spline1dbuildcubic ( AG_fluxGrid, AG_fpol, AG_s );

		// Calculate fpol on the 2D r,z mesh
		fpolzr.resize(boost::extents[nRow_][nCol_]);

		for (int j=0;j<nCol_;j++) {
			for(int i=0;i<nRow_;i++) {

				fpolzr[i][j] = alglib::spline1dcalc(AG_s,psizr[i][j]);
				bp[i][j] = fpolzr[i][j] / r[j];

			}
		}

		// Magnitude of b
		for (int j=0;j<nCol_;j++) {
			for (int i=0;i<nRow_;i++) {
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

	NcDim *rDim = dataFile.add_dim("nr",nCol_);
	NcDim *zDim = dataFile.add_dim("nz",nRow_);

	NcVar *nc_r = dataFile.add_var ("r", ncFloat, rDim );
	NcVar *nc_z = dataFile.add_var ("z", ncFloat, zDim );

	nc_r->put(&r[0],nCol_);
	nc_z->put(&z[0],nRow_);

	NcVar *nc_br = dataFile.add_var ("br", ncFloat, zDim, rDim );
	NcVar *nc_bp = dataFile.add_var ("bp", ncFloat, zDim, rDim );
	NcVar *nc_bz = dataFile.add_var ("bz", ncFloat, zDim, rDim );
	NcVar *nc_bmag = dataFile.add_var ("bmag", ncFloat, zDim, rDim );
	NcVar *nc_psizr = dataFile.add_var ("psizr", ncFloat, rDim, zDim );
	NcVar *nc_fpolzr = dataFile.add_var ("fpolzr", ncFloat, rDim, zDim );

	nc_psizr->put(&psizr[0][0],nRow_,nCol_);
	nc_fpolzr->put(&fpolzr[0][0],nRow_,nCol_);
	nc_br->put(&br[0][0],nRow_,nCol_);
	nc_bp->put(&bp[0][0],nRow_,nCol_);
	nc_bz->put(&bz[0][0],nRow_,nCol_);
	nc_bmag->put(&bmag[0][0],nRow_,nCol_);

	NcVar *nc_fpol = dataFile.add_var ("fpol", ncFloat, rDim );
	NcVar *nc_fluxGrid = dataFile.add_var ("fluxGrid", ncFloat, rDim );

	nc_fpol->put(&fpol[0],nCol_);
	nc_fluxGrid->put(&fluxGrid[0],nCol_);

	cout << "DONE." << endl;

	return 0;
}


int Ceqdsk::bForceTerms () {

	cout << "Calculating the B force terms ..." << endl;

	int Z = 1;
	arr2D_ wc(boost::extents[nRow_][nCol_]);
	arr2D_ br_B(boost::extents[nRow_][nCol_]);
	arr2D_ bp_B(boost::extents[nRow_][nCol_]);
	arr2D_ bz_B(boost::extents[nRow_][nCol_]);

	for(int j=0;j<nCol_;j++){
		for(int i=0;i<nRow_;i++){

			wc[i][j] = Z * _e * bmag[i][j] / _mi;
			br_B[i][j] = br[i][j] / bmag[i][j];
			bp_B[i][j] = bp[i][j] / bmag[i][j];
			bz_B[i][j] = bz[i][j] / bmag[i][j];

		}
	}

	arr2D_ br_B_dr(boost::extents[nRow_][nCol_]);
	arr2D_ br_B_dz(boost::extents[nRow_][nCol_]);
	arr2D_ bp_B_dr(boost::extents[nRow_][nCol_]);
	arr2D_ bp_B_dz(boost::extents[nRow_][nCol_]);
	arr2D_ bz_B_dr(boost::extents[nRow_][nCol_]);
	arr2D_ bz_B_dz(boost::extents[nRow_][nCol_]);

	gradB_r.resize(boost::extents[nRow_][nCol_]);
	gradB_z.resize(boost::extents[nRow_][nCol_]); 

	arr2D_ lnB(boost::extents[nRow_][nCol_]);
	arr2D_ lnB_dr(boost::extents[nRow_][nCol_]);
	arr2D_ lnB_dz(boost::extents[nRow_][nCol_]);

	for(int j=0;j<nCol_;j++) {
		for(int i=0;i<nRow_;i++)
			lnB[i][j] = log ( bmag[i][j] );
	}	

	// do the dr derivatives ...

	std::vector<REAL> tmpIn (nCol_);
	std::vector<REAL> tmpOut (nCol_);
	
	for(int i=0;i<nRow_;i++) {

		for(int j=0;j<nCol_;j++) tmpIn[j] = br_B[i][j];
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol_;j++) br_B_dr[i][j] = tmpOut[j];

		for(int j=0;j<nCol_;j++) tmpIn[j] = bp_B[i][j];
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol_;j++) bp_B_dr[i][j] = tmpOut[j];

		for(int j=0;j<nCol_;j++) tmpIn[j] = bz_B[i][j];
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol_;j++) bz_B_dr[i][j] = tmpOut[j];

		for(int j=0;j<nCol_;j++) tmpIn[j] = bmag[i][j];
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol_;j++) gradB_r[i][j] = tmpOut[j];

		for(int j=0;j<nCol_;j++) tmpIn[j] = lnB[i][j];
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol_;j++) lnB_dr[i][j] = tmpOut[j];
	}

	// do the dz derivatives ...

	tmpIn.resize(nRow_);
	tmpOut.resize(nRow_);
	
	for(int j=0;j<nCol_;j++) {

		for(int i=0;i<nRow_;i++) tmpIn[i] = br_B[i][j];
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow_;i++) br_B_dz[i][j] = tmpOut[i];
		
		for(int i=0;i<nRow_;i++) tmpIn[i] = bp_B[i][j];
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow_;i++) bp_B_dz[i][j] = tmpOut[i];

		for(int i=0;i<nRow_;i++) tmpIn[i] = bz_B[i][j];
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow_;i++) bz_B_dz[i][j] = tmpOut[i];

		for(int i=0;i<nRow_;i++) tmpIn[i] = bmag[i][j];
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow_;i++) gradB_z[i][j] = tmpOut[i];

		for(int i=0;i<nRow_;i++) tmpIn[i] = lnB[i][j];
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow_;i++) lnB_dz[i][j] = tmpOut[i];
	}

	arr2D_ bDotGradB_r(boost::extents[nRow_][nCol_]);
	arr2D_ bDotGradB_p(boost::extents[nRow_][nCol_]);
	arr2D_ bDotGradB_z(boost::extents[nRow_][nCol_]);

	bCurvature_r.resize(boost::extents[nRow_][nCol_]);
	bCurvature_p.resize(boost::extents[nRow_][nCol_]);
	bCurvature_z.resize(boost::extents[nRow_][nCol_]);

	for(int j=0;j<nCol_;j++) {
		for(int i=0;i<nRow_;i++) {

        	bDotGradB_r[i][j] = br_B[i][j]  * br_B_dr[i][j] 
					+ bz_B[i][j]  * br_B_dz[i][j] 
				   	- pow ( bp_B[i][j], 2 ) / r[j];

        	bDotGradB_p[i][j] = bp_B[i][j] * br_B[i][j] / r[j] 
					+ br_B[i][j] * bp_B_dr[i][j] 
					+ bz_B[i][j] * bp_B_dz[i][j];

        	bDotGradB_z[i][j] = br_B[i][j] * bz_B_dr[i][j] 
					+ bz_B[i][j] * bz_B_dz[i][j];

        	bCurvature_r[i][j] = ( 
							bp_B[i][j] * bDotGradB_z[i][j] 
							- bz_B[i][j] * bDotGradB_p[i][j] 
							) / wc[i][j];
        	bCurvature_p[i][j]   = -1.0 * ( 
							br_B[i][j] * bDotGradB_z[i][j] 
							- bz_B[i][j] * bDotGradB_r[i][j] 
							) / wc[i][j];
        	bCurvature_z[i][j] = ( 
							br_B[i][j] * bDotGradB_p[i][j] 
							- bp_B[i][j] * bDotGradB_r[i][j] 
							) / wc[i][j];
		}
	}

	// Also in here is the b.Grad b for calculating vPar

	bDotGradB.resize(boost::extents[nRow_][nCol_]);

	for(int j=0;j<nCol_;j++) {
		for(int i=0;i<nRow_;i++) {
			bDotGradB[i][j] = br_B[i][j] * gradB_r[i][j] 
					+ bz_B[i][j] * gradB_z[i][j];
		}
	}	

	bGradient_r.resize(boost::extents[nRow_][nCol_]);
	bGradient_p.resize(boost::extents[nRow_][nCol_]);
	bGradient_z.resize(boost::extents[nRow_][nCol_]);

	for(int j=0;j<nCol_;j++) {
		for(int i=0;i<nRow_;i++) {
	
			bGradient_r[i][j] = bp_B[i][j] * lnB_dz[i][j] / ( 2.0 * wc[i][j] );
			bGradient_p[i][j] = -1.0 * ( br_B[i][j] * lnB_dz[i][j] 
							- bz_B[i][j]  * lnB_dr[i][j] ) / ( 2.0 * wc[i][j] );
			bGradient_z[i][j] = -1.0 * bp_B[i][j] * lnB_dr[i][j] / ( 2.0 * wc[i][j] );
		}
	}

	cout << "DONE." << endl;

	return 0;
}
