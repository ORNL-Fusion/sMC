#include <iostream>
#include <fstream>
#include <vector>
#include "eqdsk.hpp"
#include "deriv.hpp"
#include "interpolation.h"
#include <cmath>
#include <netcdfcpp.h>
//#include "boost/multi_array.hpp"
#include "constants.hpp"
#include "array2D.hpp"

using namespace std;
using namespace constants;

//#define nRow,nCol boost::extents[nRow][nCol]

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

		fpol_.resize(nCol_);
		pres.resize(nCol_);
		ffprim.resize(nCol_);
		pprime.resize(nCol_);

		for (int j=0;j<nCol_;j++)
			inFile >> fpol_[j];

		for (int j=0;j<nCol_;j++)
			inFile >> pres[j];

		for (int j=0;j<nCol_;j++)
			inFile >> ffprim[j];

		for (int j=0;j<nCol_;j++)
			inFile >> pprime[j];

		psizr_.resize(nRow_,nCol_);

		for (int i=0;i<nRow_;i++)
		{
			for (int j=0;j<nCol_;j++)
			{
				inFile >> psizr_(i,j);
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

        ascending_flux = BCHECK;
        if ( sibry > simag )
            ascending_flux = false;

        REAL fStep = ( sibry - simag ) / ( nCol_ - 1 );

        fluxGrid_.resize(nCol_);
        for (int j=0;j<nCol_;j++)
            fluxGrid_[j] = j * fStep + simag;

    	dr_   = rdim / ( nCol_ - 1 );
    	dz_   = zdim / ( nRow_ - 1 );

    	r_.resize(nCol_);
    	z_.resize(nRow_);

    	for (int j=0;j<nCol_;j++)
    	    r_[j] = j * dr_ + rleft;

    	for (int i=0;i<nRow_;i++)
    	    z_[i] = i * dz_ + zmid - zdim / 2.0;


	} // end if inFile.good()
	else {

		cout << "ERROR: file '" << fName << "' does not exist?" << endl;
		return 1;

	}

	cout << "DONE." << endl;

	return 0;

}

// Get the eqdsk b field on a desired grid of nrow x ncol
int Ceqdsk::calc_b ( const unsigned int nrow, const unsigned int ncol ) {
		
	cout << "Calculating b field from eqdsk data ..." << endl;

	// Set the public variables
	
	nRow = nrow;
	nCol = ncol;

    // Calculate other quantities from read data

    dr   = rdim / ( nCol - 1 );
    dz   = zdim / ( nRow - 1 );

    r.resize(nCol);
    z.resize(nRow);

    for (int j=0;j<nCol;j++)
        r[j] = j * dr + rleft;

    for (int i=0;i<nRow;i++)
        z[i] = i * dz + zmid - zdim / 2.0;

    br.resize(nRow,nCol);
    bz.resize(nRow,nCol);
    bp.resize(nRow,nCol);
    bmag.resize(nRow,nCol);
	psizr.resize(nRow,nCol);
	fpolzr.resize(nRow,nCol);

	// Interpolate psizr_ from eqdsk grid onto desired r,z grid
	// using ALGLIB bicubic interpolation
	
	// Initialize the AGLIB data arrays
	alglib::real_2d_array AG_psizr;
	alglib::real_1d_array AG_r, AG_z;
	alglib::spline2dinterpolant AG_s2D;
	alglib::ae_int_t AG_M, AG_N;

	// AGLIB is double only, so copy REAL vectors to double
	std::vector<double> r_dbl(r_.begin(),r_.end());
	std::vector<double> z_dbl(z_.begin(),z_.end());
	array2D<double,BCHECK> psizr_dbl(psizr_);
	//psizr_dbl = psizr_;

	// Population the ALGLIB arrays
	AG_r.setcontent(nCol_,&r_dbl[0]);
	AG_z.setcontent(nRow_,&z_dbl[0]);
	AG_psizr.setcontent(nRow_,nCol_,&psizr_dbl(0,0));
	
	// Build the spline
	AG_M = nRow_;
	AG_N = nCol_;
	alglib::spline2dbuildbicubic ( AG_r, AG_z, AG_psizr, AG_M, AG_N, AG_s2D );

	// Interpolate psizr on the new 2D r,z mesh

	for (int j=0;j<nCol;j++) {
		for(int i=0;i<nRow;i++) {

			psizr(i,j) = alglib::spline2dcalc(AG_s2D,r[j],z[i]);

		}
	}

	// br = -dpsi/dz * 1/r
	for (int j=0;j<nCol;j++)
	{
		vector<REAL> tmpData (nRow);
		vector<REAL> tmpRes (nCol);
		for (int i=0;i<nRow;i++) 
			tmpData[i] = psizr(i,j);

    	tmpRes = deriv ( tmpData, dz, 4 );

		for (int i=0;i<nRow;i++)
			br(i,j) = -tmpRes[i] / r[j];

	}

	cout << "\t dr: " << dr << endl;
	cout << "\t dz: " << dz << endl;

	// bz = dpsi/dr * 1/r
	for (int i=0;i<nRow;i++)
	{
		vector<REAL> tmpData (nCol);
		vector<REAL> tmpRes (nRow);
		for (int j=0;j<nCol;j++) 
			tmpData[j] = psizr(i,j);

    	tmpRes = deriv ( tmpData, dr, 4 );
		for (int j=0;j<nCol;j++)
			bz(i,j) = tmpRes[j] / r[j];

	}

	// Interpolate fpol_ from fluxGrid to r,z 2D space
	// Using ALGLIB

	// Initialize the AGLIB data arrays
	alglib::real_1d_array AG_fluxGrid;
	alglib::real_1d_array AG_fpol;
	alglib::spline1dinterpolant AG_s;

	// AGLIB is double only, so copy REAL vectors to double
	std::vector<double> fluxGrid_dbl(fluxGrid_.begin(),fluxGrid_.end());
	std::vector<double> fpol_dbl(fpol_.begin(),fpol_.end());

	// Population the ALGLIB arrays
	AG_fluxGrid.setcontent(nCol_,&fluxGrid_dbl[0]);
	AG_fpol.setcontent(nCol_,&fpol_dbl[0]);

	// Build the spline
	alglib::spline1dbuildcubic ( AG_fluxGrid, AG_fpol, AG_s );

	// Calculate fpol on the 2D r,z mesh

	for (int j=0;j<nCol;j++) {
		for(int i=0;i<nRow;i++) {

			fpolzr(i,j) = alglib::spline1dcalc(AG_s,psizr(i,j));
			bp(i,j) = fpolzr(i,j) / r[j];

		}
	}

	// Magnitude of b
	for (int j=0;j<nCol;j++) {
		for (int i=0;i<nRow;i++) {
			bmag(i,j) = 
				sqrt ( pow(br(i,j),2)+pow(bp(i,j),2)+pow(bz(i,j),2) );
		}
	}

	cout << "DONE." << endl;

	return 0;

} // end calc_b

int Ceqdsk::write_ncfile ( const string fName ) {

	cout << "Writing file " << fName << endl;

	NcFile dataFile ( &fName[0], NcFile::Replace );
	if (!dataFile.is_valid()) {
		cout << "ERROR: Could not open nc file for writing." << endl;
		return 1;
	}	

	NcDim *rDim = dataFile.add_dim("nr",nCol);
	NcDim *zDim = dataFile.add_dim("nz",nRow);

	NcDim *rDim_ = dataFile.add_dim("nr_",nCol_);
	//NcDim *zDim_ = dataFile.add_dim("nz_",nRow_);

	NcVar *nc_r = dataFile.add_var ("r", ncFloat, rDim );
	NcVar *nc_z = dataFile.add_var ("z", ncFloat, zDim );

	nc_r->put(&r[0],nCol);
	nc_z->put(&z[0],nRow);

	NcVar *nc_br = dataFile.add_var ("br", ncFloat, zDim, rDim );
	NcVar *nc_bp = dataFile.add_var ("bp", ncFloat, zDim, rDim );
	NcVar *nc_bz = dataFile.add_var ("bz", ncFloat, zDim, rDim );
	NcVar *nc_bmag = dataFile.add_var ("bmag", ncFloat, zDim, rDim );
	NcVar *nc_psizr = dataFile.add_var ("psizr", ncFloat, zDim, rDim );
	NcVar *nc_fpolzr = dataFile.add_var ("fpolzr", ncFloat, zDim, rDim );

	nc_psizr->put(&psizr(0,0),nRow,nCol);
	nc_fpolzr->put(&fpolzr(0,0),nRow,nCol);

	nc_br->put(&br(0,0),nRow,nCol);
	nc_bp->put(&bp(0,0),nRow,nCol);
	nc_bz->put(&bz(0,0),nRow,nCol);
	nc_bmag->put(&bmag(0,0),nRow,nCol);

	NcVar *nc_fpol = dataFile.add_var ("fpol", ncFloat, rDim_ );
	NcVar *nc_fluxGrid = dataFile.add_var ("fluxGrid", ncFloat, rDim_ );

	nc_fpol->put(&fpol_[0],nCol_);
	nc_fluxGrid->put(&fluxGrid_[0],nCol_);

	cout << "DONE." << endl;

	return 0;
}


int Ceqdsk::bForceTerms ( const int _Z, const int amu ) {

	cout << "Calculating the B force terms ..." << endl;

    array2D<REAL,BCHECK> wc(nRow,nCol),
	    br_B(nRow,nCol),
	    bp_B(nRow,nCol),
	    bz_B(nRow,nCol);

	for(int j=0;j<nCol;j++){
		for(int i=0;i<nRow;i++){

			wc(i,j) = _Z * _e * bmag(i,j) / (amu * _mi);
			br_B(i,j) = br(i,j) / bmag(i,j);
			bp_B(i,j) = bp(i,j) / bmag(i,j);
			bz_B(i,j) = bz(i,j) / bmag(i,j);
		}
	}

    array2D<REAL,BCHECK>
	    br_B_dr(nRow,nCol),
	    br_B_dz(nRow,nCol),
	    bp_B_dr(nRow,nCol),
	    bp_B_dz(nRow,nCol),
	    bz_B_dr(nRow,nCol),
	    bz_B_dz(nRow,nCol);

	gradB_r.resize(nRow,nCol);
	gradB_z.resize(nRow,nCol); 

    array2D<REAL,BCHECK>
	    lnB(nRow,nCol),
	    lnB_dr(nRow,nCol),
	    lnB_dz(nRow,nCol);

	for(int j=0;j<nCol;j++) {
		for(int i=0;i<nRow;i++)
			lnB(i,j) = log ( bmag(i,j) );
	}	

	// do the dr derivatives ...

	std::vector<REAL> tmpIn (nCol);
	std::vector<REAL> tmpOut (nCol);
	
	for(int i=0;i<nRow;i++) {

		for(int j=0;j<nCol;j++) tmpIn[j] = br_B(i,j);
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol;j++) br_B_dr(i,j) = tmpOut[j];

		for(int j=0;j<nCol;j++) tmpIn[j] = bp_B(i,j);
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol;j++) bp_B_dr(i,j) = tmpOut[j];

		for(int j=0;j<nCol;j++) tmpIn[j] = bz_B(i,j);
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol;j++) bz_B_dr(i,j) = tmpOut[j];

		for(int j=0;j<nCol;j++) tmpIn[j] = bmag(i,j);
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol;j++) gradB_r(i,j) = tmpOut[j];

		for(int j=0;j<nCol;j++) tmpIn[j] = lnB(i,j);
    	tmpOut = deriv ( tmpIn, dr, 4 );
		for(int j=0;j<nCol;j++) lnB_dr(i,j) = tmpOut[j];
	}

	// do the dz derivatives ...

	tmpIn.resize(nRow);
	tmpOut.resize(nRow);
	
	for(int j=0;j<nCol;j++) {

		for(int i=0;i<nRow;i++) tmpIn[i] = br_B(i,j);
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow;i++) br_B_dz(i,j) = tmpOut[i];
		
		for(int i=0;i<nRow;i++) tmpIn[i] = bp_B(i,j);
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow;i++) bp_B_dz(i,j) = tmpOut[i];

		for(int i=0;i<nRow;i++) tmpIn[i] = bz_B(i,j);
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow;i++) bz_B_dz(i,j) = tmpOut[i];

		for(int i=0;i<nRow;i++) tmpIn[i] = bmag(i,j);
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow;i++) gradB_z(i,j) = tmpOut[i];

		for(int i=0;i<nRow;i++) tmpIn[i] = lnB(i,j);
    	tmpOut = deriv ( tmpIn, dz, 4 );
		for(int i=0;i<nRow;i++) lnB_dz(i,j) = tmpOut[i];
	}

    array2D<REAL,BCHECK>
	    bDotGradB_r(nRow,nCol),
	    bDotGradB_p(nRow,nCol),
	    bDotGradB_z(nRow,nCol);

	bCurvature_r.resize(nRow,nCol);
	bCurvature_p.resize(nRow,nCol);
	bCurvature_z.resize(nRow,nCol);

	for(int j=0;j<nCol;j++) {
		for(int i=0;i<nRow;i++) {

        	bDotGradB_r(i,j) = br_B(i,j)  * br_B_dr(i,j) 
					+ bz_B(i,j)  * br_B_dz(i,j) 
				   	- pow ( bp_B(i,j), 2 ) / r[j];

        	bDotGradB_p(i,j) = bp_B(i,j) * br_B(i,j) / r[j] 
					+ br_B(i,j) * bp_B_dr(i,j) 
					+ bz_B(i,j) * bp_B_dz(i,j);

        	bDotGradB_z(i,j) = br_B(i,j) * bz_B_dr(i,j) 
					+ bz_B(i,j) * bz_B_dz(i,j);

        	bCurvature_r(i,j) = ( 
							bp_B(i,j) * bDotGradB_z(i,j) 
							- bz_B(i,j) * bDotGradB_p(i,j) 
							) / wc(i,j);
        	bCurvature_p(i,j)   = -1.0 * ( 
							br_B(i,j) * bDotGradB_z(i,j) 
							- bz_B(i,j) * bDotGradB_r(i,j) 
							) / wc(i,j);
        	bCurvature_z(i,j) = ( 
							br_B(i,j) * bDotGradB_p(i,j) 
							- bp_B(i,j) * bDotGradB_r(i,j) 
							) / wc(i,j);
		}
	}

	// Also in here is the b.Grad b for calculating vPar

	bDotGradB.resize(nRow,nCol);

	for(int j=0;j<nCol;j++) {
		for(int i=0;i<nRow;i++) {
			bDotGradB(i,j) = br_B(i,j) * gradB_r(i,j) 
					+ bz_B(i,j) * gradB_z(i,j);
		}
	}	

	bGradient_r.resize(nRow,nCol);
	bGradient_p.resize(nRow,nCol);
	bGradient_z.resize(nRow,nCol);

	for(int j=0;j<nCol;j++) {
		for(int i=0;i<nRow;i++) {
	
			bGradient_r(i,j) = bp_B(i,j) * lnB_dz(i,j) / ( 2.0 * wc(i,j) );
			bGradient_p(i,j) = -1.0 * ( br_B(i,j) * lnB_dz(i,j) 
							- bz_B(i,j)  * lnB_dr(i,j) ) / ( 2.0 * wc(i,j) );
			bGradient_z(i,j) = -1.0 * bp_B(i,j) * lnB_dr(i,j) / ( 2.0 * wc(i,j) );
		}
	}

	cout << "DONE." << endl;

	return 0;
}
