pro plot_power, $
		particleRatio = particleRatio

	if strCmp ( getEnv ( 'MACHINE' ), 'jaguar' ) then scratchDir = getEnv ( 'SYSTEM_USERDIR' )
	if strCmp ( getEnv ( 'MACHINE' ), 'franklin' ) then scratchDir = getEnv ( 'SCRATCH' )
	if strCmp ( getEnv ( 'MACHINE' ), 'dlghp' ) then scratchDir = '/home/dg6/scratch' 
	if strCmp ( getEnv ( 'HOSTNAME' ), 'benten.gat.com' ) then scratchDir = '/u/greendl/scratch' 

	eqdsk_fileName   = 'data/g129x129_1051206002.01120.cmod'
	;eqdsk_fileName	= 'data/eqdsk.122993'
	;eqdsk_fileName	= 'data/g122080.03100'

	amu	= 1d0

	e_   = 1.60217646d-19
	q   = 1d0 * e_
	k   = 1.3806504d-23
	mi  = amu * 1.67262158d-27
	c   = 3.0d8

	print, '*** Using ', eqdsk_fileName
	eqdsk   = readGEQDSK ( eqdsk_fileName )

	; Read netcdf data into variables

	fileName	= 'data/pList.dav.nc.001'
	cdfId	= ncdf_open ( fileName, /noWrite )
	glob	= ncdf_inquire ( cdfId )
	ncdf_varGet, cdfId, 'power', power 
	ncdf_varGet, cdfId, 'powerPar', powerPar
	ncdf_varGet, cdfId, 'power_R_binCenters', power_R
	ncdf_varGet, cdfId, 'power_z_binCenters', power_z
	ncdf_varGet, cdfId, 'time', time 
	ncdf_close, cdfId

	fileName	= 'data/mchoi_dlg.nc'
	cdfId	= ncdf_open ( fileName, /noWrite )
	glob	= ncdf_inquire ( cdfId )
	ncdf_varGet, cdfId, 'ePlus', ePlusReal
	ncdf_varGet, cdfId, 'ePlus_img', ePlusImg
	ncdf_varGet, cdfId, 'eMinu', eMinuReal
	ncdf_varGet, cdfId, 'eMinu_img', eMinuImg
	ncdf_varGet, cdfId, 'kPer_cold', kPerReal 
	ncdf_varGet, cdfId, 'kPer_img_cold', kPerImg 
	ncdf_varGet, cdfId, 'R', eField_R 
	ncdf_varGet, cdfId, 'z', eField_z 
	ncdf_close, cdfId


	R_nBins	= n_elements ( power_R )
	z_nBins	= n_elements ( power_z )

	;	Convert into MW/m^3

	R_binCenters2D	= rebin ( power_R, R_nBins, z_nBins )
	z_binCenters2D	= transpose ( rebin ( power_z, z_nBins, R_nBins ) )	

	R_binSize	= abs ( power_R[0]-power_R[1] )
	z_binSize	= abs ( power_z[0]-power_z[1] )

	;	Power from file is in ev * e_ remember

	power	= power / time / ( R_binSize * z_binSize * R_binCenters2D * 2.0 * !pi )
	powerPar	= powerPar / time / ( R_binSize * z_binSize * R_binCenters2D * 2.0 * !pi )

	print, 'Time: ', time
;	Create flux surface averaged quantities

;--------------------------------------------
;	Create a density profile

	power_all	= 0d0
	powerPar_all	= 0d0
	
	R_all	= 0.0
	z_all	= 0.0

	for i = 0, R_nBins-1 do begin
		for j = 0, z_nBins - 1 do begin

	    q1_  = where ( ( power_R[i] - eqdsk.rbbbs gt 0 ) and ( power_z[j] - eqdsk.zbbbs gt 0 ), q1 )
	    q2_  = where ( ( power_R[i] - eqdsk.rbbbs gt 0 ) and ( power_z[j] - eqdsk.zbbbs le 0 ), q2 )
	    q3_  = where ( ( power_R[i] - eqdsk.rbbbs le 0 ) and ( power_z[j] - eqdsk.zbbbs gt 0 ), q3 )
	    q4_  = where ( ( power_R[i] - eqdsk.rbbbs le 0 ) and ( power_z[j] - eqdsk.zbbbs le 0 ), q4 )

	    if ( q1 gt 0 ) and ( q2 gt 0 ) and ( q3 gt 0 ) and ( q4 gt 0 ) then begin

			power_all	= [ power_all, double ( power[i,j] ) ]
			powerPar_all	= [ powerPar_all, double ( powerPar[i,j] ) ]

			R_all	= [ R_all, power_R[i] ]
			z_all	= [ z_all, power_z[j] ]

		endif

	    endfor 
	
	endfor
	
	power_all	= power_all[1:*]
	powerPar_all	= powerPar_all[1:*]
	R_all	= R_all[1:*]
	z_all	= Z_all[1:*]
 
	psi_all = interpolate ( eqdsk.psizr, ( R_all - eqdsk.rleft ) / eqdsk.rdim * eqdsk.nW, $
    			( z_all - min ( eqdsk.z ) ) / eqdsk.zdim * eqdsk.nH )

	psiRange	= abs ( eqdsk.siMag - eqdsk.siBry )
	psiNorm_all	= ( psi_all - eqdsk.siMag ) / psiRange			

	rho_all	= sqrt ( psiNorm_all )	


	psi_2D = interpolate ( eqdsk.psizr, ( R_binCenters2D - eqdsk.rleft ) / eqdsk.rdim * eqdsk.nW, $
    			( z_binCenters2D - min ( eqdsk.z ) ) / eqdsk.zdim * eqdsk.nH )

	psi_2D	= ( psi_2D - eqdsk.siMag ) / psiRange

	;	Bin by rho coord.

	rho_nBins	= 40
	rho_binEdges	= fIndGen ( rho_nBins+1 ) / rho_nBins 
	dRho	= abs(rho_binEdges[1]-rho_binEdges[2])
	rho_binCenters	= rho_binEdges[1:*] - dRho/2.0

	power_rho	= dblArr ( n_elements ( rho_binCenters ) )
	power_rho_median	= dblArr ( n_elements ( rho_binCenters ) )
	powerPar_rho	= dblArr ( n_elements ( rho_binCenters ) )
	powerPar_rho_median	= dblArr ( n_elements ( rho_binCenters ) )


	for i = 0, n_elements ( rho_binCenters ) - 1 do begin

		iiDx	= where ( rho_all ge rho_binCenters[i] - dRho $
							AND rho_all lt rho_binCenters[i] + dRho $
							AND psiNorm_all le 1, cnt)
		print, cnt
		if cnt gt 0 then begin
			power_rho[i]	= total ( power_all[ iiDx ] ) / cnt
			power_rho_median[i]	= median ( power_all[iiDx] )
			powerPar_rho[i]	= total ( powerPar_all[ iiDx ] ) / cnt
			powerPar_rho_median[i]	= median ( powerPar_all[iiDx] )
		endif else begin
			print, rho_binCenters[i], ' has no data'
		endelse

		
	endfor


	if not keyword_set ( particleRatio ) then particleRatio = 1.0
	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'data/power.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=12, ysize=18,xoffset=.1, yoffset=.1, /encapsul
	
	!p.multi = [0,2,3]
	!p.charSize = 1.0 
	loadct, 3
	levels	= (fIndGen(10)+1)/10 * 4e5
	colors	= reverse(bytScl ( levels, top =253 ) + 1)
	
	contour, power*particleRatio, power_r, power_z, $
			levels = levels, $
   			c_colors = colors, $
		   	color = 0, $
			/fill, $
			xTitle = 'R[m]', $
			yTitle = 'z[m]', $
			xrange = [ min(eqdsk.rbbbs),max(eqdsk.rbbbs) ], $
			yRange = [ min(eqdsk.zbbbs),max(eqdsk.zbbbs) ], $
			xStyle = 1, $
			yStyle = 1

	;contour, power*particleRatio, power_r, power_z, $
	;		levels = levels, $
   	;		c_colors = colors, $
	;	   	color = 0, $
	;		/over	
	plots, eqdsk.rbbbs, eqdsk.zbbbs, $
		   	color = 0

	plot, rho_all, power_all * particleRatio*1d-6, $
			psym = 4, $
			color = 0, $
			yRange = [0, 10.0], $
			xRange = [0,1.0], $
			yStyle = 1, $
			xTitle = 'rho', $
			yTitle = 'power [MW/m^3]', $
			symSize = 0.5, /noData
	loadct, 12
	oplot, rho_binCenters, power_rho*particleRatio*1d-6, $
		  	color = 12*16-1, $
		   	thick = 4
	;oplot, rho_binCenters, power_rho_median*particleRatio*1d-6, $
	;	  	color = 8*16-1, $
	;		thick = 4

	contour, powerPar*particleRatio, power_r, power_z, $
			levels = levels, $
   			c_colors = colors, $
		   	color = 0, $
			/fill

	;contour, powerPar*particleRatio, power_r, power_z, $
	;		levels = levels, $
   	;		c_colors = colors, $
	;	   	color = 0, $
	;		/over	
	plots, eqdsk.rbbbs, eqdsk.zbbbs, $
		   	color = 0

	plot, rho_all, powerPar_all * particleRatio, $
			psym = 4, $
			color = 0, $
			yRange = [0, 2e7], ystyle = 1, /noData


	loadct, 12
	oplot, rho_binCenters, powerPar_rho*particleRatio, $
		  	color = 12*16-1, $
		   	thick = 2
	oplot, rho_binCenters, powerPar_rho_median*particleRatio, $
		  	color = 8*16-1, $
			thick = 2


	loadct, 13, file = 'data/davect.tbl'
	levels	= ( fIndGen ( 31 ) - 15 ) / 15.0 * 6e3
	colors	= bytScl ( levels , top = 253 ) + 1
	
	contour, (ePlusReal<max(levels))>min(levels), eField_R, eField_z, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			/fill, $
			title = 'Re E+', $
			xTitle = 'R[m]', $
			yTitle = 'z[m]', $
			xrange = [ min(eqdsk.rbbbs),max(eqdsk.rbbbs) ], $
			yRange = [ min(eqdsk.zbbbs),max(eqdsk.zbbbs) ], $
			xStyle = 1, $
			yStyle = 1



	loadct, 0
	plots, eqdsk.rbbbs, eqdsk.zbbbs, $
		   	color = 0

	loadct, 13, file = 'data/davect.tbl'
	
	contour, (eMinuReal<max(levels))>min(levels), eField_R, eField_z, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			/fill, $
			title = 'Re E-'

	loadct, 0
	plots, eqdsk.rbbbs, eqdsk.zbbbs, $
		   	color = 0

	device, /close

	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'data/ePlus.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=8, ysize=8,xoffset=.1, yoffset=.1, /encapsul
	loadct, 13, file = 'data/davect.tbl'
	!p.multi = 0
	contour, (ePlusReal<max(levels))>min(levels), eField_R, eField_z, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			/fill, $
			title = 'Re E+', charSize = 1.0, $
			xRange = [min(eqdsk.rbbbs),max(eqdsk.rbbbs)], $
			yRange = [min(eqdsk.zbbbs),max(eqdsk.zbbbs)], $
			xTitle = 'R [m]', yTitle = 'z [m]', $
			xStyle = 1, yStyle = 1

	loadct, 0
	plots, eqdsk.rbbbs, eqdsk.zbbbs, $
		   	color = 0
	device, /close

	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'data/eMinu.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=8, ysize=8,xoffset=.1, yoffset=.1, /encapsul
	loadct, 13, file = 'data/davect.tbl'
	!p.multi = 0
	contour, (eMinuReal<max(levels))>min(levels), eField_R, eField_z, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			/fill, $
			title = 'Re E+', charSize = 1.0, $
			xRange = [min(eqdsk.rbbbs),max(eqdsk.rbbbs)], $
			yRange = [min(eqdsk.zbbbs),max(eqdsk.zbbbs)], $
			xTitle = 'R [m]', yTitle = 'z [m]', $
			xStyle = 1, yStyle = 1

	loadct, 0
	plots, eqdsk.rbbbs, eqdsk.zbbbs, $
		   	color = 0
	device, /close


stop
end
