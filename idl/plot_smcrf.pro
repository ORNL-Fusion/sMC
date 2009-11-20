pro plot_smcrf, $
	COMPARE = compare

	if strCmp ( getEnv ( 'MACHINE' ), 'jaguar' ) then scratchDir = getEnv ( 'SYSTEM_USERDIR' )
	if strCmp ( getEnv ( 'MACHINE' ), 'franklin' ) then scratchDir = getEnv ( 'SCRATCH' )
	if strCmp ( getEnv ( 'MACHINE' ), 'dlghp' ) then scratchDir = '/home/dg6/scratch' 
	if strCmp ( getEnv ( 'HOSTNAME' ), 'benten.gat.com' ) then scratchDir = '/u/greendl/scratch' 

	;fileName2	= scratchDir + '/p2f/d3d_D_Heidbrink_delta/data/fdis.dav.nc'
	fileName2	= scratchDir + '/p2f/cmod_t0_dlgMod_delta/data/fdis.dav.nc'
	eqdsk_fileName   = 'data/g129x129_1051206002.01120.cmod'
	;eqdsk_fileName	= 'data/eqdsk.122993'

	amu	= 1d0

	e_   = 1.60217646d-19
	q   = 1d0 * e_
	k   = 1.3806504d-23
	mi  = amu * 1.67262158d-27
	c   = 3.0d8

	print, '*** Using ', eqdsk_fileName
	eqdsk   = readGEQDSK ( eqdsk_fileName )

if keyword_set ( compare ) then begin

	cdfId	= ncdf_open ( fileName2, /noWrite )
	glob	= ncdf_inquire ( cdfId )
	ncdf_varGet, cdfId, 'f_rzvv', f_rzvv2
	ncdf_varGet, cdfId, 'vPerp_binEdges', vPerp_binEdges2
	ncdf_varGet, cdfId, 'vPar_binEdges', vPar_binEdges2
	ncdf_varGet, cdfId, 'vPerp_binCenters', vPerp_binCenters2
	ncdf_varGet, cdfId, 'vPar_binCenters', vPar_binCenters2
	nCdf_varGet, cdfId, 'nP', nP2

	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'R_nBins' ), name, R_nBins2
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'z_nBins' ), name, z_nBins2
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'vPerp_nBins' ), name, vPerp_nBins2
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'vPar_nBins' ), name, vPar_nBins2

	ncdf_close, cdfId

	dVPer2	= abs ( vPerp_binEdges2[0] - vPerp_binEdges2[1] )
	dVPar2	= abs ( vPar_binEdges2[0] - vPar_binEdges2[1] )

endif
	
	; Read netcdf data into variables

	fileName	= 'fdis.dav.nc'

	cdfId	= ncdf_open ( 'data/'+fileName, /noWrite )
	glob	= ncdf_inquire ( cdfId )

	; Read netcdf data into variables

	ncdf_varGet, cdfId, 'f_rzvv', f_rzvv
	ncdf_varGet, cdfId, 'R_binCenters', R_binCenters
	ncdf_varGet, cdfId, 'z_binCenters', z_binCenters
	ncdf_varGet, cdfId, 'vPerp_binCenters', vPerp_binCenters
	ncdf_varGet, cdfId, 'vPar_binCenters', vPar_binCenters
	ncdf_varGet, cdfId, 'R_binEdges', R_binEdges
	ncdf_varGet, cdfId, 'z_binEdges', z_binEdges
	ncdf_varGet, cdfId, 'vPerp_binEdges', vPerp_binEdges
	ncdf_varGet, cdfId, 'vPar_binEdges', vPar_binEdges
	ncdf_varGet, cdfId, 'nP', nP
	ncdf_varGet, cdfId, 'R_binSize', R_binSize
	ncdf_varGet, cdfId, 'z_binSize', z_binSize
	ncdf_varGet, cdfId, 'vPerp_binSize', vPerp_binSize
	ncdf_varGet, cdfId, 'vPar_binSize', vPar_binSize
	ncdf_varGet, cdfId, 'density', density_file

	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'R_nBins' ), name, R_nBins
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'z_nBins' ), name, z_nBins
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'vPerp_nBins' ), name, vPerp_nBins
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'vPar_nBins' ), name, vPar_nBins

	ncdf_close, cdfId

	R_nBins	= n_elements ( R_binCenters )
	z_nBins	= n_elements ( z_binCenters )	
	vPerp_nBins	= n_elements ( vPerp_binCenters )
	vPar_nBins	= n_elements ( vPar_binCenters )	

if  keyword_set ( compare ) then begin

	R_nBins2	= n_elements ( R_binCenters2 )
	z_nBins2	= n_elements ( z_binCenters2 )	
	vPerp_nBins2	= n_elements ( vPerp_binCenters2 )
	vPar_nBins2	= n_elements ( vPar_binCenters2 )	

endif

;	Create a density map

	vPerp_3D	= rebin ( vPerp_binCenters, vPerp_nBins, R_nBins, z_nBins )
	vPerp_3D	= transpose ( vPerp_3D, [2,1,0] )

	R_2D	= rebin ( R_binCenters, R_nBins, z_nBins )

	density	= total ( total ( f_rzvv, 4 ) * vPerp_3D , 3 ) * $
		vPerp_binSize * vPar_binSize * 2.0 * !pi

;	Create nP map

	nParticle_map	= density * z_binSize * R_binSize * R_2D * 2.0 * !pi

;	Create wPerp map ( perp energy PER PARTICLE )

	wPerp	= total ( 0.5 * mi * vPerp_3D^2 $
						* total ( f_rzvv, 4 ) * vPerp_3D , 3 ) $
		* vPerp_binSize * vPar_binSize * 2.0 * !pi

	wPerp	= wPerp / e_ * 1d-3 / density
	iiBad	= where ( density eq 0, iiBadCnt )
	if iiBadCnt gt 0 then wPerp[iiBad] = 0


;	Create flux surface averaged quantities

;--------------------------------------------
;	Create a density profile

	density_all	= 0.0
	wperp_all	= 0.0

	R_all	= 0.0
	z_all	= 0.0

	print, 'Building flux surface average f ...'
	for i = 0, R_nBins-1 do begin
		for j = 0, z_nBins - 1 do begin

	    q1_  = where ( ( R_binCenters[i] - eqdsk.rbbbs gt 0 ) and ( z_binCenters[j] - eqdsk.zbbbs gt 0 ), q1 )
	    q2_  = where ( ( R_binCenters[i] - eqdsk.rbbbs gt 0 ) and ( z_binCenters[j] - eqdsk.zbbbs le 0 ), q2 )
	    q3_  = where ( ( R_binCenters[i] - eqdsk.rbbbs le 0 ) and ( z_binCenters[j] - eqdsk.zbbbs gt 0 ), q3 )
	    q4_  = where ( ( R_binCenters[i] - eqdsk.rbbbs le 0 ) and ( z_binCenters[j] - eqdsk.zbbbs le 0 ), q4 )

	    if ( q1 gt 0 ) and ( q2 gt 0 ) and ( q3 gt 0 ) and ( q4 gt 0 ) then begin

			if (size(f_vv_all,/dim))[0] eq 0 then f_vv_all = reform(f_rzvv[i,j,*,*]) $
				else f_vv_all = [ [[ f_vv_all ]], [[ reform ( f_rzvv[i,j,*,*] ) ]] ]
			
			density_all	= [ density_all, density[i,j] ]
			wperp_all	= [ wperp_all, wperp[i,j] ]
			R_all	= [ R_all, R_binCenters[i] ]
			z_all	= [ z_all, z_binCenters[j] ]

		endif

	    endfor 
	
	endfor

	density_all	= density_all[1:*]
	wperp_all	= wperp_all[1:*]
	R_all	= R_all[1:*]
	z_all	= Z_all[1:*]
 
	psi_all = interpolate ( eqdsk.psizr, ( R_all - eqdsk.rleft ) / eqdsk.rdim * eqdsk.nW, $
    			( z_all - min ( eqdsk.z ) ) / eqdsk.zdim * eqdsk.nH )

	psiRange	= abs ( eqdsk.siMag - eqdsk.siBry )
	psiNorm_all	= ( psi_all - eqdsk.siMag ) / psiRange			

	rho_all	= sqrt ( psiNorm_all )	

	R_binCenters2D	= rebin ( R_binCenters, R_nBins, z_nBins )
	z_binCenters2D	= transpose ( rebin ( z_binCenters, z_nBins, R_nBins ) )	
	psi_2D = interpolate ( eqdsk.psizr, ( R_binCenters2D - eqdsk.rleft ) / eqdsk.rdim * eqdsk.nW, $
    			( z_binCenters2D - min ( eqdsk.z ) ) / eqdsk.zdim * eqdsk.nH )

	psi_2D	= ( psi_2D - eqdsk.siMag ) / psiRange

	;	Bin by rho coord.

	rho_binEdges	= fIndGen ( R_nBins/2+1 ) / (R_nBins/2)
	dRho	= abs(rho_binEdges[1]-rho_binEdges[2])
	rho_binCenters	= rho_binEdges[1:*] - dRho/2.0

	density_rho	= fltArr ( n_elements ( rho_binCenters ) )
	wperp_rho	= fltArr ( n_elements ( rho_binCenters ) )
	f_vv_rho	= fltArr ( n_elements ( rho_binCenters ), vPerp_nBins, vPar_nBins )

	for i = 0, n_elements ( rho_binCenters ) - 1 do begin

		iiDx	= where ( rho_all ge rho_binCenters[i] - dRho $
							AND rho_all lt rho_binCenters[i] + dRho $
							AND psiNorm_all le 1, cnt)
		if cnt gt 0 then begin
			density_rho[i]	= total ( density_all[ iiDx ] ) / cnt
			wperp_rho[i]	= total ( wperp_all[ iiDx ] ) / cnt
			if cnt eq 1 then $
				f_vv_rho[i,*,*]	= f_vv_all[ *, *, iiDx ] $
			else $
				f_vv_rho[i,*,*]	= total ( f_vv_all[ *, *, iiDx ], 3 ) / cnt
			print, cnt
		endif
		
	endfor

	;	Create the analytical profile AORSA uses	

	n0	= 1.022e19
	nLim	= 2.733e18
	alpha	= 0.6
	beta_	= 1.4

	;n0	= 1.522e18
	;nLim	= 0.933e17
	;alpha	= 5.0
	;beta_	= 3.0
	
	density_rho_aorsa	= nLim + ( n0 - nLim )*(1d0-rho_binCenters^beta_)^alpha

	device, decomposed = 0
	!p.background = 255
	!p.multi	= [0,1,3]
	!p.charSize = 3
	window, 1, xSize = 400, ySize = 750
	plot, sqrt ( psi_2D ), density/1e19, $
		psym = 4, $
		xtitle = 'rho', $
		ytitle = 'density [m^-3] x1e19', $
		xRange = [0, 1], $
		color = 0, $
		yRange = [0,2], $
		yStyle = 1, $
		title = 'ORBIT-rf'
	loadct, 12
	oPlot, rho_binCenters, density_rho/1e19, $
		psym = -4, $
		color = 8*16-1, $
		thick = 2.0
	oPlot, rho_binCenters, density_rho_aorsa/1e19,$
		thick = 2.0, $
		color = 12*16-1

	plot, sqrt ( psi_2D ), wperp, $
		psym = 4, $
		xtitle = 'rho', $
		ytitle = 'wPerp [keV] per particle', $
		xRange = [0, 1], $
		color = 0, $
		yRange = [0,8], $
		yStyle = 1
	oPlot, rho_binCenters, wperp_rho, $
		psym = -4, $
		color = 8*16-1, $
		thick = 2.0
	
;
;--------------------------------------------


;--------------------------------------------
;	Create a temperature profile by fitting 2D gaussian
;	functions at each pt in space

	tempProfile	= fltArr ( R_nBins, z_nBins )
	psiTemp2D	= fltArr ( R_nBins, z_nBins )
	rhoTemp2D	= fltArr ( R_nBins, z_nBins )

	for i = 0, R_nBins-1 do begin
		print, i
		for j = 0, z_nBins-1 do begin
		
			psiTemp2D[i,j] = interpolate ( eqdsk.psizr, $
				( R_binCenters[i] - eqdsk.rleft ) / eqdsk.rdim * eqdsk.nW, $
    			( z_binCenters[j] - min ( eqdsk.z ) ) / eqdsk.zdim * eqdsk.nH )

			rhoTemp2D[i,j]	= sqrt ( ( psiTemp2D[i,j] - eqdsk.siMag ) / psiRange )	

			tmp	= [reverse(reform(f_rzvv[i,j,*,*]),1),reform(f_rzvv[i,j,*,*])]
			if mean ( tmp ) gt 0 then begin
				fit	= gauss2dFit ( tmp, A )
				tempProfile[i,j]	= mi * ( $
					( A[2] * vPerp_binsize )^2 + ( A[3] * vPar_binSize )^2 $
						) / ( 2.0 * 1e3 * e_ )  
			endif

		endfor
	endfor
	
	temp_rho	= fltArr ( n_elements ( rho_binCenters ) )

	for i = 0, n_elements ( rho_binCenters ) - 1 do begin

		iiDx	= where ( rhoTemp2D[*] ge rho_binCenters[i] - dRho $
			AND rhoTemp2D[*] lt rho_binCenters[i] + dRho, cnt)
		if cnt gt 0 then begin
			temp_rho[i]	= total ( (tempProfile[*])[ iiDx ] )
			temp_rho[i]	= temp_rho[i] / cnt
			print, cnt
		endif
		
	endfor

	temp_rho	= temp_rho * 1e3

	;	Create the analytical profile AORSA uses	

	t0	= 2.0e3 
	tLim	= 0.052e3
	alpha	= 1.3 
	beta_	= 1.9

	;t0	= 24.0e3 
	;tLim	= 2.0e3
	;alpha	= 1.0 
	;beta_	= 2.8


	temp_rho_aorsa	= tLim + ( t0 - tLim )*(1d0-rho_binCenters^beta_)^alpha

	plot, rhoTemp2D, tempProfile, $
		yTitle = 'temperature [keV]', $
		xTitle = 'rho', $
		thick = 1, $
		color = 0, $
		psym = 4, $
		xRange = [0,1.0], $
		xStyle = 1, $
		yStyle = 1, $
		yRange = [0, 5]

	oPlot,rho_binCenters, temp_rho_aorsa/1e3, $
		color = 12*16-1, $
		thick = 2

	oPlot, rho_binCenters, temp_rho/1e3, $
		psym = -4, $
	 	color = 8*16-1, $
		thick = 2

;
;--------------------------------------------


	dVPer	= abs ( vPerp_binEdges[0] - vPerp_binEdges[1] )
	dVPar	= abs ( vPar_binEdges[0] - vPar_binEdges[1] )
	dR	= abs ( R_binEdges[0] - R_binEdges[1] )
	dz	= abs ( z_binEdges[0] - z_binEdges[1] )

	
	; Plot

	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'data/'+fileName+'.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=25, ysize=37.5,xoffset=.1, yoffset=.1, /encapsul

	nLevs	= 10 
	levScale	= 1d-7
	levels	= 10.0^fIndGen(nLevs)*levScale
	colors	= reverse ( bytScl ( fIndGen(nLevs), top = 253 ) + 1 )
	loadct, 0

	!p.multi = [0,R_nBins,z_nBins]
	
	for j=0,z_nBins-1 do begin
		for i=0,R_nBins-1 do begin
	
			f_vv_smooth	= f_rzvv[i,j,*,*]
		
			contour, transpose ( f_vv_smooth )>min(levels), $
				vPar_binCenters / 3.0e6, vPerp_binCenters / 3.0e6, $
				levels = levels, $
				c_colors = colors, $
				color = 0, $
				charSize = 0.01, $
				xRange = [min(vPar_binCenters),max(vPar_binCenters)] / 3.0e6,$
				yRange = [0.0,max(vPerp_binCenters)] / 3.0e6, $
				title = 'R: '+string ( r_binCenters[i], for='(f5.2)' ) $
					+ ' z: '+string ( z_binCenters[j], for='(f5.2)' ), $
				xTicks = 1, $
				yTicks = 1, $
			 	xStyle = 9, $
				yStyle = 9, $
				thick = 0.5

			print, i, j
		endfor	
	endfor	

	!p.position = 0
	!p.charSize = 2	
	loadct, 0
	xyOuts, 0.98, 0.975, 'fileName: '+fileName, color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.95, 'nP: '+string(nP,format='(e7.1)'), color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.925, 'R_nBins: '+string(R_nBins,format='(i3.3)'), color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.9, 'z_nBins: '+string(z_nBins,format='(i3.3)'), color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.875, 'vPerp_nBins: '+string(vPerp_nBins,format='(i3.3)'), color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.85, 'vPar_nBins: '+string(vPar_nBins,format='(i3.3)'), color = 0, /norm, align = 1.0

	device, /close

if keyword_set ( compare ) then begin

	set_plot, 'X'
	!p.multi = [0,3,3]
	!p.charsize = 3.0
	!p.background = 255
	window, 0, xSize = 1200, ySize = 600
	device, decomposed = 0
	loadct, 12

	pt	= where ( density eq max ( density ) )
	i	= (array_indices(density,pt))[0]
	j	= (array_indices(density,pt))[1]

	densityHere	= density[i,j]
	tempHere	= tempProfile[i,j]	

	f_vv	= reform(f_rzvv[i,j,*,*])
	f_vv2	= reform(f_rzvv2[i,j,*,*])

	vPer2D	= rebin(vPerp_binCenters,vPerp_nBins,vPar_nBins)
	vPer2D2	= rebin(vPerp_binCenters2,vPerp_nBins2,vPar_nBins2)

	;	Create analytical maxwellian


	T_keV    = tempHere ; [keV]  
	T	= T_keV * 1d3 * e_ / k

	f_vv_a	= fltArr ( size ( f_vv, /dim ) )
	
	for i = 0, n_elements ( vPerp_binCenters ) - 1 do begin
		for j = 0, n_elements ( vPar_binCenters ) - 1 do begin
			v	= sqrt ( (vPerp_binCenters[i]^2) $
					+ (vPar_binCenters[j]^2) )
			meanV	= 0.0
			f_vv_a[i,j]	= $
				 ( mi / ( 2.0 * !pi * k * T ))^1.5 * $
					exp ( -  mi * (v-meanV)^2 / (2.0*k*T) )
		endfor
	endfor
	
	f_vv_a	= densityHere * f_vv_a $
		/ total ( f_vv_a * rebin(vPerp_binCenters,vPerp_nBins,vPar_nBins) * dVpar * dVPer * 2.0 * !Pi )

	dfdupar	= dlg_pderiv ( f_vv, 2, dVpar )
	dfdupar2	= dlg_pderiv ( f_vv2, 2, dVpar2 )
	dfdupar_a	= dlg_pderiv ( f_vv_a, 2, dVpar )

	dfduper	= dlg_pderiv ( f_vv, 1, dVper )
	dfduper2	= dlg_pderiv ( f_vv2, 1, dVper2 )
	dfduper_a	= dlg_pderiv ( f_vv_a, 1, dVper )

	;	Interpoate to match grids
	
	nGridPer	= (vPerp_binCenters2 - min ( vPerp_binCenters )) $
		/ ( max ( vPerp_binCenters ) - min ( vPerp_binCenters ) ) * (vPerp_nBins-1)

	nGridPar	= (vPar_binCenters2 - min ( vPar_binCenters )) $
		/ ( max ( vPar_binCenters ) - min ( vPar_binCenters ) ) * (vPar_nBins-1)

	f_vv_a	= interpolate ( f_vv_a, nGridPer, nGridPar, /grid )
	f_vv	= interpolate ( f_vv, nGridPer, nGridPar, /grid )
	dfduper_a	= interpolate ( dfduper_a, nGridPer, nGridPar, /grid )
	dfdupar_a	= interpolate ( dfdupar_a, nGridPer, nGridPar, /grid )
	dfduper	= interpolate ( dfduper, nGridPer, nGridPar, /grid )
	dfdupar	= interpolate ( dfdupar, nGridPer, nGridPar, /grid )

	plot, vPar_binCenters2, f_vv2[1,*], color = 0 
	oplot, vPar_binCenters2, f_vv[1,*], color = 12*16-1, thick = 2.0
	oplot, vPar_binCenters2, f_vv_a[1,*], color = 2*16-1, thick = 2.0

	plot, vPar_binCenters2, f_vv2[9,*],color = 0 
	oplot, vPar_binCenters2, f_vv[9,*], color = 12*16-1, thick = 2.0
	oplot, vPar_binCenters2, f_vv_a[9,*], color = 2*16-1, thick = 2.0

	plot, vPerp_binCenters2, f_vv2[*,vPar_nBins2/2], color = 0, psym = 4
	oplot, vPerp_binCenters2, f_vv[*,vPar_nBins2/2], color = 12*16-1, thick = 1.0
	oplot, vPerp_binCenters2, f_vv_a[*,vPar_nBins2/2], color = 2*16-1, thick = 2.0

	plot, vPar_binCenters2, dfdupar2[1,*], color = 0, title = 'dfduPar'
	oplot, vPar_binCenters2, dfdupar[1,*], color = 12*16-1, thick = 2.0
	oplot, vPar_binCenters2, dfdupar_a[1,*], color = 2*16-1, thick = 2.0

	plot, vPar_binCenters2, dfdupar2[9,*], color = 0
	oplot, vPar_binCenters2, dfdupar[9,*], color = 12*16-1, thick = 2.0
	oplot, vPar_binCenters2, dfdupar_a[9,*], color = 2*16-1, thick = 2.0

	plot, vPerp_binCenters2, dfdupar_a[*,vPar_nBins2/2.0], color = 2*16-1, thick = 2.0
	oplot, vPerp_binCenters2, dfdupar2[*,vPar_nBins2/2.0], color = 0
	oplot, vPerp_binCenters2, dfdupar[*,vPar_nBins2/2.0], color = 12*16-1, thick = 2.0

	plot, vPar_binCenters2, dfduper2[1,*], color = 0, title='dfduPer'
	oplot, vPar_binCenters2, dfduper[1,*], color = 12*16-1, thick = 2.0
	oplot, vPar_binCenters2, dfduper_a[1,*], color = 2*16-1, thick = 2.0

	plot, vPar_binCenters2, dfduper2[9,*], color = 0
	oplot, vPar_binCenters2, dfduper[9,*], color = 12*16-1, thick = 2.0
	oplot, vPar_binCenters2, dfduper_a[9,*], color = 2*16-1, thick = 2.0

	plot, vPerp_binCenters2, dfduper_a[*,vPar_nBins2/2], color = 2*16-1, thick = 2.0
	oplot, vPerp_binCenters2, dfduper2[*,vPar_nBins2/2], color = 0
	oplot, vPerp_binCenters2, dfduper[*,vPar_nBins2/2], color = 12*16-1, thick = 2.0

endif
set_plot, 'X'

;;	Contour f_rho_vv
;
;for i=0,n_elements(rho_binCenters)-1 do begin
;
;	contour, transpose ( f_vv_rho[i,*,*] ), $
;		levels = 10.0^fIndgen(12)*1d-2, $
;		path_xy = path, path_info = info
;	for j=0,n_elements(info.offset)-1 do begin
;		iSurface, path[0,(info.offset)[j]:(info.offset)[j]+(info.n)[j]-1]*0+j, $
;					path[0,(info.offset)[j]:(info.offset)[j]+(info.n)[j]-1], $
;					path[1,(info.offset)[j]:(info.offset)[j]+(info.n)[j]-1], $
;					/overPlot, $
;					transparency = j*10
;	endfor
;
;endfor

stop
end

