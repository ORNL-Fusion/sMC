pro cpp_plot_orbit

	fileList = file_search ( 'output/*orbit.nc' )

	eqdsk = readgeqdsk ( 'data/eqdsk' )	

	lim_p	= fIndGen(360) * !dtor
	lim_r	= fltArr(360) + max(eqdsk.rlim)

	lim_x	= lim_r * cos ( lim_p )
	lim_y	= lim_r * sin ( lim_p )

	window, 0
	window, 1, xSize = 400, ySize = 400
	window, 2, xSize = 400, ySize = 400

	for f=0,n_elements(fileList)-1 do begin

		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'rOrb', rOrb
			ncdf_varget, cdfId, 'pOrb', pOrb
			ncdf_varget, cdfId, 'zOrb', zOrb
			ncdf_varget, cdfId, 'stat', stat

		nCdf_close,	cdfId 

		if stat eq 0 then begin

			wset, 0
			fsc_contour, eqdsk.psizr, eqdsk.r, eqdsk.z, /iso
			fsc_plot, rOrb, zOrb, /over, thick = 2
			;p = plot( rOrb, zOrb, /over, thick = 2 )

			x = rOrb * cos ( pOrb )
			y = rOrb * sin ( pOrb )
			wset, 1
			fsc_plot, lim_x, lim_y, aspect = 1.0
			fsc_plot, x, y, /over, thick = 2

			p = plot3d ( x, y, zOrb, aspect = 1.0 )

		endif

		r=get_kbrd()
	endfor
	stop

end	
