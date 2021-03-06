pro cpp_plot_orbit, color = color

	if not keyword_set(color) then color = 'red'

	fileList = file_search ( 'output/*orbit.nc' )

	eqdsk = readgeqdsk ( 'data/eqdsk' )	

	lim_p	= fIndGen(360) * !dtor
	lim_r	= fltArr(360) + max(eqdsk.rlim)

	lim_x	= lim_r * cos ( lim_p )
	lim_y	= lim_r * sin ( lim_p )

	window, 0
	window, 1, xSize = 400, ySize = 800
	window, 2, xSize = 400, ySize = 400

	ii = 1
	f = 0
	while ii lt 11 do begin

		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'rOrb', rOrb
			ncdf_varget, cdfId, 'pOrb', pOrb
			ncdf_varget, cdfId, 'zOrb', zOrb
			ncdf_varget, cdfId, 'stat', stat
			ncdf_varget, cdfId, 'vPar', vPar 

		nCdf_close,	cdfId 

			tleStr = string(f) + ' , ' + string (n_elements(rorb)) + ' ' + string(vPar)
			c = contour ( eqdsk.psizr, eqdsk.r, eqdsk.z, $
					aspect = 1.0, $
					layout = [5,2,ii], $
		   			/current, dim = [1400,800], $
				   	title = tleStr	)
			p = plot ( rOrb, zOrb, /over, thick = 1, color = color )

			x = rOrb * cos ( pOrb )
			y = rOrb * sin ( pOrb )

			ii++

		f++

	endwhile

	!p.multi = 0
	stop

end	
