pro cpp_plot_orbit

	eqdsk = readgeqdsk ( 'data/eqdsk' )	

	cdfId = ncdf_open('output/orbit.nc')

		ncdf_varget, cdfId, 'rOrb', rOrb
		ncdf_varget, cdfId, 'pOrb', pOrb
		ncdf_varget, cdfId, 'zOrb', zOrb

	nCdf_close,	cdfId 

	c = contour( eqdsk.psizr, eqdsk.r, eqdsk.z, /iso )
	p = plot( rOrb, zOrb, /over, thick = 2 )

	x = rOrb * cos ( pOrb )
	y = rOrb * sin ( pOrb )

	p = plot( x, y, /iso)

	stop

end	
