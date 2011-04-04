pro Cpp_plot_bfield

	eqdsk = readgeqdsk ( 'data/eqdsk' )	

	cdfId = ncdf_open('output/bField.nc')

		ncdf_varget, cdfId, 'r', r
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 'br',br
		ncdf_varget, cdfId, 'bp',bp
		ncdf_varget, cdfId, 'bz',bz
		ncdf_varget, cdfId, 'bmag',bmag
		ncdf_varget, cdfId, 'psizr', psizr
		ncdf_varget, cdfId, 'fpolzr',fpolzr

	nCdf_close,	cdfId 

	!p.charSize = 2
	c = contour( psizr, r, z, title = 'psizr', aspect = 1.0, layout = [3,2,1] )
	c = contour( bmag, r, z, title = 'bmag', aspect = 1.0, layout = [3,2,2], /current )
	c = contour( br, r, z, title = 'br', aspect = 1.0, layout = [3,2,3], /current )
	c = contour( bp, r, z, title = 'bp', aspect = 1.0, layout = [3,2,4], /current )
	c = contour( bz, r, z, title = 'bz', aspect = 1.0, layout = [3,2,5], /current )

	c = contour( eqdsk.psizr, 	eqdsk.r, eqdsk.z, title = 'psizr', aspect = 1.0, layout = [3,2,1] )
	c = contour( eqdsk.bmag, 	eqdsk.r, eqdsk.z, title = 'bmag', aspect = 1.0, layout = [3,2,2], /current )
	c = contour( eqdsk.br, 		eqdsk.r, eqdsk.z, title = 'br', aspect = 1.0, layout = [3,2,3], /current )
	c = contour( eqdsk.bphi,	eqdsk.r, eqdsk.z, title = 'bp', aspect = 1.0, layout = [3,2,4], /current )
	c = contour( eqdsk.bz, 		eqdsk.r, eqdsk.z, title = 'bz', aspect = 1.0, layout = [3,2,5], /current )
	
	!p.charSize = 1

stop
end
