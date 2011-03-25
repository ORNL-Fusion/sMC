pro Cpp_plot_bfield

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


	window, 0, xSize = 900, ySize = 900
	!p.multi = [0,3,2]
	!p.charSize = 2
		fsc_contour, psizr, r, z, title = 'psizr', /iso
		fsc_contour, fpolzr, r, z, title = 'fpolzr', /iso
		fsc_contour, bmag, r, z, title = 'bmag', /iso
		fsc_contour, br, r, z, title = 'br', /iso
		fsc_contour, bp, r, z, title = 'bp', /iso
		fsc_contour, bz, r, z, title = 'bz', /iso
	!p.charSize = 1
	!p.multi =0

end
