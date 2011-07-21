pro cpp_plot_f

	fileName = 'output/fdis_D_40keV_D3D_out.nc'
	;fileName = 'data/fdis_D_40keV_D3D.nc'

	cdfId = ncdf_open(fileName)

		ncdf_varget, cdfId, 'r', r
		ncdf_varget, cdfId, 'p', p
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 'status', status
		ncdf_varget, cdfId, 'vPer', vPer 
		ncdf_varget, cdfId, 'vPar', vPar 
		ncdf_varget, cdfId, 'amu', amu
		ncdf_varget, cdfId, 'weight', weight

	nCdf_close,	cdfId 

	; constants
	
	k  = 1.3806504e-23
	e_ = 1.60217646e-19
	mi = 1.67262158e-27
	c  = 3.0e8
	
	E_keV = 100.0
	
	T_joule = 2.0/3.0 * E_keV * 1e3 * e_
	vTh = sqrt ( 2.0 * T_joule / (mi*amu[0]) )

	nThermal = 5
	
	nPar = 21
	parRange = vTh * nThermal
	vPar_grid = (fIndGen(nPar)/(nPar-1)-0.5)*2*parRange
	parSize = vPar_grid[1]-vPar_grid[0]
	
	nPer = 11 
	perRange = vTh * nThermal 
	vPer_grid = fIndGen(nPer)/(nPer-1)*perRange
	perSize = vPer_grid[1]-vPer_grid[0]
	
	v_hist = hist_2d(vPar, vPer, $
			bin1=parSize, bin2=perSize, $
			min1=vPar_grid[0],max1=vPar_grid[-1],$
			min2=vPer_grid[0],max2=vPer_grid[-1])
	
	; apply Jacobian
	
	v_hist = v_hist / (2*!pi*transpose(rebin ( vPer_grid+perSize/2.+perSize/2.0, nPer, nPar )))
	
	; apply constant weight
	
	v_hist = v_hist * weight[0]
	
	contour, v_hist, (vPar_grid[0:-2]+parSize/2.0)/c*1e2, vPer_grid/c*1e2, $
			xTitle = 'vPar [%c]', yTitle='vPer [%c]', $
			levels = 10.0^fIndGen(30)*1e-15, /iso

stop
end
