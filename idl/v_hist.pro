pro v_hist, vPer, vPar, weight, vTh, $
		nThermal = nThermal, $
		nBinsPer = nBinsPer, $
		nBinsPar = nBinsPar, $
		hist_unweighted = ret_unweighted, $
		hist_weighted = ret_weighted, $
		hist_KED = ret_KED, $
		vPer_Grid = vPer_grid, $
		vPar_Grid = vPar_grid, $
		vPer2D = vPer_grid_2D, $
		vPar2D = vPar_grid_2D


	if not keyword_set(nThermal) then nThermal = 5 
	if not keyword_set(nBinsPer) then nBinsPer = 100 
	if not keyword_set(nBinsPar) then nBinsPar = 51 

	parRange = vTh * nThermal * 2
	parMin	= -parRange / 2d0
	parSize = parRange / nBinsPar 
	vPar_grid = dIndGen(nBinsPar)*parSize+parMin+parSize/2d0

	perRange = vTh * nThermal 
	perMin	= 0d0
	perSize = perRange / nBinsPer 
	vPer_grid = dIndGen(nBinsPer)*perSize+perMin+perSize/2d0
	
	ret_unweighted = hist_nd(transpose([[vPar],[vPer]]),$
			min = [parMin,perMin],$
			max = [parMin+parRange,perMin+perRange],$
			nBins = [ nBinsPar, nBinsPer ], $
			reverse_ind = ri)

	vPer_grid_2D = transpose(rebin ( vPer_grid, nBinsPer, nBinsPar ))
	vPar_grid_2D = rebin ( vPar_grid, nBinsPar, nBinsPer )

	ret_unweighted = ret_unweighted / (2*!pi*vPer_grid_2D*perSize*parSize)

	ret_weighted = fltArr(size(ret_unweighted,/dim))
	ret_KED = fltArr(size(ret_unweighted,/dim))
	
	nx = n_elements(ret_weighted[*,0])
	ny = n_elements(ret_weighted[0,*])

	; Need to find the area for a single kernel for normalization to unity
	sigma = 1.5e5
	single_kernel = fltArr(size(ret_unweighted,/dim))
	for i=0,nx-1 do begin
		for j=0,ny-1 do begin

			expArg = -1*((vPar_grid[i])^2+vPer_grid[j]^2)/(2*sigma^2)
			fHerePerParticle = 2*!pi*exp(expArg)
			single_kernel[i,j] = total(fHerePerParticle)

		endfor
	endfor

	singleKernelArea = total(single_kernel*2*!pi*vPer_grid_2D*perSize*parSize)

	for i=0,nx-1 do begin
		for j=0,ny-1 do begin

			; Standard histogram weighting
			ind = i+nx*j
			if ri[ind] lt ri[ind+1] then begin
				dV = 2*!pi*vPer_grid[j]*parSize*perSize
				ret_weighted[i,j] = total( weight[ri[ri[ind]:ri[ind+1]-1]] ) / dV
			endif

			; KED
			bArg = vPer*vPer_grid[j]/sigma^2
			expArg = -1*(vPer^2+(vPar-vPar_grid[i])^2+vPer_grid[j]^2)/(2*sigma^2)
			fHerePerParticle = 2*!pi*beselI(bArg,0,/double)*exp(expArg)*weight/singleKernelArea
			ret_KED[i,j] = total(fHerePerParticle)

		endfor
	endfor

end
