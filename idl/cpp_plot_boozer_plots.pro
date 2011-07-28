pro cpp_plot_boozer_plots

	fileList = file_search ( 'output/*orbit.nc' )

	for f=0,(n_elements(fileList)-1)<4 do begin

		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'rOrb', r
			ncdf_varget, cdfId, 'pOrb', p
			ncdf_varget, cdfId, 'zOrb', z
			ncdf_varget, cdfId, 'stat', status
			ncdf_varget, cdfId, 'E_eV', E_eV 
			ncdf_varget, cdfId, 'pitch', pitch 
			ncdf_varget, cdfId, 'nu_B', nu_B 
			ncdf_varget, cdfId, 'dt', dt 

		nCdf_close,	cdfId 

		iiBad = where ( pitch lt -1, cnt )
		if cnt gt 0 then stop
		iiBad = where ( pitch gt 1, cnt )
		if cnt gt 0 then stop
		iiBad = where ( E_eV le 0, cnt )
		if cnt gt 0 then stop

		collisionTime = 1.0 / nu_B
		time = total ( dt, /cumu )

		pitchTimes = [0.1, 2.2, 18.2, 200.0]*collisionTime
		ETimes = [5.0, 50.0, 200.0]*collisionTime
		print, ETimes

		nBins = 101
		binMin = -1.0
		binMax = 1.0
		binSize = (binMax-binMin)/nBins
		binEdges = fIndGen(nBins+1)*binSize+binMin
		binCenters = binEdges[0:-2]+binSize/2.0
		pitchHist = fltArr(nBins)

		for i=0,n_elements(pitchTimes)-1 do begin
			iiUse = where ( time le pitchTimes[i], nPts )
			print, i, nPts
			for n=0,nBins-1 do begin
				iiThisBin = where(pitch[iiUse] ge binEdges[n] and pitch[iiUse] lt binEdges[n+1], cnt)	
				if cnt gt 0 then pitchHist[n] = total((dt[iiUse])[iiThisBin])
			endfor
			pitchHist = pitchHist / total(pitchHist)
			if i eq 0 then begin
				p1 = plot ( binCenters, pitchHist, $
					thick=i+1, $
					xRange=[-1.0,1.0], $
				   	layout=[5,2,f+1], /current, dimensions=[1500,600])
			endif else begin
				p1 = plot ( binCenters, pitchHist, $
					thick=i+1, /over)
			endelse
		endfor

		binMin = 0.0
		binMax = 8.0
		binSize = (binMax-binMin)/nBins
		binEdges = fIndGen(nBins+1)*binSize+binMin
		binCenters = binEdges[0:-2]+binSize/2.0
		EHist = fltArr(nBins)
		T_eV = 2e3

		for i=0,n_elements(ETimes)-1 do begin
			iiUse = where ( time le ETimes[i], nPts )
			print, i, nPts
			for n=0,nBins-1 do begin
				iiThisBin = where((E_eV/T_eV)[iiUse] ge binEdges[n] and (E_eV/T_eV)[iiUse] lt binEdges[n+1], cnt)	
				if cnt gt 0 then EHist[n] = total((dt[iiUse])[iiThisBin]) / sqrt(binCenters[n])
			endfor
			p2 = plot ( binCenters, EHist, $
					thick=i+1, $
					xRange=[0,8.0],$
					yRange=[1e-9,1e5],$
					/yLog,$
				   	layout=[5,2,f+1+5], /current)
		endfor

	endfor

end
