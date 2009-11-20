pro plot_sMC

	amu	= 1.0
	e_	= 1.6021e-19
	mi	= amu * 1.6726e-27

	cdfId	= nCdf_open ( 'data/p_ql.nc' )
	nCdf_varGet, cdfId, 'ql_b', ql_b
	nCdf_varGet, cdfId, 'R', ql_R
	nCdf_varGet, cdfId, 'z', ql_z
	nCdf_varGet, cdfId, 'vPer', ql_vPer
	nCdf_varGet, cdfId, 'vPar', ql_vPar 
	nCdf_varGet, cdfId, 'vc_mks', vc_mks
	nCdf_varGet, cdfId, 'pscale', pscale 
	nCdf_varGet, cdfId, 'eNorm', eNorm 
	nCdf_close,	cdfId 

	ql_vPer_binSize    = ql_vPer(2) - ql_vPer(1)
    ql_vPer_min    = min ( ql_vPer ) - ql_vPer_binSize / 2.0
    ql_vPer_max    = max ( ql_vPer ) + ql_vPer_binSize / 2.0
    ql_vPer_range  = ql_vPer_max - ql_vPer_min

	fileList	= file_search ( 'data/pList.dav.nc.*' )

	;	Create my own histograms

	nBins	= 50
	pMinVal	= -1.0
	pRange	= 2.0
	pBinSize	= pRange / nBins
	pBinEdges	= fIndGen ( nBins + 1 ) * pBinSize + pMinVal 
	pBinCenters	= pBinEdges[1:*] - pBinSize / 2.0
	pHist	= pBinCenters * 0.0
	
	eMinVal	= 0.0
	eRange	= 25.5
	eBinSize	= eRange / nBins
	eBinEdges	= fIndGen ( nBins + 1 ) * eBinSize + eMinVal
	eBinCenters	= eBinEdges[1:*] - eBinSize / 2.0
	eHist	= eBinCenters * 0.0
	
	vPerMinVal	= 0.0
	vPerMaxVal	= sqrt ( 2.0 * 10d3 * e_ / mi )
	vPerRange	= vPerMaxVal - vPerMinVal 
	vPerBinSize	= vPerRange / nBins
	vPerBinEdges	= fIndGen ( nBins + 1 ) * vPerBinSize + vPerMinVal
	vPerBinCenters	= vPerBinEdges[1:*] - vPerBinSize / 2.0
	vPerHist	= vPerBinCenters * 0.0
	
	vParMinVal	= -vPerMaxVal
	vParMaxVal	= vPerMaxVal 
	vParRange	= vParMaxVal - vParMinVal 
	vParBinSize	= vParRange / nBins
	vParBinEdges	= fIndGen ( nBins + 1 ) * vParBinSize + vParMinVal
	vParBinCenters	= vParBinEdges[1:*] - vParBinSize / 2.0
	vParHist	= vParBinCenters * 0.0
	
	for i = 0, n_elements ( fileList ) - 1 do begin

		print, i, n_elements(fileList)-1

		ncId	= nCdf_open ( fileList[i] )
		nCdf_varGet, ncId, 'R', R
		nCdf_varGet, ncId, 'z', z
		nCdf_varGet, ncId, 'vPer', vPer
		nCdf_varGet, ncId, 'vPar', vPar
		nCdf_varGet, ncId, 'weight', weight
		nCdf_varGet, ncId, 'status', status 
		nCdf_close,	ncId 
	
		v	= sqrt ( vPer^2 + vPar^2 )
		pitch	= vPar / v
		E	=  mi * v^2 / ( 2.0 * e_ ) / 1e3

		;useMask	= where ( z gt -0.02 and z lt 0.02 and r gt 0.71 and r lt 0.72, $
		;	   comp = discardII ) 
	    ;weight[discardII]	= 0

		if (size ( pHist_all ))[0] eq 0 then begin 

			eIndex	= ( E - eMinVal ) / eRange * nBins 
			eHist[*]	= 0.0
			for ii = 0L, n_elements ( eIndex ) - 1 do begin
				;if eIndex[ii] le n_elements ( eHist ) - 1 then eHist[eIndex[ii]] ++ 
				if eIndex[ii] le n_elements ( eHist ) - 1 then eHist[eIndex[ii]] = eHist[eIndex[ii]] + weight[ii] 
			endfor
			eHist_all	= eHist

			iinan = where(pitch ne pitch, nanCnt)
			if nanCnt gt 0 then pitch[iiNan] = 0
			pIndex	= (( pitch - pMinVal ) / pRange * nBins)<(n_elements(pHist)-1)
			pHist[*]	= 0.0
			for ii = 0L, n_elements ( pIndex ) - 1 do pHist[pIndex[ii]] = pHist[pIndex[ii]] + 1;weight[ii] 
			pHist_all	= pHist

		endif else begin

			eIndex	= ( E - eMinVal ) / eRange * nBins 
			eHist[*]	= 0.0
			for ii = 0L, n_elements ( eIndex ) - 1 do begin
				;if eIndex[ii] le n_elements ( eHist ) - 1 then eHist[eIndex[ii]] ++ 
				if eIndex[ii] le n_elements ( eHist ) - 1 then eHist[eIndex[ii]] = eHist[eIndex[ii]] + weight[ii] 
			endfor
		    eHist_all	= [ [ eHist_all ] , [ eHist ] ]	
			iinan = where(pitch ne pitch, nanCnt)
			if nanCnt gt 0 then pitch[iiNan] = 0
			pIndex	= (( pitch - pMinVal ) / pRange * nBins)<(n_elements(pHist)-1)
			pHist[*]	= 0.0
			for ii = 0L, n_elements ( pIndex ) - 1 do pHist[pIndex[ii]]  = pHist[pIndex[ii]] + 1;weight[ii]
		    pHist_all	= [ [ pHist_all ] , [ pHist ] ]	

		endelse

		vPerZeroII	= where ( vPer eq 0, iiZCnt )
		print, iiZCnt
		weightZeroII	= where ( weight eq 0, iiWCnt )
		print, iiWCnt, total ( pHist )
		device, decomposed = 0
		plot, vPar, vPer, psym = 4, color = 255, $
				   yRange = [1d2,1d7], $
				   yStyle = 1, $
				   charSize = 2

		for ii=0,n_elements(ql_vper)-1 do begin
			oPlot, ql_vPar, ql_vPar*0+ql_vPer(ii), $
				   color = 255
		endfor
stop	
	endfor

	density	= 5e14 
	ETh	= 2.0 ; [keV]

	loadct, 1 
	device, decomposed = 0
	colorStp	= 200 / n_elements ( fileList )
	!p.background = 255

	maxwellian	= density * sqrt ( eBinCenters ) * exp ( - eBinCenters / ETh )

	window, 0, ySize = 1000
	!p.charSize = 2.0
	!p.multi = [0,1,3]
	plot, pBinCenters, pHist_all[*,0], $
		 xRange = [-1,1], color = 0, yRange = [min(pHist_all[*,1:*]), max(pHist_all[*,1:*])]
	
	oPlot, pBinCenters, pBinCenters*0+100, thick = 12, color = 15*16-1

		if n_elements ( pHist_all[*,0] ) gt 1 then $
		for i = 1, n_elements(pHist_all[0,*])-1 do $
			oPlot,pBinCenters, pHist_all[*,i], $
			color = i * colorStp

	oPlot,pBinCenters, pHist_all[*,n_elements(pHist_all[0,*])-1], $
			color = 0, thick = 4


	plot, eBinCenters, eHist_all[*,0],$
		   xRange = [0, eRange],$
		   color = 0, $
		   yRange = [0,max(eHist_all[*,1:*])], thick = 2
	oPlot, eBinCenters, maxwellian, thick = 12, color = 14*16-1

		if n_elements ( eHist_all[*,0] ) gt 1 then $
		for i = 1, n_elements(eHist_all[0,*])-1 do $
			oPlot, eBinCenters, eHist_all[*,i], color = i * colorStp

	oPlot,eBinCenters, eHist_all[*,n_elements(eHist_all[0,*])-1], $
			color = 0, thick = 4

	plot, eBinCenters, eHist_all[*,0] / sqrt ( eBinCenters ),$
		   xRange = [0, eRange], color = 0,$
		   yRange = [max(max(eHist_all[*,1:*],dim=2)/sqrt( eBinCenters ) )*1e-3,max(max(eHist_all[*,1:*],dim=2)/sqrt( eBinCenters ) )], /ylog
	oPlot, eBinCenters, maxwellian / sqrt ( eBinCenters ), thick = 12, color = 14*16-1

		if n_elements ( eHist_all[*,0] ) gt 1 then $
		for i = 1, n_elements(eHist_all[0,*])-1 do $
			oPlot, eBinCenters, eHist_all[*,i]/ sqrt ( eBinCenters ), color = i * colorStp

		oPlot, eBinCenters, eHist_all[*,n_elements(eHist_all[0,*])-1]/ sqrt ( eBinCenters ),$
			   color = 0, thick = 4 

	stop
end
