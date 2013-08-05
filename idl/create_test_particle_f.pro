pro create_test_particle_f, $
		single_energy = single_energy, $ ; all particles same energy
		standard_maxwellian = standard_maxwellian, $ ; more particles near vTh
		standard_maxwellian_1d = standard_maxwellian_1d, $ ; more particles near vTh
		weighted_maxwellian = weighted_maxwellian, $ ; uniform v particle distrubution, weighting higher near vTh
		weighted_maxwellian_XYZ = weighted_maxwellian_XYZ, $ ; for kineticj
		cql3d = cql3d, $
		per_offset = per_offset, $ ; % c
		par_offset = par_offset, $ ; % c
		eqdskFName = eqdskFName, $
		plotf = plotf, $
		rsfwc_1d = rsfwc_1d, $ ; set this to the rsfwc_1d.nc output file
		energy_keV = energy_keV, $
        density_m3 = density_m3, $
        n_particles =  n_particles, $
        OutputFileName = OutputFileName, $
		nVTh = nVTh

@constants

e_ = e

; species

amu	= _me_mi 
Z   = -1.0
;amu = 1 
;Z = 1.0
q   = Z * e_
m   = amu * mi

; create f

if keyword_set(n_particles) then nP = n_particles else nP = 500L
if keyword_set(density_m3) then n_m_3 = density_m3 else n_m_3 = 1.1d14
if keyword_set(energy_keV) then E_keV = energy_keV else E_keV = 0.5

nP_test = nP

kT_joule = E_keV * 1d3 * e_
vTh = sqrt ( 2.0*kT_joule / m )

print, 'vTh: ', vTh

if keyword_set(OutputFileName) then fName = OutputFilename else fName = 'f.nc'

; Sanity checking ...

if(E_keV gt 11) then begin
		print, 'ERROR: Do you really want E [keV] to be > 10 keV?'
		stop
endif



if keyword_set(weighted_maxwellian_XYZ) then begin

	; Create a grid to sample the pdf

	nDim = 1
	if keyword_set(nVTh) then nThermal = nVTh else nThermal = 3 

	if nDim eq 3 then begin

		nPts_vx = 100
		vxRange = vTh * nThermal * 2
		vxMin	= -vxRange / 2d0
		vxSize = vxRange / nPts_vx
		vx_grid = dIndGen(nPts_vx)*vxSize+vxMin+vxSize/2d0

		nPts_vy =3 
		vyRange = vTh * nThermal * 2
		vyMin	= -vyRange / 2d0
		vySize = vyRange / nPts_vy
		vy_grid = dIndGen(nPts_vy)*vySize+vyMin+vySize/2d0

		nPts_vz =1
		vzRange = vTh * nThermal * 2
		vzMin	= -vzRange / 2d0
		vzSize = vzRange / nPts_vz
		vz_grid = dIndGen(nPts_vz)*vzSize+vzMin+vzSize/2d0

		vx_grid_3D = rebin ( vx_grid, nPts_vx, nPts_vy, nPts_vz )
		vy_grid_3D = transpose(rebin ( vy_grid, nPts_vy, nPts_vz, nPts_vx ),[2,0,1])
		vz_grid_3D = transpose(rebin ( vz_grid, nPts_vz, nPts_vx, nPts_vy ),[1,2,0])
		nP = nPts_vx * nPts_vy * nPts_vz

		v = sqrt ( vx_grid_3D^2 + vy_grid_3D^2 + vz_grid_3D^2 )
		;f_m_3_analytic = n_m_3 / (sqrt(2*!pi)*vTh)^3 * exp ( -v^2 / (2*vTh^2) )
		f_m_3_analytic = exp ( -v^2 / vTh^2 )

		; Check the density at this point

		dV = n_m_3 / total ( f_m_3_analytic )

		print, "Density: ", total(f_m_3_analytic)*dV

		weight = f_m_3_analytic[*]*dV
		v_x = vx_grid_3D[*]
		v_y = vy_grid_3D[*]
		v_z = vz_grid_3D[*]
		vPer = v_x*0
		vPar = v_y*0
		vMag = sqrt ( v_x^2 + v_y^2 + v_z^2 )
	endif

	if nDim eq 1 then begin

		dx = 1
		dy = 1
		dz = 1

		dvx = 1
		dvy = 1
		dvz = 1

		nPts_vx = 500 
		vxRange = vTh * nThermal * 2
		vxMin	= -vxRange / 2d0
		vxSize = vxRange / nPts_vx
		vx_grid = dIndGen(nPts_vx)*vxSize+vxMin+vxSize/2d0

		dvx = vx_grid[1]-vx_grid[0]

		vy = 0
		vz = 0

		nP = nPts_vx

		vx_grid_1D = vx_grid
		vy_grid_1D = fltArr(nP) + vy
		vz_grid_1D = fltArr(nP) + vz

		v = sqrt ( vx_grid_1d^2 + vy_grid_1d^2 + vz_grid_1d^2 )
		;f_m_3_analytic = n_m_3 / (sqrt(2*!pi)*vTh)^3 * exp ( -v^2 / (2*vTh^2) )
		f_m_3_analytic = n_m_3 / (vTh*sqrt(!pi)) * exp ( -v^2 / (vTh^2) )

		tmp_n_m_3 = 0d0
		for i=1,nP-1 do begin
			tmp_n_m_3 = tmp_n_m_3 + dvx * (f_m_3_analytic[i-1]+f_m_3_analytic[i]) / 2
		endfor

		;scalef = n_m_3 / tmp_n_m_3
		;f_m_3_analytic = f_m_3_analytic * scalef
		;print, "Density: ", total(f_m_3_analytic)*dvx*dvy*dvz
		print, "Density [trap]: ", tmp_n_m_3 

		weight = f_m_3_analytic[*];*dvx;
		v_x = vx_grid_1D[*]
		v_y = vy_grid_1D[*]
		v_z = vz_grid_1D[*]

		; Placeholders
		vPer = v_y*0 
		vPar = v_y*0
		vMag = sqrt ( v_x^2 + v_y^2 + v_z^2 )

		;; Test the weighted f with the un-weighted f.	

		;seed1 = 3.0

		;randX	= randomN ( seed1, nP_test )
		;randY	= randomN ( seed1, nP_test )
		;randZ	= randomN ( seed1, nP_test )

		;totalParticleCnt = n_m_3*dx*dy*dz
		;weight_test = totalParticleCnt / nP_test

		;v_x_test = randX * vTh;/sqrt(2.0)

		;dvx_test = vTh*0.05 
		;dvy_test = 1
		;dvz_test = 1
		;f_test = histogram ( v_x_test/vTh, binSize = dvx_test/vTh, locations = locations )*weight_test/(dvx_test*dvy_test*dvz_test)

		;print, 'Density test: ', total(f_test)*dvx_test*dvy_test*dvz_test

		;p = barplot ( locations+(locations[1]-locations[0])/2.0, f_test,title ='1D f [s/m^2]' )	
		;!null=plot(v_x/vTh,f_m_3_analytic,xtitle='v/vTh',/over)


	endif
	
endif

; Spatial point

x_r = fltArr(nP) + 10.0;(max(eqdsk.rbbbs)-eqdsk.rmaxis)/2 + eqdsk.rmaxis
x_z = fltArr(nP)*0.0
x_p = fltArr(nP)*0.0

x_x = x_r * cos ( x_p )
x_y = x_r * sin ( x_p )

print, 'R: ', x_r[0]

;   Interpolate B to particle locations

if ( NOT keyword_set ( rsfwc_1d ) ) then begin

	if not keyword_set(eqdskFName) then eqdskFName	= 'data/g129x129_1051206002.01120.cmod'
	eqdsk   = readGEQDSK ( eqdskFName )

	r_b0 = eqdsk.r
	z_b0 = eqdsk.z
	b0_r = eqdsk.br
	b0_p = eqdsk.bphi
	b0_z = eqdsk.bz

	b_r  = interpolate ( b0_r, ( x_r - min(r_b0) ) / (max(r_b0)-min(r_b0)) * (n_elements(r_b0)-1.0), $
	    ( x_z - min (z_b0) ) / (max(z_b0)-min(z_b0)) * (n_elements(z_b0)-1.0) )
	b_p  = interpolate ( b0_p, ( x_r - min(r_b0) ) / (max(r_b0)-min(r_b0)) * (n_elements(r_b0)-1.0), $
	    ( x_z - min (z_b0) ) / (max(z_b0)-min(z_b0)) * (n_elements(z_b0)-1.0) )
	b_z  = interpolate ( b0_z, ( x_r - min(r_b0) ) / (max(r_b0)-min(r_b0)) * (n_elements(r_b0)-1.0), $
	    ( x_z - min (z_b0) ) / (max(z_b0)-min(z_b0)) * (n_elements(z_b0)-1.0) )

endif else begin

	cdfId = ncdf_open(rsfwc_1d)

		ncdf_varget, cdfId, 'B0_r', b0_r 
		ncdf_varget, cdfId, 'B0_p', b0_p
		ncdf_varget, cdfId, 'B0_z', b0_z
		ncdf_varget, cdfId, 'r', r_b0

	nCdf_close,	cdfId 

	z_b0 = r_b0 * 0

	b_r  = interpolate ( b0_r, ( x_r - min(r_b0) ) / (max(r_b0)-min(r_b0)) * (n_elements(r_b0)-1.0) )
	b_p  = interpolate ( b0_p, ( x_r - min(r_b0) ) / (max(r_b0)-min(r_b0)) * (n_elements(r_b0)-1.0) )
	b_z  = interpolate ( b0_z, ( x_r - min(r_b0) ) / (max(r_b0)-min(r_b0)) * (n_elements(r_b0)-1.0) )

endelse

bMag    = sqrt ( b_r^2 + b_p^2 + b_z^2 )

bu_r = b_r / bMag
bu_p = b_p / bMag
bu_z = b_z / bMag
	
if keyword_set(weighted_maxwellian) then begin

		print, 'WEIGHTED_MAXWELLIAN'
	; Create a grid to sample the pdf

	nThermal = 5 
	
	nPtsPar = round((-1+sqrt(1+8*nP))/2)
	parRange = vTh * nThermal * 2
	parMin	= -parRange / 2d0
	parSize = parRange / nPtsPar 
	vPar_grid = dIndGen(nPtsPar)*parSize+parMin+parSize/2d0

	nPtsPer = nPtsPar / 2
	perRange = vTh * nThermal 
	perMin	= 0d0
	perSize = perRange / nPtsPer 
	vPer_grid = dIndGen(nPtsPer)*perSize+perMin+perSize/2d0
	
	vPer_grid_2D = transpose(rebin ( vPer_grid, nPtsPer, nPtsPar ))
	vPar_grid_2D = rebin ( vPar_grid, nPtsPar, nPtsPer )

	v = sqrt ( vPer_grid_2D^2 + vPar_grid_2D^2 )
	f_m_3_analytic = n_m_3 / (sqrt(2*!pi)*vTh)^3 * exp ( -v^2 / (2*vTh^2) )

	dV = 2 * !pi * vPer_grid_2D * perSize * parSize

	weight = f_m_3_analytic[*] * dV[*]
	vPar = vPar_grid_2D[*]
	vPer = vPer_grid_2D[*]
	vMag = sqrt ( vPer^2 + vPar^2 )

endif

if keyword_set(standard_maxwellian) then begin

	; Create a template Maxwellian for a single spatial point
	
	randX	= randomN ( undefined, nP )
	randY	= randomN ( undefined, nP )
	randZ	= randomN ( undefined, nP )

	v_x = randX * vTh
	v_y = randY * vTh
	v_z = randZ * vTh
	
	;	Convert velocity vector to cylindrical
	;	pg. 39, Cheng.
	
	v_r	=  cos ( x_p ) * v_x + sin ( x_p ) * v_y
	v_p	= -sin ( x_p ) * v_x + cos ( x_p ) * v_y
	
	vMag = sqrt ( v_x^2 + v_y^2 + v_z^2 )
	vPar = v_r * bu_r + v_z * bu_z + v_p * bu_p
	vPer = sqrt ( vMag^2 - vPar^2 )

	weight = fltArr(nP) + n_m_3 / nP

	if keyword_set(per_offset) then $
			vPer = vPer + per_offset * c
	if keyword_set(par_offset) then $
			vPar = vPar + par_offset * c
	if keyword_set(par_offset) OR keyword_set(per_offset) then $
			vMag = sqrt ( vPer^2 + vPar^2 )

	print, 'Density: ', total (weight)

endif ; standard_mawellian

if keyword_set(standard_maxwellian_1d) then begin

	; Create a template Maxwellian for a single spatial point
	
	seed1 = 3.0
	randX	= randomN ( seed1, nP )
	seed2 	= randX[-1]
	randY	= randomN ( seed2, nP )
	seed3 	= randY[-1]
	randZ	= randomN ( seed3, nP )

	v_x = randX * vTh
	v_y = randY * vTh 
	v_z = randZ * vTh

	;	Convert velocity vector to cylindrical
	;	pg. 39, Cheng.
	
	v_r	=  cos ( x_p ) * v_x + sin ( x_p ) * v_y
	v_p	= -sin ( x_p ) * v_x + cos ( x_p ) * v_y
	
	vMag = sqrt ( v_x^2 + v_y^2 + v_z^2 )
	vPar = v_r * bu_r + v_z * bu_z + v_p * bu_p
	vPer = sqrt ( vMag^2 - vPar^2 )

	weight = fltArr(nP) + n_m_3 / nP

	if keyword_set(per_offset) then $
			vPer = vPer + per_offset * c
	if keyword_set(par_offset) then $
			vPar = vPar + par_offset * c
	if keyword_set(par_offset) OR keyword_set(per_offset) then $
			vMag = sqrt ( vPer^2 + vPar^2 )

	print, 'Density: ', total (weight)

	;h = histogram ( v_x/vTh, binSize = 0.05, locations = locations )
	;p = barplot ( locations, h )	

endif ; standard_mawellian



if keyword_set(single_energy) then begin

	vMag = fltArr(nP) + vTh
	vPar = vMag * 0
	vPer = sqrt ( vMag^2 - vPar^2 )
	weight = fltArr(nP) + n_m_3 / nP

endif ; single energy

if keyword_set(cql3d) then begin

	cdfId = ncdf_open(cql3d)

		ncdf_varGet, cdfId, 'vnorm', vnorm ; velocity (momentum-per-mass) norm cms/sec
		ncdf_varGet, cdfId, 'enorm', enorm
		ncdf_varGet, cdfId, 'f', f ; vnorm**3/(cm**3*(cm/sec)**3)
		ncdf_varGet, cdfId, 'rya', rya
		ncdf_varGet, cdfId, 'x', x ; normalized momentum-per-mass
		ncdf_varGet, cdfId, 'y', y
		ncdf_varGet, cdfId, 'iy_', iy_

	nCdf_close,	cdfId 

	; Create a grid to sample the pdf

	nThermal = 5 
	
	nPtsPar = round((-1+sqrt(1+8*nP))/2)
	parRange = vTh * nThermal * 2
	parMin	= -parRange / 2d0
	parSize = parRange / nPtsPar 
	vPar_grid = dIndGen(nPtsPar)*parSize+parMin+parSize/2d0

	nPtsPer = nPtsPar / 2
	perRange = vTh * nThermal 
	perMin	= 0d0
	perSize = perRange / nPtsPer 
	vPer_grid = dIndGen(nPtsPer)*perSize+perMin+perSize/2d0
	
	vPer_grid_2D = transpose(rebin ( vPer_grid, nPtsPer, nPtsPar ))
	vPar_grid_2D = rebin ( vPar_grid, nPtsPar, nPtsPer )

	vMag_ms_grid = sqrt(vPer_grid_2D^2+vPar_grid_2D^2)
	pitch_rad_grid = acos(vPar_grid_2D / vMag_ms_grid)

	i_rya = 4 
	;i_time = 0
	f_cql = f[0:iy_[i_rya]-1,*,i_rya]*1e6/vnorm^3*1d6

	vMag_ms_cql = x * vnorm * 1d-2
	pitch_rad_cql = y[0:iy_[i_rya]-1,i_rya]

	; Get f_ values at grid locations
	f_interp = interpolate ( f_cql, $
		(pitch_rad_grid-min(pitch_rad_cql))/(max(pitch_rad_cql)-min(pitch_rad_cql))*(n_elements(pitch_rad_cql)-1), $
		(vMag_ms_grid-min(vMag_ms_cql))/(max(vMag_ms_cql)-min(vMag_ms_cql))*(n_elements(vMag_ms_cql)-1))

	vMag_ms_2D_cql = transpose(rebin(vMag_ms_cql,n_elements(vMag_ms_cql),n_elements(pitch_rad_cql)))
	pitch_rad_2D_cql = rebin(pitch_rad_cql,n_elements(pitch_rad_cql),n_elements(vMag_ms_cql))

	vPar_ms_cql = cos(pitch_rad_2D_cql[*]) * vMag_ms_2D_cql[*]
	vPer_ms_cql = sqrt(vMag_ms_2D_cql[*]^2-vPar_ms_cql^2)

	seed = 1.2
	rand = randomN(seed,n_elements(vpar_ms_cql))*1e-4

	vPar_ms_cql += rand
	vPer_ms_cql += abs(rand)
	levels = 10.0^fIndGen(20)*1e-18
	contour,  f_cql[*], vPar_ms_cql/vTh, vPer_ms_cql/vTh, /irreg, levels = levels, /iso

	dV = 2 * !pi * vPer_grid_2D * perSize * parSize

	weight = f_interp[*] * dV[*]
	vPar = vPar_grid_2D[*]
	vPer = vPer_grid_2D[*]
	vMag = sqrt ( vPer^2 + vPar^2 )

endif

pitch = vPar / vMag

J_to_eV = 1/e_

E_eV = m * vMag^2 / 2.0 * J_to_eV

mu = ( amu * mi ) * vPer^2 / ( 2.0 * bMag )

status = intArr(nP)

if keyword_set(plotf) then begin

	v_hist, vPer, vPar, weight, vTh, $
		hist_unweighted = f_unweighted, $
		hist_weighted = f_weighted, $
		hist_KED = f_KED, $
		vPer_grid = vPer_grid, $
		vPar_grid = vPar_grid, $
		vPer2D = vPer_grid_2D, $
		vPar2D = vPar_grid_2D, $
		nBinsPer = 50, $
		nBinsPar = 101

	levels = 10.0^fIndGen(25)*1e-15
	c3 = contour ( f_weighted, vPar_grid/vTh, vPer_grid/vTh, aspect=1.0, c_value = levels ) 
	c4 = contour ( f_KED  , vPar_grid/vTh, vPer_grid/vTh, aspect=1.0, c_value = levels )

	v = sqrt ( vPer_grid_2D^2 + vPar_grid_2D^2 )
	f_m_3_analytic = n_m_3 / (sqrt(2*!pi)*vTh)^3 * exp ( -v^2 / (2*vTh^2) )

	c4 = contour ( f_m_3_analytic, vPar_grid/vTh, vPer_grid/vTh, aspect=1.0, c_value = levels, $
		   /over, c_thick=levels*0+6.0, c_color = strArr(n_elements(levels))+'blue', $
		  transpa = 80 )

	; Calculate the total number of particles to check

	nP_check1 = 0.0
	nP_check2 = 0.0
	nP_check3 = 0.0

	parSize = vPar_grid[1]-vPar_grid[0]
	perSize = vPer_grid[1]-vPer_grid[0]
	nx = n_elements(f_weighted[*,0])
	ny = n_elements(f_weighted[0,*])

	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			dV = 2 * !pi * vPer_grid[j] * parSize * perSize
			nP_check1 += f_KED[i,j] * dV
			nP_check2 += f_weighted[i,j] * dV
			nP_check3 += f_m_3_analytic[i,j] * dV
		endfor
	endfor

	print, nP_check1, nP_check2, nP_check3, n_m_3 

	p = plot(f_m_3_analytic[nx/2,*]*vPer_grid, thick=6.0, transp = 80)
	p = plot(f_weighted[nx/2,*]*vPer_grid, /over)
	p = plot(f_KED[nx/2,*]*vPer_grid, thick = 2.0, /over)

endif

; Write test netCDF file for reading into AORSA

	nc_id = nCdf_create ( fName, /clobber )

	nCdf_control, nc_id, /fill
	
	np_id = nCdf_dimDef ( nc_id, 'nP', nP )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	vPer_id = nCdf_varDef ( nc_id, 'vPer', np_id, /float )
	vPar_id = nCdf_varDef ( nc_id, 'vPar', np_id, /float )
	if (keyword_set(standard_maxwellian) OR keyword_set(standard_maxwellian_1d) OR keyword_set(weighted_maxwellian_XYZ)) $
		AND ( NOT keyword_set(per_offset) ) then begin
		vx_id = nCdf_varDef ( nc_id, 'vx', np_id, /float )
		vy_id = nCdf_varDef ( nc_id, 'vy', np_id, /float )
		vz_id = nCdf_varDef ( nc_id, 'vz', np_id, /float )
		nThermal_id = nCdf_varDef ( nc_id, 'nThermal', scalar_id, /short )
		vTh_id = nCdf_varDef ( nc_id, 'vTh', scalar_id, /float )
	endif
	E_eV_id = nCdf_varDef ( nc_id, 'E_eV', np_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', np_id, /float )
	p_id = nCdf_varDef ( nc_id, 'p', np_id, /float )
	z_id = nCdf_varDef ( nc_id, 'z', np_id, /float )
	x_id = nCdf_varDef ( nc_id, 'x', np_id, /float )
	y_id = nCdf_varDef ( nc_id, 'y', np_id, /float )
	weight_id = nCdf_varDef ( nc_id, 'weight', np_id, /float )
	status_id = nCdf_varDef ( nc_id, 'status', np_id, /short )
	amu_id = nCdf_varDef ( nc_id, 'amu', np_id, /float )
	_Z_id = nCdf_varDef ( nc_id, 'Z', np_id, /short )
	mu_id = nCdf_varDef ( nc_id, 'mu', np_id, /float )

	nCdf_control, nc_id, /enDef
	
	nCdf_varPut, nc_id, r_id, x_r
	nCdf_varPut, nc_id, p_id, x_p 
	nCdf_varPut, nc_id, z_id, x_z 
	nCdf_varPut, nc_id, x_id, x_x 
	nCdf_varPut, nc_id, y_id, x_y 
	nCdf_varPut, nc_id, vPer_id, vPer 
	nCdf_varPut, nc_id, vPar_id, vPar 
	nCdf_varPut, nc_id, weight_id, weight
	nCdf_varPut, nc_id, status_id, status
	nCdf_varPut, nc_id, E_eV_id, E_eV
	nCdf_varPut, nc_id, amu_id, amu + intArr(nP)
	nCdf_varPut, nc_id, _Z_id, Z + intArr(nP)
	nCdf_varPut, nc_id, mu_id, mu
	if (keyword_set(standard_maxwellian) or keyword_set(standard_maxwellian_1d) OR keyword_set(weighted_maxwellian_XYZ)) $
		   	AND ( NOT keyword_set(per_offset) ) then begin
		nCdf_varPut, nc_id, vx_id, v_x 
		nCdf_varPut, nc_id, vy_id, v_y 
		nCdf_varPut, nc_id, vz_id, v_z 
		if(keyword_set(weighted_maxwellian_XYZ)) then $
			nCdf_varPut, nc_id, nThermal_id, nThermal
		nCdf_varPut, nc_id, vTh_id, vTh
	endif

nCdf_close, nc_id

 
end
