pro create_test_particle_f, $
		single_energy = single_energy, $ ; all particles same energy
		standard_maxwellian = standard_maxwellian, $ ; more particles near vTh
		weighted_maxwellian = weighted_maxwellian, $ ; uniform v particle distrubution, weighting higher near vTh
		per_offset = per_offset, $ ; % c
		par_offset = par_offset, $ ; % c
		plotf = plotf

; constants

k  = 1.3806504d-23
e_ = 1.60217646d-19
mi = 1.67262158d-27
c  = 3.0d8

; species

amu	= 1.0
Z   = 1.0
q   = Z * e_
m   = amu * mi

; fileNames

fName = 'fdis_D_40keV_D3D'
eqdskFName	= 'data/g129x129_1051206002.01120.cmod'
eqdsk   = readGEQDSK ( eqdskFName )

; create f

nP    = 1000L
n_m_3 = 4d19
E_keV = 2d0
T_joule = 2d0/3d0 * E_keV * 1d3 * e_
vTh = sqrt ( T_joule / (mi*amu) )

; Spatial point

x_r = fltArr(nP) + (max(eqdsk.rbbbs)-eqdsk.rmaxis)/2 + eqdsk.rmaxis
x_z = fltArr(nP)*0.0
x_p = fltArr(nP)*0.0

print, 'R: ', x_r[0]

;   Interpolate B to particle locations

b_r  = interpolate ( eqdsk.bR, ( x_r - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    ( x_z - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0) )
b_p  = interpolate ( eqdsk.bPhi, ( x_r - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    ( x_z - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0) )
b_z  = interpolate ( eqdsk.bz, ( x_r - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    ( x_z - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0) )

bMag    = sqrt ( b_r^2 + b_p^2 + b_z^2 )

bu_r = b_r / bMag
bu_p = b_p / bMag
bu_z = b_z / bMag
	
if keyword_set(weighted_maxwellian) then begin

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

endif ; standard_mawellian


if keyword_set(single_energy) then begin

	vMag = fltArr(nP) + vTh
	vPar = vMag * 0
	vPer = sqrt ( vMag^2 - vPar^2 )
	weight = fltArr(nP) + n_m_3 / nP

endif ; single energy

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

	levels = 10.0^fIndGen(10)*1e-5
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

	nc_id = nCdf_create ( 'data/'+fName+'.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	np_id = nCdf_dimDef ( nc_id, 'nP', nP )
	
	vPer_id = nCdf_varDef ( nc_id, 'vPer', np_id, /float )
	vPar_id = nCdf_varDef ( nc_id, 'vPar', np_id, /float )
	E_eV_id = nCdf_varDef ( nc_id, 'E_eV', np_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', np_id, /float )
	p_id = nCdf_varDef ( nc_id, 'p', np_id, /float )
	z_id = nCdf_varDef ( nc_id, 'z', np_id, /float )
	weight_id = nCdf_varDef ( nc_id, 'weight', np_id, /float )
	status_id = nCdf_varDef ( nc_id, 'status', np_id, /short )
	amu_id = nCdf_varDef ( nc_id, 'amu', np_id, /short )
	_Z_id = nCdf_varDef ( nc_id, 'Z', np_id, /short )
	mu_id = nCdf_varDef ( nc_id, 'mu', np_id, /float )

	nCdf_control, nc_id, /enDef
	
	nCdf_varPut, nc_id, r_id, x_r
	nCdf_varPut, nc_id, p_id, x_p 
	nCdf_varPut, nc_id, z_id, x_z 
	nCdf_varPut, nc_id, vPer_id, vPer 
	nCdf_varPut, nc_id, vPar_id, vPar 
	nCdf_varPut, nc_id, weight_id, weight
	nCdf_varPut, nc_id, status_id, status
	nCdf_varPut, nc_id, E_eV_id, E_eV
	nCdf_varPut, nc_id, amu_id, amu + intArr(nP)
	nCdf_varPut, nc_id, _Z_id, Z + intArr(nP)
	nCdf_varPut, nc_id, mu_id, mu

nCdf_close, nc_id

stop 
end
