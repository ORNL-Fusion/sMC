pro create_test_particle_f, $
		plotf = plotf

; constants

k  = 1.3806504e-23
e_ = 1.60217646e-19
mi = 1.67262158e-27
c  = 3.0e8

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

nP    = 10000L
n_m_3 = 4.0e19
E_keV = 2.0

; Create a template Maxwellian for a single spatial point

randX	= randomN ( undefined, nP )
randY	= randomN ( undefined, nP )
randZ	= randomN ( undefined, nP )

T_joule = 2.0/3.0 * E_keV * 1e3 * e_
vTh = sqrt ( 2.0 * T_joule / (mi*amu) )

v_x = randX * vTh
v_y = randY * vTh
v_z = randZ * vTh

; Spatial point

x_r = fltArr(nP) + (max(eqdsk.rbbbs)-eqdsk.rmaxis)/2 + eqdsk.rmaxis
x_z = fltArr(nP)*0.0
x_p = fltArr(nP)*0.0

print, 'R: ', x_r[0]

;	Convert velocity vector to cylindrical
;	pg. 39, Cheng.

v_r	=  cos ( x_p ) * v_x + sin ( x_p ) * v_y
v_p	= -sin ( x_p ) * v_x + cos ( x_p ) * v_y

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

vMag = v_x*0+vTh;sqrt ( v_x^2 + v_y^2 + v_z^2 )
vPar = vMag*0;v_r * bu_r + v_z * bu_z + v_p * bu_p
vPer = sqrt ( vMag^2 - vPar^2 )

;; Offset in vPer/vPar
;; -------------------
;
;vParOffSet = 0.02 * c
;vPerOffSet = 0.02 * c
;
;vPar = vPar + vParOffSet
;vPer = vPer + vPerOffSet
;vMag = sqrt ( vPar^2 + vPer^2 )
;
;; --------------------

pitch = vPar / vMag

J_to_eV = 1/e_

E_eV = m * vMag^2 / 2.0 * J_to_eV

mu = ( amu * mi ) * vPer^2 / ( 2.0 * bMag )

weight = fltArr(nP) + n_m_3 / nP
status = intArr(nP)

if keyword_set(plotf) then begin
	
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
