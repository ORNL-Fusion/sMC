module read_namelist

implicit none
integer, parameter :: LONG = selected_int_kind(9)
integer(kind=LONG) :: nP = 1000 ! number of particles to trace
character(len=103) :: pl_fileName = 'data/fdis_25keV_D_6_flat_D3D.dav'
character(len=100) :: eqdsk_fileName = 'data/g122080.03100'
character(len=100) :: mchoi_fileName = 'data/mchoi_dlg.nc'
integer :: sR = 6 ! number of grid points to calculate particle gaussian over.
logical :: plotOrbit = .false.
logical :: NFO = .false. ! Compute orbits without finite orbits

!   Background species
integer :: amu_bgnd = 1
integer :: atomicZ_bgnd = 1

!   Background temperature rho profile
real :: E_bgnd_0 = 2d3 ! [eV] Background species maxwellian temp
real :: E_bgnd_lim = 0.052d3
real :: E_bgnd_alpha = 1.3
real :: E_bgnd_beta = 1.9

!   Background density rho profile
real :: density_bgnd_0 = 1.28d20 ! [m^-3] Background species density 
real :: density_bgnd_lim = 3.42d19
real :: density_bgnd_alpha = 0.6
real :: density_bgnd_beta = 1.4

!   Simulating species density rho profile (only for the collision frequency)
real :: density_0 = 1.022d19 ! [m^-3] Background species density 
real :: density_lim = 2.733d18
real :: density_alpha = 0.6
real :: density_beta = 1.4

real :: runTime = 1d-3 ! [s] run time
integer :: nSteps = 1 ! number of steps of length runTime to do
integer :: amu = 1 ! Atomic mass number, i.e., number of protons + neutrons 

logical :: use_QL = .false.
logical :: use_pitchScattering = .true.
logical :: use_kick_diffusion = .false.
integer :: atomicZ = 1 ! Atomic number, i.e., number of protons

integer :: startn = 0 ! time step to start from, useful for continuing.
 
logical :: pList_is_nCDF = .false.
logical :: use_QL_serial = .false.

namelist / sMCIN / nP,&
    pl_fileName, &
    eqdsk_fileName, &
    amu, &
    amu_bgnd, &
    plotOrbit, &
    NFO, &
    runTime, &
    nSteps, &
    atomicZ_bgnd, &
    use_QL, &
    use_pitchScattering, &
    atomicZ, &
    mchoi_fileName, &
    use_kick_diffusion, &
    E_bgnd_0, &
    E_bgnd_lim, &
    E_bgnd_alpha, &
    E_bgnd_beta, &
    density_bgnd_0, &
    density_bgnd_lim, &
    density_bgnd_alpha, &
    density_bgnd_beta, &
    density_0, &
    density_lim, &
    density_alpha, &
    density_beta, &
    pList_is_nCDF, &
    startn, &
    use_QL_serial


contains
subroutine init_namelist

    implicit none
    character(len=100) :: nml_fileName

    nml_fileName    = 'sMC+rf.nml'
    
    open ( 7, file = nml_fileName, delim = 'APOSTROPHE' )
    read ( unit = 7, nml = sMCIN )
    close ( 7 )

end subroutine init_namelist

end module read_namelist
