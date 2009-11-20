module constants
    use read_namelist
    implicit none
    save
   
    integer, parameter :: DBL = selected_real_kind ( 13, 300 )
    real(kind=DBL), parameter :: e_ = 1.602176d-19
    real(kind=DBL), parameter :: e_cgs = 4.803d-10 
    real(kind=DBL) :: mi, mi_bgnd
    real, parameter :: c = 3.0e8
    real(kind=dbl), parameter :: pi = 3.14159265
    real(kind=DBL), parameter :: k_ = 1.3806503d-23
    real(kind=DBL), parameter :: e0 = 8.854187d-12 
    real(kind=dbl), parameter :: protonMass = 1.67262d-27
    real(kind=dbl), parameter :: electronMass = 9.109382d-31 

    real(kind=dbl) :: q
   
contains
    subroutine set_constants 
        implicit none 
 
        mi = amu * protonMass
        mi_bgnd = amu_bgnd * protonMass 
        q = atomicZ * e_

    end subroutine set_constants

end module constants
