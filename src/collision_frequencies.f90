module collision_frequencies

contains

subroutine get_coll_freqs ( amu_a, &
                            amu_b, &
                            T_a, &
                            T_b, &
                            Z_a, &
                            Z_b, &
                            v, &
                            n_a, &
                            n_b, &
                            nu_D_ab, &
                            nu_S_ab, &
                            nu_par_ab, &
                            dnuPardE, &
                            electrons )

    use constants
    use erf_external

    implicit none

    integer, intent(in) :: amu_a, amu_b, Z_a, Z_b
    real, intent(in) :: T_a, T_b, n_a, n_b
    real(kind=dbl) :: vTh_a, vTh_b
    real(kind=dbl), intent(in) :: v
    real(kind=dbl), intent(out) :: nu_D_ab, nu_S_ab, nu_par_ab
    real(kind=dbl), intent(out) :: dnuPardE
    real(kind=dbl) :: nu_ab, e_a, e_b
    real(kind=dbl) :: phi_x, G_x, xVar_a, xVar_b
    real(kind=dbl) :: cLog_, m_a, m_b
    logical, optional, intent(in) :: electrons
    real :: dE
    real(kind=dbl) :: G_x_n, G_x_p, nu_par_ab_n, nu_par_ab_p, &
        v_n, v_p, xVar_n, xVar_p

    m_a = protonMass * amu_a
    m_b = protonMass * amu_b

    if ( present ( electrons ) ) then
      
        m_b     = electronMass 
        cLog_   = 24d0 - log ( sqrt ( n_b / 1d6 ) / T_b )

    else

        cLog_  = 23d0 - log ( &
            Z_a * Z_b * ( amu_a + amu_b ) / ( amu_a * T_b + amu_b * T_a ) &
            * sqrt ( n_a / 1d6 * Z_a**2 / T_a &
                    + n_b / 1d6 * Z_b**2 / T_b ) )
    endif
   
    vTh_a = sqrt ( 2d0 * T_a * e_ / m_a )
    vTh_b = sqrt ( 2d0 * T_b * e_ / m_b )

    xVar_a  = v / vTh_a
    xVar_b  = v / vTh_b

    e_a = Z_a * e_
    e_b = Z_b * e_

    phi_x   = erf_external_fn ( xVar_b )
    G_x   = ( -2d0 * exp ( -xVar_b**2 ) * xVar_b / sqrt ( pi ) + erf_external_fn ( xVar_b ) ) &
        / ( 2d0 * xVar_b**2 )

    !   From Sigmar pg. 38

    nu_ab   = n_b * e_a**2 * e_b**2 * cLog_ / &
        ( 4d0 * pi * e0**2 * m_a**2 * vTh_a**3  )
    
    nu_D_ab = nu_ab * ( phi_x - G_x ) / xVar_a**3
    nu_S_ab = nu_ab * 2d0 * T_a / T_b * (1d0 * m_b / m_a ) * G_x / xVar_a
    nu_par_ab   = 2d0 * nu_ab * G_x / xVar_a

    !   Utilising the energy and slowing down operator presneted by O.A.
    !   Shyshkin @ 35th EPS Conference on Plasma Phys. 9-13 June 2008 requires
    !   the derivative of the nu_par_ab term wrt energy
    
    dE  = 3d0/2d0 * T_a * 0.1
    vTh_a = sqrt ( 2d0 * T_a * e_ / m_a )

    v_n    = sqrt ( 2.0 * (3d0/2d0*T_a-dE) * e_ / m_a )
    v_p    = sqrt ( 2.0 * (3d0/2d0*T_a+dE) * e_ / m_a )
    
    xVar_n  = ( v_n ) / vTh_b 
    xVar_p  = ( v_p ) / vTh_b 
    
    G_x_n   = ( -2d0 * exp ( -xVar_n**2 ) * xVar_n &
        / sqrt ( pi ) + erf_external_fn ( xVar_n ) ) &
        / ( 2d0 * xVar_n**2 )
    G_x_p   = ( -2d0 * exp ( -xVar_p**2 ) * xVar_p &
        / sqrt ( pi ) + erf_external_fn ( xVar_p ) ) &
        / ( 2d0 * xVar_p**2 )

    nu_par_ab_n   = 2d0 * nu_ab * G_x_n / xVar_n
    nu_par_ab_p   = 2d0 * nu_ab * G_x_p / xVar_p

    dnuPardE  = ( nu_par_ab_p - nu_par_ab_n ) / ( 2d0 * dE )  

end subroutine get_coll_freqs

end module collision_frequencies
