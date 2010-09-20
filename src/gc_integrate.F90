!   ToDo List:
!       Overhaul the nanCatch crap in the collision operator. Or at the very
!       least add a nanCatch counter to see how important an issue it is.
!       Perhaps collect statistics or error codes about where the nans are
!       caught. Also, put in gotos following a nanCatch so it only stores the
!       first catch.

module gc_integrate

contains
    function dlg_gc_velocity ( pos, vPerp, vPar, NFO )
        use gc_terms
        use interp
        use eqdsk
        implicit none

        real :: dlg_gc_velocity(3)
        real, intent(IN) :: pos(3), vPerp, vPar
        real :: grad_R, grad_phi, grad_z, &
            curv_R, curv_phi, curv_z, &
            unitb_R, unitb_phi, unitb_z
        real :: surf2, bMagHere, bHere(3), &
            vgc_R, vgc_phi, vgc_z
        logical, intent(in) :: NFO

        grad_R  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bGradient_R, nw, zp_bGrad_R, sigma )
        grad_phi  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bGradient_phi, nw, zp_bGrad_phi, sigma )
        grad_z  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bGradient_z, nw, zp_bGrad_z, sigma )

        curv_R  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bCurvature_R, nw, zp_bCurv_R, sigma )
        curv_phi  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bCurvature_phi, nw, zp_bCurv_phi, sigma )
        curv_z  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bCurvature_z, nw, zp_bCurv_z, sigma )

        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )
        unitb_R = bHere(1) / bMagHere
        unitb_phi   = bHere(2) / bMagHere
        unitb_z = bHere(3) / bMagHere

        if ( NFO ) then
 
            vgc_R   = vPar * unitb_R 
            vgc_phi   = vPar * unitb_phi 
            vgc_z   = vPar * unitb_z 

        else
            
            vgc_R   = vPar * unitb_R + vPerp**2 * grad_R + vPar**2 * curv_R 
            vgc_phi   = vPar * unitb_phi + vPerp**2 * grad_phi + vPar**2 * curv_phi
            vgc_z   = vPar * unitb_z + vPerp**2 * grad_z + vPar**2 * curv_z

        endif

        dlg_gc_velocity = (/ vgc_R, vgc_phi, vgc_z /)

    end function dlg_gc_velocity

    function dlg_vPerp ( pos, u )
        use constants
        implicit none

        real :: bHere(3),bMagHere
        real, intent(IN) :: pos(3)
        real(kind=DBL), intent(in) :: u
        real :: dlg_vPerp

        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )
        dlg_vPerp   = sqrt ( 2d0 * u * bMagHere / mi )

    end function dlg_vPerp

    function dlg_vPar ( pos, u )
        use constants
        use gc_terms
        use interp
        use eqdsk
        implicit none
        
        real :: dlg_vPar, bDotGradB_here, surf2
        real, intent(IN) :: pos(3)
        real(kind=DBL), intent(IN) :: u

        bDotGradB_here  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bDotGradB, nw, zp_bDotGradB, sigma )

        dlg_vPar    = -u / mi * bDotGradB_here 

    end function dlg_vPar

    function dlg_interpB ( pos, bMagHere )
        use eqdsk
        use interp
        implicit none
        
        real :: bR_here, bPhi_here, bz_here
        real, intent(IN) :: pos(3)
        real :: dlg_interpB(3)
        real, optional, intent(OUT) :: bMagHere
        real :: surf2

        bR_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bR, nw, zp_bR, sigma )
        bPhi_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bPhi, nw, zp_bPhi, sigma )
        bz_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bz, nw, zp_bz, sigma )

        if ( present (bMagHere) ) &
            bMagHere    = sqrt ( bR_here**2 + bPhi_here**2 + bz_here**2 )

        !if ( bR_here .ne. bR_here &
        !    .or. bPhi_here .ne. bPhi_here &
        !    .or. bz_here .ne. bz_here ) then 

        !    write(*,*) 'dlg: ERROR - dlg_intepB NaN'
        !    write(*,*) 'pos: ', pos
        !    write(*,*) 'bHere: ', bR_here, bPhi_here, bz_here

        !endif

        dlg_interpB(1)  = bR_here 
        dlg_interpB(2)  = bPhi_here
        dlg_interpB(3)  = bz_here 
    
    end function dlg_interpB

    subroutine gc_orbit ( start_R, start_z, start_vPerp, &
         start_vPar, start_weight, start_status, final_R, final_z, &
         final_vPerp, final_vPar, weight_out, status_out, plot )
        use eqdsk
        use gc_terms
        use constants
        use luxury
        use communications

#if USE_DISLIN
        use dislin 
#endif
        use interp
        !use rzvv_grid
        use init_mpi
        use read_namelist
        use read_ql
        use read_mchoi
        use powerAbsGrid
        use erf_external
        use collision_frequencies
        use bessJ_mod
        use airy_functions_real_single

        implicit none
       
        real, intent(IN) :: start_R, start_z,&
            start_vPerp, start_vPar, start_weight
        integer, intent(in) :: start_status
            
        real, intent(out) :: final_R, final_z, &
            final_vPerp, final_vPar, weight_out
        integer, intent(out) :: status_out

        integer(kind=LONG), parameter :: maxSteps = 50000
        integer(kind=LONG) :: stepCnt
        integer :: var_dt, i 
        real(kind=DBL) :: u, u_n
        real :: tau, dTau, pos(3), pos_old(3), &
            bMagHere, bHere(3), vPerp, vPar, &
            vgc(3), dt, dtMin, dtMax, vPerp_n, vPar_n, &
            vPerp_old, vPar_old, bHere_n(3), bMagHere_n, &
            pos_extra_old(3), pos_extra_old_n(3)
        logical :: stillIn, firstOrbit
        logical, optional, intent(IN) :: plot
        real :: k1_vPar, k1_vgc(3), k2_vPar, k2_vgc(3), &
            k3_vPar, k3_vgc(3), k4_vPar, k4_vgc(3), &
            psi_here, surf2, vgc_n(3), k1_vgc_n(3), k2_vgc_n(3), &
            k3_vgc_n(3), k4_vgc_n(3), pos_(3), k1_vPar_n, &
            k2_vPar_n, k3_vPar_n, k4_vPar_n, EStep, EStep_NFO
        
        real :: T, T_bgnd, coulombLog, vTh, start_T, rho_here
        real :: E_bgnd, density_bgnd, density
        
#if USE_DISLIN
 
        !   Plotting variables
    
        integer :: nLevs, nxPag, nyPag
        real :: levStep
        real, allocatable :: levels(:)
        real, allocatable :: rTrack(:), zTrack(:), &
            vPerpTrack(:), vParTrack(:), &
            rTrack_nfo(:), zTrack_nfo(:), vPerpTrack_nfo(:), &
            vParTrack_nfo(:)

#endif
        
        character(len=15) :: randPhrase
        real(kind=DBL) :: xVar, vB, vB_, vD, v0, pitch0, pitch1
        real(kind=DBL) :: energy0, energy1, phi_x, psi_x
        real(kind=DBL) :: vPerp_update, vPar_update, v0_update
        real :: genUnf, rand_sgn, rand_sgnE, &
            rand_lux_sgn(50000), rand_lux_sgnE(50000), &
            rand_sgnE_total
        integer :: seed1, seed2
        real(kind=DBL) :: dE, v0_n, v0_p, xVar_n, xVar_p, &
            psi_x_n, psi_x_p, vE_n, vE_p, dvE_dE, vE, xVarE, &
            psi_xE, dLnvE_dLnE
        integer :: values(8)
        character(len=8) :: date
        character(len=10) :: time
        character(len=5) :: zone

        !   QL variables

        integer :: ql_rIndex, ql_zIndex, ql_vPerIndex, ql_vParIndex
        integer :: mpi_coords_hasData(2), mpi_pId_hasData
        real :: ql_recv(4), ql_bHere, ql_cHere, ql_eHere, ql_fHere
        real :: ql_bHereCyl, ql_cHereCyl, ql_eHereCyl, ql_fHereCyl
        real :: ql_thHere
        real :: ql_bHereP, ql_cHereP, ql_eHereP, ql_fHereP
        real :: ql_bHereCylP, ql_cHereCylP, ql_eHereCylP, ql_fHereCylP
        real :: ql_thHereP, divD_per1
 
        integer, dimension ( MPI_STATUS_SIZE ) :: mpi_stat_getData

        real :: kPar, antFreq_rads, omegaHere, resonanceCondition
        integer :: harmonicNumber, toroidalModeNumber
        real :: antFreq 
        real(kind=dbl) :: bArg, dvPerp, dvPar
        real :: powerAdd, powerParAdd
        real(kind=DBL) :: bF
        real :: dS, bF1, bF2, bHere_plus(3), bHere_minu(3), &
            bMagHere_plus, bMagHere_minu, dBdS, &
            dBdS_, bOld(3), bMagOld, &
            dvPardS
        real(kind=dbl) :: interactionTime, interactionTimeB,&
            interactionTimeC, interactionTimeMin
        real(kind=dbl) :: bfn_m2, bFn_m1, bFn, bFn_p1, bFn_p2
        real(kind=dbl) :: D_vPerp, mean_vPerp
        complex :: ePlusHere, eMinuHere, kPerHere
        real :: ePlusRealHere, eMinuRealHere, kPerRealHere
        real :: ePlusImgHere, eMinuImgHere, kPerImgHere
        integer :: power_R_index, power_z_index
        logical :: wasResonant
        real :: rcA, rcB
        real :: bHereLeft(3), bHereRight(3), bMagLeft, bMagRight, &
            omegaLeft, omegaRight, rcLeft, rcRight
        real :: distance, distance_NFO
        real(kind=dbl) :: distance3D
        real :: xx, yy, zz, xxOld, yyOld, zzOld
        real :: posStep1(3)
        logical :: nanCatch
        integer :: nanErrorMsg
        real(kind=dbl) :: lVar, Cf, lVarV, G_lVarV, vE_, vD_
        real(kind=dbl) :: nu_D_ii, nu_S_ii, nu_par_ii, dnuPardE_ii
        real(kind=dbl) :: nu_D_ei, nu_S_ei, nu_par_ei, dnuPardE_ei
        real(kind=dbl) :: nu_D, nu_S, nu_par, dnuPardE, dE_ii, dE_ei, dE_1
        real(kind=dbl) :: vPerpK, vPerpK_old
        real :: posK(3), vParK, pos_old_n(3), posK_old(3), &
            vPar_old_n, vPerp_old_n, vParK_old, &
            bHereK(3), bMagHereK, uK, phase, randNo

        real(kind=dbl) :: c1, c2, d2BdS2, res_accel
        real :: bExtraOld(3), bMagExtraOld, posK_extra_old(3)
        real :: ai, dai, airyArg
        integer :: airyErr
        logical :: kickNaN, collNaN
        integer :: kickErr, collErr
        real :: dtR, dtZ, dtvPer, dtvPar, dt_rfSerial
        real :: time_off_QLgrid, time_bad_coll

        !   Initialize variables
        
        stepCnt = 1
        firstOrbit  = .true.
        stillIn = .true.
        dTau    = 0.0
        tau = 0.0
        var_dt  = 1
        weight_out = start_weight 
        i = 0
        wasResonant = .false.
        rcA = 0
        rcB = 0
        nanCatch    = .false.
        nanErrorMsg = -9999 
        time_off_QLgrid = 0
        time_bad_coll = 0
        status_out  = start_status

        !   Initialize the random variables

        call date_and_time ( date, time, zone, values )
        call rLuxGo ( 223, abs ( values(8)*values(7)*(mpi_pId+1) ) + 1, 0, 0)
        call ranLux ( rand_lux_sgn, 50000 )
        call rLuxGo ( 223, abs ( values(8)*values(6)*(mpi_pId+1) ) + 1, 0, 0)
        call ranLux ( rand_lux_sgnE, 50000 )

        rand_sgnE_total = 0.0

        pos(1)  = start_R
        pos(2)  = 0.0
        pos(3)  = start_z

        vPerp   = start_vPerp
        vPar    = start_vPar

        if ( start_status < 0 ) go to 643

        start_T   = mi * ( vPar**2 + vPerp**2 ) / ( 3.0 * e_ )
 
        !if ( NFO ) then 

        !    pos_    = pos
        !    vPar_n  = vPar
        !    vPerp_n = vPerp

        !endif

        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )

        u   = mi * vPerp**2 / ( 2.0 * bMagHere ) 
        !if ( NFO ) u_n   = mi * vPerp_n**2 / ( 2.0 * bMagHere ) 

        dtMin   = 1e-10 
        dtMax   = 5e-6
        
        dt  = dtMin

#if USE_DISLIN

        if ( present ( plot ) .and. mpi_pId == 0 ) then 
            if ( plot ) then  
 
                allocate ( rTrack(maxSteps), zTrack(maxSteps), &
                    vPerpTrack(maxSteps), vParTrack(maxSteps) )
           
                call scrMod ( 'REVERS' )
                call setPag ( 'DA4P' )! nxPag = 2100, nyPag = 2970
                call metaFl ( 'XWIN' )
                call disIni ()
                call winMod ( 'NONE' )
                call erase ()
                call noChek ()
                call unit ( 0 )
                call getPag ( nxPag, nyPag )
                call axsPos ( 300, 1700 )
                call axsLen ( 1500, 1400 )
                call graf ( 0.5, 1.0, 0.5, 0.25, -0.25, 0.25, -0.25, 0.25 ) 
                call noClip ()
                call color ('BLUE')         
                
                nLevs   = 10 
                levStep    = (sibry-simag)/nLevs 
                if ( .not. allocated ( levels ) ) &
                    allocate ( levels(nLevs) )
                levels = (/ (i*levStep,i=0,nLevs) /)+simag
                do i=1,nLevs
                    call contur ( r, nw, z, nh, psizr, levels(i) ) 
                end do
                call color ( 'MAGENTA' )
                call curve ( rbbbs, zbbbs, nbbbs )

                call endGrf ()
                call color ('RED' )
                !call setGrf ( 'NONE', 'NONE', 'NONE', 'NONE' ) 

            end if 
        end if 
#endif

        main_loop: &
        do 

            if ( vPar .ne. vPar .or. vPerp .ne. vPerp &
                .or. abs(pos(1)) .gt. 1e4 .or. abs(pos(3)) .gt. 1e4 ) then 
                nanCatch = .true.
                nanErrorMsg = -332
                go to 72
            end if
  
            vPerp_old   = vPerp
            vPar_old    = vPar
            pos_extra_old   = pos_old
            pos_old     = pos

            !if (NFO) then

            !        if ( vPar_n .ne. vPar_n .or. vPerp_n .ne. vPerp_n &
            !        .or. abs(pos_(1)) .gt. 1e4 .or. abs(pos_(3)) .gt. 1e4 ) then 
            !            nanCatch = .true.
            !            nanErrorMsg = -332
            !            go to 72
            !        end if
 
            !        pos_extra_old_n   = pos_old_n
            !        pos_old_n = pos_
            !        vPar_old_n  = vPar_n
            !        vPerp_old_n = vPerp_n
            !endif

            ! Include finite orbits

            vPerp   = dlg_vPerp ( pos, u ) 
            vgc = dlg_gc_velocity ( pos, vPerp, vPar, nfo = nfo )
            k1_vPar   = dt * dlg_vPar ( pos, u ) 
            k1_vgc  = dt * vgc
            
            vPerp   = dlg_vPerp ( pos + k1_vgc / 2, u ) 
            vgc = dlg_gc_velocity ( pos + k1_vgc / 2, vPerp, vPar + k1_vPar / 2, nfo = nfo )
            k2_vPar   = dt * dlg_vPar ( pos + k1_vgc / 2, u ) 
            k2_vgc  = dt * vgc
           
            vPerp   = dlg_vPerp ( pos + k2_vgc / 2, u ) 
            vgc = dlg_gc_velocity ( pos + k2_vgc / 2, vPerp, vPar + k2_vPar / 2, nfo = nfo )
            k3_vPar   = dt * dlg_vPar ( pos + k2_vgc / 2, u ) 
            k3_vgc  = dt * vgc
            
            vPerp   = dlg_vPerp ( pos + k3_vgc, u ) 
            vgc = dlg_gc_velocity ( pos + k3_vgc, vPerp, vPar + k3_vPar, nfo = nfo )
            k4_vPar   = dt * dlg_vPar ( pos + k3_vgc, u ) 
            k4_vgc  = dt * vgc
           
            vPar    = vPar + ( k1_vPar + 2 * k2_vPar + 2 * k3_vPar &
                + k4_vPar ) / 6
           
            pos   = pos + ( k1_vgc + 2 * k2_vgc + 2 * k3_vgc + k4_vgc ) / 6
            if ( pos(2) .ge. 2d0 * pi ) pos(2) = pos(2) - 2d0 * pi
            if ( pos(2) < 0 ) pos(2) = pos(2) + 2d0 * pi

                !   Adjust time step prior to collision operators so we
                !   can base it on conservation of energy between steps.
                !   Yeah, I know not the greatest way since the error will 
                !   accumulate, but OK for now.
 
            if ( vPar .ne. vPar .or. vPerp .ne. vPerp &
            .or. abs(pos(1)) .gt. 1e4 .or. abs(pos(3)) .gt. 1e4 ) then 
                nanCatch = .true.
                nanErrorMsg = -377
                go to 72
            end if
                

            EStep = abs ( sqrt ( vPerp**2 + vPar**2 ) &
                    - sqrt ( vPerp_old**2 + vPar_old**2 ) ) &
                    / sqrt ( start_vPerp**2 + start_vPar**2 ) 
  
            !pos = pos_old
            !vPerp   = vPerp_old
            !vPar    = vPar_old 

            !no_finite_orbits: &
            !if ( NFO ) then
           
            !    vPerp_n   = dlg_vPerp ( pos_, u_n ) 
            !    vgc_n = dlg_gc_velocity ( pos_, vPerp_n, vPar_n, NFO = 1 )
            !    k1_vPar_n   = dt * dlg_vPar ( pos_, u_n ) 
            !    k1_vgc_n  = dt * vgc_n
            !    
            !    vPerp_n   = dlg_vPerp ( pos_ + k1_vgc_n / 2.0, u_n ) 
            !    vgc_n = dlg_gc_velocity ( pos_ + k1_vgc_n / 2.0, vPerp_n, vPar_n + k1_vPar_n / 2.0, NFO = 1 )
            !    k2_vPar_n   = dt * dlg_vPar ( pos_ + k1_vgc_n / 2.0, u_n ) 
            !    k2_vgc_n  = dt * vgc_n
           
            !    vPerp_n   = dlg_vPerp ( pos_ + k2_vgc_n / 2.0, u_n ) 
            !    vgc_n = dlg_gc_velocity ( pos_ + k2_vgc_n / 2.0, vPerp_n, vPar_n + k2_vPar_n / 2.0, NFO = 1 )
            !    k3_vPar_n   = dt * dlg_vPar ( pos_ + k2_vgc_n / 2.0, u_n ) 
            !    k3_vgc_n  = dt * vgc_n
            !    
            !    vPerp_n   = dlg_vPerp ( pos_ + k3_vgc_n, u_n ) 
            !    vgc_n = dlg_gc_velocity ( pos_ + k3_vgc_n, vPerp_n, vPar_n + k3_vPar_n, NFO = 1 )
            !    k4_vPar_n   = dt * dlg_vPar ( pos_ + k3_vgc_n, u_n ) 
            !    k4_vgc_n  = dt * vgc_n
           
            !    vPar_n    = vPar_n + ( k1_vPar_n + 2.0 * k2_vPar_n + 2.0 * k3_vPar_n &
            !        + k4_vPar_n ) / 6.0
           
            !    pos_   = pos_ + ( k1_vgc_n + 2.0 * k2_vgc_n + 2.0 * k3_vgc_n + k4_vgc_n ) / 6.0
            !    if ( pos_(2) .ge. 2d0 * pi ) pos_(2) = pos_(2) - 2d0 * pi
            !    if ( pos_(2) < 0 ) pos_(2) = pos_(2) + 2d0 * pi

            !    EStep_NFO = abs ( sqrt ( vPerp_n**2 + vPar_n**2 ) &
            !            - sqrt ( vPerp_old_n**2 + vPar_old_n**2 ) ) &
            !            / sqrt ( start_vPerp**2 + start_vPar**2 ) 

            !    EStep   = minVal ( (/ EStep, EStep_NFO /) ) 

            !    if ( vPar_n .ne. vPar_n .or. vPerp_n .ne. vPerp_n &
            !    .or. abs(pos_(1)) .gt. 1e4 .or. abs(pos_(3)) .gt. 1e4) then 
            !        nanCatch = .true.
            !        nanErrorMsg = -488
            !        go to 72
            !    end if
 
            !endif no_finite_orbits
            
            !   Interpolate to get psi at new pos
            !   and if psi(pos) is outside the last 
            !   closed flux surface then dump the particle
            !   The 0.99 is a factor suggested by EFJ to keep
            !   away boundary of the flux grid.

            !if ( mpi_pId .eq. 0 ) then

            !    write(*,*) '------------------------------'
            !    write(*,'(a20,1x,3f8.3)') 'pos: ', pos
            !    write(*,'(a20,1x,3f8.3)') 'pos_old: ', pos_old
            !    write(*,'(a20,1x,3f8.3)') 'diff: ', pos-pos_old
            !    write(*,'(a20,1x,3f8.3)') 'distance: ', (pos-pos_old)*(/1e0,pos(1),1e0/)
            !    write(*,'(a20,1x,3e12.2)') 'v: ', vPerp, vPar
            !    write(*,'(a20,1x,3e12.2)') 'v old: ', vPerp_old, vPar_old
            !    write(*,'(a20,1x,3f8.3)') 'total d: ', sqrt(sum(((pos-pos_old)*(/1e0,pos(1),1e0/))**2))
            !    write(*,'(a20,1x,3f8.3)') 'total d R,phi,z: ', sqrt( &
            !        (pos(3)-pos_old(3))**2 &
            !        +pos(1)**2+pos_old(1)**2 &
            !        -2.0*pos(1)*pos_old(1)*cos( (pos(2)-pos_old(2)) ) )

            !endif

            !if (NFO) then

            !    psi_here = surf2 ( pos_(1), pos_(3), nw, nh, r, z, &
            !        psizr, nw, zp_psi, sigma )
            !    rho_here    = sqrt ( ( psi_here - siMag ) / abs ( siMag - siBry ) )
 
            !else

                psi_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
                    psizr, nw, zp_psi, sigma )
                rho_here    = sqrt ( ( psi_here - siMag ) / abs ( siMag - siBry ) )

            !endif

            if ( psi_here > sibry * 0.99 ) then 
                if ( mpi_pId == 0 ) write(*,*) 'Wall'
                nP_wall = nP_wall + 1
                weight_out = 0.0
                status_out  = -8000
                stillIn = .false.
            end if

            tau = tau + dt
            if ( tau > runTime ) firstOrbit = .false.

            if ( stepCnt+1 >= maxSteps .and. plot ) firstOrbit = .false.

            v0  = sqrt ( vPar**2 + vPerp**2 )
            !if  ( NFO ) then 
            !    v0  = sqrt ( vPar_n**2 + vPerp_n**2 )
            !endif

            T   = mi * v0**2 / ( 3.0 * e_ )

            if ( stepCnt+1 >= maxSteps*10 ) then 
                
                write(*,"('Too slow:', f8.2,' eV ', i10.10, 3(2x,f12.6))") &
                    0.5 * mi * v0**2 / e_ , stepCnt, tau/runTime, tau, runTime
                nP_bad  = nP_bad + 1
                weight_out = 0.0
                firstOrbit = .false. 
                status_out  = -4000
            end if

72          continue 

            if ( nanCatch ) then

                write(*,*) 'dlg: ERROR - trace - ', nanErrorMsg
                write(*,*) '---------------------------------'
                write(*,*) 'stepCnt: ', stepCnt
                write(*,*) 'dt: ', dt
                write(*,*) 'vPer: ', vPerp
                write(*,*) 'vPar: ', vPar
                write(*,*) 'vPer_old: ', vPerp_old
                write(*,*) 'vPar_old: ', vPar_old
                write(*,*) 'pos: ', pos
                write(*,*) 'pos_old: ', pos_old
                write(*,*) 'k?_vPar: ', k1_vPar, k2_vPar, k3_vPar, k4_vPar
                write(*,*) 'k1_vgc: ', k1_vgc
                write(*,*) 'k2_vgc: ', k2_vgc
                write(*,*) 'k3_vgc: ', k3_vgc
                write(*,*) 'k4_vgc: ', k4_vgc
                write(*,*) 'u:', u
                !if ( NFO ) then
                !write(*,*) 'vPer_n: ', vPerp_n
                !write(*,*) 'vPar_n: ', vPar_n
                !write(*,*) 'vPer_old_n: ', vPerp_old_n
                !write(*,*) 'vPar_old_n: ', vPar_old_n
                !write(*,*) 'pos_: ', pos_
                !write(*,*) 'pos_old_n: ', pos_old_n
                !write(*,*) 'k?_vPar_n: ', k1_vPar_n, k2_vPar_n, k3_vPar_n, k4_vPar_n
                !write(*,*) 'k1_vgc_n: ', k1_vgc_n
                !write(*,*) 'k2_vgc_n: ', k2_vgc_n
                !write(*,*) 'k3_vgc_n: ', k3_vgc_n
                !write(*,*) 'k4_vgc_n: ', k4_vgc_n
                !write(*,*) 'u_n:', u_n
                !endif
                
                nP_bad  = nP_bad + 1
                weight_out  = 0.0
                status_out  = nanErrorMsg
                 
            endif

            exitLoop: &
            if ( (.not. stillIn) .or. (.not. firstOrbit) .or. nanCatch) then
               
                if ( mpi_pId == 0 ) &
                    write(*,"(1x,'StepCnt: ',i9,3x,'t: ',e8.2,3x,'vPer: ',"// & 
                                "f10.1,3x,'vPar: ',f10.1,3x,'E[eV]: ',f9.1,3x,'dE: ',f8.2,3x,L4)") &
                        stepCnt, tau, vPerp, vPar, 3.0/2.0 * T, (T-start_T) / start_T * 100.0, wasResonant 
                exit

            endif exitLoop


        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )
       ! if ( NFO ) bHere_n   = dlg_interpB ( pos_, bMagHere = bMagHere_n )

        distance    = sqrt ( ( pos_old(1) - pos(1) )**2 + &
            ( pos_old(3) - pos(3) )**2 )
        !if ( NFO ) then

        !    distance_nfo   = sqrt ( ( pos_old_n(1) - pos_(1) )**2 + &
        !        ( pos_old_n(3) - pos_(3) )**2 )
        !    distance = maxVal ( (/ distance, distance_NFO /) )

        !endif


    !   Update to new COM magnetic moment
            
        u   = mi * vPerp**2 / ( 2.0 * bMagHere ) 
        !if ( NFO ) u_n   = mi * vPerp_n**2 / ( 2.0 * bMagHere_n ) 

    !   +------------------------------------
    !   Fixed k diffusion 
    !

        kick_diffusion: &
        if ( use_kick_diffusion ) then 


                posK = pos
                posK_old    = pos_old
                posK_extra_old  = pos_extra_old
                vPerpK  = vPerp
                vParK   = vPar
                vParK_old   = vPar_old

            !   Did the particle cross a resonant surface?

            bHereK  = dlg_interpB ( posK, bMagHere = bMagHereK )

            toroidalModeNumber = 10
            kPar  = toroidalModeNumber / posK(1) 
            antFreq = 80e6
            antFreq_rads = antFreq * 2.0 * pi
            harmonicNumber  = 1
            omegaHere   = q * bMagHereK / mi
            resonanceCondition  = &
                antFreq_rads - harmonicNumber * omegaHere - kPar * vParK

            power_R_index = ( posK(1) - power_R_min ) / power_R_range * power_R_nBins + 1
            power_z_index = ( posK(3) - power_z_min ) / power_z_range * power_z_nBins + 1

            rcA = rcB
            rcB = resonanceCondition

            if_resonant: &
            if ( rcA * rcB < 0 .and. stepCnt .gt. 2 ) then 
                
                kickNaN = .false.
                kickErr = -7777 
                wasResonant = .true.

            !   Get dB/dS here - both along field and along orbit

            !   Get E+, E- and kPer here
            
                ePlusRealHere  = surf2 ( posK(1), posK(3), mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, &
                    mchoi_ePlus_real, mchoi_nR, zp_ePlus_real, sigma )
                eMinuRealHere  = surf2 ( posK(1), posK(3), mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, &
                    mchoi_eMinu_real, mchoi_nR, zp_eMinu_real, sigma )
                kPerRealHere  = surf2 ( posK(1), posK(3), mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, &
                    mchoi_kPer_real, mchoi_nR, zp_kPer_real, sigma )
 
                ePlusImgHere  = surf2 ( posK(1), posK(3), mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, &
                    mchoi_ePlus_img, mchoi_nR, zp_ePlus_img, sigma )
                eMinuImgHere  = surf2 ( posK(1), posK(3), mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, &
                    mchoi_eMinu_img, mchoi_nR, zp_eMinu_img, sigma )
                kPerImgHere  = surf2 ( posK(1), posK(3), mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, &
                    mchoi_kPer_img, mchoi_nR, zp_kPer_img, sigma )

                ePlusHere  = cmplx ( ePlusRealHere, ePlusImgHere )
                eMinuHere  = cmplx ( eMinuRealHere, eMinuImgHere )
                kPerHere   = cmplx ( kPerRealHere, kPerImgHere )

            !   Calculate the interaction time

                bOld  = dlg_interpB ( posK_old, bMagHere = bMagOld )
                bExtraOld  = dlg_interpB ( posK_extra_old, bMagHere = bMagExtraOld)

                distance3D = sqrt( &
                    (1d0*posK(3)-1d0*posK_old(3))**2 &
                    +1d0*posK(1)**2+1d0*posK_old(1)**2 &
                    -2d0*posK(1)*posK_old(1)*cos( (1d0*posK(2)-1d0*posK_old(2)) ) )

                dBdS_   = ( bMagHereK - bMagOld ) / distance3D 
                d2BdS2  = (bMagExtraOld - 2d0 * bMagOld + bMagHereK)/distance3D**2

                if (dBdS_ .ne. dBdS_) then
                     kickErr = 638
                     !kickNaN   = .true.
                     go to 321
                endif 

                if (d2BdS2 .ne. d2BdS2 ) then
                     kickErr = 644
                     !kickNan = .true.
                     go to 321
                endif

                interactionTime = sqrt ( 2d0 * pi / ( q / mi * abs ( vParK * dBdS_ ) ) )

                if ( interactionTime .ne. interactionTime ) then 
                        kickErr = 674
                        kickNaN = .true.
                        go to 321
                endif

                res_accel   = ( vParK - vParK_old ) / dt

                if ( res_accel .ne. res_accel ) then
                        kickNaN = .true.
                        kickErr = 654
                        go to 321
                endif

                c1  = q*dBdS_/mi/2d0*vParK
                c2  = -1d0/6d0*( q*dBdS_/mi*res_accel + d2BdS2 * vParK**2)

                airyArg = -c1**2/(3d0*abs(c2))**(4d0/3d0)

                if ( airyArg .ne. airyArg ) then 
                    kickNaN = .true.
                    kickErr = 650
                    go to 321 
                endif

                call airy_ai ( airyArg, ai, dai, airyErr )
                interactionTimeC = abs ( 2.0 * pi / ( 2d0 * abs ( c2 ) )**(1d0/3d0) * ai )
               
                if ( interactionTimeC .ne. interactionTimeC ) then 
                        kickErr = 674
                        kickNaN = .true.
                        go to 321
                endif

                interactionTimeMin = minVal ( (/ interactionTime, interactionTimeC /) ) 

                if ( interactionTimeMin > 1e-4 ) then

                    write(*,*) 'dlg: - large interaction time - correcting to 1e-4 from ', &
                        interactionTimeMin
                    interactionTimeMin = 1e-4

                endif

            !   Apply kick to vPer
               
                bArg    = real ( kPerHere ) * vPerpK / omegaHere 
                !
                !bFn_m2 = bessJ ( harmonicNumber-2, bArg) 
                bFn_m1 = bessJ ( harmonicNumber-1, bArg) 
                !bFn = bessJ ( harmonicNumber, bArg) 
                bFn_p1 = bessJ ( harmonicNumber+1, bArg) 
                !bFn_p2 = bessJ ( harmonicNumber+2, bArg) 
    
                !D_vPerp =  ( q / ( 2d0 * mi ) * abs (ePlusHere))**2 * interactionTimeMin
                D_vPerp =  ( q / ( 2d0 * mi ) &
                    * ( abs (ePlusHere) * bFn_m1 &
                        + abs (eMinuHere) * bFn_p1 ))**2 * interactionTimeMin

                mean_vPerp  = 0!D_vPerp / vPerpK * interactionTimeMin

                randNo  = rand_lux_sgnE( int(mod(stepCnt,50000)+1) )-0.5
                !randNo  = rand_lux_sgnE( stepCnt )-0.5

                randNo  = randNo / abs ( randNo )

                if ( randNo .ne. randNo ) randNo = 1

                dvPerp  = mean_vPerp &
                    +  randNo * sqrt ( 2d0 *  D_vPerp * interactionTimeMin )

                dvPar   = kPar * vPerpK * dvPerp / ( harmonicNumber * omegaHere )

                vPerpK = vPerpK + dvPerp
                if ( vPerpK .lt. 0 ) vPerpK = abs ( vPerpK )
                vParK  = vParK + dvPar 

                powerAdd    = mi * dvPerp**2 * real(start_weight,kind=dbl) / 2d0
                powerParAdd    = mi * dvPar**2 * real(start_weight,kind=dbl) / 2d0

                if ( power_R_index > 1 .and. power_R_index .le. power_R_nBins &
                        .and. power_z_index > 1 .and. power_z_index .le.  power_z_nBins ) then

                    power(power_R_index,power_z_index) = power(power_R_index,power_z_index) + powerAdd 
                    powerPar(power_R_index,power_z_index) = powerPar(power_R_index,power_z_index) + powerParAdd

                else
                    kickErr = -766 
                    kickNaN = .true.
                    go to 321
                endif

                if ( vPerpK .eq. vPerpK .and. vParK .eq. vParK .and. vPerpK > 0 ) then

                        vPerp  = vPerpK
                        vPar   = vParK

                else
                    kickNan = .true.
                    kickErr = 743
                endif

321         continue

            if ( kickNaN ) then 

                    write(*,*) 'dlg: ERROR - KICK - ', kickErr       
                    write(*,*) '--------------------------------'
                    write(*,*) 'mi: ', mi
                    write(*,*) 'dvPerpK: ', dvPerp
                    write(*,*) 'dvPerpK^2: ', dvPerp**2
                    write(*,*) 'dvPar: ', dvPar
                    write(*,*) 'dvPar^2: ', dvPar**2
                    write(*,*) 'dvParK: ', kPar * vPerpK * dvPerp / ( harmonicNumber * omegaHere )
                    write(*,*) 'start_weight: ', start_weight
                    write(*,*) 'powerAdd: ', powerAdd
                    write(*,*) 'powerParAdd: ', powerParAdd
                    write(*,*) 'interactionTime: ', interactionTimeMin
                    write(*,*) 'bF1: ', bF1
                    write(*,*) 'bF2: ', bF2
                    write(*,*) 'cAbs(ePlusHere): ', cAbs(ePlusHere)
                    write(*,*) 'cAbs(eMinuHere): ', cAbs(eMinuHere)
                    write(*,*) 'vParK: ', vParK
                    write(*,*) 'vParK_old: ', vParK_old
                    write(*,*) 'vPerpK: ', vPerpK
                    write(*,*) 'omegaHere: ', omegaHere
                    write(*,*) 'real(kPerHere): ', real(kPerHere)
                    write(*,*) 'posK: ', posK
                    write(*,*) 'posK_old: ', posK_old
                    write(*,*) 'bArg: ', bArg 
                    write(*,*) 'dBdS_: ', dBdS_
                    write(*,*) 'd2Bds2: ',d2BdS2 
                    write(*,*) 'bMagHereK: ', bMagHereK
                    write(*,*) 'bMagOld: ', bMagOld
                    write(*,*) 'bMagExtraOld: ', bMagExtraOld
                    write(*,*) 'distance3D: ', distance3D
                    write(*,*) 'C: ', interactionTimeC
                    write(*,*) 'B: ', interactionTimeB
                    write(*,*) 'A: ', interactionTime
                    write(*,*) 'ai: ', ai
                    write(*,*) 'airyArg: ', airyArg
                    write(*,*) 'q:', q
                    write(*,*) 'pi: ', pi
                    write(*,*) 'intA sqrt arg: ', 2d0 * pi / ( q / mi * abs ( vParK * dBdS_ ) )
                    write(*,*) 'c1: ', c1
                    write(*,*) 'c2: ', c2
                    write(*,*) 'res_accel: ', res_accel
                    write(*,*) 'stepCnt: ', stepCnt
            
                endif

            endif if_resonant

        endif kick_diffusion

    !
    !   End fixed k diffusion  
    !   +------------------------------------


    !   +------------------------------------
    !   Serial QL diffusion 
    !

        ql_serial_diffusion: &
        if ( use_ql_serial ) then

            !if ( NFO ) then
            !    posK = pos_
            !    posK_old    = pos_old_n
            !    posK_extra_old  = pos_extra_old_n
            !    vPerpK  = vPerp_n
            !    vParK   = vPar_n
            !    vParK_old   = vPar_old_n
            !    vPerpK_old   = vPerp_old_n
            !else
                posK = pos
                posK_old    = pos_old
                posK_extra_old  = pos_extra_old
                vPerpK  = vPerp
                vParK   = vPar
                vParK_old   = vPar_old
                vPerpK_old   = vPerp_old
            !endif 
 
            ql_rIndex   = int ( ( posK(1)- ql_R_min ) / ql_R_range * ql_nR_serial ) + 1
            ql_zIndex   = int ( ( posK(3)- ql_z_min ) / ql_z_range * ql_nz_serial ) + 1
            ql_vPerIndex   = int ( ( vPerpK - ql_vPer_min ) / ql_vPer_range * ql_nvPer_serial ) + 1
            ql_vParIndex   = int ( ( vParK - ql_vPar_min ) / ql_vPar_range * ql_nvPar_serial ) + 1

            if ( ql_vPerIndex < 1 ) ql_vPerIndex = 1
            if ( ql_vPerIndex > ql_nvPer_serial-1 ) ql_vPerIndex = ql_nvPer_serial-1

            if_inQL_grid: &
            if( ql_rIndex > 0 .and. ql_rIndex .le. ql_nR_serial &
                .and. ql_zIndex > 0 .and. ql_zIndex .le. ql_nz_serial &
                .and. ql_vPerIndex > 0 .and. ql_vPerIndex .le. ql_nvPer_serial-1 &
                .and. ql_vParIndex > 0 .and. ql_vParIndex .le. ql_nvPar_serial) then

                ql_bHereCyl = ql_b(ql_rIndex,ql_zIndex,ql_vPerIndex,ql_vParIndex) 
                ql_bHereCylP = ql_b(ql_rIndex,ql_zIndex,ql_vPerIndex+1,ql_vParIndex) 
 
                divD_per1   = ( ql_bHereCylP - ql_bHereCyl ) &
                    / ( ql_vPer(ql_vPerIndex+1) - ql_vPer(ql_vPerIndex) )

                randNo  = rand_lux_sgnE( int(mod(stepCnt,50000)+1) )-0.5
                randNo  = randNo / abs ( randNo )
                if ( randNo .ne. randNo ) randNo = 1

                !divD_per1 = 0

                if ( ql_bHereCyl > 0 ) then 

                    dvPerp  = divD_per1 * dt &
                        + randNo * sqrt ( 2d0 * ql_bHereCyl * dt )

                else
                        
                    dvPerp = 0

                endif

                !   Ensure the particle hits every QL box at least once
                !   by adjusting the time step so the step in 4D phase
                !   space is smaller than the QL resolution

                dtR = dt * ql_R_binSize / abs(posK(1) -posK_old(1))
                dtz = dt * ql_z_binSize / abs(posK(3) -posK_old(3))  
                dtvPer = dt * ql_vPer_binSize / abs(vPerpK - vPerpK_old) 
                dtvPar = dt * ql_vPar_binSize / abs(vParK -vParK_old)

                dt_rfSerial = minVal ( (/ dtR, dtz, dtvPer, dtvPar /) )

                if ( dt_rfSerial < 0 ) read(*,*)

                if ( dvPerp .ne. dvPerp ) then
                        write(*,*) 'THIS SHOUDL NOT HAPPEN'
                        dvPerp = 0
                endif
                !if ( vPerp + dvPerp > 0 ) dvPerp = -1d0 * dvPerp   

                                    
                    vPerpK = vPerpK + dvPerp
                    if ( vPerpK < 0 ) vPerpK = abs (vPerpK)
                    power_R_index = ( posK(1) - power_R_min ) / power_R_range * power_R_nBins + 1
                    power_z_index = ( posK(3) - power_z_min ) / power_z_range * power_z_nBins + 1

                    powerAdd    = mi * dvPerp**2 * real(start_weight,kind=dbl) / 2d0

                    if ( power_R_index > 1 .and. power_R_index .le. power_R_nBins &
                            .and. power_z_index > 1 .and. power_z_index .le.  power_z_nBins ) then

                        power(power_R_index,power_z_index) = power(power_R_index,power_z_index) + powerAdd 

                    endif

                !endif

                !if ( NFO ) then
                !        if (vPerpK .eq. vPerpK .and. vParK .eq. vParK) then
                !    vPerp_n  = vPerpK
                !    vPar_n   = vParK
                !        endif
                !else
                        if (vPerpK .eq. vPerpK .and. vParK .eq. vParK) then
                    vPerp  = vPerpK
                    vPar   = vParK
                        endif
                !endif 
               
            else ! particle off the QL grid

                if ( time_off_QLgrid .ne. time_off_QLgrid ) then 
                    write(*,*) time_off_QLgrid, dt
                    read(*,*)
                else
                    time_off_QLgrid   = time_off_QLgrid + dt
                endif

            endif if_inQL_grid
            
        endif ql_serial_diffusion
 
    !
    !   End serial QL diffusion  
    !   +------------------------------------


    !   +------------------------------------
    !   Pitch angle scattering operator 
    !

            pitch_scattering: & 
            if ( use_pitchScattering ) then

                  collErr  = -8888 
                  collNaN   = .false. 

            !   Background temperature and density profiles

            E_bgnd  = (E_bgnd_lim &
                + ( E_bgnd_0 - E_bgnd_lim ) &
                * ( 1d0 - rho_here**E_bgnd_beta)**E_bgnd_alpha)
            density_bgnd  = (density_bgnd_lim &
                + ( density_bgnd_0 - density_bgnd_lim ) &
                * ( 1d0 - rho_here**density_bgnd_beta)**density_bgnd_alpha)

            !   Also, the density of the simulating species is also required.
            !   However, I am avoiding a self consistent calculation of this as
            !   it will be slow. Fortunately the density here is only used in
            !   the collision frequency and the density profile does not change
            !   too much, if at all really. So this should be OK for the time
            !   being.

            density  = (density_lim &
                + ( density_0 - density_lim ) &
                * ( 1d0 - rho_here**density_beta)**density_alpha)

            T_bgnd  = 2.0 / 3.0 * E_bgnd

            if ( vPar**2 + vPerp**2 < 0 ) then 
                collNaN = .true.
                collErr = -465
                go to 765 
            end if 

            if ( vPar .ne. vPar .or. vPerp .ne. vPerp ) then 
                collNaN = .true.
                collErr = -474
                go to 765
            end if
 
            v0  = sqrt ( vPar**2 + vPerp**2 )
            vTh = sqrt ( 2d0 * T_bgnd * e_ / mi )

            T   = mi * v0**2 / ( 3.0 * e_ )

                coulombLog  = 23.0 &
                    - log ( &
                    atomicZ * atomicZ_bgnd * ( amu + amu_bgnd ) / ( amu * T_bgnd + amu_bgnd * T ) &
                    * sqrt ( density / 1d6 * atomicZ**2 / T &
                                + density_bgnd / 1d6 * atomicZ_bgnd**2 / T_bgnd ) )

                if ( coulombLog .ne. coulombLog ) then 
                    collNaN = .true.
                    collErr = -498
                    go to 765
                end if

                vB  = coulombLog / 10d0 / 3d6 * sqrt ( 2d0 / amu ) * density_bgnd / 1d6 &
                    / T_bgnd**1.5

                if ( vB .ne. vB .or. vB < 0 ) then 
                    collNaN = .true.
                    collErr = -510
                    go to 765
                end if

                xVar    = v0 / vTh

                phi_x   = erf_external_fn ( xVar )
                psi_x   = ( -2d0 * exp ( -xVar**2 ) * xVar / sqrt ( pi ) + erf_external_fn ( xVar ) ) &
                    / ( 2d0 * xVar**2 )

                if ( psi_x < 0 ) then
                    collNan = .true.
                    collErr = -893
                    go to 765
                endif

                !   Get the dvE/dE term

                energy0 = 3.0 / 2.0 * T 
                dE = energy0 * 0.1

                if ( 2.0 * (energy0-dE) * e_ / mi < 0 ) then 
                    collNaN = .true.
                    collErr = -542
                    go to 765
                end if

                if ( 2.0 * (energy0-dE) * e_ / mi < 0 ) then 
                    collNaN = .true.
                    collErr = -549
                    go to 765
                end if

                v0_n    = sqrt ( 2.0 * (energy0-dE) * e_ / mi )
                v0_p    = sqrt ( 2.0 * (energy0+dE) * e_ / mi )

                xVar_n  = ( v0_n ) / vTh 
                xVar_p  = ( v0_p ) / vTh 

                psi_x_n   = ( -2d0 * exp ( -xVar_n**2 ) * xVar_n &
                    / sqrt ( pi ) + erf_external_fn ( xVar_n ) ) &
                    / ( 2d0 * xVar_n**2 )
                psi_x_p   = ( -2d0 * exp ( -xVar_p**2 ) * xVar_p &
                    / sqrt ( pi ) + erf_external_fn ( xVar_p ) ) &
                    / ( 2d0 * xVar_p**2 )

                vE_n  = 3d0 * sqrt ( pi / 2d0 ) * vB * psi_x_n / xVar_n
                vE_p  = 3d0 * sqrt ( pi / 2d0 ) * vB * psi_x_p / xVar_p
                !vE_n = Cf * psi_x_n * xVar**2 / (v0_n*1d3)**3
                !vE_p = Cf * psi_x_n * xVar**2 / (v0_p*1d3)**3

                if ( vE_n < 0 .or. vE_p < 0 ) then
                        collNan = .true.
                        collErr = -940
                        go to 765
                endif


                dvE_dE  = ( vE_p - vE_n ) / ( 2d0 *dE )  
                dLnvE_dLnE  =  ( log(vE_p) - log(vE_n) ) / ( log(energy0+dE) - log(energy0-dE) ) 

                if ( dvE_dE .ne. dvE_dE ) then 
                    collNaN = .true.
                    collErr = -957
                    go to 765
                end if
 

                vD  = 3d0 / 2d0 * sqrt ( pi / 2d0 ) * vB &
                    * ( phi_x - psi_x ) / xVar**3 
                vE  = 3d0 * sqrt ( pi / 2d0 ) * vB * psi_x / xVar

                Cf  = 8.0 * pi * density_bgnd / 1d6 * atomicZ_bgnd * atomicZ * &
                        e_cgs**4 * coulombLog / (mi*1d3)**2
                vE_ = Cf * psi_x * xVar**2 / (v0*1d3)**3*1d3
                vD_ = Cf * (phi_x-psi_x) / (2.0*(v0*1d3)**3)*1d3

                call get_coll_freqs ( amu, amu_bgnd, T, T_bgnd, atomicZ, &
                    atomicZ_bgnd, v0, density, density_bgnd, &     
                    nu_D_ii, nu_S_ii, nu_par_ii, dnuPardE_ii )

                call get_coll_freqs ( amu, amu_bgnd, T, T_bgnd, atomicZ, &
                    atomicZ_bgnd, v0, density, density_bgnd, &     
                    nu_D_ei, nu_S_ei, nu_par_ei, dnuPardE_ei, electrons = .true. )


                nu_D    = nu_D_ii + nu_D_ei
                nu_S    = nu_S_ii + nu_S_ei
                nu_par  = nu_par_ii + nu_par_ei
                dnuPardE    = dnuPardE_ii + dnuPardE_ei

                !if ( mpi_pId .eq. 0 ) write(*,*) vE_, vE, vD_, vD, T * 2.0 / 3.0

                if ( vD .ne. vD ) then 
                    collNaN = .true.
                    collErr = -527
                    go to 765
                end if

                if ( vE < 0 ) then
                    collNaN = .true.
                    collErr = -994
                    go to 765    
                endif

                rand_sgn    = (rand_lux_sgn(int(mod(stepCnt,50000)+1))-0.5)*2.0
                rand_sgn    = rand_sgn / abs ( rand_sgn )
                rand_sgnE    = (rand_lux_sgnE(int(mod(stepCnt,50000)+1))-0.5)*2.0
                rand_sgnE    = rand_sgnE / abs ( rand_sgnE )

                !   Catch when a 0 random number is produced

                if ( rand_sgn .ne. rand_sgn ) rand_sgn = 1.0
                if ( rand_sgnE .ne. rand_sgnE ) rand_sgnE = 1.0

                rand_sgnE_total = rand_sgnE_total + rand_sgnE

                pitch0  = vPar / v0

                if ( ( 1d0 - pitch0**2 ) * vD * dt < 0 ) then 
                    collNaN = .true.
                    collErr = -783
                    go to 765
                end if

                pitch1  = pitch0 * ( 1d0 - nu_D * dt ) &
                    + rand_sgn * sqrt ( ( 1d0 - pitch0**2 ) * nu_D * dt )

                if ( pitch1 > 1 .or. pitch1 < -1 ) pitch1 = pitch0

                !pitch1  = pitch0 * ( 1d0 - vD * dt ) &
                !    + rand_sgn * sqrt ( ( 1d0 - pitch0**2 ) * vD * dt )

                !energy1 = energy0 - ( 2.0 * nu_S_ii * dt ) * &
                !    ( energy0 - ( 3.0 / 2.0 + energy0 / nu_S_ii * dnuPardE_ii ) &
                !        * T_bgnd ) &
                !        + rand_sgnE * 2d0 * sqrt ( T_bgnd * energy0 * ( nu_S_ii * dt ) )  

                energy1 = energy0 - ( 2.0 * nu_S * dt ) * &
                    ( energy0 - ( 3.0 / 2.0 + energy0 / vE * dvE_dE ) &
                        * T_bgnd ) &
                        + rand_sgnE * 2d0 * sqrt ( T_bgnd * energy0 * ( nu_S * dt ) )  

                dE_1    = - ( 2.0 * vE * dt ) * &
                    ( energy0 - ( 3.0 / 2.0 + energy0 / vE * dvE_dE ) &
                        * T_bgnd )

                dE_ii   = -2d0 * energy0 * dt &
                    * ( mi / ( mi + mi_bgnd ) * nu_S_ii - 5d0/2d0 * nu_par_ii &
                            - energy0 * dnuPardE_ii )
                dE_ei   = -2d0 * energy0 * dt &
                    * ( mi / ( mi + electronMass ) * nu_S_ei - 5d0/2d0 * nu_par_ei &
                            - energy0 * dnuPardE_ei )
                    
                !energy1 = energy0 + ( dE_ii ) & 
                !        + rand_sgnE * 2d0 * energy0 * sqrt ( nu_par_ii * dt )  
                !if ( xVar <0.2 ) then 
                !if ( mpi_pId == 0 ) then 
                !    write(*,*) '------------------'
               !endif


                if ( rand_sgnE .ne. rand_sgnE .or. rand_sgn .ne. rand_sgn ) then
                    collNan = .true.
                    collErr = -799
                    go to 765
                end if

                if (  2.0 * energy1 * e_ / mi < 0 ) then
                    collNaN = .true.
                    collErr = -807
                    go to 765
                end if
                
                v0_update   = sqrt ( 2.0 * energy1 * e_ / mi )
                vPar_update = pitch1 * v0_update 

                if ( v0_update**2 - vPar_update**2 < 0 ) then
                    collNaN = .true. 
                    collErr = -818
                    go to 765
                end if

                vPerp_update    = sqrt ( v0_update**2 - vPar_update**2 )

                if ( vPar_update .ne. vPar_update &
                        .or. vPerp_update .ne. vPerp_update ) then
                    collNaN = .true. 
                    collErr = -1278
                    go to 765
                end if

                !if ( mpi_pId .eq. 0 ) &
                !    write(*,*) energy0, energy1, vE, dvE_dE, energy0 / vE * dvE_dE, stepCnt,tau

                !if ( vD * dt  > 0.1 ) write(*,*) 'WARNING: vD * dt ~ 1', vD * dt, &
                !    0.5 * mi * v0**2 / e_ / 1d3

765             continue

                if ( collNan ) then

                    !write(*,*) 'dlg: ERROR - coll - ', collErr
                    !write(*,*) '--------------------------------'
                    !write(*,*) 'lVar: ', lVar
                    !write(*,*) 'Cf: ', Cf
                    !write(*,*) 'lVar*v: ', lVarV
                    !write(*,*) 'xVar: ', xVar
                    !write(*,*) 'G(lVar*v): ', G_lVarV
                    !write(*,*) 'v: ', v0
                    !write(*,*) 'vTh: ', vTh
                    !write(*,*) 'kT: ', k_ * T_bgnd
                    !write(*,*) 'vB: ', vB
                    !write(*,*) 'vE: ', vE
                    !write(*,*) 'vE_: ', vE_
                    !write(*,*) 'vD: ', vD
                    !write(*,*) 'vD_: ', vD_ 
                    !write(*,*) 'Cf / v^3: ', Cf / (v0*1d3)**3
                    !write(*,'(a15,3x,e9.2)') 'nu_D_ii: ', nu_D_ii
                    !write(*,'(a15,3x,e9.2)') 'nu_S_ii: ', nu_S_ii
                    !write(*,'(a15,3x,e9.2)') 'nu_par_ii: ', nu_par_ii
                    !write(*,'(a15,3x,e9.2)') 'nu_D_ei: ', nu_D_ei
                    !write(*,'(a15,3x,e9.2)') 'nu_S_ei: ', nu_S_ei
                    !write(*,'(a15,3x,e9.2)') 'nu_par_ei: ', nu_par_ei
                    !write(*,'(a15,3x,e9.2)') 'dnuPardE_ii: ', dnuPardE_ii
                    !write(*,'(a15,3x,e9.2)') 'dnuPardE_ei: ', dnuPardE_ei
                    !write(*,'(a15,3x,e9.2)') 'dvE_dE: ', dvE_dE 
                    !write(*,'(a15,3x,e9.2)') 'dE_ii: ', dE_ii
                    !write(*,'(a15,3x,e9.2)') 'dE_ei: ', dE_ei 
                    !write(*,'(a15,3x,e9.2)') 'dE_1: ', dE_1
                    !write(*,'(a15,3x,e9.2)') 'E/vE *dvE_dE: ', energy0 / vE * dvE_dE 
                    !write(*,'(a15,3x,e9.2)') 'dLnvE_dLnE: ', dLnvE_dLnE 

                    time_bad_coll   = time_bad_coll + dt

                else

                    vPar    = vPar_update
                    vPerp   = vPerp_update

                endif

            endif pitch_scattering

    !
    !   End pitch angle scattering  
    !   +------------------------------------


    !   +------------------------------------
    !   QL operator stuff
    !
            ql_diffusion: &    
            if ( use_QL ) then

                ql_rIndex   = int ( ( pos(1)- ql_R_min ) / ql_R_range * ql_nR ) + 1
                ql_zIndex   = int ( ( pos(3)- ql_z_min ) / ql_z_range * ql_nz ) + 1

                !if ( vPar .ge. ql_vMin(ql_rIndex,ql_zIndex) .and. &
                !    vPar .le. ql_vMax(ql_rIndex,ql_zIndex) ) then

                    ql_vPerIndex   = int ( ( vPerp - ql_vPer_min ) / ql_vPer_range * ql_nvPer ) + 1
                    ql_vParIndex   = int ( ( vPar - ql_vPar_min ) / ql_vPar_range * ql_nvPar ) + 1
                
                    !   Which processor has the data?

                    mpi_coords_hasData  = & ! mpi indices start at 0 remember
                        (/ (ql_rIndex-1) / chunk_nR, (ql_zIndex-1) / chunk_nz /) 

                    call mpi_cart_rank ( mpi_comm_cart, mpi_coords_hasData, mpi_pId_hasData, mpi_iErr )

                    call dlg_sendData ()

                    if ( mpi_pId_hasData .eq. mpi_pId ) then

                        write(*,*) 'Had data:', ql_rIndex, ql_zIndex, &
                            ql_vPerIndex, ql_vParIndex, mpi_pId_hasData, &
                            ql_b(ql_rIndex,ql_zIndex,ql_vPerIndex,ql_vParIndex), &
                            ql_c(ql_rIndex,ql_zIndex,ql_vPerIndex,ql_vParIndex), &
                            ql_e(ql_rIndex,ql_zIndex,ql_vPerIndex,ql_vParIndex), &
                            ql_f(ql_rIndex,ql_zIndex,ql_vPerIndex,ql_vParIndex)

                    else

                        call mpi_send ( (/ mpi_pId,ql_rIndex,ql_zIndex,ql_vPerIndex,ql_vParIndex /), &
                            5, MPI_INTEGER, &
                            mpi_pId_hasData, mpi_tag_wantsData, mpi_comm_cart, mpi_iErr)

                        write(*,*) mpi_pId, 'is wating on data from ', mpi_pId_hasData
                        call mpi_irecv ( ql_recv, 4, MPI_REAL, &
                            mpi_pId_hasData, mpi_tag_getData, mpi_comm_cart, mpi_req_getData, mpi_iErr ) 

                        do !    Deal with any data requests while waiting 

                            call mpi_test ( mpi_req_getData, mpi_flag_getData, &
                                mpi_stat_getData, mpi_iErr )

                            if ( mpi_flag_getData ) exit

                            call dlg_sendData ()
                        end do
 
                        write(*,*) 'Received data:', ql_rIndex, ql_zIndex, &
                            ql_vPerIndex, ql_vParIndex, mpi_pId_hasData, &
                            ql_recv

                    end if
 
                !end if

            end if ql_diffusion

    !
    !   End QL operator
    !   +------------------------------------

            
            !   Update to new COM magnetic moment

            u   = mi * vPerp**2 / ( 2.0 * bMagHere ) 
            !if (NFO) u_n   = mi * vPerp_n**2 / ( 2.0 * bMagHere_n ) 
          
            !   Update the time step
      
            if ( EStep < 5.0e-7 ) dt = dt * 1.2
            if ( EStep > 2.5e-6 ) dt = dt * 0.8

            if ( use_pitchScattering ) then
                if ( nu_D * dt > 0.01 ) dt = dt * 0.5
                if ( nu_par * dt > 0.01 ) dt = dt * 0.5
            endif

            if ( distance > 0.005 ) dt = dt * 0.005 / distance

            if ( dt_rfSerial < dt .and. use_ql_serial .and. stepCnt > 2 ) &
                dt = dt_rfSerial

            if ( dt < dtMin .or. dt .ne. dt .or. dt*0 .ne. 0 ) then
                dt = dtMin
            endif
            if ( dt > dtMax ) dt = dtMax
 
#if USE_DISLIN
           
            !   Plot track

            if ( present ( plot ) .and. mpi_pId .eq. 0 ) then
                if ( plot ) then  

                    rTrack(stepCnt) = pos(1)
                    zTrack(stepCnt) = pos(3)
                    vPerpTrack(stepCnt) = vPerp
                    vParTrack(stepCnt)  = vPar
        
                    !if ( NFO ) then 
                    !    rTrack_nfo(stepCnt) = pos_(1)
                    !    zTrack_nfo(stepCnt) = pos_(3)
                    !    vPerpTrack_nfo(stepCnt) = vPerp_n
                    !    vParTrack_nfo(stepCnt)  = vPar_n
                    !endif

                    if ( stepCnt > 2 .and. mpi_pId == 0 ) then

                        call axsPos ( 300, 1700 )
                        call axsLen ( 1500, 1400 )
                        call color ( 'WHITE' )
                        call graf ( 0.5, 1.0, 0.5, 0.25, -0.25, 0.25, -0.25, 0.25 ) 
                        call color ( 'RED' )
                        call linWid ( 5 )
                        call curve ( rTrack(1:stepCnt-1), zTrack(1:stepCnt-1), stepCnt-1 )
                        !if ( NFO ) then 
                        !    call color ( 'GREEN' )
                        !    call curve ( rTrack_nfo(1:stepCnt-1), zTrack_nfo(1:stepCnt-1), stepCnt-1 )
                        !endif
                        call linWid ( 1 )
                        call endGrf ()
                        call axsPos ( 300, 2600 )
                        call axsLen ( 1500, 700 )
                        call graf ( start_R-0.01, start_R+0.01, start_R-0.01,&
                            0.005, start_z-0.01, start_z+0.01,start_z-0.01, 0.005 ) 
                        call color ( 'RED' )
                        call curve ( rTrack(1:stepCnt-1), zTrack(1:stepCnt-1), stepCnt-1 )
                        call endGrf ()

                    end if

                endif
            end if
#endif
 
            stepCnt = stepCnt + 1
 
            if ( vPerp .ne. vPerp .or. vPar .ne. vPar &
            .or. abs(pos(1)) .gt. 1e4 .or. abs(pos(3)) .gt. 1e4) then 
                nanCatch = .true.
                nanErrorMsg = -1473
                go to 72
            end if

            !if (NFO) then
            !if ( vPerp_n .ne. vPerp_n .or. vPar_n .ne. vPar_n &
            !.or. abs(pos_(1)) .gt. 1e4 .or. abs(pos_(3)) .gt. 1e4) then 
            !    nanCatch = .true.
            !    nanErrorMsg = -1479
            !    go to 72
            !end if
            !endif

        enddo main_loop


        if ( vPerp .ne. vPerp .or. vPar .ne. vPar) then 
            nanCatch    = .true.
            nanErrorMsg = -1284
            status_out  = nanErrorMsg
        end if

        !if ( NFO ) then 

        !    final_R = pos_(1)
        !    final_z = pos_(3)
        !    final_vPerp = vPerp_n
        !    final_vPar  = vPar_n        

        !else

            final_R = pos(1)
            final_z = pos(3)
            final_vPerp = vPerp
            final_vPar  = vPar        

        !endif

643     continue 

        if ( tau > 0 ) then
            nP_off_QLgrid   = nP_off_QLgrid + time_off_QLgrid / tau
            nP_bad_coll     = nP_bad_coll + time_bad_coll / tau
        endif

        if ( nanCatch .or. start_status < 0 ) then

            final_R = start_R
            final_z = start_z 
            final_vPerp = start_vPerp
            final_vPar  = start_vPar        
            weight_out  = start_weight

            !write(*,*) 'dlg: ERROR - loop nan - ', nanErrorMsg

        endif
    
   end subroutine gc_orbit 

end module gc_integrate
