module read_particle_list
    implicit none
    save

    real, allocatable :: p_psi(:), p_theta(:), &
        p_R(:), p_z(:), p_E(:), p_lambda(:), p_weight(:), &
        p_vPerp(:), p_vPar(:), p_R_final(:), &
        p_z_final(:), p_vPerp_final(:), p_vPar_final(:), &
        p_weight_final(:) 
    integer, allocatable :: p_status(:), p_status_final(:)
   
contains    
    subroutine read_pl ()
        use constants
        use init_mpi
        use mpi
        use read_namelist
        use netcdf
        use dlg
        use powerAbsGrid

        integer(kind=LONG) :: i
        real, dimension(nP) :: p_v
        integer :: nc_nP, nc_id, R_id, z_id, vPer_id, vPar_id, weight_id, R_dim_ids(1), &
            power_id, powerPar_id, power_dim_ids(2), nc_power_nR, nc_power_nz, &
            status_id
        real, allocatable :: nc_R(:), nc_z(:), nc_vPer(:), nc_vPar(:), nc_weight(:), &
            nc_power(:,:), nc_powerPar(:,:)
        integer, allocatable :: nc_status(:)

        if ( .not. allocated ( p_R ) ) &
        allocate ( p_psi(nP), p_theta(nP), &
            p_R(nP), p_z(nP), p_E(nP), p_lambda(nP), p_weight(nP), &
            p_vPerp(nP), p_vPar(nP), p_R_final(nP), p_z_final(nP), &
            p_vPerp_final(nP), p_vPar_final(nP), p_weight_final(nP), &
            p_status(nP), p_status_final(nP) )

        p_R_final   = 0.0
        p_z_final   = 0.0
        p_vPerp_final   = 0.0
        p_vPar_final    = 0.0
        p_weight_final  = 0.0
        p_status_final  = 0

        p_vPerp   = 0.0
        p_vPar    = 0.0
        p_weight  = 0.0
        p_status    = 0
        p_R = 0
        p_z = 0

        p_psi   = 0.0
        p_theta   = 0.0
        p_E = 0
        p_lambda = 0

        powerPrev   = 0
        powerParPrev    = 0

        !   Read in variables from .dav particle list file

        if ( pList_is_nCDF ) then 

            !   Read in variables from .nc particle list file

            call dlg_check ( nf90_open ( path = pl_fileName, mode = nf90_nowrite, ncid = nc_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'R', R_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'z', z_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'vPer', vPer_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'vPar', vPar_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'weight', weight_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'status', status_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'power', power_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'powerPar', powerPar_id ) )
         
            call dlg_check ( nf90_inquire_variable ( &
                nc_id, R_id, dimIds = R_dim_ids ) )
 
            call dlg_check ( nf90_inquire_variable ( &
                nc_id, power_id, dimIds = power_dim_ids ) )
            
            call dlg_check ( nf90_inquire_dimension ( nc_id, R_dim_ids(1), &
                len = nc_nP ) ) 
            call dlg_check ( nf90_inquire_dimension ( nc_id, power_dim_ids(1), &
                len = nc_power_nR ) ) 
            call dlg_check ( nf90_inquire_dimension ( nc_id, power_dim_ids(2), &
                len = nc_power_nz ) ) 
          
            allocate ( nc_R ( nc_nP ), &
                nc_z ( nc_nP ), &
                nc_vPer ( nc_nP ), &
                nc_vPar ( nc_nP ), &
                nc_weight ( nc_nP ), &
                nc_power ( nc_power_nR, nc_power_nz ), &
                nc_powerPar ( nc_power_nR, nc_power_nz ), &
                nc_status ( nc_nP ) )

            call dlg_check ( nf90_get_var ( nc_id, R_id, nc_R ) )
            call dlg_check ( nf90_get_var ( nc_id, z_id, nc_z ) )
            call dlg_check ( nf90_get_var ( nc_id, vPer_id, nc_vPer ) )
            call dlg_check ( nf90_get_var ( nc_id, vPar_id, nc_vPar ) )
            call dlg_check ( nf90_get_var ( nc_id, weight_id, nc_weight ) )
            call dlg_check ( nf90_get_var ( nc_id, status_id, nc_status ) )
            call dlg_check ( nf90_get_var ( nc_id, power_id, nc_power ) )
            call dlg_check ( nf90_get_var ( nc_id, powerPar_id, nc_powerPar ) )
           
            call dlg_check ( nf90_close ( nc_id ) )

            p_R = nc_R(1:nP)
            p_z = nc_z(1:nP)
            p_vPerp = nc_vPer(1:nP)
            p_vPar  = nc_vPar(1:nP)
            p_weight    = nc_weight(1:nP)
            p_status    = nc_status(1:nP)
            powerPrev   = nc_power
            powerParPrev    = nc_powerPar
 
            deallocate ( nc_R, &
                nc_z, &
                nc_vPer, &
                nc_vPar, &
                nc_weight, &
                nc_power, &
                nc_powerPar, &
                nc_status )

        else 
        
            open ( unit = 8, file = pl_fileName, status = 'OLD', action = 'READ' )

            do i=1,nP
               
                read (8,2000) p_psi(i), p_theta(i), &
                    p_R(i), p_z(i), p_E(i), p_lambda(i), p_weight(i)

            end do 

            2000 format ( 7e16.6 )

            close ( unit = 8 )

            !   Convert R,z to [m] from [cm]

            p_R = p_R * 1e-2
            p_z = p_z * 1e-2

            !   Calculate particle vPerp and vPar

            p_v = sqrt ( 2.0 * p_E * 1.0e3 * e_ / mi ) 
            p_vPar  = p_v * p_lambda
            p_vPerp = sqrt ( p_v**2 - p_vPar**2 )
            p_status = 0

        endif

    end subroutine read_pl


    subroutine erase_pl ()

        use constants
        use init_mpi
        use mpi
        use read_namelist

        implicit none
        deallocate ( p_psi, p_theta, &
            p_R, p_z, p_E, p_lambda, p_weight, &
            p_vPerp, p_vPar, p_R_final, p_z_final, &
            p_vPerp_final, p_vPar_final, p_weight_final )

    end subroutine erase_pl

end module read_particle_list
