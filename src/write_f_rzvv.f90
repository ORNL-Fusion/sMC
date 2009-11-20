module write_f_rzvv
    
contains
    subroutine check ( status )
        use netcdf
        integer, intent ( in) :: status
      
        if(status /= nf90_noerr) then 
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if

    end subroutine check 
    
    subroutine write_f ( n, ncFileName )
        !use rzvv_grid
        use netcdf
        use init_mpi
        use read_particle_list
        use read_namelist
        use powerAbsGrid

        implicit none

        integer, intent(in) :: n
    
        character ( len = 103 ), intent(in) :: ncFileName
        integer :: nc_id, &
            nP_dim_id, &
            scalar_id, &
            R_id, &
            z_id, &
            vPer_id, &
            vPar_id, &
            nP_id, &
            weight_id, &
            power_R_nBins_id, &
            power_z_nBins_id, &
            power_R_id, &
            power_z_id, &
            power_id, &
            powerPar_id, &
            time_id, &
            status_id

        call check ( nf90_create ( ncFileName, nf90_clobber, nc_id ) )
        call check ( nf90_def_dim ( nc_id, "nP", nP, nP_dim_id ) )
        call check ( nf90_def_dim ( nc_id, "scalar", 1, scalar_id ) )
        call check ( nf90_def_dim ( nc_id, "power_R_nBins", power_R_nBins, power_R_nBins_id ) )
        call check ( nf90_def_dim ( nc_id, "power_z_nBins", power_z_nBins, power_z_nBins_id ) )

        call check ( nf90_def_var ( nc_id, "R", nf90_float, &
            (/ nP_dim_id /), R_id ) )
        call check ( nf90_def_var ( nc_id, "z", nf90_float, &
            (/ nP_dim_id /), z_id ) )
        call check ( nf90_def_var ( nc_id, "vPer", nf90_float, &
            (/ nP_dim_id /), vPer_id ) )
        call check ( nf90_def_var ( nc_id, "vPar", nf90_float, &
            (/ nP_dim_id /), vPar_id ) )
        call check ( nf90_def_var ( nc_id, "weight", nf90_float, &
            (/ nP_dim_id /), weight_id ) )
        call check ( nf90_def_var ( nc_id, "status", NF90_INT, &
            (/ nP_dim_id /), status_id ) )

        call check ( nf90_def_var ( nc_id, "power", nf90_float, &
            (/ power_R_nBins_id, power_z_nBins_id /), power_id ) )
        call check ( nf90_def_var ( nc_id, "powerPar", nf90_float, &
            (/ power_R_nBins_id, power_z_nBins_id /), powerPar_id ) )
 
        call check ( nf90_def_var ( nc_id, "power_R_binCenters", nf90_float, &
            (/ power_R_nBins_id /), power_R_id ) )
        call check ( nf90_def_var ( nc_id, "power_z_binCenters", nf90_float, &
            (/ power_z_nBins_id /), power_z_id ) )
        call check ( nf90_def_var ( nc_id, "time", nf90_float, &
            scalar_id, time_id ) )
   

        call check ( nf90_def_var ( nc_id, "nP", NF90_INT, &
            (/ scalar_id /), nP_id ) )
     
        call check ( nf90_enddef ( nc_id ) )

        call check ( nf90_put_var ( nc_id, nP_id, nP ) )
        call check ( nf90_put_var ( nc_id, R_id, p_R_final ) )
        call check ( nf90_put_var ( nc_id, z_id, p_z_final ) )
        call check ( nf90_put_var ( nc_id, vPer_id, p_vPerp_final ) )
        call check ( nf90_put_var ( nc_id, vPar_id, p_vPar_final ) )
        call check ( nf90_put_var ( nc_id, weight_id, p_weight_final ) )
        call check ( nf90_put_var ( nc_id, power_id, power+powerPrev ) )
        call check ( nf90_put_var ( nc_id, powerPar_id, powerPar+powerParPrev ) )
        call check ( nf90_put_var ( nc_id, power_R_id, power_R_binCenters ) )
        call check ( nf90_put_var ( nc_id, power_z_id, power_z_binCenters ) )
        call check ( nf90_put_var ( nc_id, time_id, n * runTime ) )
        call check ( nf90_put_var ( nc_id, status_id, p_status_final ) )

        call check ( nf90_close ( nc_id ) )

    end subroutine write_f

end module write_f_rzvv
