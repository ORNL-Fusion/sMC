program sMCrf
    use eqdsk
    use gc_terms
    use gc_integrate 
    use read_particle_list
    !use rzvv_grid
    use init_mpi
    use constants
    !use plot
    use write_f_rzvv
    use read_namelist
    use read_ql
    use communications
    use read_mchoi
    use powerAbsGrid
    use interp

    implicit none
    
    integer(kind=LONG) :: i 
    real :: T1, T2
    integer :: mpi_count, n
    integer :: nP_wall_total, nP_bad_total, nP_off_vGrid_total
    integer :: nP_lowV
    real :: nP_off_QLgrid_total, nP_bad_coll_total
    integer :: tmp_done
    real, allocatable :: power_total(:,:), powerPar_total(:,:)
    character ( len = 103 ) :: ncFileName

    call init_namelist ()
    call start_mpi ()
    if ( mpi_pId .eq. 0 ) write(*,*) 'Load constants ...' 
    call set_constants ()
    if ( mpi_pId .eq. 0 ) write(*,*) 'Load eqdsk ...' 
    call read_geqdsk ( eqdsk_fileName, plot = .false. )
    if ( mpi_pId .eq. 0 ) write(*,*) 'Calculate curvature terms ...' 
    call bCurvature ()
    if ( mpi_pId .eq. 0 ) write(*,*) 'Calculate gradient terms ...' 
    call bGradient ()
    if ( mpi_pId .eq. 0 ) write(*,*) 'Init interpolation arrays ...' 
    call init_interp ()
    if ( mpi_pId .eq. 0 ) write(*,*) 'Init power abs grid  ...' 
    call init_powerAbs_grid ()

    call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )

    if ( mpi_pId .eq. 0 ) write(*,*) 'Load AORSA fields ...' 
    if ( use_QL ) call open_ql ()
    if ( use_kick_diffusion ) call read_aorsa_fields ()
    if ( use_QL_serial ) then
            if ( mpi_pId .eq. 0 ) write(*,*) &
                'Allocating memory and reading serial QL ...'
            call open_ql_serial ()
            if ( mpi_pId .eq. 0 ) write(*,*) &
                'DONE'
    endif
 
    allocate ( power_total(power_R_nBins,power_z_nBins), &
        powerPar_total(power_R_nBins,power_z_nBins) )

    power_total  = 0
    powerPar_total = 0

    !call init_rzvv_grid ()

    call cpu_time (T1)

    allocate ( nP_done ( mpi_nProcs ) )

    call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )

    do n = startn, startn+nSteps

        call read_pl ()    
        call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )

        power   = 0
        powerPar    = 0

        nP_done = 0

        if ( n .gt. 0 ) then 

            !   Initialize the data request handles

            if ( use_QL ) then
            if ( n .eq. 1 ) then  

                call mpi_irecv ( mpi_pId_and_indices, 5, MPI_INTEGER, &
                    MPI_ANY_SOURCE, mpi_tag_wantsData, mpi_comm_cart, &
                    mpi_req_wantsData, mpi_iErr)
     
                call mpi_irecv ( mpi_finished_pId, 1, MPI_INTEGER, &
                    MPI_ANY_SOURCE, mpi_tag_finished, mpi_comm_cart, &
                    mpi_req_finished, mpi_iErr)
    
            else

                call dlg_sendData ()
                call dlg_check_finished ()
                nP_done = 0

            end if
            end if
 
            do i = mpi_start_, mpi_end_

                if ( mod(i,100) .eq. 0 .and. mpi_pId .eq. 0 ) &
                    write(*,*) mpi_pId, n, nSteps, i-mpi_start_, mpi_end_ - mpi_start_

                call gc_orbit ( p_R(i), p_z(i), p_vPerp(i),&
                    p_vPar(i), p_weight(i), p_status(i), p_R_final(i), p_z_final(i), &
                    p_vPerp_final(i), p_vPar_final(i), p_weight_final(i), &
                    p_status_final(i), plot = plotOrbit )

                if ( use_QL ) call dlg_check_finished ()
 
            end do

            nP_done ( mpi_pId + 1 ) = 1
            write(*,*) mpi_pId, 'Finished my work, waiting for others ...' 

            !   Tell all the other processors I'm done

            if ( use_QL ) then 
            do i = 0, mpi_nProcs -1 

                write (*,*) 'Telling everyone im done', i

                call dlg_check_finished ()

                call dlg_wakeUp (i)

                if ( mpi_pId .ne. i ) call mpi_send ( mpi_pId, &
                    1, MPI_INTEGER, &
                    i, mpi_tag_finished, mpi_comm_cart, mpi_iErr)

            end do

            write(*,*) 'mpi_pID: ', mpi_pId, nP_done

            do  !    While waiting deal with send requests 

                tmp_done    = sum ( nP_done )

                call dlg_check_finished ()
                
                if ( sum ( nP_done ) > tmp_done ) &
                    write(*,*) 'mpi_pID: ', mpi_pId, nP_done
 
                if ( sum ( nP_done ) .eq. mpi_nProcs ) exit

                !   If not then send data while waiting
                call dlg_sendData_wait ()

            end do

            !   Wake up any left over waiting procs
            do i = 0, mpi_nProcs -1 

                call dlg_wakeUp (i)

            end do
            end if 


        else

            do i = mpi_start_, mpi_end_

                p_R_final(i) = p_R(i)
                p_z_final(i) = p_z(i)
                p_vPerp_final(i) = p_vPerp(i)
                p_vPar_final(i) = p_vPar(i)
                p_weight_final(i) = p_weight(i) 
                p_status_final(i) = p_status(i) 

            end do

        end if

        if ( mpi_pId == 0 ) write(*,*) 'Waiting for all processors to finish ...'
        call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )

        mpi_count   = mpi_end_ - mpi_start_ + 1

        !   Gather the new particle locations from all processors.

        if ( mpi_pId == 0 ) write(*,*) 'Gathering data for file write'
        if ( mpi_pId == 0 ) then 

            call mpi_gather ( MPI_IN_PLACE, mpi_count, MPI_REAL, &
                p_R_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( MPI_IN_PLACE, mpi_count, MPI_REAL, &
                p_z_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( MPI_IN_PLACE, mpi_count, MPI_REAL, &
                p_vPerp_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( MPI_IN_PLACE, mpi_count, MPI_REAL, &
                p_vPar_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( MPI_IN_PLACE, mpi_count, MPI_REAL, &
                p_weight_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( MPI_IN_PLACE, mpi_count, MPI_INTEGER, &
                p_status_final(mpi_start_:mpi_end_), mpi_count, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_iErr )

        else

            call mpi_gather ( p_R_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, &
                p_R_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( p_z_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, &
                p_z_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( p_vPerp_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, &
                p_vPerp_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( p_vPar_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, &
                p_vPar_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( p_weight_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, &
                p_weight_final(mpi_start_:mpi_end_), mpi_count, MPI_REAL, 0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_gather ( p_status_final(mpi_start_:mpi_end_), mpi_count, MPI_INTEGER, &
                p_status_final(mpi_start_:mpi_end_), mpi_count, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_iErr )

        end if 

        !   Collect power absorption data

        write(*,*) mpi_pId, sum ( power ), sum(powerPrev)
        if ( mpi_pId == 0 ) write(*,*) 'Collecting power data ...'
        call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )
        if ( mpi_pId == 0 ) then 

            call mpi_reduce ( MPI_IN_PLACE, power, power_R_nBins * power_z_nBins, MPI_REAL, MPI_SUM, &
                0, MPI_COMM_WORLD, mpi_iErr )
            call mpi_reduce ( MPI_IN_PLACE, powerPar, power_R_nBins * power_z_nBins, MPI_REAL, MPI_SUM, &
                0, MPI_COMM_WORLD, mpi_iErr )

        else

           call mpi_reduce ( power, 0, power_R_nBins * power_z_nBins, MPI_REAL, MPI_SUM, &
               0, MPI_COMM_WORLD, mpi_iErr )
           call mpi_reduce ( powerPar, 0, power_R_nBins * power_z_nBins, MPI_REAL, MPI_SUM, &
               0, MPI_COMM_WORLD, mpi_iErr )
    
        endif
 
        call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )
        write(*,*) mpi_pId, sum ( power ), sum(powerPrev)
       
        !   Sum the particle diagnositics from each processor.

        call mpi_reduce ( nP_wall, nP_wall_total, 1, MPI_INTEGER, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_iErr )
        call mpi_reduce ( nP_bad, nP_bad_total, 1, MPI_INTEGER, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_iErr )
        call mpi_reduce ( nP_off_vGrid, nP_off_vGrid_total, 1, MPI_INTEGER, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_iErr )
        call mpi_reduce ( nP_off_QLgrid, nP_off_QLgrid_total, 1, MPI_REAL, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_iErr )
        call mpi_reduce ( nP_bad_coll, nP_bad_coll_total, 1, MPI_REAL, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_iErr )


        call cpu_time (T2)

        write(ncFileName,"( 'data/pList.dav.nc.',i3.3)") n 
        pl_fileName = ncFileName
        pList_is_nCdf = .true.

        if ( mpi_pId == 0 ) then

            nP_lowV = 0
            do i=1,nP
                if ( p_vPerp_final(i) < 2e4 ) nP_lowV = nP_lowV + 1
            enddo

            write (*,*) 'Time taken: ', T2-T1 
            write (*,*) '%Wall: ', real ( nP_wall_total ) / real ( nP ) * 100.0, nP_wall_total
            write (*,*) '%Bad: ', real ( nP_bad_total ) / real ( nP ) * 100.0, nP_bad_total
            write (*,*) '%off_vGrid: ', real ( nP_off_vGrid_total ) / real ( nP ) * 100.0, nP_off_vGrid_total
            write (*,*) '%off_QLGrid: ', real ( nP_off_QLgrid_total ) / real ( nP ) * 100.0, nP_off_QLgrid_total
            write (*,*) '%bad_coll: ', real ( nP_bad_coll_total ) / real ( nP ) * 100.0, nP_bad_coll_total
            write (*,*) '%lowV: ', real ( nP_lowV ) / real ( nP ) * 100.0, nP_lowV

            call write_f ( n, ncFileName )

        end if

        call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )
        !call erase_pl ()

    end do

    call stop_mpi () 

end program sMCrf 
