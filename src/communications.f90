module communications
use mpi
use init_mpi
use read_ql

integer :: mpi_tag_getData = 4, mpi_tag_wantsData = 3, mpi_tag_finished = 5
integer :: mpi_pId_wantsData
logical :: mpi_flag_wantsData, mpi_flag_getData, mpi_flag_finished, mpi_flag_finSend
integer :: mpi_pId_and_indices(5)
integer :: mpi_req_getData, mpi_req_finished, mpi_req_wantsData 
integer :: mpi_req_finSend
integer, dimension ( MPI_STATUS_SIZE ) :: mpi_stat_wantsData
integer, allocatable :: nP_done(:)
integer :: mpi_finished_pId
integer, dimension ( MPI_STATUS_SIZE ) :: mpi_stat_finished
integer, dimension ( MPI_STATUS_SIZE ) :: mpi_stat_finSend

contains

subroutine dlg_sendData ()
    use read_ql
    implicit none

    real :: ql_send(4)

    !   Check if there is a request waiting 
    call mpi_test ( mpi_req_wantsData, mpi_flag_wantsData, &
        mpi_stat_wantsData, mpi_iErr )

    !   If there is send data

    if ( mpi_flag_wantsData ) then 

        mpi_pId_wantsData   = mpi_pId_and_indices(1)

        if ( mpi_pId_wantsData .lt. mpi_nProcs ) then 

            ql_send = (/ ql_b( mpi_pId_and_indices(2), &
                                  mpi_pId_and_indices(3), &
                                  mpi_pId_and_indices(4), &
                                  mpi_pId_and_indices(5)), &
                         ql_c( mpi_pId_and_indices(2), &
                                  mpi_pId_and_indices(3), &
                                  mpi_pId_and_indices(4), &
                                  mpi_pId_and_indices(5)), &
                         ql_e( mpi_pId_and_indices(2), &
                                  mpi_pId_and_indices(3), &
                                  mpi_pId_and_indices(4), &
                                  mpi_pId_and_indices(5)), &
                         ql_f( mpi_pId_and_indices(2), &
                                  mpi_pId_and_indices(3), &
                                  mpi_pId_and_indices(4), &
                                  mpi_pId_and_indices(5) ) /)

            call mpi_send ( ql_send, 4, MPI_REAL, &
                mpi_pId_wantsData, mpi_tag_getData, mpi_comm_cart, mpi_iErr)

        end if

        !   Reset request

        call mpi_irecv ( mpi_pId_and_indices, 5, MPI_INTEGER, &
            MPI_ANY_SOURCE, mpi_tag_wantsData, mpi_comm_cart, &
            mpi_req_wantsData, mpi_iErr)

    end if

end subroutine dlg_sendData


subroutine dlg_sendData_wait ()
    use read_ql
    implicit none

    real :: ql_send(4)

    !   Wait for a request 
    call mpi_wait ( mpi_req_wantsData, &
        mpi_stat_wantsData, mpi_iErr )

    mpi_pId_wantsData   = mpi_pId_and_indices(1)

    !   Check if this is a legit request or an 
    !   I'm finished you can stop waiting for me to ask
    !   for data request ;-)

    if ( mpi_pId_wantsData .lt. mpi_nProcs ) then 

        ql_send = (/ ql_b( mpi_pId_and_indices(2), &
                              mpi_pId_and_indices(3), &
                              mpi_pId_and_indices(4), &
                              mpi_pId_and_indices(5)), &
                     ql_c( mpi_pId_and_indices(2), &
                              mpi_pId_and_indices(3), &
                              mpi_pId_and_indices(4), &
                              mpi_pId_and_indices(5)), &
                     ql_e( mpi_pId_and_indices(2), &
                              mpi_pId_and_indices(3), &
                              mpi_pId_and_indices(4), &
                              mpi_pId_and_indices(5)), &
                     ql_f( mpi_pId_and_indices(2), &
                              mpi_pId_and_indices(3), &
                              mpi_pId_and_indices(4), &
                              mpi_pId_and_indices(5) ) /)

        call mpi_send ( ql_send, 4, MPI_REAL, &
            mpi_pId_wantsData, mpi_tag_getData, mpi_comm_cart, mpi_iErr)

    end if

    !   Reset request

    call mpi_irecv ( mpi_pId_and_indices, 5, MPI_INTEGER, &
        MPI_ANY_SOURCE, mpi_tag_wantsData, mpi_comm_cart, &
        mpi_req_wantsData, mpi_iErr)

end subroutine dlg_sendData_wait


subroutine dlg_check_finished ()
    implicit none

    call mpi_test ( mpi_req_finished, mpi_flag_finished, &
        mpi_stat_finished, mpi_iErr )
    
    if ( mpi_flag_finished ) then 
    
        nP_done ( mpi_finished_pId+1 ) = 1
    
        if ( sum ( nP_done ) .ne. mpi_nProcs ) &
        call mpi_irecv ( mpi_finished_pId, 1, MPI_INTEGER, &
            MPI_ANY_SOURCE, mpi_tag_finished, mpi_comm_cart, &
            mpi_req_finished, mpi_iErr)
    
    end if

end subroutine dlg_check_finished 

subroutine dlg_wakeUp (i)
    implicit none
    integer, intent(in) :: i

    call mpi_send ( (/ mpi_nProcs+1,0,0,0,0 /), &
        5, MPI_INTEGER, &
        i, mpi_tag_wantsData, mpi_comm_cart, mpi_iErr)

end subroutine dlg_wakeUp 

end module communications
