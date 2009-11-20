module read_ql
    use init_mpi
    use mpi
#   include <pnetcdf.inc>

    integer :: mpi_comm_cart
    real, allocatable :: ql_b(:,:,:,:),ql_c(:,:,:,:),ql_e(:,:,:,:),ql_f(:,:,:,:)
    real, allocatable :: ql_vMin(:,:), ql_vMax(:,:)
    real, allocatable :: ql_R(:), ql_z(:), ql_vPer(:), ql_vPar(:)
    integer :: chunk_nR, chunk_nz, chunk_nvPer, chunk_nvPar
    integer(KIND=MPI_OFFSET_KIND) :: ql_nR, ql_nz, ql_nvPer, ql_nvPar
    real :: ql_R_min, ql_R_max, ql_R_binSize, ql_R_range
    real :: ql_z_min, ql_z_max, ql_z_binSize, ql_z_range
    real :: ql_vPer_min, ql_vPer_max, ql_vPer_binSize, ql_vPer_range
    real :: ql_vPar_min, ql_vPar_max, ql_vPar_binSize, ql_vPar_range
    integer :: ql_nR_serial, ql_nz_serial, ql_nvPer_serial, ql_nvPar_serial
    real :: ql_vc_mks, ql_pScale

contains

subroutine dlg_p_check ( status )
    integer, intent ( in) :: status

    if(status /= nf_noerr) then 

        print *, trim(nfmpi_strerror(status))

        stop "p_check: Stopped"
    end if

end subroutine dlg_p_check 
 
subroutine open_ql_serial ()

    use netcdf
    use dlg

    implicit none
    integer :: ncId, qlB_id, ql_dim_ids(4), qlC_id, qlE_id, qlF_id
    character(len=100) :: ql_fileName
    integer :: ql_R_id, ql_z_id, ql_pScale_id
    integer :: ql_vPer_id, ql_vPar_id, ql_vc_mks_id
    integer :: aStat


    ql_fileName = 'data/p_ql.nc'

    !   Read in variables from .nc particle list file

    call dlg_check ( nf90_open ( path = ql_fileName, mode = nf90_nowrite, ncid = ncId ) )

    call dlg_check ( nf90_inq_varId ( ncId, 'ql_b', qlB_id ) )
    !call dlg_check ( nf90_inq_varId ( ncId, 'ql_c', qlC_id ) )
    !call dlg_check ( nf90_inq_varId ( ncId, 'ql_e', qlE_id ) )
    !call dlg_check ( nf90_inq_varId ( ncId, 'ql_f', qlF_id ) )

    call dlg_check ( nf90_inq_varId ( ncId, 'R', ql_R_id ) )
    call dlg_check ( nf90_inq_varId ( ncId, 'z', ql_z_id ) )
    call dlg_check ( nf90_inq_varId ( ncId, 'vPer', ql_vPer_id ) )
    call dlg_check ( nf90_inq_varId ( ncId, 'vPar', ql_vPar_id ) )
    call dlg_check ( nf90_inq_varId ( ncId, 'vc_mks', ql_vc_mks_id ) )
    call dlg_check ( nf90_inq_varId ( ncId, 'pscale', ql_pScale_id ) )

    call dlg_check ( nf90_inquire_variable ( &
        ncId, qlB_id, dimIds = ql_dim_ids ) )
 
    call dlg_check ( nf90_inquire_dimension ( ncId, ql_dim_ids(1), &
        len = ql_nR_serial ) ) 
    call dlg_check ( nf90_inquire_dimension ( ncId, ql_dim_ids(2), &
        len = ql_nz_serial ) ) 
    call dlg_check ( nf90_inquire_dimension ( ncId, ql_dim_ids(3), &
        len = ql_nvPer_serial ) ) 
    call dlg_check ( nf90_inquire_dimension ( ncId, ql_dim_ids(4), &
        len = ql_nvPar_serial ) ) 

    allocate ( ql_R ( ql_nR_serial ), &
        ql_z ( ql_nz_serial ), &
        ql_vPer ( ql_nvPer_serial ), &
        ql_vPar ( ql_nvPar_serial ) )

    call dlg_check ( nf90_get_var ( ncId, ql_R_id, ql_R ) )
    call dlg_check ( nf90_get_var ( ncId, ql_z_id, ql_z ) )
    call dlg_check ( nf90_get_var ( ncId, ql_vPer_id, ql_vPer ) )
    call dlg_check ( nf90_get_var ( ncId, ql_vPar_id, ql_vPar ) )
    call dlg_check ( nf90_get_var ( ncId, ql_vc_mks_id, ql_vc_mks ) )
    call dlg_check ( nf90_get_var ( ncId, ql_pScale_id, ql_pScale ) )

    ql_R_binSize    = ql_R(2) - ql_R(1)
    ql_R_min    = minVal ( ql_R ) - ql_R_binSize / 2.0
    ql_R_max    = maxVal ( ql_R ) + ql_R_binSize / 2.0
    ql_R_range  = ql_R_max - ql_R_min

    ql_z_binSize    = abs ( ql_z(2) - ql_z(1) )
    ql_z_min    = minVal ( ql_z ) - ql_z_binSize / 2.0
    ql_z_max    = maxVal ( ql_z ) + ql_z_binSize / 2.0
    ql_z_range  = ql_z_max - ql_z_min

    ql_vPer_binSize    = ql_vPer(2) - ql_vPer(1)
    ql_vPer_min    = minVal ( ql_vPer ) 
    ql_vPer_max    = maxVal ( ql_vPer ) 
    ql_vPer_range  = ql_vPer_max - ql_vPer_min

    ql_vPar_binSize    = ql_vPar(2) - ql_vPar(1)
    ql_vPar_min    = minVal ( ql_vPar ) - ql_vPar_binSize / 2.0
    ql_vPar_max    = maxVal ( ql_vPar ) + ql_vPar_binSize / 2.0
    ql_vPar_range  = ql_vPar_max - ql_vPar_min

    if ( ql_R_binSize < 0 ) read(*,*)
    if ( ql_z_binSize < 0 ) read(*,*)
    if ( ql_vPer_binSize < 0 ) read(*,*)
    if ( ql_vPar_binSize < 0 ) read(*,*)


    allocate ( ql_b ( ql_nR_serial, ql_nz_serial, ql_nvPer_serial, ql_nvPar_serial ), stat = aStat )
        if ( aStat .ne. 0 ) write(*,*) 'DLG: ERROR - allocate failed for QL B, 113'
    !allocate ( ql_c ( ql_nR_serial, ql_nz_serial, ql_nvPer_serial, ql_nvPar_serial ), stat = aStat )
    !    if ( aStat .ne. 0 ) write(*,*) 'DLG: ERROR - allocate failed for QL C, 115'
    !allocate ( ql_e ( ql_nR_serial, ql_nz_serial, ql_nvPer_serial, ql_nvPar_serial ), stat = aStat )
    !    if ( aStat .ne. 0 ) write(*,*) 'DLG: ERROR - allocate failed for QL E, 117'
    !allocate ( ql_f ( ql_nR_serial, ql_nz_serial, ql_nvPer_serial, ql_nvPar_serial ), stat = aStat )
    !    if ( aStat .ne. 0 ) write(*,*) 'DLG: ERROR - allocate failed for QL F, 118'

    call dlg_p_check ( nf90_get_var ( &
        ncId, qlB_id, ql_b ) )
    !call dlg_p_check ( nf90_get_var ( &
    !    ncId, qlC_id, ql_c ) )
    !call dlg_p_check ( nf90_get_var ( &
    !    ncId, qlE_id, ql_e ) )
    !call dlg_p_check ( nf90_get_var ( &
    !    ncId, qlF_id, ql_f ) )

    ql_b    = ql_b * ql_pScale
    !ql_c    = ql_c * ql_pScale
    !ql_e    = ql_e * ql_pScale
    !ql_f    = ql_f * ql_pScale

    call dlg_check ( nf90_close ( ncId ) )

    call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )

end subroutine open_ql_serial


subroutine open_ql ()

    implicit none
    integer :: ncId, qlB_id, ql_dim_ids(4), qlC_id, qlE_id, qlF_id
    character(len=100) :: ql_fileName
    integer(KIND=MPI_OFFSET_KIND) :: start(4), cnt(4)
    integer :: mpi_nDims, mpi_dimSize(2), mpi_iErr
    logical :: mpi_dimPeriodicity(2), mpi_reOrder 
    integer :: mpi_coords(2)
    integer :: ql_vMin_id, ql_vMax_id, ql_R_id, ql_z_id
    integer :: ql_vPer_id, ql_vPar_id

    !   Test vars
    integer :: iiR, iiz, iivPer, iivPar
    integer :: mpi_coords_hasData(2), mpi_pId_hasData, &
        mpi_pId_wantsData
    integer :: mpi_tag_getData, mpi_tag_whoWants
    integer, dimension ( MPI_STATUS_SIZE ) :: mpi_status
    real :: ql_b_recv

    ql_fileName = 'data/p_ql.nc'

    !   Read in variables from .nc particle list file

    call dlg_p_check ( nfmpi_open ( MPI_COMM_WORLD, &
        ql_fileName, nf_nowrite, MPI_INFO_NULL, ncId ) )

    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'ql_b', qlB_id ) )
    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'ql_c', qlC_id ) )
    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'ql_e', qlE_id ) )
    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'ql_f', qlF_id ) )

    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'maxV', ql_vMax_id ) )
    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'minV', ql_vMin_id ) )
    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'R', ql_R_id ) )
    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'z', ql_z_id ) )
    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'vPer', ql_vPer_id ) )
    call dlg_p_check ( nfmpi_inq_varId ( ncId, 'vPar', ql_vPar_id ) )

    call dlg_p_check ( nfmpi_inq_varDimId ( &
        ncId, qlB_id, ql_dim_ids ) )
    
    call dlg_p_check ( nfmpi_inq_dimLen ( ncId, ql_dim_ids(1), ql_nR ) ) 
    call dlg_p_check ( nfmpi_inq_dimLen ( ncId, ql_dim_ids(2), ql_nz ) ) 
    call dlg_p_check ( nfmpi_inq_dimLen ( ncId, ql_dim_ids(3), ql_nvPer ) ) 
    call dlg_p_check ( nfmpi_inq_dimLen ( ncId, ql_dim_ids(4), ql_nvPar ) ) 

    allocate ( ql_vMin ( ql_nR, ql_nz ), ql_vMax ( ql_nR, ql_nz ) )
    allocate ( ql_R ( ql_nR ), ql_z ( ql_nz ), ql_vPer ( ql_nvPer ), ql_vPar ( ql_nvPar ) )

    !   Give every processor a copy of the vMax/vMin arrays

    call dlg_p_check ( nfmpi_get_var_real_all ( ncId, ql_vMax_id, ql_vMax ) )
    call dlg_p_check ( nfmpi_get_var_real_all ( ncId, ql_vMin_id, ql_vMin ) )
    call dlg_p_check ( nfmpi_get_var_real_all ( ncId, ql_R_id, ql_R ) )
    call dlg_p_check ( nfmpi_get_var_real_all ( ncId, ql_z_id, ql_z ) )
    call dlg_p_check ( nfmpi_get_var_real_all ( ncId, ql_vPer_id, ql_vPer ) )
    call dlg_p_check ( nfmpi_get_var_real_all ( ncId, ql_vPar_id, ql_vPar ) )

    ql_R_binSize    = ql_R(2) - ql_R(1)
    ql_R_min    = minVal ( ql_R ) - ql_R_binSize / 2.0
    ql_R_max    = maxVal ( ql_R ) + ql_R_binSize / 2.0
    ql_R_range  = ql_R_max - ql_R_min

    ql_z_binSize    = abs ( ql_z(2) - ql_z(1) )
    ql_z_min    = minVal ( ql_z ) - ql_z_binSize / 2.0
    ql_z_max    = maxVal ( ql_z ) + ql_z_binSize / 2.0
    ql_z_range  = ql_z_max - ql_z_min

    ql_vPer_binSize    = ql_vPer(2) - ql_vPer(1)
    ql_vPer_min    = minVal ( ql_vPer ) - ql_vPer_binSize / 2.0
    ql_vPer_max    = maxVal ( ql_vPer ) + ql_vPer_binSize / 2.0
    ql_vPer_range  = ql_vPer_max - ql_vPer_min

    ql_vPar_binSize    = ql_vPar(2) - ql_vPar(1)
    ql_vPar_min    = minVal ( ql_vPar ) - ql_vPar_binSize / 2.0
    ql_vPar_max    = maxVal ( ql_vPar ) + ql_vPar_binSize / 2.0
    ql_vPar_range  = ql_vPar_max - ql_vPar_min

    !   Setup mpi grid

    mpi_nDims   = 2
    mpi_dimSize(1)  = int ( sqrt ( real ( mpi_nProcs ) ) )
    mpi_dimSize(2)  = int ( sqrt ( real ( mpi_nProcs ) ) )
    mpi_dimPeriodicity(1)   = .false.
    mpi_dimPeriodicity(2)   = .false.
    mpi_reOrder = .true.

    call mpi_cart_create ( MPI_COMM_WORLD, mpi_nDims, mpi_dimSize, &
        mpi_dimPeriodicity, mpi_reOrder, mpi_comm_cart, mpi_iErr ) 
    
    !   Get chunk processor location on grid

    call mpi_cart_coords ( mpi_comm_cart, mpi_pId, 2, mpi_coords, mpi_iErr )

    !   Setup ql chunk sizes

    if ( mod ( ql_nR, mpi_dimSize(1) ) .ne. 0 ) then 

        write(*,*) 'ERROR: sqrt(nProcs) must equal nR & nz'
        read(*,*)

    end if

    chunk_nR = ql_nR / mpi_dimSize(1) 
    chunk_nz = ql_nz / mpi_dimSize(2) 
    chunk_nvPer  = ql_nvPer
    chunk_nvPar  = ql_nvPar

    start(1)    = mpi_coords(1)*chunk_nR+1
    start(2)    = mpi_coords(2)*chunk_nz+1
    start(3)    = 1
    start(4)    = 1

    cnt(1)  = chunk_nR
    cnt(2)  = chunk_nz
    cnt(3)  = chunk_nvPer
    cnt(4)  = chunk_nvPar

    write(*,*) 'Allocating memory for QL operator'
 
    allocate ( ql_b ( start(1):start(1)+chunk_nR-1,&
                     start(2):start(2)+chunk_nz-1,&
                     chunk_nvPer, chunk_nvPar ) )
    allocate ( ql_c ( start(1):start(1)+chunk_nR-1,&
                     start(2):start(2)+chunk_nz-1,&
                     chunk_nvPer, chunk_nvPar ) )
    allocate ( ql_e ( start(1):start(1)+chunk_nR-1,&
                     start(2):start(2)+chunk_nz-1,&
                     chunk_nvPer, chunk_nvPar ) )
    allocate ( ql_f ( start(1):start(1)+chunk_nR-1,&
                     start(2):start(2)+chunk_nz-1,&
                     chunk_nvPer, chunk_nvPar ) )

    write (*,*) mpi_pId, start(1:2), cnt(1:2), mpi_coords
    
    call dlg_p_check ( nfmpi_get_vara_real_all ( &
        ncId, qlB_id, start, cnt, ql_b ) )
    call dlg_p_check ( nfmpi_get_vara_real_all ( &
        ncId, qlC_id, start, cnt, ql_c ) )
    call dlg_p_check ( nfmpi_get_vara_real_all ( &
        ncId, qlE_id, start, cnt, ql_e ) )
    call dlg_p_check ( nfmpi_get_vara_real_all ( &
        ncId, qlF_id, start, cnt, ql_f ) )
 
    call dlg_p_check ( nfmpi_close ( ncId ) )

    write(*,*) 'QL operator read ;-)'
    call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )

    !   Test read values

    iiR = 8
    iiz = 4
    iivPer  = 64
    iivPar  = 60

    !   Which processor has the data?

    mpi_coords_hasData  = (/ (iiR-1) / chunk_nR, (iiz-1) / chunk_nz /)
    call mpi_cart_rank ( mpi_comm_cart, mpi_coords_hasData, mpi_pId_hasData, mpi_iErr )
   
    if ( mpi_pId_hasData .eq. mpi_pId ) then 

        write(*,*) mpi_pId, ql_b(iiR,iiz,iivPer,iivPar), mpi_coords

    end if 

    !  Try a data request 

    mpi_tag_whoWants    = 1
    mpi_tag_getData = 2

    if ( mpi_pId .eq. mpi_pId_hasData ) &
    call mpi_recv ( mpi_pId_wantsData, 1, MPI_INTEGER, &
        MPI_ANY_SOURCE, mpi_tag_whoWants, mpi_comm_cart, mpi_status, mpi_iErr)
 
    if ( mpi_pID .eq. 2 ) then

        if ( mpi_pId .ne. mpi_pId_hasData ) &
        call mpi_send ( mpi_pId, 1, MPI_INTEGER, &
            mpi_pId_hasData, mpi_tag_whoWants, mpi_comm_cart, mpi_iErr)
    
        if ( mpi_pId .ne. mpi_pId_hasData ) &
        call mpi_recv ( ql_b_recv, 1, MPI_REAL, &
            mpi_pId_hasData, mpi_tag_getData, mpi_comm_cart, mpi_status, mpi_iErr ) 
     
    end if

    if ( mpi_pId .eq. mpi_pId_hasData ) &
    call mpi_send ( ql_b(iiR,iiz,iivPer,iivPar), 1, MPI_REAL, &
        mpi_pId_wantsData, mpi_tag_getData, mpi_comm_cart, mpi_iErr)
    
    write(*,*) mpi_pId, ql_b_recv, mpi_pId_hasData, mpi_coords_hasData

    call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )

end subroutine open_ql

end module read_ql
