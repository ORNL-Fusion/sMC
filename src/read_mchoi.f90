module read_mchoi
    implicit none

    real, allocatable :: mchoi_rho(:,:), &
        mchoi_theta(:,:), &
        mchoi_R(:), &
        mchoi_z(:)

    complex, allocatable :: &
        mchoi_ePlus(:,:), &
        mchoi_eMinu(:,:), &
        mchoi_kPer(:,:)

    real, allocatable :: mChoi_ePlus_real(:,:), &
            mChoi_ePlus_img(:,:), &
            mChoi_eMinu_real(:,:), &
            mChoi_eMinu_img(:,:), &
            mChoi_kPer_real(:,:), &
            mChoi_kPer_img(:,:)

    real, allocatable :: zp_ePlus_real(:), zp_eMinu_real(:), zp_kPer_real(:)
    real, allocatable :: zp_ePlus_img(:), zp_eMinu_img(:), zp_kPer_img(:)

    integer :: mchoi_nR, mchoi_nz

contains

subroutine read_aorsa_fields ()

    use read_namelist
    use interp
    use netcdf
    use netcdf_check

    implicit none

        integer :: i, j, nN, nM, aStat

        !   Interpolation vars

        integer :: islpsw, iErr
        real, allocatable :: temp(:)
        real, allocatable :: zx1(:), zxm(:), zy1(:), zyn(:)
        real :: zxy11, zxym1, zxy1n, zxymn

        !   netCDF vars

        integer :: nc_id, ePlus_id, dim_ids(2)
        integer :: eMinu_id, ePlus_img_id, eMinu_img_id, &
            kPer_cold_id, kPer_img_cold_id, R_id, z_id

        write(*,*) 'Reading mchoi data file'

        !   Read in variables from .nc particle list file

        call check ( nf90_open ( path = mchoi_fileName, mode = nf90_nowrite, ncid = nc_id ) )

        write(*,*) 'File opened'

        call check ( nf90_inq_varId ( nc_id, 'ePlus', ePlus_id ) )
        call check ( nf90_inq_varId ( nc_id, 'eMinu', eMinu_id ) )
        call check ( nf90_inq_varId ( nc_id, 'ePlus_img', ePlus_img_id ) )
        call check ( nf90_inq_varId ( nc_id, 'eMinu_img',eMinu_img_id ) )
        call check ( nf90_inq_varId ( nc_id, 'kPer_cold',kPer_cold_id ) )
        call check ( nf90_inq_varId ( nc_id, 'kPer_img_cold',kPer_img_cold_id ) )
        call check ( nf90_inq_varId ( nc_id, 'R', R_id ) )
        call check ( nf90_inq_varId ( nc_id, 'z', z_id ) )

        call check ( nf90_inquire_variable ( &
            nc_id, ePlus_id, dimIds = dim_ids ) )
        
        call check ( nf90_inquire_dimension ( nc_id, dim_ids(1), &
            len = mchoi_nR ) ) 
         call check ( nf90_inquire_dimension ( nc_id, dim_ids(2), &
            len = mchoi_nz ) ) 

        write(*,*) 'nR, nZ: ', mchoi_nR, mchoi_nz

         allocate ( mChoi_R(mchoi_nR), &
                    mChoi_z(mchoi_nz), &
                    mChoi_ePlus_real(mchoi_nR,mchoi_nz), &
                    mChoi_eMinu_real(mchoi_nR,mchoi_nz), &
                    mChoi_kPer_real(mchoi_nR,mchoi_nz), &
                    mChoi_ePlus_img(mchoi_nR,mchoi_nz), &
                    mChoi_eMinu_img(mchoi_nR,mchoi_nz), &
                    mChoi_kPer_img(mchoi_nR,mchoi_nz), &
                    mChoi_ePlus(mchoi_nR,mchoi_nz), &
                    mChoi_eMinu(mchoi_nR,mchoi_nz), &
                    mChoi_kPer(mchoi_nR,mchoi_nz), stat = aStat )
        if ( aStat .ne. 0 ) &
            write(*,*) 'dlg: ERROR - allocation failed in read_mchoi.f90'


        call check ( nf90_get_var ( nc_id, ePlus_id, mchoi_ePlus_real ) )
        call check ( nf90_get_var ( nc_id, ePlus_img_id, mchoi_ePlus_img ) )
        call check ( nf90_get_var ( nc_id, eMinu_id, mchoi_eMinu_real ) )
        call check ( nf90_get_var ( nc_id, eMinu_img_id, mchoi_eMinu_img ) )
        call check ( nf90_get_var ( nc_id, kPer_cold_id, mchoi_kPer_real ) )
        call check ( nf90_get_var ( nc_id, kPer_img_cold_id, mchoi_kPer_img ) )
        call check ( nf90_get_var ( nc_id, R_id, mchoi_R ) )
        call check ( nf90_get_var ( nc_id, z_id, mchoi_z ) )
        
        call check ( nf90_close ( nc_id ) )

        write(*,*) 'Mchoi file closed'

        mchoi_ePlus = cmplx ( mChoi_ePlus_real, mChoi_ePlus_img )
        mchoi_eMinu = cmplx ( mChoi_eMinu_real, mChoi_eMinu_img )
        mchoi_kPer = cmplx ( mChoi_kPer_real, mChoi_kPer_img )

        !mchoi_ePlus = mChoi_ePlus_real
        !mchoi_eMinu =  mChoi_eMinu_real
        !mchoi_kPer =  mChoi_kPer_real

        !deallocate ( mChoi_ePlus_real, mChoi_ePlus_img, &
        !                mChoi_eMinu_real, mChoi_eMinu_img, &
        !                mChoi_kPer_real, mChoi_kPer_img )

!        open ( unit = 8, file = mchoi_fileName, status = 'OLD', action = 'READ' )
!
!        read (8,2000 ) nN, nM
!
!        allocate ( mChoi_rho(nN,nM), &
!                    mChoi_theta(nN,nM), &
!                    mChoi_R(nN,nM), &
!                    mChoi_z(nN,nM), &
!                    mChoi_ePlus(nN,nM), &
!                    mChoi_eMinu(nN,nM), &
!                    mChoi_kPer(nN,nM), &
!                    mChoi_n(nN,nM), &
!                    mChoi_m(nN,nM) )
!
!        do i=1,nN
!        do j=1,nM
!          
!            write(*,*) i, j 
!            read (8,1313) mchoi_n(i,j), &
!                mchoi_m(i,j), &
!                mchoi_rho(i,j), &
!                mchoi_theta(i,j), &
!                mchoi_R(i,j), &
!                mchoi_z(i,j), &
!                mchoi_ePlus(i,j), &
!                mchoi_eMinu(i,j), &
!                mchoi_kPer(i,j)
! 
!        end do
!        end do 
!
!        2000 format ( 2i10 )
!        1313 format ( 2i10,7e12.4 )
!
!        close ( unit = 8 )
       
        !   Initialise the ePlus & eMinu interpolations

        islpsw  = 255 

        allocate ( zx1(mchoi_nR), zxm(mchoi_nR), zy1(mchoi_nz), zyn(mchoi_nz), stat = aStat )
        if ( aStat .ne. 0 ) &
            write(*,*) 'dlg: ERROR - allocation failed in read_mchoi.f90'

        allocate ( temp(mchoi_nz+mchoi_nz+mchoi_nR), &
            zp_ePlus_real(3*mchoi_nR*mchoi_nz), zp_eMinu_real(3*mchoi_nR*mchoi_nz), &
            zp_kPer_real(3*mchoi_nR*mchoi_nz), &
            zp_ePlus_img(3*mchoi_nR*mchoi_nz), zp_eMinu_img(3*mchoi_nR*mchoi_nz), &
            zp_kPer_img(3*mchoi_nR*mchoi_nz), stat = aStat )
        if ( aStat .ne. 0 ) &
            write(*,*) 'dlg: ERROR - allocation failed in read_mchoi.f90'

        temp    = 0
        zx1 = 0
        zxm = 0
        zy1 = 0
        zyn = 0
        zp_ePlus_real   = 0
        zp_eMinu_real   = 0
        zp_kPer_real    = 0
        zp_ePlus_img    = 0
        zp_eMinu_img    = 0
        zp_kPer_img = 0

        write(*,*) 'About to initialise the mchoi interpolations'

        call surf1 ( mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, mchoi_ePlus_real, mchoi_nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_ePlus_real, temp, sigma, iErr)

        call surf1 ( mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, mchoi_eMinu_real, mchoi_nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_eMinu_real, temp, sigma, iErr)

        call surf1 ( mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, mchoi_kPer_real, mchoi_nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_kPer_real, temp, sigma, iErr)

        call surf1 ( mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, mchoi_ePlus_img, mchoi_nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_ePlus_img, temp, sigma, iErr)
 
        call surf1 ( mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, mchoi_eMinu_img, mchoi_nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_eMinu_img, temp, sigma, iErr)

        call surf1 ( mchoi_nR, mchoi_nz, mchoi_R, mchoi_z, mchoi_kPer_img, mchoi_nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_kPer_img, temp, sigma, iErr)

end subroutine read_aorsa_fields

end module read_mchoi
