module powerAbsGrid

    use constants
     
    implicit none
    
    real :: power_R_range, power_z_Range
    real :: power_R_binSize, power_z_binSize
    real, allocatable :: power_R_binEdges (:), &
            power_z_binEdges (:)
    real, allocatable :: power_R_binCenters (:), &
            power_z_binCenters (:)
    real, allocatable :: power(:,:), powerPar(:,:)
    real, allocatable :: powerPrev(:,:), powerParPrev(:,:)
    real :: power_R_min, power_z_min
    integer :: power_R_nBins, power_z_nBins

contains
    subroutine init_powerAbs_grid ()
        use eqdsk
        use read_namelist
        implicit none

        integer :: i
        i = 0
        power_R_nBins   = 128 
        power_z_nBins   = 128  

        allocate ( power_R_binEdges ( power_R_nBins + 1 ), &
                   power_z_binEdges ( power_z_nBins + 1 ), &
                   power_R_binCenters ( power_R_nBins ), &
                   power_z_binCenters ( power_z_nBins ) )
        allocate ( power ( power_R_nBins, power_z_nBins ) )
        allocate ( powerPar ( power_R_nBins, power_z_nBins ) )
        allocate ( powerPrev ( power_R_nBins, power_z_nBins ) )
        allocate ( powerParPrev ( power_R_nBins, power_z_nBins ) )
   

        power_R_min = minVal ( rbbbs ) 
        power_R_range = maxVal ( rbbbs ) - minVal ( rbbbs )
        power_R_binSize   = power_R_range / power_R_nBins
        power_R_binEdges  = (/ (i*power_R_binSize,i=0,power_R_nBins) /) + power_R_min 
        power_R_binCenters    = power_R_binEdges(1:power_R_nBins) + power_R_binSize / 2.0

        power_z_min = minVal ( zbbbs )
        power_z_range = maxVal ( zbbbs ) - minVal ( zbbbs )
        power_z_binSize   = power_z_range / power_z_nBins
        power_z_binEdges  = (/ (i*power_z_binSize,i=0,power_z_nBins) /) + power_z_min
        power_z_binCenters    = power_z_binEdges(1:power_z_nBins) + power_z_binSize / 2.0

        power    = 0d0
        powerPar    = 0d0

    end subroutine init_powerAbs_grid

end module powerAbsGrid
