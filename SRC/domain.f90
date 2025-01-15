! Set the domain dimensionality, size and number of subdomains.

module domain

!  VOCALS; Case 4 ICMW 2016
!       integer, parameter :: YES3D = 1  ! Domain dimensionality: 1 - 3D, 0 - 2D
!       integer, parameter :: nx_gl = 200 ! Number of grid points in X
!       integer, parameter :: ny_gl = 200 ! Number of grid points in Y
!       integer, parameter :: nz_gl = 160 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 4 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 8 ! No of subdomains in y

!       integer, parameter :: YES3D = 1  ! Domain dimensionality: 1 - 3D, 0 - 2D
!       integer, parameter :: nx_gl = 10 ! Number of grid points in X
!       integer, parameter :: ny_gl = 10 ! Number of grid points in Y
!       integer, parameter :: nz_gl = 160 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 1 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 1 ! No of subdomains in y

! ARM LES case sgp20150627

       ! integer, parameter :: YES3D = 1  ! Domain dimensionality: 1 - 3D, 0 - 2D
       ! integer, parameter :: nx_gl = 10 ! Number of grid points in X
       ! integer, parameter :: ny_gl = 10 ! Number of grid points in Y
       ! integer, parameter :: nz_gl = 337 ! Number of pressure (scalar) levels
       ! integer, parameter :: nsubdomains_x  = 1 ! No of subdomains in x
       ! integer, parameter :: nsubdomains_y  = 1 ! No of subdomains in y

! LASSO ENA test cases

       integer, parameter :: YES3D = 1  ! Domain dimensionality: 1 - 3D, 0 - 2D
       ! integer, parameter :: nx_gl = 256 ! Number of grid points in Y
       ! integer, parameter :: ny_gl = 256 ! Number of grid points in Y
       ! integer, parameter :: nz_gl = 250 ! Number of pressure (scalar) levels
       ! integer, parameter :: nsubdomains_x  = 16 ! No of subdomains in x
       ! integer, parameter :: nsubdomains_y  = 16 ! No of subdomains in y
       integer, parameter :: nx_gl = 1024 ! Number of grid points in Y
       integer, parameter :: ny_gl = 1024 ! Number of grid points in Y
       integer, parameter :: nz_gl = 250 ! Number of pressure (scalar) levels
       integer, parameter :: nsubdomains_x  = 32 ! No of subdomains in x
       integer, parameter :: nsubdomains_y  = 32 ! No of subdomains in y

       ! define # of points in x and y direction to average for 
       !   output relating to statistical moments.
       ! For example, navgmom_x = 8 means the output will be   
       !  8 times coarser grid than the original.
       ! If don't wanna such output, just set them to -1 in both directions. 
       ! See Changes_log/README.UUmods for more details.
       integer, parameter :: navgmom_x = -1 
       integer, parameter :: navgmom_y = -1 

       integer, parameter :: ntracers = 0 ! number of transported tracers (dotracers=.true.)
       
! Note:
!  * nx_gl and ny_gl should be a factor of 2,3, or 5 (see User's Guide)
!  * if 2D case, ny_gl = nsubdomains_y = 1 ;
!  * nsubdomains_x*nsubdomains_y = total number of processors
!  * if one processor is used, than  nsubdomains_x = nsubdomains_y = 1;
!  * if ntracers is > 0, don't forget to set dotracers to .true. in namelist 

end module domain
