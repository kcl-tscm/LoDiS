MODULE SUBSTRATE
 !===================================
 ! Global variables for MgO substrate
 !===================================
 IMPLICIT NONE
 INTEGER                        :: sub_geom          
!values default = 0 no support; allowed = 1 for a double square (MgO) ; =2 for a hexagon ~ graphene  

 CHARACTER(LEN=2), DIMENSION(2) :: mgo_met           ! Names of 2 metal interacting with MgO      
 REAL(8)                        :: alpha_sub, rc_sub
 REAL(8)                        :: beta_sub, beta_sub_angstrom
 REAL(8)                        :: amgo, mgo_alat    ! latt constant MgO (O-O dist) [A]
 REAL(8), DIMENSION(2,3,3,3)    :: ccc               ! Parameters c(i,j,1:3) for atom A and B
 REAL(8), DIMENSION(2,2,3)      :: ccc4              ! extra c parameters for MOT effect
 REAL(8)                        :: ener_sub          ! potential energy substrate (?)

 REAL(8), ALLOCATABLE :: mgo_cn(:), mgo_s(:,:), mgo_dS(:,:) 
 REAL(8), ALLOCATABLE :: mgo_dcn_dx(:), mgo_dcn_dy(:), mgo_dcn_dz(:)  

 REAL(8), ALLOCATABLE :: mgo_cn_mot(:), mgo_s_mot(:,:), mgo_dS_mot(:,:) 
 REAL(8), ALLOCATABLE :: mgo_dcn_mot_dx(:), mgo_dcn_mot_dy(:), mgo_dcn_mot_dz(:)  

 ! Normalised positions
 REAL(8), ALLOCATABLE :: mgo_x(:), mgo_y(:), mgo_z(:)

 REAL(8), ALLOCATABLE :: mgo_fx(:), mgo_fy(:), mgo_fz(:)

END MODULE SUBSTRATE
