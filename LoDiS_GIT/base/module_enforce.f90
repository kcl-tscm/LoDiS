MODULE ENFORCE
USE PARACLUSTER
!
IMPLICIT NONE
 !
 ! CALCULATED VALUES
 !Position
 REAL*8 :: x0(nsiz),y0(nsiz),z0(nsiz) !!initial coordinates for atom i
 REAL*8, ALLOCATABLE :: x(:),y(:),z(:) !!position array for atoms
 REAL*8, ALLOCATABLE :: u(:),v(:),w(:) ! displacement at time t                          
 REAL*8, ALLOCATABLE :: du(:),dv(:),dw(:) ! disp at time t-dt
 REAL*8, ALLOCATABLE :: v_acf(:,:,:),vel_act_est(:) !velocity and acf for calculation 
 !real*8, allocatable :: vel_zero_x(:),vel_zero_y(:),vel_zero_z(:)
 !real*8 :: vel_acf_zero
 !Voisin and Bigvoi
 INTEGER, ALLOCATABLE :: nv4(:),iv4(:,:)
 INTEGER, ALLOCATABLE :: nvois(:),ivois(:,:)
 ! Energy/force
 REAL*8 :: ener,etot,ecin,temp
 REAL*8, ALLOCATABLE ::  potener(:) ! potential energy of each atom
 REAL*8, ALLOCATABLE ::  fx(:),fy(:),fz(:)             
 REAL*8, ALLOCATABLE ::  dfx(:),dfy(:),dfz(:)
 REAL*8, ALLOCATABLE ::  vx(:),vy(:),vz(:)
! REAL*8, ALLOCATABLE ::  vel_zero_x(:),vel_zero_y(:),vel_zero_z(:)
! REAL*8 :: vel_acf_zero
 !
 INTEGER :: nfile    !for writing output-files
 REAL*8 :: edelta,tpar !averaged energy and T
 CHARACTER (LEN=15), ALLOCATABLE :: filename(:)
 !
 ! More outputs
 REAL*8 :: rshell
 !
END MODULE ENFORCE
