SUBROUTINE force_geometry

  USE PARACLUSTER
  USE CLUSTER 
  USE POTENTIAL
  USE ENFORCE
  USE DISTANCE
  USE SUBSTRATE
  USE module_function_var

  IMPLICIT NONE
  
  INTEGER :: n,i,j
  INTEGER :: mgo_type
  REAL (8):: aa,ener_mgo,s3
  
CALL mgo_coord

! Now I have cn(i) and and its derivative w.r.t. x, y, z
 ! They are called 
 ! mgo_cn(i)                                    adimensional
 ! mgo_dcn_dx(i), mgo_dcn_dy(i), mgo_dcn_dz(i)  [1/L]
 ! In this moment the unit is [1/Angstrom], and can be modified in mgo_coord_no

 ! Normalised positions  (in Angstrom and then *mgo_alat)

 ! Inizialisation of local variables
ener_mgo = 0.d0
 !------------------------------------------------------
 ! Inizialisation of global variables
ener_sub = 0.d0 ! This is global  
mgo_fx   = 0.d0
mgo_fy   = 0.d0
mgo_fz   = 0.d0
  

mgo_x = (x + u) *arete(1) ! [A]
mgo_y = (y + v) *arete(1) ! [A]
mgo_z = (z + w) *arete(1) ! [A]

aa = mgo_alat             ! [1/A]
s3 = sqrt(3.0d0)


IF (geometry_type == 1) THEN
 DO n = 1, natom
    mgo_type = itype(n)
    !-------------------------
    funct1 = (cos(aa*mgo_x(n))+cos(aa*mgo_y(n)))
    funct2 = (cos(aa*(mgo_x(n)+mgo_y(n)))+cos(aa*(mgo_x(n)-mgo_y(n))))
    funct1_dx = sin(aa*mgo_x(n))
    funct2_dx = (-sin(aa*(mgo_x(n)+mgo_y(n)))*aa-sin(aa*(mgo_x(n)-mgo_y(n)))*aa)
    funct1_dy = sin(aa*mgo_y(n))
    funct2_dy = (-sin(aa*(mgo_x(n)+mgo_y(n)))*aa+sin(aa*(mgo_x(n)-mgo_y(n)))*aa)
    !-------------------------
    DO j = 1,3      
       b(:,j)    = ccc(mgo_type,:,j,1) +ccc(mgo_type,:,j,2)*funct1 &
                & +ccc(mgo_type,:,j,3)*funct2
       db_dx(:,j)=-ccc(mgo_type,:,j,2)*funct1_dx*aa+ccc(mgo_type,:,j,3)* funct2_dx
       db_dy(:,j)=-ccc(mgo_type,:,j,2)*funct1_dy*aa+ccc(mgo_type,:,j,3)* funct2_dy
    END DO
    !--------------------------
    exp_term = exp(-mgo_cn(n)/b(:,3))
    !--------------------------
    a_   = b(:,1) +b(:,2)*exp_term
    da_dx= db_dx(:,1) +  db_dx(:,2)*exp_term +  &
            & b(:,2)*exp_term *((-mgo_dcn_dx(n)*b(:,3) +mgo_cn(n)*db_dx(:,3))/b(:,3)**2)
    da_dy= db_dy(:,1) +  db_dy(:,2)*exp_term +  &
            & b(:,2)*exp_term *((-mgo_dcn_dy(n)*b(:,3) +mgo_cn(n)*db_dy(:,3))/b(:,3)**2)
    da_dz=    b(:,2)*exp_term *(-mgo_dcn_dz(n) /b(:,3))
    !--------------------------
    exp_term3= (mgo_z(n)-a_(3))
    exp_term1 = exp(-a_(2)*exp_term3)
     !exp_term2 = exp(-2.d0*a_(2)*exp_term3)
    exp_term2 = exp_term1*exp_term1
    !--------------------------
    ener_mgo = a_(1)     *(exp_term2-2.d0*exp_term1)  
    mgo_fx(n)= da_dx(1) *(exp_term2-2.d0*exp_term1) &
                & +a_(1) *(exp_term2 *(-2.d0*da_dx(2)*exp_term3+2.d0*a_(2)*da_dx(3)) &
                &         -2.d0*exp_term1*(-da_dx(2)*exp_term3 +a_(2)*da_dx(3)) ) 
    mgo_fy(n)= da_dy(1) *(exp_term2-2.d0*exp_term1) &
                & +a_(1) *(exp_term2 *(-2.d0*da_dy(2)*exp_term3+2.d0*a_(2)*da_dy(3)) &
                &         -2.d0*exp_term1*(-da_dy(2)*exp_term3 +a_(2)*da_dy(3)) ) 
    mgo_fz(n)= da_dz(1) *(exp_term2-2.d0*exp_term1) &
                & +a_(1) *(exp_term2 *(-2.d0*da_dz(2)*exp_term3-2.d0*a_(2)*(1.d0-da_dz(3))) &
                &         -2.d0*exp_term1*(-da_dz(2)*exp_term3 -a_(2)*(1.d0-da_dz(3))) ) 

    ! Converting to arete(1) units
    mgo_fx(n) = -mgo_fx(n) *arete(1)
    mgo_fy(n) = -mgo_fy(n) *arete(1)
    mgo_fz(n) = -mgo_fz(n) *arete(1)
    
    ! Energy from the substrate
    ener_sub = ener_sub + ener_mgo
    !--------------------------
 END DO
 ELSE IF (geometry_type == 2) THEN
  DO n = 1, natom
    mgo_type = itype(n)
    !-------------------------
    funct1 = (cos(aa*(mgo_x(n)+mgo_y(n)/s3)) + cos(aa*(mgo_x(n)-mgo_y(n)/s3)) + cos(2*aa*(mgo_y(n))/s3))
    funct2 = (cos(aa*(mgo_x(n)+mgo_y(n)*s3)) + cos(aa*(mgo_x(n)-mgo_y(n)*s3)) + cos(2*aa*(mgo_x(n))))
    funct1_dx = (-aa*(sin(aa*(mgo_x(n)+mgo_y(n)/s3))) - aa*(sin(aa*(mgo_x(n)-mgo_y(n)/s3))))
    funct2_dx=(-aa*sin(aa*(mgo_x(n)+mgo_y(n)*s3))-aa*sin(aa*(mgo_x(n)-mgo_y(n)*s3)) -2*aa*sin(2*aa*(mgo_x(n))))
    funct1_dy=(-aa/s3*(sin(aa*(mgo_x(n)+mgo_y(n)/s3)))+aa/s3*(sin(aa*(mgo_x(n)-mgo_y(n)/s3)))+aa*2/s3*sin((2*aa*(mgo_y(n))/s3)))
    funct2_dy = (-aa*s3*sin(aa*(mgo_x(n)+mgo_y(n)*s3)) + aa/s3*sin(aa*(mgo_x(n)-mgo_y(n)*s3)))
    !-------------------------
    DO j = 1,3      
       b(:,j)    = ccc(mgo_type,:,j,1) +ccc(mgo_type,:,j,2)*funct1 &
                & +ccc(mgo_type,:,j,3)*funct2
       db_dx(:,j)=-ccc(mgo_type,:,j,2)*funct1_dx*aa+ccc(mgo_type,:,j,3)* funct2_dx
       db_dy(:,j)=-ccc(mgo_type,:,j,2)*funct1_dy*aa+ccc(mgo_type,:,j,3)* funct2_dy
    END DO
    !--------------------------
    exp_term = exp(-mgo_cn(n)/b(:,3))
    !--------------------------
    a_   = b(:,1) +b(:,2)*exp_term
    da_dx= db_dx(:,1) +  db_dx(:,2)*exp_term +  &
            & b(:,2)*exp_term *((-mgo_dcn_dx(n)*b(:,3) +mgo_cn(n)*db_dx(:,3))/b(:,3)**2)
    da_dy= db_dy(:,1) +  db_dy(:,2)*exp_term +  &
            & b(:,2)*exp_term *((-mgo_dcn_dy(n)*b(:,3) +mgo_cn(n)*db_dy(:,3))/b(:,3)**2)
    da_dz=    b(:,2)*exp_term *(-mgo_dcn_dz(n) /b(:,3))
    !--------------------------
    exp_term3= (mgo_z(n)-a_(3))
    exp_term1 = exp(-a_(2)*exp_term3)
     !exp_term2 = exp(-2.d0*a_(2)*exp_term3)
    exp_term2 = exp_term1*exp_term1
    !--------------------------
    ener_mgo = a_(1)     *(exp_term2-2.d0*exp_term1)  
    mgo_fx(n)= da_dx(1) *(exp_term2-2.d0*exp_term1) &
                & +a_(1) *(exp_term2 *(-2.d0*da_dx(2)*exp_term3+2.d0*a_(2)*da_dx(3)) &
                &         -2.d0*exp_term1*(-da_dx(2)*exp_term3 +a_(2)*da_dx(3)) ) 
    mgo_fy(n)= da_dy(1) *(exp_term2-2.d0*exp_term1) &
                & +a_(1) *(exp_term2 *(-2.d0*da_dy(2)*exp_term3+2.d0*a_(2)*da_dy(3)) &
                &         -2.d0*exp_term1*(-da_dy(2)*exp_term3 +a_(2)*da_dy(3)) ) 
    mgo_fz(n)= da_dz(1) *(exp_term2-2.d0*exp_term1) &
                & +a_(1) *(exp_term2 *(-2.d0*da_dz(2)*exp_term3-2.d0*a_(2)*(1.d0-da_dz(3))) &
                &         -2.d0*exp_term1*(-da_dz(2)*exp_term3 -a_(2)*(1.d0-da_dz(3))) ) 

    ! Converting to arete(1) units
    mgo_fx(n) = -mgo_fx(n) *arete(1)
    mgo_fy(n) = -mgo_fy(n) *arete(1)
    mgo_fz(n) = -mgo_fz(n) *arete(1)
    
    ! Energy from the substrate
    ener_sub = ener_sub + ener_mgo
    !--------------------------
  END DO
! ELSE IF (geometry_type == 3) THEN
 !
 !
 !
 !
 !
 !
 !
 !
! ELSE WRITE(*,*) '######################################'
!      WRITE(*,*) 'Error:'
!      WRITE(*,*) 'Geometry type not supported'
!      WRITE(*,*) '######################################'
!      WRITE(*,*) 'Check implemented geometries'
!      WRITE(*,*) '######################################'
 END IF

END SUBROUTINE
