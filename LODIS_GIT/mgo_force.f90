SUBROUTINE mgo_force

  USE PARACLUSTER
  USE CLUSTER 
  USE POTENTIAL
  USE ENFORCE
  USE DISTANCE
  USE SUBSTRATE

  IMPLICIT NONE
  
  INTEGER :: n,i,j
  INTEGER :: mgo_type
  REAL(8) :: ener_mgo, aa 
  REAL(8) :: t1, t2, t3, t4, t5, t6, t8, t9, t10
  REAL(8), DIMENSION(3)   :: t7

  REAL(8), DIMENSION(3)   :: a_, da_dx, da_dy, da_dz
  REAL(8), DIMENSION(3,3) :: b, db_dx, db_dy
  
 !-----------------------------------------------------
 ! Inizialisation of local variables
 ener_mgo = 0.d0
 !------------------------------------------------------
 ! Inizialisation of global variables
 ener_sub = 0.d0 ! This is global  
 mgo_fx   = 0.d0
 mgo_fy   = 0.d0
 mgo_fz   = 0.d0

 ! Calculating coordination number
 CALL mgo_coord
 ! Now I have cn(i) and and its derivative w.r.t. x, y, z
 ! They are called 
 ! mgo_cn(i)                                    adimensional
 ! mgo_dcn_dx(i), mgo_dcn_dy(i), mgo_dcn_dz(i)  [1/L]
 ! In this moment the unit is [1/Angstrom], and can be modified in mgo_coord_no

 ! Normalised positions  (in Angstrom and then *mgo_alat)  
 mgo_x = (x + u) *arete(1) ! [A]
 mgo_y = (y + v) *arete(1) ! [A]
 mgo_z = (z + w) *arete(1) ! [A]

 aa = mgo_alat             ! [1/A]

 DO n = 1, natom
    mgo_type = itype(n)
    !-------------------------
    t1 = (cos(aa*mgo_x(n))+cos(aa*mgo_y(n)))
    t2 = (cos(aa*(mgo_x(n)+mgo_y(n)))+cos(aa*(mgo_x(n)-mgo_y(n))))
    t3 = sin(aa*mgo_x(n))
    t4 = (-sin(aa*(mgo_x(n)+mgo_y(n)))*aa-sin(aa*(mgo_x(n)-mgo_y(n)))*aa)
    t5 = sin(aa*mgo_y(n))
    t6 = (-sin(aa*(mgo_x(n)+mgo_y(n)))*aa+sin(aa*(mgo_x(n)-mgo_y(n)))*aa)
    !-------------------------
    DO j = 1,3      
       b(:,j)    = ccc(mgo_type,:,j,1) +ccc(mgo_type,:,j,2)*t1 &
                & +ccc(mgo_type,:,j,3)*t2
       db_dx(:,j)=-ccc(mgo_type,:,j,2)*t3*aa+ccc(mgo_type,:,j,3)* t4
       db_dy(:,j)=-ccc(mgo_type,:,j,2)*t5*aa+ccc(mgo_type,:,j,3)* t6
    ENDDO
    !--------------------------
    t7 = exp(-mgo_cn(n)/b(:,3))
    !--------------------------
    a_   = b(:,1) +b(:,2)*t7
    da_dx= db_dx(:,1) +  db_dx(:,2)*t7 +  &
            & b(:,2)*t7 *((-mgo_dcn_dx(n)*b(:,3) +mgo_cn(n)*db_dx(:,3))/b(:,3)**2)
    da_dy= db_dy(:,1) +  db_dy(:,2)*t7 +  &
            & b(:,2)*t7 *((-mgo_dcn_dy(n)*b(:,3) +mgo_cn(n)*db_dy(:,3))/b(:,3)**2)
    da_dz=    b(:,2)*t7 *(-mgo_dcn_dz(n) /b(:,3))
    !--------------------------
    t10= (mgo_z(n)-a_(3))
    t8 = exp(-a_(2)*t10)
    !t9 = exp(-2.d0*a_(2)*t10)
    t9 = t8*t8 
    !--------------------------
    ener_mgo = a_(1)     *(t9-2.d0*t8)  
    mgo_fx(n)= da_dx(1) *(t9-2.d0*t8) &
                & +a_(1) *(t9 *(-2.d0*da_dx(2)*t10+2.d0*a_(2)*da_dx(3)) &
                &         -2.d0*t8*(-da_dx(2)*t10 +a_(2)*da_dx(3)) ) 
    mgo_fy(n)= da_dy(1) *(t9-2.d0*t8) &
                & +a_(1) *(t9 *(-2.d0*da_dy(2)*t10+2.d0*a_(2)*da_dy(3)) &
                &         -2.d0*t8*(-da_dy(2)*t10 +a_(2)*da_dy(3)) ) 
    mgo_fz(n)= da_dz(1) *(t9-2.d0*t8) &
                & +a_(1) *(t9 *(-2.d0*da_dz(2)*t10-2.d0*a_(2)*(1.d0-da_dz(3))) &
                &         -2.d0*t8*(-da_dz(2)*t10 -a_(2)*(1.d0-da_dz(3))) ) 

    ! Converting to arete(1) units
    mgo_fx(n) = -mgo_fx(n) *arete(1)
    mgo_fy(n) = -mgo_fy(n) *arete(1)
    mgo_fz(n) = -mgo_fz(n) *arete(1)
    
    ! Energy from the substrate
    ener_sub = ener_sub + ener_mgo   

    !WRITE(900,*) ipas, n, mgo_fx(n), mgo_fy(n), mgo_fz(n)
    !WRITE(900,*) ipas, n, mgo_cn(n), mgo_dcn_dx(n), mgo_dcn_dy(n), mgo_dcn_dz(n)
 ENDDO
END SUBROUTINE mgo_force
