SUBROUTINE mgo_force_mot

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
  REAL(8), DIMENSION(4)   :: t7

  REAL(8), DIMENSION(4)   :: a_, da_dx, da_dy, da_dz
  REAL(8), DIMENSION(4,3) :: b, db_dx, db_dy
  
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
 CALL mgo_coord_no_mot
 ! Now I have cn(i) and and its derivative w.r.t. x, y, z
 ! They are called 
 ! mgo_cn(i)                                    adimensional
 ! mgo_dcn_dx(i), mgo_dcn_dy(i), mgo_dcn_dz(i)  [1/L]
 ! In this moment the unit is [1/Angstrom], and can be modified in mgo_coord_no

 ! And also: mgo_cn_mot(i)
 ! mgo_dcn_mot_dx(i), mgo_dcn_mot_dy(i), mgo_dcn_mot_dz(i)

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
       b(1:3,j)    = ccc(mgo_type,1:3,j,1) +ccc(mgo_type,1:3,j,2)*t1 &
                & +ccc(mgo_type,1:3,j,3)*t2
       db_dx(1:3,j)=-ccc(mgo_type,1:3,j,2)*t3*aa+ccc(mgo_type,1:3,j,3)* t4
       db_dy(1:3,j)=-ccc(mgo_type,1:3,j,2)*t5*aa+ccc(mgo_type,1:3,j,3)* t6 
    ENDDO  
    DO j = 1,2      
       b(4,j)    = ccc4(mgo_type,j,1) +ccc4(mgo_type,j,2)*t1 &
                & +ccc4(mgo_type,j,3)*t2
       db_dx(4,j)=-ccc4(mgo_type,j,2)*t3*aa+ccc4(mgo_type,j,3)* t4
       db_dy(4,j)=-ccc4(mgo_type,j,2)*t5*aa+ccc4(mgo_type,j,3)* t6 
    ENDDO        
    !--------------------------
    t7(1:3) =  exp(-mgo_cn(n)/b(1:3,3))
    t7(4)   =  exp( mgo_cn_mot(n)*b(4,2))
    !--------------------------
    a_(1:3)   = b(1:3,1) +b(1:3,2)*t7(1:3)
    a_(4)     = 1.d0     +b(4  ,1)*t7(4) 
    da_dx(1:3)= db_dx(1:3,1) +  db_dx(1:3,2)*t7(1:3) +  &
            & b(1:3,2)*t7(1:3) *((-mgo_dcn_dx(n)*b(1:3,3) +mgo_cn(n)*db_dx(1:3,3))/b(1:3,3)**2)
    da_dy(1:3)= db_dy(1:3,1) +  db_dy(1:3,2)*t7(1:3) +  &
            & b(1:3,2)*t7(1:3) *((-mgo_dcn_dy(n)*b(1:3,3) +mgo_cn(n)*db_dy(1:3,3))/b(1:3,3)**2)
    da_dz(1:3)= b(1:3,2)*t7(1:3) *(-mgo_dcn_dz(n) /b(1:3,3))
    da_dx(4) = db_dx(4,1)*t7(4) + b(4,1)*t7(4)*(db_dx(4,2)*mgo_cn_mot(n) + b(4,2)* mgo_dcn_mot_dx(n))
    da_dy(4) = db_dy(4,1)*t7(4) + b(4,1)*t7(4)*(db_dy(4,2)*mgo_cn_mot(n) + b(4,2)* mgo_dcn_mot_dy(n))
    da_dz(4) = b(4,1)*t7(4)*(b(4,2)* mgo_dcn_mot_dz(n))
    !--------------------------
    t10= (mgo_z(n)-a_(3))
    t8 = exp(-a_(2)*t10)
    !t9 = exp(-2.d0*a_(2)*t10)
    t9 = t8*t8 
    !--------------------------
    ener_mgo = a_(1)     *(t9-2.d0*a_(4)*t8)  
    mgo_fx(n)= da_dx(1) *(t9-2.d0*a_(4)*t8) &
                & +a_(1) *(t9 *(-2.d0*da_dx(2)*t10+2.d0*a_(2)*da_dx(3)) &
                &         -2.d0*a_(4)*t8*(-da_dx(2)*t10 +a_(2)*da_dx(3))&
                &         -2.d0*da_dx(4)*t8 ) 
    mgo_fy(n)= da_dy(1) *(t9-2.d0*a_(4)*t8) &
                & +a_(1) *(t9 *(-2.d0*da_dy(2)*t10+2.d0*a_(2)*da_dy(3)) &
                &         -2.d0*a_(4)*t8*(-da_dy(2)*t10 +a_(2)*da_dy(3))&
                &         -2.d0*da_dy(4)*t8 ) 
    mgo_fz(n)= da_dz(1) *(t9-2.d0*a_(4)*t8) &
                & +a_(1) *(t9 *(-2.d0*da_dz(2)*t10-2.d0*a_(2)*(1.d0-da_dz(3))) &
                &         -2.d0*a_(4)*t8*(-da_dz(2)*t10 -a_(2)*(1.d0-da_dz(3)))&
                &         -2.d0*da_dz(4)*t8 ) 

    ! Converting to arete(1) units
    mgo_fx(n) = -mgo_fx(n) *arete(1)
    mgo_fy(n) = -mgo_fy(n) *arete(1)
    mgo_fz(n) = -mgo_fz(n) *arete(1)
    
    ! Energy from the substrate
    ener_sub = ener_sub + ener_mgo   

    !WRITE(900,*) ipas, n, mgo_fx(n), mgo_fy(n), mgo_fz(n)
    !WRITE(901,*) ipas, n, mgo_cn(n), mgo_dcn_dx(n), mgo_dcn_dy(n), mgo_dcn_dz(n)
    !WRITE(902,*) ipas, n, mgo_cn_mot(n), mgo_dcn_mot_dx(n), mgo_dcn_mot_dy(n), mgo_dcn_mot_dz(n)
 ENDDO
END SUBROUTINE mgo_force_mot
