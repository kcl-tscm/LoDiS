SUBROUTINE sub_dsq_force

  USE PARACLUSTER
  USE CLUSTER 
  USE POTENTIAL
  USE ENFORCE
  USE DISTANCE
  USE SUBSTRATE

  IMPLICIT NONE
  
  INTEGER :: n,i,j
  INTEGER :: mgo_type
  REAL(8) :: ener_mgo, aa , xx, yy, zz
  REAL(8) :: t1, t2, t3, t1_dx, t1_dy, t2_dx, t2_dy, t3_dx, t3_dy, t3_dz
  REAL(8) :: t4, t8, t4_dx,t4_dy,t4_dz
  !!, t6, t5, t9, t10
  REAL(8), DIMENSION(3)   :: a_, da_dx, da_dy, da_dz
  REAL(8), DIMENSION(3,3) :: b, db_dx, db_dy
  
 !-----------------------------------------------------
 ! Inizialisation of local variables
 ener_mgo = 0.d0
 !------------------------------------------------------
 ! Inizialisation of global variables
 ener_sub = 0.d0 ! This is global  

 aa = mgo_alat ! mgo_alat = 2.*pi*sqrt(2) / amgo in [1/A]  as in mgo_read.f90 
  if (ipas == 1) write(*,*) "substrate_dsq_force:> mgo_alat [1/A]", aa,"aretebim[A]",aretebim
  ! Calculating coordination number
  if (ipas == 1) write(*,*) "substrate_dsq_force:> calling mgo_coord"
 CALL mgo_coord
 if (ipas == 1) write(*,*) "substrate_dsq_force:> after calling mgo_coord"
 ! Now I have cn(i) and and its derivative w.r.t. x, y, z
 ! They are called 
 ! mgo_cn(i)                   adimensional
 ! mgo_dcn_dx(i), mgo_dcn_dy(i), mgo_dcn_dz(i)  [1/L]
 ! In this moment the unit is [1/Angstrom], and can be modified in mgo_coord_no
!
if (ipas==1) WRITE(900,*)'#step ','atom ','xyz [A] ',' CN(n) ','t4= ', 'xyz-forces [eV/arete]','en_sub [eV]'
!if (ipas==1) WRITE(899,*)'#step ','atom ','xyz [A] ', 'a_(1)=  ', 'a_(2)=  ','a_(3)=   '

 DO n = 1, natom
     b(:,:) =0.d0
     ener_mgo = 0.d0
     mgo_fx(n)   = 0.d0
     mgo_fy(n)   = 0.d0
     mgo_fz(n)   = 0.d0
	 !
     mgo_type = itype(n)  !!this is to identify the chemical species of atom n 
     mgo_x(n) = (x(n) + u(n))*aretebim !in A
     mgo_y(n) = (y(n) + v(n))*aretebim !in A
     mgo_z(n) = (z(n) + w(n))*aretebim !in A 
	 !
     xx = aa*mgo_x(n)
     yy = aa*mgo_y(n)
     zz = mgo_z(n) 
    !-------------------------
     t1 = cos(xx) + cos(yy)
     t2 = cos(xx+yy) + cos(xx-yy)
	!!t1 differentiated w.r.t. x and y
    t1_dx = -sin(xx)
    t2_dx = -sin(yy)
	!!t2 differentiated w.r.t x and y
    t1_dy = -sin(xx+yy) - sin(xx-yy)
    t2_dy = -sin(xx+yy) + sin(xx-yy)
    !-------------------------
! formula for the VMG potential
!b(i,j) = c_ij1 + c_ij2*t1 + c_ij3*t2  	
!c_ij1,c_ij2,c_ij3 read in mgo_read
!
DO i=1,3
    DO j = 1,3      
       b(i,j)    = ccc(mgo_type,i,j,1) +ccc(mgo_type,i,j,2)*t1 &
                & +ccc(mgo_type,i,j,3)*t2
				!in 1/A all of them?
       db_dx(i,j)=(ccc(mgo_type,i,j,2)*t1_dx + ccc(mgo_type,i,j,3)*t2_dx)*aa
       db_dy(i,j)=(ccc(mgo_type,i,j,2)*t1_dy + ccc(mgo_type,i,j,3)*t2_dy)*aa
    ENDDO  !!su j
    !--------------------------
	!the a_i coefficients units table
	! a_i(x,y) = b_i1(x,y) + b_i2(x,y)*t3
	! only a_1 is in eV => b(1,1) and b(1,2) are in eV
	! a_2 is in 1/A, as b(2,1) and b(2,2)
	! a_3 is in A, as b(3,1) and b(3,2)
	!b(1,3), b(2,3)and b(3,3) are adim
	!da_dx(1) is in eV/A
	!da_dx(2) is in 1/A^2   = [ccc(212)+ccc(223)]/A
	!da_dx(3) is adim       = [ccc(312)+ccc(323)]/A
	!
    t3 = exp(-mgo_cn(n)/b(i,3)) !!this term is by definition adimension
	!
    t3_dx = ( -mgo_dcn_dx(n)*b(i,3) +mgo_cn(n)*db_dx(i,3) ) /b(i,3)**2
    t3_dy = ( -mgo_dcn_dy(n)*b(i,3) +mgo_cn(n)*db_dy(i,3) ) /b(i,3)**2
    t3_dz = -mgo_dcn_dz(n)/b(i,3) !!from mgo_coord units of the force in 1/A
    !--------------------------
    a_(i)   = b(i,1) +b(i,2)*t3
	!--------------------------
    da_dx(i)= db_dx(i,1) +  db_dx(i,2)*t3 +  &
            			& b(i,2)*t3_dx*t3
			!
    da_dy(i)= db_dy(i,1) +  db_dy(i,2)*t3 +  &
            			& b(i,2)*t3_dy*t3
			!
    da_dz(i)=  b(i,2)*t3_dz*t3
    !--------------------------

!write(910,'(2i6,3f9.4,i5,3f12.5)') ipas, n, xx, yy, zz, i, b(i,1), b(i,2), b(i,3)
!write(911,'(2i6,3f9.4,i5,3f12.5)') ipas, n, xx, yy, zz, i, db_dx(i,1), db_dx(i,2), db_dx(i,3)
!write(912,'(2i6,3f9.4,i5,3f12.5)') ipas, n, xx, yy, zz, i, db_dy(i,1), db_dy(i,2), db_dy(i,3)
!write(913,'(2i6,3f9.4,i5,4f12.5)') ipas, n, xx, yy, zz, i, t3, t3_dx, t3_dy, t3_dz

ENDDO !!su i
!write(899,'(2i6,3f9.4,3f13.6)') ipas, n, xx, yy, zz, a_(1), a_(2), a_(3)
!write(901,'(2i6,3f9.4,a5,4f12.5)') ipas, n, xx, yy, zz, 'CN=', mgo_cn(n), mgo_dcn_dx(n), mgo_dcn_dy(n), mgo_dcn_dz(n)
!write(902,'(2i6,3f9.4,3f13.6)') ipas, n, xx, yy, zz, da_dz(1), da_dz(2), da_dz(3)
!
    ! slave quantities
    t4 = exp(a_(2)*(a_(3)-zz)) !!units a_(3) in A, a_(2) in 1/A, zz is in AA 
    t8 = t4*t4
	!
	!esub = a_(1) * [exp(2a_(2)*(a_(3)-z)) - 2exp(a_(2)*(a(3)-z)) ]	=
	!     =a_(1) * [t8 - 2*exp(t4)]
    ener_mgo = a_(1)*(t8-2.d0*t4)  ! in eV
    ! Energy from the substrate
    ener_sub = ener_sub + ener_mgo
	!
	!calculating the force on atom n-th
    t4_dx = da_dx(2)*(a_(3)-zz) + a_(2)*da_dx(3)
    t4_dy = da_dy(2)*(a_(3)-zz) + a_(2)*da_dy(3)
    t4_dz = da_dz(2)*(a_(3)-zz) + a_(2)*(da_dz(3)-1.d0)
    !
    t8= t4*t4
	!t8_dx = 2*t4*t4_dx and similarly on y, and z
    !--------------------------
	!force due to the substrate
	!differentiate -grad(esub)/dx
	!and similarly on y  
	!
    mgo_fx(n) = da_dx(1)*(t8-2.d0*t4) &
		     & + 2.d0*a_(1)*t4_dx*(t8-t4)
	!mgo_fx(n) = da_dx(1)*(t8-2.d0*t4) + 2.d0 
	!
    mgo_fy(n) = da_dy(1)*(t8-2.d0*t4) &
		     & + 2.d0*a_(1)*t4_dy*(t8-t4)
	!
    mgo_fz(n) = da_dz(1)*(t8-2.d0*t4) &
	         & + 2.d0*a_(1)*t4_dz*(t8-t4)
	!
	!LP's implementation	  
 !               & a_(1) *(t9 *(-2.d0*da_dx(2)*t10+2.d0*a_(2)*da_dx(3)) &
 !               &         -2.d0*t8*(-da_dx(2)*t10 +a_(2)*da_dx(3)) ) 
 !   mgo_fy(n)= da_dy(1) *(t9-2.d0*t8) &
 !               & +a_(1) *(t9 *(-2.d0*da_dy(2)*t10+2.d0*a_(2)*da_dy(3)) &
 !               &         -2.d0*t8*(-da_dy(2)*t10 +a_(2)*da_dy(3)) ) 
 !   mgo_fz(n)= da_dz(1) *(t9-2.d0*t8) &
 !               & +a_(1) *(t9 *(-2.d0*da_dz(2)*t10-2.d0*a_(2)*(1.d0-da_dz(3)))&
 !               &         -2.d0*t8*(-da_dz(2)*t10 -a_(2)*(1.d0-da_dz(3))) ) 

    ! force due to the substrate are in eV/arete 
	!adding the minus from the force to potential definition
	!to be sure that the differentiate of the CN doesn't take it into account
    mgo_fx(n) = -mgo_fx(n)
    mgo_fy(n) = -mgo_fy(n)
    mgo_fz(n) = -mgo_fz(n)
    !
    !
    mgo_fx(n) = mgo_fx(n) * arete(1)
    mgo_fy(n) = mgo_fy(n) * arete(1)
    mgo_fz(n) = mgo_fz(n) * arete(1)
	!
    if(MOD(ipas,scrivo)==0) WRITE(900,'(2i6,5f9.4,3(f12.5,1x),f12.6)') &
	& ipas, n, xx, yy, zz, mgo_cn(n), t4, mgo_fx(n), mgo_fy(n), mgo_fz(n), ener_mgo 
	!
 ENDDO
 
 !!if(MOD(ipas,scrivo)==0) write(*,*) "substrate_dsq_force:> at step",ipas,'substrate energy;=',ener_sub,"[eV]"
END SUBROUTINE sub_dsq_force
