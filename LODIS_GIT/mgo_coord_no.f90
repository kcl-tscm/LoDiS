SUBROUTINE mgo_coord
 ! Luca
 ! Calculating the analytic coordination number for MgO

 !USE OMP_LIB                                                 !>>> CPU_TIME
 USE PARACLUSTER
 USE CLUSTER
 USE POTENTIAL
 USE ENFORCE
 USE DISTANCE
 USE SUBSTRATE

 IMPLICIT NONE

 integer               :: i, j, k
 REAL(8)      		   :: comfactor, gamma_sub, gamma_sub_1, gamma_sub_2
 REAL(8), ALLOCATABLE  :: dij(:,:), df_dx(:,:), df_dy(:,:), df_dz(:,:)
 !REAL(8) :: t_start, t_end                                   !>>> CPU_TIME
 !t_start = OMP_GET_WTIME()                                   !>>> CPU_TIME

IF (ipas.EQ.1) WRITE(*,*) 'mgo_coord:> first calculation of CN for MgO substrate'
IF (ipas==1) WRITE(*,*) 'mgo_coord:> aretebim [A]',aretebim
ALLOCATE(dij(natom,natom), df_dx(natom,natom),df_dy(natom,natom),df_dz(natom,natom))
!
df_dx(:,:) = 0.d0
df_dy(:,:) = 0.d0
df_dz(:,:) = 0.d0 
!
!coordination number using a smooth Fermi function as
!CN(i) = sum_{j\=i} 1 / [1+exp^alpha(rij/rc-1)]
!beta_sub = alpha_sub/ (rc_sub/aretebim) 
!cn(i) = sum_{j/=i} 1 / 1+gamma_sub
!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(i, j), SCHEDULE(STATIC)
   DO i = 1, natom
	   dij(:,i) = 0.d0
      DO j = 1, natom
         ! alpha/rc*r(i,j) -alpha
         dij(j,i) = (beta_sub * pair_dist(j,i) ) - alpha_sub
      END DO
   END DO
!$OMP END DO
if(ipas==1) write(*,*) "mgo_coord:>alpha_sub and beta_sub",alpha_sub, beta_sub
!!write(*,*) "mgo_coord:> dist, dij(1,natom)",ipas,pair_dist(1,natom),dij(1,natom)

!$OMP DO PRIVATE(i), SCHEDULE(STATIC)
   ! Elements on the diagonal
   DO i = 1, natom
      mgo_s(i,i) = 0.d0
	  mgo_dS(i,j) = 0.d0
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j,gamma_sub), SCHEDULE(DYNAMIC,20)
   DO i = 1, natom-1
      DO j = i+1, natom
            gamma_sub = exp(dij(j,i)) 
            gamma_sub_1 = gamma_sub + 1.d0
			gamma_sub_2 = gamma_sub /(gamma_sub_1 * gamma_sub_1)
			!! for large cluster the pair distance may grow up and become very
			!! large (too large) for large value th eaddition of 1 change little or nothing, gamma_sub_1 ~ gamma_sub -> 
			!gamma_sub_2 ~ 1/ gamma_sub
		    if(dij(j,i)>=80.d0) gamma_sub_2 = 1.d0 / gamma_sub 
           ! The sum of these s(i,j) is the analytic coord no.
            mgo_s(j,i) = 1.d0/(gamma_sub + 1.d0)
            mgo_s(i,j) = mgo_s(j,i)

            ! Derivative of s w.r.t r(i,j)
            ! Use beta_sub or beta_sub_angstrom to have arete(1) or Angstrom as units
            mgo_dS(j,i) = -beta_sub_angstrom * gamma_sub_2
            mgo_dS(i,j) = mgo_dS(j,i)

      END DO ! on j
	  !!write(500,'(a25,2i4,4e13.5)') 'mgo_coord:> atom j=60',ipas, i, dij(60,i), gamma_sub, gamma_sub_1, mgo_dS(60,i)
   END DO ! on i
!$OMP END DO
!!write(*,*) "mgo_coord:> mgo_s(1,natom) mgo_s(1,2)",mgo_s(1,natom),mgo_s(1,2)

!$OMP DO PRIVATE(i), SCHEDULE(STATIC)
   ! LP: Calculating the sum of s(i,j), which is mgo_CN(i), global.
   ! This is the coord num of any single atom i
   ! No problem because s(i,i) is always zero
   DO i = 1, natom
      mgo_cn(i) = SUM(mgo_s(1:natom,i))
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j,comfactor), SCHEDULE(STATIC)
   !-------------------------------------------------------------------------
   ! Calculating the sum on j of derivatives of s
   ! with respect to coordinates of every atoms
   DO i = 1, natom
      DO j = 1, natom
         IF (j .NE. i) THEN
            ! Calculating the derivative of coordination number to atom i
            ! multiplied for the derivative of distance between atoms i and j
            ! with respect to position atom i
            !
			!xyz_dist(1,j,i) = x(i)-x(j)
			!pair_dist(j,i) = sqrt((x(i)-x(j))^2 + on y + on z)
			!comfactor doesn't depend on i,j order
			!
            comfactor  = mgo_dS(j,i)/ pair_dist(j,i)
            df_dx(j,i) = comfactor * xyz_dist(1,j,i) ! order of indexes depends on the definition of xyz_dist
            df_dy(j,i) = comfactor * xyz_dist(2,j,i) ! order of indexes depends on the definition of xyz_dist
            df_dz(j,i) = comfactor * xyz_dist(3,j,i) ! order of indexes depends on the definition of xyz_dist
         ELSE
            df_dx(j,i) = 0.d0
            df_dy(j,i) = 0.d0
            df_dz(j,i) = 0.d0
         END IF
      END DO
	  !!write(501,*) 'mgo_coord:> force atom j=60',ipas, i, mgo_dS(60,i), df_dx(60,i)
   END DO
!$OMP END DO

!$OMP SECTIONS PRIVATE(i)
!$OMP SECTION
Do k = 1,natom
    ! No factor 2 here: this is derivative of cn(i) w.r.o. x(i), not total cn w.r.o. x(i)
  !sum over all j/=k and subtract df_dx(k,i)
	mgo_dcn_dx(k) = SUM(df_dx(1:natom, k))
   DO i = 1, natom 
	   if(i/=k)  mgo_dcn_dx(k) = mgo_dcn_dx(k) - df_dx(k,i) 
   END DO     
 ENDDO     
!$OMP SECTION
   DO k = 1, natom 
     mgo_dcn_dy(k) = SUM(df_dy(1:natom, k)) 
    DO i = 1, natom 
       if(i/=k)  mgo_dcn_dy(k) = mgo_dcn_dy(k) - df_dy(k,i) 
     END DO     
   END DO 
!$OMP SECTION
DO k = 1, natom 
    mgo_dcn_dz(k) = SUM(df_dz(1:natom, k)) 
   DO i = 1, natom 
     if(i/=k)  mgo_dcn_dz(k) = mgo_dcn_dz(k) - df_dz(k,i) 
   END DO  
 END DO 
!$OMP END SECTIONS
!$OMP END PARALLEL
!
!LP's implementation
!   DO i = 1, natom
      ! No factor 2 here: this is derivative of cn(i) w.r.o. x(i), not total cn w.r.o. x(i)
!      mgo_dcn_dx(i) = SUM(df_dx(1:natom, i))
!   END DO
!!$OMP SECTION
!!same on y and z
!DO i = 1, natom
!   ! No factor 2 here: this is derivative of cn(i) w.r.o. x(i), not total cn w.r.o. x(i)
!   mgo_dcn_dy(i) = SUM(df_dy(1:natom, i))
!END DO
!!$OMP SECTION
!DO i = 1, natom
   ! No factor 2 here: this is derivative of cn(i) w.r.o. x(i), not total cn w.r.o. x(i)
 !  mgo_dcn_dz(i) = SUM(df_dz(1:natom, i))
 !END DO
!!$OMP SECTION
!
!!$OMP END SECTIONS
!!$OMP END PARALLEL
   !---------------------------------------------------------------------------
   !t_end = OMP_GET_WTIME()                                     !>>> CPU_TIME
   !! Writing the cpu time for the subroutine                   !>>> CPU_TIME
   !WRITE(500,'(1I10, 4X, 1F10.8)') ipas, t_end -t_start        !>>> CPU_TIME
END SUBROUTINE mgo_coord