SUBROUTINE mgo_coord_no_mot
 ! Luca 9Aug14 - 1Nov14
 ! Calculating the anisotropic analytic coordination number for MgO
 ! as well as the isotropic coordination number for MgO
     
 !USE OMP_LIB                                                 !>>> CPU_TIME
 USE PARACLUSTER  
 USE CLUSTER   
 USE POTENTIAL
 USE ENFORCE
 USE DISTANCE
 USE SUBSTRATE
   
 IMPLICIT NONE
   
 integer               :: i, j
 double precision      :: comfactor, gamma_sub
 REAL(8), DIMENSION(nsiz,nsiz)   :: dij, df_dx, df_dy, df_dz
 REAL(8), DIMENSION(nsiz,nsiz)   :: cos_theta, df_mot_dx, df_mot_dy, df_mot_dz
 !REAL(8) :: t_start, t_end                                   !>>> CPU_TIME
 !t_start = OMP_GET_WTIME()                                   !>>> CPU_TIME
   
IF (ipas.EQ.1) WRITE(*,*) 'mgo_coord_no_mot> first calculation of isotropic and anisotropic CN for MgO substrate'

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(i), SCHEDULE(STATIC) 
   ! Elements on the diagonal
   DO i = 1, natom
      cos_theta(i,i) = 0.d0
   END DO
!$OMP END DO 

!$OMP DO PRIVATE(i, j), SCHEDULE(STATIC) 
   DO i = 1, natom-1
      DO j = i+1, natom 
         !====================================================
         ! Note that here the cosine of theta is always from
         ! 0 to 1. Matrix cos_theta is antisymmetric
         !==================================================== 
         cos_theta(j,i) = xyz_dist(3,j,i)/ pair_dist(j,i)
         cos_theta(i,j) = -cos_theta(j,i)
      END DO
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i, j), SCHEDULE(STATIC) 
   DO i = 1, natom
      DO j = 1, natom                                       
         ! alpha/rc*r(i,j) -alpha
         dij(j,i) = (beta_sub * pair_dist(j,i) ) - alpha_sub
      END DO
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i), SCHEDULE(STATIC) 
   ! Elements on the diagonal
   DO i = 1, natom
      mgo_s(i,i) = 0.d0
   END DO
!$OMP END DO  

!$OMP DO PRIVATE(i,j,gamma_sub), SCHEDULE(DYNAMIC,20)
   DO i = 1, natom-1
      DO j = i+1, natom      
         
            gamma_sub = exp(dij(j,i))
          
            ! The sum of these s(i,j) is the analytic coord no.
            mgo_s(j,i) = 1.d0/(1.d0 + gamma_sub)  
            mgo_s(i,j) = mgo_s(j,i)
    
            ! Derivative of s w.r.t r(i,j)  
            ! Use beta_sub or beta_sub_angstrom to have arete(1) or Angstrom as units
            mgo_dS(j,i) = -beta_sub_angstrom * gamma_sub / (gamma_sub+1.d0)**2
            mgo_dS(i,j) = mgo_dS(j,i)
         
      END DO ! on j 
   END DO ! on i
!$OMP END DO

!$OMP DO PRIVATE(i), SCHEDULE(STATIC)
  ! mgo_s_mot is antisymmetric
   DO i = 1, natom     
      mgo_s_mot(1:natom,i) = mgo_s(1:natom,i) * cos_theta(1:natom,i)
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i), SCHEDULE(STATIC)
   ! LP: Calculating the sum of s(i,j), which is mgo_cn(i), global.
   ! This is the coord num of any single atom i
   ! No problem because s(i,i) is always zero
   ! The sum can be done on any index because mgo_s is a symmetric matrix
   DO i = 1, natom 
      mgo_cn(i) = SUM(mgo_s(1:natom,i))
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i), SCHEDULE(STATIC)
   ! LP: Calculating the sum of s(i,j), which is mgo_CN(i), global.
   ! This is the coord num of any single atom i
   ! No problem because s(i,i) is always zero
   ! This time the sum should be done on the second index,
   ! because mgo_s_mot is antisymmetric, but if I
   ! multiply by -1 I can sum on the first index
   DO i = 1, natom 
      mgo_cn_mot(i) = -SUM(mgo_s_mot(1:natom,i))
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j,comfactor), SCHEDULE(STATIC)
   !-------------------------------------------------------------------------
   ! ISOTROPIC CN 
   ! Calculating the sum on j of derivatives of s
   ! with respect to coordinates of every atoms   
   DO i = 1, natom 
      DO j = 1, natom 
         IF (j .NE. i) THEN
            ! Calculating the derivative of coordination number with
            ! respect to distance between atoms i and j, multiplied for 
            ! the derivative of distance between atoms i and j
            ! with respect to position of atom i
            comfactor  = mgo_dS(j,i)/ pair_dist(j,i)
            df_dx(j,i) = comfactor * xyz_dist(1,j,i) ! Here there should be a minus, but later on [*] I will omit another one
            df_dy(j,i) = comfactor * xyz_dist(2,j,i) ! Here there should be a minus, but later on [*] I will omit another one
            df_dz(j,i) = comfactor * xyz_dist(3,j,i) ! Here there should be a minus, but later on [*] I will omit another one
         ELSE
            df_dx(j,i) = 0.d0
            df_dy(j,i) = 0.d0
            df_dz(j,i) = 0.d0 
         END IF
      END DO
   END DO     
!$OMP END DO

!$OMP DO PRIVATE(i,j,comfactor), SCHEDULE(STATIC)
   !-------------------------------------------------------------------------
   ! ANISOTROPIC CN 
   ! Calculating the sum on j of derivatives of s
   ! with respect to coordinates of every atoms   
   DO i = 1, natom 
      DO j = 1, natom 
         IF (j .NE. i) THEN
            ! Calculating the derivative of coordination number with
            ! respect to distance between atoms i and j, multiplied for 
            ! the derivative of distance between atoms i and j
            ! with respect to position of atom i
            comfactor = xyz_dist(3,j,i) /pair_dist(j,i)**3  

            df_mot_dx(j,i) = df_dx(j,i) *cos_theta(j,i) + mgo_s(j,i)*(comfactor *xyz_dist(1,j,i)) 
            df_mot_dy(j,i) = df_dy(j,i) *cos_theta(j,i) + mgo_s(j,i)*(comfactor *xyz_dist(2,j,i)) 
            df_mot_dz(j,i) = df_dz(j,i) *cos_theta(j,i) + mgo_s(j,i)*(comfactor *xyz_dist(3,j,i) -1.d0/pair_dist(j,i))
         ELSE
            df_mot_dx(j,i) = 0.d0
            df_mot_dy(j,i) = 0.d0
            df_mot_dz(j,i) = 0.d0 
         END IF
      END DO
   END DO     
!$OMP END DO
  
!$OMP SECTIONS PRIVATE(i)
!$OMP SECTION 
   DO i = 1, natom 
      ! No factor 2 here: this is derivative of cn(i) w.r.o. x(i), not total cn w.r.o. x(i).
      ! df_dx is antisymmetric. The sum should be done on the second index, but a minus has
      ! been introduced previously (in [*]).
      mgo_dcn_dx(i) = SUM(df_dx(1:natom, i)) 
   END DO        
!$OMP SECTION
   DO i = 1, natom 
      mgo_dcn_dy(i) = SUM(df_dy(1:natom, i)) 
   END DO 
!$OMP SECTION
   DO i = 1, natom 
      mgo_dcn_dz(i) = SUM(df_dz(1:natom, i)) 
   END DO 
!-------------------------------------------
!$OMP SECTION 
   DO i = 1, natom 
      ! No factor 2 here: this is derivative of cn(i) w.r.o. x(i), not total cn w.r.o. x(i)
      ! df_mot_dx is symmetric, no matter what is the index of the sum.
      mgo_dcn_mot_dx(i) = SUM(df_mot_dx(1:natom, i)) 
   END DO        
!$OMP SECTION
   DO i = 1, natom 
      mgo_dcn_mot_dy(i) = SUM(df_mot_dy(1:natom, i)) 
   END DO 
!$OMP SECTION
   DO i = 1, natom 
      mgo_dcn_mot_dz(i) = SUM(df_mot_dz(1:natom, i)) 
   END DO 
!$OMP END SECTIONS
!$OMP END PARALLEL
   !---------------------------------------------------------------------------   
   !t_end = OMP_GET_WTIME()                                     !>>> CPU_TIME
   !! Writing the cpu time for the subroutine                   !>>> CPU_TIME 
   !WRITE(500,'(1I10, 4X, 1F10.8)') ipas, t_end -t_start        !>>> CPU_TIME
END SUBROUTINE mgo_coord_no_mot

