SUBROUTINE ENVCN2
   
 USE PARACLUSTER  
 USE CLUSTER   
 USE POTENTIAL
 USE ENFORCE
 USE DISTANCE
 USE ENVIRONMENT

  IMPLICIT NONE

  integer               :: i, j
  double precision      :: gamma_env
  REAL(8), DIMENSION(nsiz,nsiz)   :: dij, df_dx, df_dy, df_dz

 !write(203,*) ipas, beta_env_angstrom, rc_env, beta_env

   DO i = 1, natom
      DO j = 1, natom 
         dij(j,i) = (beta_env * pair_dist(j,i) ) - alpha_env
         !write(202,*) ipas, pair_dist(j,i), dij(j,i)                                      
      END DO
   END DO

   ! Elements on the diagonal
   DO i = 1, natom
      env_s(i,i) = 0.d0
   END DO

   DO i = 1, natom-1
      DO j = i+1, natom      
            if (abs(dij(j,i)) .le. 16.d0) then
            !write(185,*) ipas, i, j, dij(j,i)
            gamma_env = exp(dij(j,i))
            env_s(j,i) = 1.d0/(1.d0 + gamma_env)  
            env_s(i,j) = env_s(j,i)
            env_ds(j,i) = -beta_env_angstrom * gamma_env / (gamma_env+1)**2
            env_ds(i,j) = env_ds(j,i)
            else

            gamma_env  = 10000
            env_s(j,i) = 0.d0
            env_s(i,j) = env_s(j,i)
            env_ds(j,i) = 0.d0
            env_ds(i,j) = env_ds(j,i)
            end if

!      write(183,*) ipas, i, j, dij(j,i), gamma_env, env_s(j,i)

      END DO ! on j 
   END DO ! on i

   DO i = 1, natom 
      env_cn(i) = SUM(env_s(1:natom,i))
      !WRITE(184,*) ipas, i, env_cn(i)
      IF (env_cn(i) .le. 5.d0) THEN
         env_cn(i) = 5.d0
      END IF
      IF (env_cn(i) .gt. 12.d0) THEN
         env_cn(i) = 12.d0
      END IF
   END DO

   DO i = 1, natom 
      DO j = 1, natom 
         IF (j .NE. i) THEN
            env_comfactor  = env_dS(j,i)/ pair_dist(j,i)
            df_dx(j,i) = env_comfactor * xyz_dist(1,j,i) ! order of indexes depends on the definition of xyz_dist
            df_dy(j,i) = env_comfactor * xyz_dist(2,j,i) ! order of indexes depends on the definition of xyz_dist
            df_dz(j,i) = env_comfactor * xyz_dist(3,j,i) ! order of indexes depends on the definition of xyz_dist
         ELSE
            df_dx(j,i) = 0.d0
            df_dy(j,i) = 0.d0
            df_dz(j,i) = 0.d0 
         END IF
      END DO
   END DO     
  
   DO i = 1, natom 
      ! No factor 2 here: this is derivative of cn(i) w.r.o. x(i), not total cn w.r.o. x(i)
      env_dcn_dx(i) = SUM(df_dx(1:natom, i)) 
      env_dcn_dy(i) = SUM(df_dy(1:natom, i))  
      env_dcn_dz(i) = SUM(df_dz(1:natom, i))
!      write(204,*) ipas, i, env_dcn_dx(i), env_dcn_dy(i), env_dcn_dz(i)
   END DO 

  
END SUBROUTINE ENVCN2

