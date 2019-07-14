SUBROUTINE ENVCN
!======================================
! Environment CNi
!======================================
! Calculating the analytic CNi for a cluster with an interacting evnvironment.
! To use the Fermi formulation of the CNi use 2nd version of the subroutine, ENVCN2.
!======================================
     
   USE PARACLUSTER  
   USE CLUSTER   
   USE POTENTIAL
   USE ENFORCE
   USE DISTANCE
   USE ENVIRONMENT
   
   IMPLICIT NONE
   
   integer               :: i, j, env_mminusn_pwr
   real(8), dimension(3) :: env_rij0 ! 3 possible cases: atoms A-A, B-B, A-B
   REAL(8), DIMENSION(nsiz,nsiz)   :: env_dij, env_df_dx, env_df_dy, env_df_dz

   
   if (ipas.eq.1) write(*,*) 'cnum> first calculation of sigmoid function (CN)'

   !---------------------------------------------------
   env_mminusn_pwr = env_m_pwr - env_n_pwr
   
   ! Evaluating the values of rij0 parameter in dependence of atom pair kind
   env_rij0 = dist * env_rzero * dsqrt(2.d0)

   !---------------------------------------------------
   ! Calculating the analytic coordination number
   ! The sum of these s(i,j) is the analytic coord no.

   ! Default Case: (dij(j,i) .LE. 0.d0)
   env_s  = 1.d0
   ! Derivative of s
   env_dS = 0.d0

   DO i = 1, natom
      DO j = 1, natom                                       
         ! Distance between atom i and atom j
         ! minus the ref distance for the couple of atoms
         env_dij(j,i) = pair_dist(j,i) - dist(pairkindmat(j,i))
      END DO
   END DO

   ! Elements on the diagonal
   DO i = 1, natom
      env_s(i,i) = 0.d0
   END DO

   ! Case (dij(j,i) .GT. 0.d0)
   DO i = 1, natom-1
      DO j = i+1, natom      
         IF (env_dij(j,i) .GT. 0.d0) THEN
            env_ratio = env_dij(j,i)/env_rij0(pairkindmat(j,i))
            ! I need to rise to (X_pwr-1) for the derivetive  
            env_ratio_ton = env_ratio**(env_n_pwr-1)    
            env_ratio_tom = env_ratio_ton*(env_ratio**(env_mminusn_pwr))
            
            ! Calculationg some pieces of derivative
            env_rnp = env_n_pwr * env_ratio_ton
            env_rmp = env_m_pwr * env_ratio_tom
            
            ! Now ratio_ton becomes really ratio to n_pwr   
            env_ratio_ton = env_ratio_ton * env_ratio
            
            ! Another piece of derivative
            env_num = env_rmp - env_rnp - env_mminusn_pwr*env_ratio_ton*env_ratio_tom
            ! same as:
            ! num = rmp - rnp + ( n_pwr - m_pwr )*ratio**(n_pwr + m_pwr - 1)

            ! Now ratio_tom becomes really ratio to m_pwr
            env_ratio_tom = env_ratio_tom * env_ratio
            
            ! The sum of these s(i,j) is the analytic coord no.
            env_s(j,i) = (1.d0-env_ratio_ton)/(1.d0-env_ratio_tom)  ! this is f(rij)
            env_s(i,j) = env_s(j,i)
    
            ! Derivative of s
            env_dS(j,i) = env_num/( env_rij0(pairkindmat(j,i))* (1.d0 - env_ratio_tom)**2)
            env_dS(i,j) = env_dS(j,i)
         END IF
      END DO ! on j 
   END DO ! on i

   ! KR: Calculating the sum of s(i,j), which is CNatom(i), global.
   ! This is the coord num of any single atom i
   ! No problem because s(i,i) is always zero
   DO i = 1, natom 
      env_cn(i) = SUM(env_s(:,i))
      WRITE(184,*) i, env_cn(i)
   END DO


   !-------------------------------------------------------------------------
   ! Calculating the sum on j of derivatives of s
   ! with respect to coordinates of every atoms   
   DO i = 1, natom 
      DO j = 1, natom 
         IF (j .NE. i) THEN
            ! Calculating the derivative of coordination number with
            ! respect to distance between atoms i and j, multiplied for 
            ! the derivative of distance between atoms i and j
            ! with respect to position of atom i
            env_comfactor  = env_dS(j,i)/ pair_dist(j,i)
            env_dS_dx(j,i) = env_comfactor * xyz_dist(1,j,i) ! order of indexes depends on the definition of xyz_dist
            env_dS_dy(j,i) = env_comfactor * xyz_dist(2,j,i) ! order of indexes depends on the definition of xyz_dist
            env_dS_dz(j,i) = env_comfactor * xyz_dist(3,j,i) ! order of indexes depends on the definition of xyz_dist
         ELSE
            env_dS_dx(j,i) = 0.d0
            env_dS_dy(j,i) = 0.d0
            env_dS_dz(j,i) = 0.d0 
         END IF
      END DO
   END DO     


   DO i = 1, natom 
      env_dcn_dx(i) = 2.d0 * SUM(env_dS_dx(1:natom, i))
      env_dcn_dy(i) = 2.d0 * SUM(env_dS_dy(1:natom, i))
      env_dcn_dz(i) = 2.d0 * SUM(env_dS_dz(1:natom, i))
      write(204,*) i, env_dcn_dx(i), env_dcn_dy(i), env_dcn_dz(i)
   END DO  

  
END SUBROUTINE envcn

