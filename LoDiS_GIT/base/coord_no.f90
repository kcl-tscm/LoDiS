SUBROUTINE cnum
!======================================
! Coordination Number
!======================================
! Calculates the Analytic Coordination Number for monometallic clusters
! Derivatives are calculated in this version
! To improve performance, switch to the cnum light version
!======================================
     
   !USE OMP_LIB                                                 !>>> CPU_TIME
   USE PARACLUSTER  
   USE CLUSTER   
   USE POTENTIAL
   USE ENFORCE
   USE DISTANCE
   USE META
   
   IMPLICIT NONE
   
   integer               :: i, j, mminusn_pwr
   double precision      :: ratio_ton, ratio_tom, rnp, rmp, num
   double precision      :: ratio, comfactor
   real(8), dimension(3) :: rij0 ! 3 possible cases: atoms A-A, B-B, A-B
   REAL(8), DIMENSION(nsiz,nsiz)   :: dij, df_dx, df_dy, df_dz
   !REAL(8) :: t_start, t_end                                   !>>> CPU_TIME
   !t_start = OMP_GET_WTIME()                                   !>>> CPU_TIME
   
   if (ipas.eq.1) write(*,*) 'cnum> first calculation of sigmoid function (CN)'

   !---------------------------------------------------
   mminusn_pwr = m_pwr - n_pwr
   
   ! Evaluating the values of rij0 parameter in dependence of atom pair kind
   rij0 = dist * rzero * dsqrt(2.d0)

   !---------------------------------------------------
   ! Calculating the analytic coordination number
   ! The sum of these s(i,j) is the analytic coord no.

   ! Default Case: (dij(j,i) .LE. 0.d0)
   s  = 1.d0
   ! Derivative of s
   dS = 0.d0

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(i, j), SCHEDULE(STATIC) 
   DO i = 1, natom
      DO j = 1, natom                                       
         ! Distance between atom i and atom j
         ! minus the ref distance for the couple of atoms
         dij(j,i) = pair_dist(j,i) - dist(pairkindmat(j,i))
      END DO
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i), SCHEDULE(STATIC) 
   ! Elements on the diagonal
   DO i = 1, natom
      s(i,i) = 0.d0
   END DO
!$OMP END DO  

!$OMP DO PRIVATE(i, j, ratio, ratio_ton, ratio_tom, rnp, rmp, num), SCHEDULE(DYNAMIC,20)
   ! Case (dij(j,i) .GT. 0.d0)
   DO i = 1, natom-1
      DO j = i+1, natom      
         IF (dij(j,i) .GT. 0.d0) THEN
            ratio = dij(j,i)/rij0(pairkindmat(j,i))
            ! I need to rise to (X_pwr-1) for the derivetive  
            ratio_ton = ratio**(n_pwr-1)    
            ratio_tom = ratio_ton*(ratio**(mminusn_pwr))
            
            ! Calculationg some pieces of derivative
            rnp = n_pwr * ratio_ton
            rmp = m_pwr * ratio_tom
            
            ! Now ratio_ton becomes really ratio to n_pwr   
            ratio_ton = ratio_ton * ratio
            
            ! Another piece of derivative
            num = rmp - rnp - mminusn_pwr*ratio_ton*ratio_tom
            ! same as:
            ! num = rmp - rnp + ( n_pwr - m_pwr )*ratio**(n_pwr + m_pwr - 1)

            ! Now ratio_tom becomes really ratio to m_pwr
            ratio_tom = ratio_tom * ratio
            
            ! The sum of these s(i,j) is the analytic coord no.
            s(j,i) = (1.d0-ratio_ton)/(1.d0-ratio_tom)  ! this is f(rij)
            s(i,j) = s(j,i)
    
            ! Derivative of s
            dS(j,i) = num/( rij0(pairkindmat(j,i))* (1.d0 - ratio_tom)**2)
            dS(i,j) = dS(j,i)
         END IF
      END DO ! on j 
   END DO ! on i
!$OMP END DO

!$OMP DO PRIVATE(i), SCHEDULE(STATIC)
   ! LP: Calculating the sum of s(i,j), which is CNatom(i), global.
   ! This is the coord num of any single atom i
   ! No problem because s(i,i) is always zero
   DO i = 1, natom 
      CNatom(i) = SUM(s(:,i))
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j,comfactor), SCHEDULE(STATIC)
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
            comfactor  = dS(j,i)/ pair_dist(j,i)
            df_dx(j,i) = comfactor * xyz_dist(1,j,i) ! order of indexes depends on the definition of xyz_dist
            df_dy(j,i) = comfactor * xyz_dist(2,j,i) ! order of indexes depends on the definition of xyz_dist
            df_dz(j,i) = comfactor * xyz_dist(3,j,i) ! order of indexes depends on the definition of xyz_dist
         ELSE
            df_dx(j,i) = 0.d0
            df_dy(j,i) = 0.d0
            df_dz(j,i) = 0.d0 
         END IF
      END DO
   END DO     
!$OMP END DO
  
!$OMP SECTIONS PRIVATE(i)
!$OMP SECTION  
   DO i = 1, natom 
      dS_dx(i) = 2.d0 * SUM(df_dx(1:natom, i))
   END DO       
!$OMP SECTION
   DO i = 1, natom 
      dS_dy(i) = 2.d0 * SUM(df_dy(1:natom, i))
   END DO  
!$OMP SECTION
   DO i = 1, natom 
      dS_dz(i) = 2.d0 * SUM(df_dz(1:natom, i))
   END DO  
!$OMP END SECTIONS
!$OMP END PARALLEL
   !---------------------------------------------------------------------------   
   !t_end = OMP_GET_WTIME()                                     !>>> CPU_TIME
   !! Writing the cpu time for the subroutine                   !>>> CPU_TIME 
   !WRITE(500,'(1I10, 4X, 1F10.8)') ipas, t_end -t_start        !>>> CPU_TIME
END SUBROUTINE cnum

