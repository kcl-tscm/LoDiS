SUBROUTINE cnfnum
!========================================
! Common Neighbour Function
!========================================
! Calculating the number of pairs that have the same number of common nearest neighbours. 
! The derivatives are also calculated.
! OLD VERSION (BUT USING LOCAL VARIABLES)
! For computationally taxing runs, use the cnfnum_light subroutine instead.
!========================================
   
   !USE OMP_LIB                                                 !>>> CPU_TIME
   USE PARACLUSTER  
   USE CLUSTER   
   USE POTENTIAL
   USE ENFORCE
   USE DISTANCE
   USE META
   
   IMPLICIT NONE
   
   INTEGER :: i, j, pp, k
   REAL(8) :: cnf_num_loc
   REAL(8), DIMENSION(3) :: sum_a, sum_b, sum_c
   REAL(8), DIMENSION(nsiz) :: friends
   REAL(8), DIMENSION(nsiz, nsiz) :: l_ij
   LOGICAL, DIMENSION(nsiz, nsiz) :: cutoffmatwin
   REAL(8), DIMENSION(nsiz, nsiz, 2) :: sigma_r, lambda_l
   !REAL(8) :: t_start, t_end                                   !>>> CPU_TIME
   !t_start = OMP_GET_WTIME()                                   !>>> CPU_TIME
  
   IF (ipas.eq.1) WRITE(*,*) 'cnfnum> first calculation of Common Neighbour Function (CNF)'
  
   !---------------------------------------------------------------------------
   ! Storing in sigma_r the values of the sigmoid function and its
   ! derivatives for each pair i,j. So they are not computed more than once.
   ! > sigma_r(i,j,1) : value of the sigmoid function
   ! > sigma_r(i,j,2) : deriv of the sigmoid function with respect of r_ij(i,j)
   !
   ! sigma_r(i,j,1:2) is a symmetric matrix. Only a part needs to be computed.
   !--------------------------------------------------------------------------- 
   sigma_r = 0.d0
   l_ij    = 0.d0
   cutoffmatwin = .FALSE.
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO PRIVATE(i, j), SCHEDULE(STATIC)
   DO i = 1, natom-1
      DO j = i+1, natom
         IF (cutoffmat(j,i)) THEN
            sigma_r(j,i,1:2) = f_sigmoid_der(pair_dist(j,i), pairkindmat(j,i))
            sigma_r(i,j,:) = sigma_r(j,i,:)
         END IF
      END DO !j
   END DO !i
!$OMP END DO
 
   !---------------------------------------------------------------------------
   ! Storing in l_ij the value of the term l_ij(i,j) in the calculation
   ! of the CV CNF. It is a triangular matrix
   !---------------------------------------------------------------------------
!$OMP DO PRIVATE(i, j, friends), SCHEDULE(STATIC)
   DO i = 1, natom-1
      DO j = i+1, natom
         IF (cutoffmat(j,i)) THEN
            WHERE (cutoffmat(:,i) .AND. cutoffmat(:,j))
               friends(1:natom) = sigma_r(j,i,1) *sigma_r(1:natom,i,1) *sigma_r(1:natom,j,1)
            ELSE WHERE
               friends(1:natom) = 0.d0
            END WHERE
            l_ij(j,i) = SUM(friends(1:natom))
         END IF
      END DO !j
   END DO !i
!$OMP END DO
 
!$OMP SINGLE PRIVATE(i,j) 
   !---------------------------------------------------------------------------
   ! Creating the matrix cutoffmatwin
   !---------------------------------------------------------------------------
   !(the values of 'inflim' and 'suplim' are given in the main) 
   DO i = 1, natom-1
      DO j = i+1, natom
         cutoffmatwin(j,i) = ((l_ij(j,i).GT. inflim).AND.(l_ij(j,i).LT. suplim))
         cutoffmatwin(i,j) = cutoffmatwin(j,i)
      END DO
   END DO

   !---------------------------------------------------------------------------
   ! Storing in lambda_l the values of the window function and its
   ! derivatives for each pair i,j. So they are not computed more than once.
   ! > lambda_l(i,j,1) : window function
   ! > lambda_l(i,j,2) : its derivative with respect of l_ij(i,j)
   ! 
   ! It is a triangular matrix, it needs to be defined square symmetric for
   ! the derivative
   !---------------------------------------------------------------------------
   lambda_l = 0.d0
   DO i = 1, natom-1
      DO j = i+1, natom 
         IF (cutoffmatwin(j,i)) THEN
            lambda_l(j,i,1:2) = f_window_der(l_ij(j,i))
            lambda_l(i,j,2) = lambda_l(j,i,2)
         END IF
      END DO !j
   END DO !i
   
   !---------------------------------------------------------------------------
   ! Return value for the CV CNF
   !---------------------------------------------------------------------------
   cnf_num_loc = 0.d0
   DO i = 1, natom-1
      DO j = i+1, natom
         cnf_num_loc = cnf_num_loc + lambda_l(j,i,1)
      END DO !j
   END DO !i
   cnf_num = cnf_num_loc 
!$OMP END SINGLE
   
   !---------------------------------------------------------------------------
   ! Calculating the values of derivatives:
   ! dCNF/dx(i) 
   ! dCNF/dy(i)
   ! dCNF/dz(i)
   !---------------------------------------------------------------------------   
!$OMP DO PRIVATE(i, j, k, pp, sum_a, sum_b, sum_c), SCHEDULE(DYNAMIC,1)
   DO pp = 1, natom
      sum_a = 0.d0
      sum_c = 0.d0
      DO i = 1, natom
         IF (cutoffmatwin(i,pp)) THEN
            sum_b = 0.d0
	    DO k = 1, natom
	       IF (cutoffmat(i,pp) .AND. cutoffmat(k,i) .AND. cutoffmat(k,pp)) THEN             
	          sum_b = sum_b +                                      & 
	          &  (sigma_r(i,pp,2) * xyz_dist(:,i,pp)/pair_dist(i,pp)       &
	          & * sigma_r(k,i,1) * sigma_r(k,pp,1))                 &
	          & +(sigma_r(i,pp,1) * sigma_r(k,i,1) * sigma_r(k,pp,2) &
	          & * xyz_dist(:,k,pp)/pair_dist(k,pp))
	       END IF
	    END DO
	    sum_a = sum_a + lambda_l(i,pp,2) * sum_b         	    
         END IF    
      END DO
      DO i = 1, natom-1
         DO j = i+1, natom
            IF (cutoffmatwin(i,pp) .AND. cutoffmat(j,i) .AND. cutoffmat(i,pp) .AND. cutoffmat(j,pp)) THEN
               sum_c = sum_c +lambda_l(j,i,2) *sigma_r(j,i,1) *              &
               &  (sigma_r(i,pp,2) *xyz_dist(:,i,pp)/pair_dist(i,pp) *sigma_r(j,pp,1) &
               &  +sigma_r(i,pp,1) *sigma_r(j,pp,2) *xyz_dist(:,j,pp)/pair_dist(j,pp) )    
            END IF
         END DO
      END DO
      dScnn_dx(pp) = sum_a(1) + sum_c(1)
      dScnn_dy(pp) = sum_a(2) + sum_c(2)
      dScnn_dz(pp) = sum_a(3) + sum_c(3)
   END DO
!$OMP END DO
!$OMP END PARALLEL
   !t_end = OMP_GET_WTIME()                                     !>>> CPU_TIME
   !! Writing the cpu time for the subroutine                   !>>> CPU_TIME 
   !WRITE(500,'(1I10, 4X, 1F10.8)') ipas, t_end -t_start        !>>> CPU_TIME

END SUBROUTINE cnfnum

