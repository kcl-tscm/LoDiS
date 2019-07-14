SUBROUTINE cnfnum_light
!=======================================
! Common Neighbour Function light
!=======================================
! Calculates the number of pairs that have the same number of common nearest neighbours. 
! The derivative are not calculated in order to improve performance.
!=======================================
   
   USE PARACLUSTER  
   USE CLUSTER   
   USE POTENTIAL
   USE ENFORCE
   USE DISTANCE
   USE META
   
   IMPLICIT NONE
   
   INTEGER :: i, j, k
   REAL(8) :: cnf_num_loc
   REAL(8), DIMENSION(nsiz) :: friends
      
   IF (ipas.eq.1) WRITE(*,*) 'cnfnum_light> first calculation of CNN.'
   
   cnf_num_loc = 0.d0
   
   DO i = 1, natom-1
      DO j = i+1, natom
         IF (cutoffmat(j,i)) THEN
            WHERE (cutoffmat(:,i) .AND. cutoffmat(:,j))
               friends(1:natom) = f_sigmoid(pair_dist(j,i), pairkindmat(j,i))* &
               & f_sigmoid(pair_dist(:,i), pairkindmat(:,i))* &
               & f_sigmoid(pair_dist(:,j), pairkindmat(:,j))
            ELSE WHERE
               friends(1:natom) = 0.d0
            END WHERE
            cnf_num_loc = cnf_num_loc + f_window(SUM(friends(1:natom)))
         END IF
      END DO !j
   END DO !i

   ! Return value
   cnf_num = cnf_num_loc 
   
END SUBROUTINE cnfnum_light

