subroutine pair_distances
   ! Calculate pair_dist(i,j) : distances for all the pairs of atoms
   !              x_dist(i,j) is xyz_dist(1,i,j)
   !              y_dist(i,j) is xyz_dist(2,i,j)
   !              z_dist(i,j) is xyz_dist(3,i,j)
   !         cutoffmat(i,j) : the cutoff matrix (logical), if needed

   use cluster, only : natom
   use enforce, only : x, y, z, u, v, w 
   use distance
   !USE OMP_LIB                                                 !>>> CPU_TIME

   implicit none
   integer :: i, j  
   !REAL(8) :: t_start, t_end                                   !>>> CPU_TIME

   !t_start = OMP_GET_WTIME()                                   !>>> CPU_TIME
!$OMP PARALLEL DEFAULT(SHARED) 
   !----------------------------------------------
   ! Creating x_dist, y_dist, z_dist and pair_dist 
   ! This should replicate the results of
   ! MolDyn M2.3; I need a way of dealing with the
   ! problem of reproducibility 
   
!$OMP SECTIONS PRIVATE(i, j)     
!$OMP SECTION
   do j=1, natom   
      do i=1, natom
         xyz_dist(1,i,j) = x(j) +u(j) -x(i) -u(i)       
      end do
   end do

!$OMP SECTION
   do j=1, natom   
      do i=1, natom
         xyz_dist(2,i,j) = y(j) +v(j) -y(i) -v(i)     
      end do
   end do

!$OMP SECTION
   do j=1, natom   
      do i=1, natom
         xyz_dist(3,i,j) = z(j) +w(j) -z(i) -w(i)       
      end do
   end do
!$OMP END SECTIONS
   
!$OMP DO PRIVATE(i, j), SCHEDULE(STATIC) 
   do i=1, natom   
      do j=1, natom    
         pair_dist(j,i) = DSQRT(SUM(xyz_dist(:,j,i)**2))               
      end do
   end do
!$OMP END DO
   !!write(*,*) "pair_distances:> pair(1,2) pair(1,natom)",pair_dist(1,2), pair_dist(1,natom)
!
   !---------------------------------------------------------
   ! Creating the cutoff matrix (logical),
   ! considering different pair kinds: (1)A-A, (2)B-B, (3)A-B
!$OMP SINGLE PRIVATE(i,j)
   IF (cutoffmat_req) THEN
      DO i = 1, natom-1
        DO j = i+1, natom     
           cutoffmat(j,i) = (pair_dist(j,i) .LE. dcutoff(pairkindmat(j,i) ) )
           cutoffmat(i,j) = cutoffmat(j,i)
        END DO
      END DO     
   END IF
!$OMP END SINGLE
!$OMP END PARALLEL
   !t_end = OMP_GET_WTIME()                                     !>>> CPU_TIME
   !! Writing the cpu time for the subroutine                   !>>> CPU_TIME 
   !WRITE(500,'(1I10, 4X, 1F10.8)') ipas, t_end -t_start        !>>> CPU_TIME

end subroutine pair_distances

