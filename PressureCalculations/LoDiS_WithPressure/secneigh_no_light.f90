SUBROUTINE c2num_light
!======================================
! Calculating the stacking fault number
! The derivative are not calculated in
! order toimprove the performance
!======================================
   
   USE PARACLUSTER  
   USE CLUSTER   
   USE POTENTIAL
   USE ENFORCE
   USE DISTANCE
   USE META
   
   IMPLICIT NONE
   
   integer               :: i, j, itypik, mminusn_pwr
   double precision      :: ratio_ton, ratio_tom
   double precision      :: dij, ratio
   real(8), dimension(3) :: rij0 ! 3 possible cases: atoms A-A, B-B, A-B
   real(8), dimension(3) :: dij0 ! 3 possible cases: atoms A-A, B-B, A-B
        

   if (ipas.eq.1) write(*,*) '> c2num_light: first calculation of C2N.'
   
   ! IMPROVEMENT: All this could become global! ++++++++++++++++++++++++++++
   mminusn_pwr = m2n_pwr - n2n_pwr
   
   ! Evaluating the values of rij0 parameter in dependence of atom pair kind
   rij0 = dist * r2n * dsqrt(2.d0) ! rij0 and dist are 3 element vectors
   
   ! Evaluating the ref distance for the couple of atoms (THIS IS DIFFERENT) 
   dij0 = dist * d2n * dsqrt(2.d0) 
   ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   

   do i = 1, natom-1
      do j = i+1, natom
         
         ! Checking what kind of atoms the couple i,j is made of
         if      ((itype(i).eq.1).and.(itype(j).eq.1)) then          
            itypik = 1  ! same metal A          
         else if ((itype(i).eq.2).and.(itype(j).eq.2)) then             
            itypik = 2  ! same metal B          
         else             
            itypik = 3  ! A-B interaction                
         end if
                
         ! Calculating the distance between atom i and atom j
         dij = pair_dist (i,j) - dij0(itypik)       
         
         ! Calculating the analytic C2N          
         ratio = dij/rij0(itypik)
         ! I need to rise to (X_pwr-1) for the derivetive  
         ratio_ton = ratio**(n2n_pwr)    
         ratio_tom = ratio_ton*(ratio**(mminusn_pwr))
            
         ! The sum of these s(i,j) is the analytic coord no.
         s2n(i,j) = (1.d0-ratio_ton)/(1.d0-ratio_tom)  ! this is f(rij)
         s2n(j,i) = s2n(i,j)
 
      end do ! on j 
   end do ! on i

   ! LP: Calculating the sum of s2n(i,j), which is C2Natom(i,j), global.
   ! This is the second neighbour number of any single atom i
   C2Natom(:) = 0.d0 ! Initializing
   do i = 1, natom
      do j = 1, natom 
         if (j.ne.i) C2Natom(i) = C2Natom(i) + s2n(i, j)
      end do
   end do

END SUBROUTINE c2num_light

