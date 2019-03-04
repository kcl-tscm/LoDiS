SUBROUTINE bim_cn_light
   !=====================
   ! CN_bim light version
   ! (No derivatives) 
   !=====================
   USE PARACLUSTER  
   USE CLUSTER   
   USE POTENTIAL
   USE ENFORCE
   USE DISTANCE
   USE META
   
   IMPLICIT NONE
   
   integer               :: i, j, itypik, mminusn_pwr
   double precision      :: ratio_ton, ratio_tom
   double precision      :: dij0, dij, ratio
   real(8), dimension(3) :: rij0 ! 3 possible cases: atoms A-A, B-B, A-B
   

   if (ipas.eq.1) write(*,*) 'bim_cn_light> first calculation of CN_bim'
   
   mminusn_pwr = cn_m_pwr - cn_n_pwr
   
   ! Evaluating the values of rij0 parameter in dependence of atom pair kind
   rij0 = dist * cn_rzero * dsqrt(2.d0)
 
   ! For 'false' pairs and elements on the diagonal 
   cn_s = 0.d0
   
   do i = 1, natom-1
      do j = i+1, natom
         IF (pair_yesno_mat(j,i)) THEN   
		 ! Checking what kind of atoms the couple i,j is made of
		 if      ((itype(i).eq.1).and.(itype(j).eq.1)) then          
		    itypik = 1  ! same metal A          
		 else if ((itype(i).eq.2).and.(itype(j).eq.2)) then             
		    itypik = 2  ! same metal B          
		 else             
		    itypik = 3  ! A-B interaction                
		 end if
		        
		 ! Evaluating the ref distance for the couple of atoms 
		 dij0 = dist(itypik) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~LP: check if it's ok

		 ! Calculating the distance between atom i and atom j
		 dij = pair_dist (i,j) - dij0       
		 
		 ! Calculating the analytic coordination number  
		 if (dij.le.0) then 
		    ! The sum of these s(i,j) is the analytic coord no.
		    cn_s(i,j) = 1.d0
		    cn_s(j,i) = cn_s(i,j)
		 else
		    ratio = dij/rij0(itypik)
		    ! I need to rise to (X_pwr-1) for the derivetive  
		    ratio_ton = ratio**(cn_n_pwr)    
		    ratio_tom = ratio_ton*(ratio**(mminusn_pwr))
		    
		    ! The sum of these s(i,j) is the analytic coord no.
		    cn_s(i,j) = (1.d0-ratio_ton)/(1.d0-ratio_tom)  ! this is f(rij)
		    cn_s(j,i) = cn_s(i,j)
	 	 end if 
         ENDIF      	     
      end do ! on j 
   end do ! on i

   ! LP: Calculating the sum of s(i,j), which is CNatom(i), global.
   ! This is the coord num of any single atom i
   cn_CNatom(:) = 0.d0 ! Initializing
   do i = 1, natom
      do  j = 1, natom 
         cn_CNatom(i) = cn_CNatom(i) + cn_s(j, i)
      end do
   end do
   
END SUBROUTINE bim_cn_light

