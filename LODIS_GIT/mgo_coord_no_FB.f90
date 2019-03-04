SUBROUTINE mgo_coord
 ! Luca
 ! Calculating the analytic coordination number for MgO
     
 USE PARACLUSTER
 USE CLUSTER 
 USE POTENTIAL
 USE ENFORCE
 USE DISTANCE
 USE SUBSTRATE

 IMPLICIT NONE
   
 INTEGER               :: i, j, k
 REAL      :: comfactor, gamma_sub
 REAL, ALLOCATABLE	   :: dij(:,:), df_dx(:,:), df_dy(:,:), df_dz(:,:)
   
IF (ipas==1) WRITE(*,*) "mgo_coord:> first calculation of CN for MgO substrate"
ALLOCATE(dij(natom,natom), df_dx(natom,natom),df_dy(natom,natom),df_dz(natom,natom))
dij(:,:) =0.d0
df_dx(:,:) = 0.d0
df_dy(:,:) = 0.d0
df_dz(:,:) = 0.d0 
!coordination number using a smooth Fermi function as
!CN(i) = sum_{j\=i} 1 / [1+exp^alpha(rij/rc-1)]
!beta_sub = alpha_sub/ (rc_sub/aretebim) 
!cn(i) = sum_{j/=i} 1 / 1+gamma_sub
!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(i, j), SCHEDULE(STATIC) 
   DO i = 1, natom-1
      DO j = i+1, natom                                       
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
          
            ! The sum of these s(i,j) is the analytic metallic coord no.
            mgo_s(j,i) = 1.d0/(1.d0 + gamma_sub)  
            mgo_s(i,j) = mgo_s(j,i)
    
            ! Derivative of s w.r.t r(i,j)  
            ! Use beta_sub to have aretebim as units
            mgo_dS(j,i) = -beta_sub * gamma_sub / (gamma_sub+1)**2
            mgo_dS(i,j) = mgo_dS(j,i)
			!	           
      END DO ! on j 
   END DO ! on i
!$OMP END DO

!$OMP DO PRIVATE(i), SCHEDULE(STATIC)
   ! LP: Calculating the sum of s(i,j), which is mgo_CN(i), global.
   ! This is the coord num of any single atom i
   ! No problem because s(i,i) is always zero
   DO i = 1, natom
        mgo_cn(i) = SUM(mgo_s(1:natom,i))
		write(898,*) 'mgo_coord:> metallic coord', i, mgo_cn(i)
   END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j,comfactor), SCHEDULE(STATIC)
   !-------------------------------------------------------------------------
   ! Calculating the sum on j of derivatives of s
   ! with respect to coordinates of every atoms   
   ! dcn(i)_dx(k)=-sum_{j/=k} beta_sub * gamma_sub /(1 + gamma_sub)^2 xij/rij 
   !  			+ beta_sub * gamma_sub /(1+gamma_sub)^2 * xij/rij
   ! forza su atom k due to the coordination of atom i
   ! mgo_dcn(i)_dx(k) = [sum_{j/=k} mgo_dS(k,j) * xkj/rkj] - mgo_dS(i,k)*xik/rik
   DO k = 1, natom    !force over kth-atom due to coordination of atom i
      DO j = 1, natom 
         IF (j .NE. k) THEN
            ! Calculating the derivative of coordination number with
            ! respect to distance between atoms i and j, multiplied for 
            ! the derivative of distance between atoms i and j
            ! with respect to position of atom i
            comfactor  = mgo_dS(k,j)/ pair_dist(k,j)
            df_dx(k,j) = comfactor * xyz_dist(1,j,k) ! order of indexes depends on the definition of xyz_dist = x(k) - x(j)
            df_dy(k,j) = comfactor * xyz_dist(2,j,k) ! order of indexes depends on the definition of xyz_dist
            df_dz(k,j) = comfactor * xyz_dist(3,j,k) ! order of indexes depends on the definition of xyz_dist
         ELSE
            df_dx(k,j) = 0.d0
            df_dy(k,j) = 0.d0
            df_dz(k,j) = 0.d0 
         END IF
      END DO
   END DO     
!$OMP END DO
  
!$OMP SECTIONS PRIVATE(i)
!$OMP SECTION 
Do k = 1,natom
	mgo_dcn_dx(k) = + SUM(df_dx(1:natom, k))
   DO i = 1, natom 
      ! No factor 2 here: this is derivative of cn(i) w.r.o. x(i), not total cn w.r.o. x(i)
	  !sum over all j/=k and subtract df_dx(k,i)
	  !force on k-th atom is a sum of contribution from all atom i and j/=k
      mgo_dcn_dx(k) = mgo_dcn_dx(k) - df_dx(k,i) 
   END DO     
 ENDDO     
!$OMP SECTION
   DO k = 1, natom 
     mgo_dcn_dy(k) = SUM(df_dy(1:natom, k)) 
     DO i = 1, natom 
         mgo_dcn_dy(k) = mgo_dcn_dy(k) - df_dy(k,i) 
     END DO     
   END DO 
!$OMP SECTION
DO k = 1, natom 
    mgo_dcn_dz(k) = SUM(df_dz(1:natom, k)) 
   DO i = 1, natom 
       mgo_dcn_dz(k) = mgo_dcn_dz(k) - df_dz(k,i) 
   END DO  
 END DO 
!$OMP END SECTIONS
!$OMP END PARALLEL
   !---------------------------------------------------------------------------   
   !t_end = OMP_GET_WTIME()                                     !>>> CPU_TIME
   !! Writing the cpu time for the subroutine                   !>>> CPU_TIME 
   !WRITE(500,'(1I10, 4X, 1F10.8)') ipas, t_end -t_start        !>>> CPU_TIME
   DEALLOCATE (dij, df_dx,df_dy,df_dz)
END SUBROUTINE mgo_coord

