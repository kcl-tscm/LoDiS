SUBROUTINE metaforce
 
   USE PARACLUSTER
   USE CLUSTER
   USE META
   USE DISTANCE
   
   implicit none
   integer :: i, j, k, l_cv
   double precision :: vgauss, expcommon
   
   !_________________________________________________________________________________
   ! Storing CV values at every metadynamic step tau_G
   !---------------------------------------------------------------------------------   
   IF (mod(ipas,metaperiod).eq.0) THEN
      NG = NG + 1  ! NG is initialized in main.f90
      
      s_of_t(:,NG) = coll_var ! Array operation
      
      !Writing in a file the position of added gaussians (named metahistory.out)
      ! (using the inner loop statement)
      WRITE (110, *) (coll_var(l_cv), l_cv=1,num_cv)  ! Add time step and height
   END IF

   SELECT CASE (num_cv) 
      !_____  _____  _____  _____  _____  _____  _____  _____  _____  _____  _____  _____  _____ 
      CASE (1) ! Only one CV

         ! Initialization of variables 
         dVg_ds = 0.d0 ! Derivative
         !_________________________________________________________________________________
         ! Computing the derivative of the history-dependent potential Vg with respect to s
         !---------------------------------------------------------------------------------
         ! Sum over the past gaussians
         do i = 1, NG
            ! Computing some quantities  
            expoc(1)  = coll_var(1)- s_of_t(1,i)
            vgauss = exp(-(expoc(1)* expoc(1))/deno(1))
            
            ! Derivative (it is a sum) 
            ! cfr. LAIO --> Improve the performance with dVg_ds global!!
            dVg_ds(1) = dVg_ds(1) -vgauss *(expoc(1)/halfdeno(1))   
         end do
         dVg_ds(1) = dVg_ds(1) * gheight

         SELECT CASE (collvar_case(1))   
            !====================================================================================
            CASE (1) ! CN              
               !_________________________________________________________________________________
               ! Force contribution to be added to the real force
               !---------------------------------------------------------------------------------                      
               fgx = - dVg_ds(1)* dS_dx
               fgy = - dVg_ds(1)* dS_dy
               fgz = - dVg_ds(1)* dS_dz
               
            !===================================================================================
            CASE (2) ! C2N          
               !_________________________________________________________________________________
               ! Force contribution to be added to the real force
               !---------------------------------------------------------------------------------
               ! Initializing force contributions for every atoms             
               fg2nx = - dVg_ds(1)* dS2n_dx
               fg2ny = - dVg_ds(1)* dS2n_dy
               fg2nz = - dVg_ds(1)* dS2n_dz
               
            !====================================================================================
            CASE (3) ! SFN              
               !_________________________________________________________________________________
               ! Force contribution to be added to the real force
               !---------------------------------------------------------------------------------             
               fgsfx = - dVg_ds(1)* dSsf_dx
               fgsfy = - dVg_ds(1)* dSsf_dy
               fgsfz = - dVg_ds(1)* dSsf_dz
               
            !====================================================================================
            CASE (4) ! CNN          
               !_________________________________________________________________________________
               ! Force contribution to be added to the real force
               !---------------------------------------------------------------------------------              
               fgcnnx = - dVg_ds(1)* dScnn_dx
               fgcnny = - dVg_ds(1)* dScnn_dy
               fgcnnz = - dVg_ds(1)* dScnn_dz
               
            !====================================================================================
            CASE (5) ! CN_bim              
               !_________________________________________________________________________________
               ! Force contribution to be added to the real force
               !---------------------------------------------------------------------------------                      
               cn_fgx = - dVg_ds(1)* cn_dS_dx
               cn_fgy = - dVg_ds(1)* cn_dS_dy
               cn_fgz = - dVg_ds(1)* cn_dS_dz
               
            !===================================================================================
            CASE (6) ! d_com              
               !_________________________________________________________________________________
               ! Force contribution to be added to the real force
               !---------------------------------------------------------------------------------                      
               d_com_fgx = - dVg_ds(1)* d_com_dcv_dx
               d_com_fgy = - dVg_ds(1)* d_com_dcv_dy
               d_com_fgz = - dVg_ds(1)* d_com_dcv_dz
               
            !===================================================================================
         END SELECT
      !_____  _____  _____  _____  _____  _____  _____  _____  _____  _____  _____  _____  _____ 
      CASE (2:)                                                                  ! 2 or more CVs
         dVg_ds(:) = 0.d0 ! Initialization of derivative
        
         ! Calculation of the derivative of the history-dependent
         ! meta-potential Vg with respect of every CV (s)
         !----------------------------------------------------
         ! Sum over the past gaussians
         DO i = 1, NG
            ! Computing some quantities
            expoc = coll_var - s_of_t(:,i) ! Array operation
            
            expcommon = 0.d0 ! Initialization
            DO l_cv = 1, num_cv
               expcommon = expcommon -((expoc(l_cv)* expoc(l_cv))/deno(l_cv))
            END DO
            
            vgauss = exp(expcommon)
            ! Derivative (it is a sum)
            DO l_cv = 1, num_cv
               dVg_ds(l_cv) = dVg_ds(l_cv) -vgauss *(expoc(l_cv)/halfdeno(l_cv))
            END DO
         END DO
         DO l_cv = 1, num_cv
            dVg_ds(l_cv) = dVg_ds(l_cv) * gheight
         END DO 

         DO l_cv = 1, num_cv
            SELECT CASE (collvar_case(l_cv))   
               !====================================================================================
               CASE (1) ! CN
                  !_________________________________________________________________________________
                  ! Force contribution to be added to the real force
                  !---------------------------------------------------------------------------------
		  fgx = - dVg_ds(l_cv)* dS_dx
		  fgy = - dVg_ds(l_cv)* dS_dy
		  fgz = - dVg_ds(l_cv)* dS_dz
                  
               !===================================================================================
               CASE (2) ! C2N          
                  !_________________________________________________________________________________
                  ! Force contribution to be added to the real force
                  !---------------------------------------------------------------------------------
		  fg2nx = - dVg_ds(l_cv)* dS2n_dx
		  fg2ny = - dVg_ds(l_cv)* dS2n_dy
		  fg2nz = - dVg_ds(l_cv)* dS2n_dz
                  
               !====================================================================================
               CASE (3) ! SFN          
                  !_________________________________________________________________________________
                  ! Force contribution to be added to the real force
                  !---------------------------------------------------------------------------------                     
		  fgsfx = - dVg_ds(l_cv)* dSsf_dx
		  fgsfy = - dVg_ds(l_cv)* dSsf_dy
		  fgsfz = - dVg_ds(l_cv)* dSsf_dz
                  
               !====================================================================================
               CASE (4) ! CNN          
                  !_________________________________________________________________________________
                  ! Force contribution to be added to the real force
                  !---------------------------------------------------------------------------------                  
		  fgcnnx = - dVg_ds(l_cv)* dScnn_dx
		  fgcnny = - dVg_ds(l_cv)* dScnn_dy
		  fgcnnz = - dVg_ds(l_cv)* dScnn_dz
                  
               !====================================================================================
               CASE (5) ! CN_bim
                  !_________________________________________________________________________________
                  ! Force contribution to be added to the real force
                  !---------------------------------------------------------------------------------
		  cn_fgx = - dVg_ds(l_cv)* cn_dS_dx
		  cn_fgy = - dVg_ds(l_cv)* cn_dS_dy
		  cn_fgz = - dVg_ds(l_cv)* cn_dS_dz
                  
               !===================================================================================
               CASE (6) ! d_com
                  !_________________________________________________________________________________
                  ! Force contribution to be added to the real force
                  !---------------------------------------------------------------------------------
		  d_com_fgx = - dVg_ds(l_cv)* d_com_dcv_dx
		  d_com_fgy = - dVg_ds(l_cv)* d_com_dcv_dy
		  d_com_fgz = - dVg_ds(l_cv)* d_com_dcv_dz
                  
               !===================================================================================
            END SELECT
         END DO ! on l_cv
      !_____  _____  _____  _____  _____  _____  _____  _____  _____  _____  _____  _____  _____  __
   END SELECT

END SUBROUTINE metaforce
