SUBROUTINE collective
!======================================
! Collective Variables Choice
!======================================
! Subroutine to calculate collective variables
! collective variables considered:
! > coordination number (CN)      - sigmoid function
! > 2nd neighbour number (C2N)    - window function
! > stacking fault number (SFN)   - window function
! > common neighbour number (CNN) - common neighbour function
! > CN_bim - sigmoid function for a subset of pairs (AA and/or BB and/or AB) 
! > d_com  - squared distance between CoM of atoms A and B
!======================================

   USE META
   USE CLUSTER
   USE PARACLUSTER
   USE ENFORCE
   USE POTENTIAL
   USE DISTANCE

   IMPLICIT NONE

   integer :: i, l_cv, j
   

   DO l_cv = 1, num_cv
      
      SELECT CASE (collvar_case(l_cv))
         !===========================
         CASE (1) ! CN   
            ! Calculating CN of single atoms (in the global variable CNatom)
            if (metadyn.eq.'ya') then 
               call cnum         ! to calculate coord number and derivative
            else
               call cnum_light   ! to calculate only coord number
            end if
         
            ! Calculating the total coordination number
            coll_var(l_cv) = SUM(CNatom)
         !===========================
         CASE (2) ! C2N   
            ! Calculating C2N of single atoms (in the global variable C2Natom)
            if (metadyn.eq.'ya') then 
               call c2num         ! to calculate C2N and derivative
            else
               call c2num_light   ! to calculate only C2N
            end if
         
            ! Calculating the total C2N
            coll_var(l_cv) = SUM(C2Natom)
         !===========================
         CASE (3) ! SFN   
            ! Calculating SFN of single atoms (in the global variable SFNatom)
            if (metadyn.eq.'ya') then 
               call sfnum         ! to calculate SFN and derivative
            else
               call sfnum_light   ! to calculate only SFN
            end if
         
            ! Calculating the total SFN
            coll_var(l_cv) = SUM(SFNatom)
         !===========================
         CASE (4) ! CNN   
            ! Calculating CNN (common neighbour function)
            IF (metadyn.eq.'ya') then 
               call cnfnum         ! to calculate CNN and derivative
            ELSE
               call cnfnum_light   ! to calculate only CNN   
            END IF
            ! total CNN is calculated in cnfnum, being function of the pairs
            coll_var(l_cv) = cnf_num
         !=========================== 
         CASE (5) ! CN_bim   
            ! Calculating CN_bim of single atoms (in the global variable cn_CNatom)
            if (metadyn.eq.'ya') then 
               call bim_cn         ! to calculate coord number and derivative
            else
               call bim_cn_light   ! to calculate only coord number
            end if
         
            ! Calculating the total coordination number bimetallic
            coll_var(l_cv) = SUM(cn_CNatom)
         !===========================   
         CASE (6) ! d_com   
            ! Calculating CN_bim of single atoms (in the global variable cn_CNatom)
            if (metadyn.eq.'ya') then 
               call d_com          ! to calculate cv and derivative
            else
               call d_com_light    ! to calculate only cv
            end if
         
            ! Calculating the d_com
            coll_var(l_cv) = d_com_cv 
         CASE DEFAULT
            WRITE(*,*) 'collective> Error: selected CV is not implemented.'
            STOP
      END SELECT 

   ENDDO ! on l_cv different collective variables

END SUBROUTINE collective
