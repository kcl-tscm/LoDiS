SUBROUTINE ev_time
 !=======================================
 ! This subroutine is called in:
 ! main_MD (3 places)
 !======================================= 
 USE PARACLUSTER  
 USE CLUSTER   
 USE ENFORCE
 USE DISTANCE
 USE META
 USE SUBSTRATE
 USE ENVIRONMENT
 USE module_sticky
 
 IMPLICIT NONE

 INTEGER :: i, np, l_cv
 REAL(8) :: secd, firstd, rvec 
 REAL(8) :: delta
 REAL(8) :: partialx(nsiz), partialy(nsiz), partialz(nsiz)
 REAL(8), DIMENSION(nsiz, num_cv) :: fxmeta, fymeta, fzmeta !Used to check the forces when collvar_wanted is .true.


 ! Metadynamics/CVs and MgO substrate stuff use their own shared
 ! parallel soubroutine, called only once per time step, when these functions are swithed on.
 ! However standard MD does not use it. 

 !%%%%%%%%%%%%%%%%%%%%%%%%%
 !=========================
 ! Main MD loop (over npas)
 !=========================
 DO ipas=1, npas      
    tempo=(ipas+(nd_proc-1)*npas)*tstep*1.d12   !!TIME in picoseconds

    partialx(1:natom)=0.d0
    partialy(1:natom)=0.d0
    partialz(1:natom)=0.d0

    IF(type_potential=='rgl') THEN
       !Verlet cage to save time in force_rgl
       firstd=0.d0
       secd=0.d0

       DO i=1,natom
          partialx(i)=partialx(i)+du(i)    !du dv dw of natom must be initialised in case of growth!        
          partialy(i)=partialy(i)+dv(i)            
          partialz(i)=partialz(i)+dw(i) 
          delta=partialx(i)**2+partialy(i)**2+partialz(i)**2

          IF ((delta.GT.secd).AND.(delta.LE.firstd)) THEN
             secd=delta
          ENDIF
          IF (delta.GT.firstd) THEN
             secd=firstd
             firstd=delta
          ENDIF
       ENDDO
       rvec=DSQRT(firstd)+DSQRT(secd)

       IF (rvec.GT.rshell) THEN
          DO i=1,natom
             partialx(i)=0.d0
             partialy(i)=0.d0
             partialz(i)=0.d0
          ENDDO
          CALL bigvoi     
       ENDIF
    ENDIF !! bigvoi is called only for RGL-type potential

    !IN TIME: call to voisin, force, and vel at time t, finally a call to therma for the new position at tiem t+dt
    !
    CALL voisin   !!calculating the number of first neighbour
    !
    CALL force_choice
    !
    !
    ! RGL force, meta-force and substrate force are then in eV/arete <-- force in eV/A are divided by aretebim
    !KEEP THE UNTS OF THE FORCES in eV/arete as the time step is done with respect that choice.
    !
    !-----------------------------------------------
    ! Start of Metadynamics/CVs related calculations
    !-----------------------------------------------
    IF ((collvar_wanted).OR.(metadyn.eq.'ya').OR.mgo_substrate .OR. impl_env) THEN
       ! Calling the subroutine that calculate pair distances
       ! (It needs to be done for all the implemented CVs)
       CALL pair_distances

       if ((collvar_wanted).or.(metadyn.eq.'ya')) then
          CALL collective
        
          if ((metadyn.eq.'ya') .and. (ipas.GE.metaperiod)) then
             CALL metaforce

             ! Adding the meta-forces to the "real" forces 
             DO l_cv = 1, num_cv
                SELECT CASE (collvar_case(l_cv))
                   !=================================================
                   CASE (1) ! CN
                      ! Adding the force contribution of metadynamics 
                      fx = fx + fgx
                      fy = fy + fgy
                      fz = fz + fgz     
                      !=================================================
                   CASE (2) ! C2N
                      ! Adding the force contribution of metadynamics 
                      fx = fx + fg2nx
                      fy = fy + fg2ny
                      fz = fz + fg2nz     
                      !=================================================
                   CASE (3) ! SFN
                      ! Adding the force contribution of metadynamics 
                      fx = fx + fgsfx
                      fy = fy + fgsfy
                      fz = fz + fgsfz     
                      !=================================================
                   CASE (4) ! CNN
                      ! Adding the force contribution of metadynamics 
                      fx = fx + fgcnnx
                      fy = fy + fgcnny
                      fz = fz + fgcnnz     
                      !=================================================
                   CASE (5) ! CN_bim
                      ! Adding the force contribution of metadynamics 
                      fx = fx + cn_fgx
                      fy = fy + cn_fgy
                      fz = fz + cn_fgz     
                      !=================================================
                   CASE (6) ! d_com
                      ! Adding the force contribution of metadynamics 
                      fx = fx + d_com_fgx
                      fy = fy + d_com_fgy
                      fz = fz + d_com_fgz     
                      !=================================================
                END SELECT
             END DO ! on l_cv 
          end if

          if (collvar_wanted) then
             !-------------------------------------------------------------------------------------------
             ! Writing the value of collective variables and forces at each step (AS A CHECK) 
             !
             if (ipas.eq.1) then
                write (101, *) '# ipas, ', trim(collvar_name(1)), ', ', trim(collvar_name(2))
                IF (metadyn.eq.'ya') THEN 
                   WRITE(102, *) '# ipas,  Average |force (real+MT)|,  Average |MT-force(1)|,  Average |MT-force(2)|'
                END IF  
             end if
             IF ((metadyn.eq.'ya').AND.(ipas.GE.metaperiod)) THEN         
                DO l_cv = 1, num_cv
                   SELECT CASE (collvar_case(l_cv))
                      !=================================================
                      CASE (1) ! CN       
                         fxmeta(1:natom,l_cv) = fgx
                         fymeta(1:natom,l_cv) = fgy
                         fzmeta(1:natom,l_cv) = fgz
                      !=================================================
                      CASE (2) ! C2N  
                         fxmeta(1:natom,l_cv) = fg2nx
                         fymeta(1:natom,l_cv) = fg2ny
                         fzmeta(1:natom,l_cv) = fg2nz     
                      !=================================================
                      CASE (3) ! SFN   
                         fxmeta(1:natom,l_cv) = fgsfx
                         fymeta(1:natom,l_cv) = fgsfy
                         fzmeta(1:natom,l_cv) = fgsfz     
                      !=================================================
                      CASE (4) ! CNN       
                         fxmeta(1:natom,l_cv) = fgcnnx
                         fymeta(1:natom,l_cv) = fgcnny
                         fzmeta(1:natom,l_cv) = fgcnnz     
                      !=================================================
                      CASE (5) ! CN_bim       
                         fxmeta(1:natom,l_cv) = cn_fgx
                         fymeta(1:natom,l_cv) = cn_fgy
                         fzmeta(1:natom,l_cv) = cn_fgz
                      !=================================================
                      CASE (6) ! d_com       
                         fxmeta(1:natom,l_cv) = d_com_fgx
                         fymeta(1:natom,l_cv) = d_com_fgy
                         fzmeta(1:natom,l_cv) = d_com_fgz
                      !=================================================
                   END SELECT
                END DO ! on l_cv 
             END IF
             if (num_cv.eq.1) then
                write (101, *) ipas, coll_var(1)
                IF ((metadyn.eq.'ya').AND.(ipas.GE.metaperiod)) THEN 
                   WRITE(102, '(1I12, 2F16.8)') ipas, SUM(SQRT(fx*fx+fy*fy+fz*fz))/natom, &
                             &SUM(SQRT(fxmeta(1:natom,1)**2 +fymeta(1:natom,1)**2 +fzmeta(1:natom,1)**2))/natom
                END IF 
             else
                write (101, *) ipas, coll_var(1), coll_var(2)
                IF ((metadyn.eq.'ya').AND.(ipas.GE.metaperiod)) THEN 
                   WRITE(102, '(1I12, 3F16.8)') ipas, SUM(SQRT(fx*fx+fy*fy+fz*fz))/natom, &
                             &SUM(SQRT(fxmeta(1:natom,1)**2 +fymeta(1:natom,1)**2 +fzmeta(1:natom,1)**2))/natom, &
                             &SUM(SQRT(fxmeta(1:natom,2)**2 +fymeta(1:natom,2)**2 +fzmeta(1:natom,2)**2))/natom
                END IF
             end if
             !-------------------------------------------------------------------------------------------
          end if ! (collvar_wanted)

       end if ! on the check related to collvar_wanted and metadyn 
    ENDIF
    !-----------------------------------------------
    ! End of Metadynamics/CVs related calculations
    !-----------------------------------------------

    !--------------------------------------------------------------
    ! Adding forces from the substrate MgO
    IF (mgo_substrate) THEN
       IF (metal_on_top) THEN
          CALL force_mgo_mot
       ELSE
          CALL force_mgo
       ENDIF

       fx = fx + mgo_fx
       fy = fy + mgo_fy
       fz = fz + mgo_fz 

       ! Summing the energy of the sub to the pot energy of the system

       ener = ener + ener_sub
    ENDIF
    !-------------------------------------------------------------

    !--------------------------------------------------------------
    ! Adding forces for the sticky atoms
    IF (sticky_atoms_wanted) THEN
       CALL sticky_force
 
       fx = fx + sticky_fx
       fy = fy + sticky_fy
       fz = fz + sticky_fz    
    ENDIF
    !-------------------------------------------------------------

    !--------------------------------------------------------------
    ! Adding forces for the atoms - environment interaction
    IF (impl_env) THEN

       CALL fenvironment

   !DO i = 1, natom
   !write(164,*) ipas, i, fx(i), fy(i), fz(i), env_fx(i), env_fy(i), env_fz(i)
   !END DO

       fx = fx + env_fx
       fy = fy + env_fy
       fz = fz + env_fz    

   !DO i = 1, natom
   !write(224,*) ipas, i, fx(i), fy(i), fz(i)
   !END DO

       ener = ener + ener_env
    ENDIF
    !-------------------------------------------------------------


    CALL vel
   
    !write(70,*) ipas*tstep , etot
    !write(71,*) ipas*tstep , temp
    
	

    IF ((canonical == 'ya') .and. (vel_af) ) THEN
	do i=1,natom
		v_acf(ipas,i,1)=vx(i)
		v_acf(ipas,i,2)=vy(i)
		v_acf(ipas,i,3)=vz(i)
		
	end do
    !write(32,*) ipas,vx(1),vy(1),vz(1)
    !CALL vel_acf_calc
    END IF
    !
    ! WRITING PHYSICAL QUANTITIES at time t
    !
    IF(quenching=='ya') CALL scrittoq   !!output file for quenching
    !
    IF((caloric=='ya').OR.(canonical=='ya')) CALL scrittoc     !!output file at time t
    !
    IF(coalescence=='ya') CALL writing_coal     !!output file at time t
    !
    IF(deposizione=='ya') CALL scrittod !!output file for natom 
    !
    if(metadyn.eq.'ya') call scrittometa !output file for metadynamics
    !
    !Calculate the new position time t+tstep
    !
    CALL therma  

    !
 ENDDO !!end loop on ipas

 !if ((canonical == 'ya') .and. (vel_af)) then
!	 call dff_test
 !end if
 !%%%%%%%%%%%%%%%%%%%%%%%%%
 !=========================
 ! End of the main MD loop

END SUBROUTINE ev_time
