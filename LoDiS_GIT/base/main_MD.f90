!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Copyright 2020 KCL-LoDiS-Baletto
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                            LoDiS vers. 1.0.0
! The program LoDiS-MD born from a classical Molecular Dynamics (MD) simulation
! initially developed by F. Baletto, R. Ferrando, C. Mottet, G. Treglia. 
! A collaboration from the Physics Department in Genova (Italy) and the CNRS in Marseille.
!
! In the last decade the Baletto's group based at the Physics Department at KCL as developed 
! different tools for study low dimensional systems (LoDiS).
! 
! KCL-Developers: 
! Luca Pavan, Kevin Rossi, Raphael Pinto-Miles, Matteo Tiberi, Laia Delgado, Francesca Baletto
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  LoDiS does Molecular Dynamics simulations for cluster, wires, surfaces (and bulk to be tested).
!  The equation of motion are solved using a velocity-Verlet algorithm
!  Different semi-empirical potential are available:
!  RGL (Gupta like) potntial      -- for late-transition metals (and their alloys)
! The RGL are also known as second-moment approximation if the tight-binding
!  Pacheco or Girifalco potential -- for C60 clusters
!  Lennard-Jones                  -- for Ar
!  Morse potential and modiefied Morse will be available soon (under test)
!!!!!!!!!!!!
!  LoDiS tools:
!  @quenching/annealing (see:   F. Baletto et al. JCP (2002), F. Baletto and R. Ferrando RMP 77 (2005) 371)
!
!  @one-by-one growth 
! of monometallic ( F. Baletto et al. Surf. Sci. (2000), F. Baletto et al. PRL (2000)) 
! of nanoalloys (both core shell or with a fixed probability, F. Baletto PRB (2002) )
! see I Parsina, F Baletto, Tailoring the structural motif of AgCo nanoalloys: core/shell versus Janus-like
! The Journal of Physical Chemistry C 114 (3), 1504-1511
!
!  @freezing/melting of nanoparticles/nanoalloys (see F. Baletto et al. CPL (2002))
!  subsequent updates 
![1] L Pavan, F Baletto, R Novakovic, Multiscale approach for studying melting transitions in CuPt nanoparticles, 
! Physical Chemistry Chemical Physics 17 (42), 28364-28371
![2] K Rossi, T Ellaby, LO Paz-Borbón, I Atanasov, L Pavan, F Baletto, Melting of large Pt@ MgO (1 0 0) icosahedra,
!Journal of Physics: Condensed Matter 29 (14), 145402
![3] K Rossi, LB Pártay, G Csanyi, F Baletto, Thermodynamics of CuPt nanoalloys
!Scientific reports 8 (1), 1-9 
![4]  L Delgado-Callico, K Rossi, R Pinto-Miles, P Salzbrenner, F Baletto, Universal signature in the melting of metallic nanoparticles
! arXiv preprint arXiv:2007.12218
!
!@ Metadynamics for studyng structural rearrangements
!detalis in L. Pavan, PhD 2014,  King's College London and in K. Rossi, PhD 2018, King's College London
! useful references:
![5] L Pavan, C Di Paola, F Baletto, Sampling the energy landscape of Pt13 with metadynamics
! The European Physical Journal D 67 (2), 24
![6] L Pavan, K Rossi, F Baletto, Metallic nanoparticles meet metadynamics, The Journal of chemical physics 143 (18), 184304
![7] K Rossi, F Baletto, The effect of chemical ordering and lattice mismatch on structural transitions in phase segregating nanoalloys
! Physical Chemistry Chemical Physics 19 (18), 11057-11063
![8] K Rossi, L Pavan, YY Soon, F Baletto, The effect of size and composition on structural transitions in monometallic nanoparticles
! The European Physical Journal B 91 (2), 1-8
!
!@ Coalescence (soon to be added in the public version) 
! MSci thesis by Matteo Tiberi
!---------------------------------------------------------------------

PROGRAM main
  USE PARACLUSTER 
  USE CLUSTER     
  USE POTENTIAL
  USE ENFORCE
  USE DISTANCE
  USE META
  USE SUBSTRATE
  USE ENVIRONMENT
  USE module_sticky
  USE OMP_LIB     ! to have the execution time in the output
  USE module_function_var

  IMPLICIT NONE
  INTEGER :: is, i, nat, j, l_cv
  INTEGER :: icolor
  REAL(8) :: xuar,yvar,zwar
  REAL(8) :: vargau
  LOGICAL :: no_cluster
  REAL(8) :: time_start, time_end, t_seconds
  INTEGER :: t_days, t_hours, t_minuts

  ! Default values
  writeheader         = .true. !  To write the header only the first time in output growth
  sticky_atoms_wanted = .false.
  yesno_wanted        = .false.
  cutoffmat_req       = .false.
  d_com_wanted        = .false.
  error_counter       = 0
  nome                = 0

  OPEN(UNIT=uniterr, FILE='error.out', STATUS='replace')

  !IN THIS SUBROUTINE DEFINE THE VARIABLES FROM INPUT FILE
  WRITE (*,*) 'main_MD> Reading input file'
  CALL read_input
  !--------------------------------------------------------
  ! Read type_process and irand
  ! check what it is read as type_process
  !

  !---------------------------------------------------------
  ! Read the parameters for MgOsubstrate if wanted
  IF (mgo_substrate) THEN
     IF (metal_on_top) THEN
        WRITE (*,*) 'main_MD> Reading MgO substrate (MOT) parameters'
        CALL mgo_read_mot
     ELSE
        WRITE (*,*) 'main_MD> Reading substrate parameters'
        CALL mgo_read
     ENDIF
  ENDIF
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Read list of sticky_atoms if any
  IF (sticky_atoms.GT.0) THEN
     IF (sticky_atoms.GT.natom) THEN
        WRITE(*,*) 'main_MD> Error: Number of sticky_atoms is too large'
        STOP
     ENDIF
     WRITE (*,*) 'main_MD> There are ', sticky_atoms, 'atoms bound to their starting positions'
     WRITE (*,*) 'main_MD> by a spring of elastic constant [eV/A^2]: ', sticky_k

     WRITE (*,*) 'main_MD> Allocating variables for sticky_atoms'
     ALLOCATE (sticky_labels(sticky_atoms))
     ALLOCATE (sticky_fx(natom), sticky_fy(natom), sticky_fz(natom))

     WRITE (*,*) 'main_MD> Reading file ./sticky.list'
     OPEN(UNIT=19, FILE='sticky.list', STATUS='old')
     READ(19,*) sticky_labels
     CLOSE(19)
     WRITE (*,*) 'main_MD> They are labelled: ', sticky_labels

     ! Conversion of elastic constant from eV/A^2 to eV/arete(1)^2
     sticky_k = sticky_k *(arete(1)*arete(1))

     sticky_atoms_wanted = .true.

     ! This is done in main_MD because always only the same atoms have theese
     ! quantities different from zero
     sticky_fx = 0.d0
     sticky_fy = 0.d0
     sticky_fz = 0.d0
  ENDIF
  !---------------------------------------------------------

  IF (canonical == 'ya') THEN
     metaperiod = npas ! LP: metadyn never start
  ENDIF
  !
  !THE CALL TO choice_pot IS DONE IN READ_INPUT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  IF(quenching=='ya')  THEN
     nome=npas/scrivo
  ENDIF
  !
  IF(canonical=='ya')  THEN
     nome=npas/scrivo
  ENDIF
  IF(coalescence=='ya')  THEN
     nome=npas/scrivo
  ENDIF
  !
  IF(deposizione=='ya') THEN
     nome=ndepmax*npas/scrivo
  ENDIF
  !
  IF(caloric=='ya')   THEN
     npassifreezing=ABS(tcaloric-tinit)/deltat ! I moved this here from down there...
     nome=npassifreezing*npas/scrivo
  ENDIF
  IF(metadyn=='ya') THEN
      nome = npas / metaframe
  END IF
  !
  !=========================================================================
  ! This has always to be natom (but in the case of deposizione)
  IF(deposizione=='ya') THEN
  nat = natom + ndepmax
  ELSE
  nat=natom 
  ENDIF
  ! nat is the max number of atoms, not natom, witch is the initial
  !=========================================================================

  !-------------------------------------------------------------------------
  ! ALLOCATION OF GLOBAL VARIABLES
  WRITE (*,*) 'main_MD> Allocating global variables'
  !
  IF ((caloric=='ya') .OR. (quenching=='ya') .OR. (deposizione=='ya') .OR. (canonical=='ya') .OR. (coalescence == 'ya')) THEN
     ALLOCATE(filename(nome))
     WRITE(*,*) 'main_MD> Number of outputs that will be created: ', nome
  ENDIF
  ALLOCATE(x(nat),y(nat),z(nat))
  ALLOCATE(fx(nat),fy(nat),fz(nat))
  ALLOCATE(dfx(nat),dfy(nat),dfz(nat))
  ALLOCATE(vx(nat),vy(nat),vz(nat))
  ALLOCATE(u(nat),v(nat),w(nat))
  ALLOCATE(du(nat),dv(nat),dw(nat))
  ALLOCATE(nvois(nat),ivois(nvsiz,nat))
  ALLOCATE(potener(nat))
  IF ((canonical == 'ya') .and. (vel_af) ) THEN
      allocate(v_acf(npas,natom,3))
      allocate(vel_act_est(npas))
  end if

  ! NEIGHBOURS FOR METALLIC SYSTEMS
  !
  IF(type_potential=='rgl') THEN
     ALLOCATE(nv4(nat),iv4(nvsiz,nat))
     WRITE(*,*)'main_MD> Size of nv and nvoiz:',SIZE(nv4),SIZE(iv4)
     WRITE(*,*)'main_MD> Size of nv and nvoiz:',SIZE(nv4),SIZE(iv4)
  ENDIF

  ! Check for correct allocation
  WRITE(*,*) 'main_MD> Number of time-steps (npas)   =', npas
  WRITE(*,*) 'main_MD> Number of CVs        (num_cv) =', num_cv
  WRITE(*,*) 'main_MD> Max allowed number of atoms   =', nsiz
  WRITE(*,*) 'main_MD> Max number of atoms           =', nat
  IF (nat.GT.nsiz) THEN
     WRITE(*,*) 'main_MD> Error: too many atoms!'
     WRITE(*,*) 'main_MD> Change parameter nsiz in module_parameters.f90 and recompile the code'
     STOP
  ENDIF
  IF (2*nat.LT.nsiz) THEN
     WRITE(*,*) 'main_MD> Consider changing parameter nsiz in module_parameters.f90 for performance'
  ENDIF
  !--------------------------------------------------------------------------------
  ! Check for the part of code that will use the subroutine pair_distances and more
  IF (((metadyn .EQ. 'ya').OR.(collvar_wanted).OR.impl_env) .AND.(clusters .NEQV. .true.)) THEN
     ! Check on the calculation of distances for periodic systems
     write (*, *) 'main_MD> Calculation for no cluster system is not implemented'
     write (*, *) 'main_MD> Subroutine pair_distances is only for clusters'
     stop
  END IF
  IF (collvar_wanted.AND.(deposizione=='ya')) THEN
     WRITE(*,*) 'main_MD> The calculation of CVs and Growth have not been tested'
     WRITE(*,*) 'main_MD> Since some errors might occour, the program is going to stop'
     STOP
  END IF
  IF (mgo_substrate.AND.(deposizione=='ya')) THEN
     WRITE(*,*) 'main_MD> Growth upon a substrate has not been tested'
     WRITE(*,*) 'main_MD> Since some errors might occour, the program is going to stop'
     STOP
  END IF
  IF ((collvar_wanted.OR.(metadyn.EQ.'ya').OR.impl_env).AND.(nat.NE.natom)) THEN
     WRITE(*,*) 'main_MD> Error in variables allocations! nat not equal to natom'
     STOP
  ENDIF
  !--------------------------------------------------------------------------------
  IF ((collvar_wanted).OR.(metadyn=='ya').OR.mgo_substrate.OR.impl_env) THEN
     WRITE (*,*) 'main_MD> Allocating global variables to calculate pair distances'
     ALLOCATE(pair_dist (nat, nat))
     ALLOCATE(xyz_dist (3, nat, nat))
  ENDIF
  IF (mgo_substrate) THEN
     WRITE (*,*) 'main_MD> Allocating global variables for substrate'
     ALLOCATE( mgo_cn(nat), mgo_s(nat, nat), mgo_dS(nat, nat) )
     ALLOCATE( mgo_dcn_dx(nat), mgo_dcn_dy(nat), mgo_dcn_dz(nat) )
     ALLOCATE( mgo_x(nat), mgo_y(nat), mgo_z(nat) )
     ALLOCATE( mgo_fx(nat), mgo_fy(nat), mgo_fz(nat) )
     IF (metal_on_top) THEN
        ALLOCATE( mgo_cn_mot(nat), mgo_s_mot(nat, nat), mgo_dS_mot(nat, nat) )
        ALLOCATE( mgo_dcn_mot_dx(nat), mgo_dcn_mot_dy(nat), mgo_dcn_mot_dz(nat) )
     ENDIF
  ENDIF

  IF (impl_env) THEN
     WRITE (*,*) 'main_MD> Allocating global variables for implicit environment'
    ALLOCATE( env_cn(nat), ener_env_atom(nat), env_s(nat,nat), env_dS(nat,nat) )
    ALLOCATE( env_dcn_dx(nat), env_dcn_dy(nat), env_dcn_dz(nat) )
    ALLOCATE( env_dS_dx(nat,nat), env_dS_dy(nat,nat), env_dS_dz(nat,nat) )
    ALLOCATE( env_fx(nat), env_fy(nat), env_fz(nat))
  END IF

  WRITE(*,*) 'main_MD> collvar_wanted =' ,collvar_wanted
  !IF ((quenching=='ya').OR.(canonical=='ya').OR.(metadyn=='ya')) THEN
  IF ((collvar_wanted).OR.(metadyn=='ya').OR.impl_env) THEN
     WRITE (*,*) 'main_MD> Allocating global variables for Metadynamics and calculation of CVs...'
     ALLOCATE( s_of_t(num_cv,(npas/metaperiod)) ) ! Centres of gaussians (1 or many-dimensional)
     s_of_t(:,:) = 0.d0        ! Initialization
     ALLOCATE(pairkindmat (nat, nat))    ! These are used only with meta and
     ALLOCATE(cutoffmat (nat,nat))       ! colvar_wanted !!
     ! Preassigning a value to all the entries of cutoffmat (useful for elements in the diagonal)
     cutoffmat = .FALSE.

     ALLOCATE(coll_var(num_cv))
     ALLOCATE(gausswidth(num_cv))
     ALLOCATE(deno(num_cv))
     ALLOCATE(halfdeno(num_cv))
     ALLOCATE(expoc(num_cv))
     ALLOCATE(dVg_ds(num_cv))
     collvar_case(:) = 0
     WRITE(*,*) 'main_MD> Gaussian height [eV]=', gheight
     WRITE(*,*) 'main_MD> Meta-period [steps] =', metaperiod
     DO l_cv= 1, num_cv
        IF (TRIM(collvar_name(l_cv)).EQ.'coord_number') THEN
           WRITE(*,*) 'main_MD> gaussian width (CN) =', gwidth
           WRITE(*,*) 'main_MD> Allocating variables for CN (CV case 1)'
           collvar_case(l_cv) = 1
           gausswidth(l_cv) = gwidth
           ALLOCATE(s(nat, nat))
           ALLOCATE(CNatom(nat))
           ALLOCATE(fgx(nat), fgy(nat), fgz(nat))
           ALLOCATE(dS_dx(nat), dS_dy(nat), dS_dz(nat))
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'C2N') THEN
           WRITE(*,*) 'main_MD> gaussian width (C2N) =',g2nwidth
           WRITE(*,*) 'main_MD> Allocating variables for C2N (CV case 2)'
           collvar_case(l_cv) = 2
           gausswidth(l_cv) = g2nwidth
           ALLOCATE(s2n(nat, nat))
           ALLOCATE(C2Natom(nat))
           ALLOCATE(fg2nx(nat), fg2ny(nat), fg2nz(nat))
           ALLOCATE(dS2n_dx(nat), dS2n_dy(nat), dS2n_dz(nat))
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'SFN') THEN
           WRITE(*,*) 'main_MD> gaussian width (SFN) =',gsfwidth
           WRITE(*,*) 'main_MD> Allocating variables for SFN (CV case 3)'
           collvar_case(l_cv) = 3
           gausswidth(l_cv) = gsfwidth
           ALLOCATE(ssf(nat, nat))
           ALLOCATE(SFNatom(nat))
           ALLOCATE(fgsfx(nat), fgsfy(nat), fgsfz(nat))
           ALLOCATE(dSsf_dx(nat), dSsf_dy(nat), dSsf_dz(nat))
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'CNN') THEN
           WRITE(*,*) 'main_MD> gaussian width (CNN) =', gcnnwidth
           WRITE(*,*) 'main_MD> Allocating variables for CNN (CV case 4)'
           collvar_case(l_cv) = 4
           gausswidth(l_cv) = gcnnwidth
           cutoffmat_req = .TRUE.
           ALLOCATE(fgcnnx(nat), fgcnny(nat), fgcnnz(nat))
           ALLOCATE(dScnn_dx(nat), dScnn_dy(nat), dScnn_dz(nat))
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'CN_bim') THEN
           WRITE(*,*) 'main_MD> gaussian width (CN_bim) =', cn_gwidth
           WRITE(*,*) 'main_MD> Allocating variables for CN_bim (CV case 5)'
           collvar_case(l_cv) = 5
           gausswidth(l_cv) = cn_gwidth
           yesno_wanted = .true.
           ALLOCATE(cn_s(nat, nat))
           ALLOCATE(cn_CNatom(nat))
           ALLOCATE(cn_fgx(nat), cn_fgy(nat), cn_fgz(nat))
           ALLOCATE(cn_dS_dx(nat), cn_dS_dy(nat), cn_dS_dz(nat))
           ALLOCATE(pair_yesno_mat(nat,nat))
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'d_com') THEN
           WRITE(*,*) 'main_MD> gaussian width (d_com) =', d_com_gwidth
           WRITE(*,*) 'main_MD> Allocating variables for d_com (CV case 6)'
           collvar_case(l_cv) = 6
           gausswidth(l_cv) = d_com_gwidth
           d_com_wanted = .true.
           ALLOCATE(d_com_dcv_dx(nat), d_com_dcv_dy(nat), d_com_dcv_dz(nat))
           ALLOCATE(d_com_fgx(nat), d_com_fgy(nat), d_com_fgz(nat))
        ELSE
           WRITE(*,*) 'main_MD> Error in the CVs names'
           STOP
        ENDIF
        !-----------------------------------------------------------
        ! Allocating dS, used in the calculation of CVs, cases 1,2,3,5 not 4 and 6
        ! but not when the onlylight version is wanted
        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  REMEMBER TO ADD NEW CASES!
        IF (metadyn .EQ. 'ya') THEN
           IF ( (collvar_case(l_cv).LE.3).OR.(collvar_case(l_cv).EQ.5) ) THEN
              IF (.NOT. ALLOCATED(dS)) ALLOCATE(dS(nat,nat))
           END IF
        END IF
        !-----------------------------------------------------------
     ENDDO
     ! Global variable 'deno', array operation
     deno = 2.d0* gausswidth**2
     halfdeno =   gausswidth**2
  ENDIF

  !===================================!
  ! Conversion from A to arete units  !
  !===================================!
  aretebim=arete(1)                   !
  DO is=1,natom                       !
     x(is)=x0(is)/aretebim            !
     y(is)=y0(is)/aretebim            !
     z(is)=z0(is)/aretebim            !
  ENDDO                               !
  !===================================!

  IF( clusters ) no_cluster = .FALSE.
  !
  IF( no_cluster ) THEN
     IF ( wires ) pbcz = (1.d0+dilat)*pbcz/aretebim
     !
     ! dilat is the coef. of linear expansion (due to thermal motion)
     !
     IF ( surface ) pbcx = (1.d0+dilat)*pbcx/aretebim
     IF ( surface ) pbcy = (1.d0+dilat)*pbcy/aretebim
     !
     IF ( bulk ) pbcx = (1.d0+dilat)*pbcx/aretebim
     IF ( bulk ) pbcy = (1.d0+dilat)*pbcy/aretebim
     IF ( bulk ) pbcz = (1.d0+dilat)*pbcz/aretebim
  ENDIF
  !
  ! FINAL TEMPERATURE AT THE BEGINNING OF THE SIMULATION
  !
  tfin=tinit
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! WRITING THE INPUT PARAMETRS/VARIABLES IN DEFINITION.OUT
  !
  OPEN(UNIT=3,file='definition.out',status='replace')
  WRITE(3,*) '########################'
  WRITE(3,*) 'Procedure  ', TRIM(type_process)
  WRITE(3,*) '########################'
  WRITE(3,*) 'PARAMETERS SIMULATIONS'
  WRITE(3,*) 'Random number seed  = ',irand
  WRITE(3,*) 'Step for each run   = ',npas
  WRITE(3,*) 'Writing is for each = ',scrivo
  !
  WRITE(3,*) 'DETAILS OF THERMALIZATION FOR ', npast, ' STEPS'
  WRITE(3,*) 'DETAILS OF THERMALIZATION AT  ', tinit, ' (k)'
  !
  !!!!!!!!!!PROCEDURE DETAILS
  !
  WRITE(3,*) 'PROCEDURE DETAILS'
  WRITE(3,*) '#####################################'
  !
  IF(quenching=='ya') WRITE(3,*) 'QUENCHING START AFTER ',itremp,' steps'
  IF(canonical=='ya') WRITE(3,*) 'FINITE TEMPERATURE RUN AT ',tfin,' (k)'
  IF(coalescence == 'ya') WRITE(3, *) 'Number of atoms in cluster two = ', natom2
  !
  IF(deposizione=='ya') THEN
     !
     WRITE(3,*) 'DEPOSITION OF ',ndepmax,' OVER A CORE OF ',natinizio, ' ATOMS '
     WRITE(3,*) 'AT TEMPERATURE ', tfin, 'T(k)'
     WRITE(3,*) 'DEPOSITION TYPE:',lcs
     !
     IF(lcs==1) THEN
        WRITE(3,*) 'Monometallic deposition of ',ndepmax,' of ',elem1
     ELSEIF(lcs==2) THEN
        WRITE(3,*) 'Bimetallic deposition of ',elem1,' AND ',elem2
        WRITE(3,*) 'Probability ',prob
        WRITE(3,*) 'IF randomperc>prob select ', elem1, ' IS DEPOSITED'
        WRITE(3,*)  'Stop to deposit ', elem2, 'after ',at_tipo2, ' atoms'
     ELSEIF(lcs==3) THEN
        WRITE(3,*) &
           &'core/shell deposition',elem2,' onto ',elem1
     ENDIF
  ENDIF
  !
  IF(caloric=='ya') THEN
     if(tcaloric > tinit)  WRITE(3,*) 'MELTING PROCESS:'
     if(tcaloric < tinit)  WRITE(3,*) 'FREEZING PROCESS:'
     WRITE(3,*) 'FROM T_init(K)=',tinit
     WRITE(3,*) 'UP TO T(K)=',tcaloric
     WRITE(3,*) 'CHANGING T OF =',deltat,' (k) EACH ',npas, ' STEPS'
  ENDIF
  !
  WRITE(3,*) '########################'
  WRITE(3,*) 'INFORMATION on the SYSTEM'
  WRITE(3,*) '########################'
  WRITE(3,*) 'INITIAL COMPOSITION: '
  WRITE(3,*) ntipo1,' ATOMS OF ',elem1,' AND ',ntipo2,' OF ', elem2
  !
  vargau=SQRT(2.d0*t2m(1)*cbol*tfin)
  IF(sys == 'bim' ) vargau=SQRT((t2m(1)+t2m(2))*cbol*tfin)
  WRITE(3,*) 'T =',tfin,' (k)', 'AND vargau=',vargau
  !
  IF( no_cluster ) THEN
     WRITE(3,*) 'WHICH PERIODICITY, in ARETE UNITS'
     WRITE(3,*) '########################'
     WRITE(3,*) 'ARETE', aretebim
     IF ( wires )   WRITE (3,*) 'NANOWIRES ALONG Z',pbcz
     IF ( surface ) WRITE (3,*) 'SURFACE ALONG x,y',pbcx,pbcy
     IF ( bulk )    WRITE (3,*) 'BULK',pbcx,pbcy,pbcz
  ENDIF
  !
  WRITE(3,*) 'In energy.out FILE_OUTPUT for energy WHERE there are:'
  WRITE(3,*) 'tempo,ener,etot,ecin,edelta,emedia,tpar'
  WRITE(3,*) 'In depo.out FILE_OUTPUT for deposited atoms only for growth'
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  WRITE(3,*) '########################'
  WRITE(3,*) 'DEFINITION OF POTENTIAL PARAMETERS:'
  WRITE(3,*) 'KIND OF POTENTIAL:', type_potential
  WRITE(3,*) '########################'
  !
  IF(type_potential =='rgl') THEN
     WRITE(3,*) 'RGL POTENTIAL PARAMETERS FOR:'
     IF(sys == 'mon' ) THEN
        WRITE(3,*) 'Metal ',leg(imet1)
        WRITE(3,'(a30,1x,f7.5)')'Lattice constant, arete =',arete(1)
        WRITE(3,'(a30,1x,f7.5)')'T2M quantity = ',t2m(1)
        WRITE(3,*)'Potential parameters: A',' ','qsi',' ','p',' ','q',' ','ecoh'
        WRITE(3,'(a2,1x,a2,5f12.7)') 'Metal A=',elem1,a(1),qsi(1),p(1),q(1),ecoh(1)
     !
     ELSEIF(sys == 'bim' ) THEN
        WRITE(3,*) 'Metal ',leg(imet1),' AND ', leg(imet2)
        WRITE(3,'(a30,1x,f7.5,a5,f7.5)')'Lattice constant, arete =',arete(1),' AND ',arete(2)
        WRITE(3,'(a30,1x,f7.5,a5,f7.5)')'T2M quantity = ',t2m(1),' AND ', arete(2)
        WRITE(3,*)'Potential parameters: A',' ','qsi',' ','p',' ','q',' ','ecoh'
        WRITE(3,'(a12,a2,1x,5f12.7)') 'Metal A = ',elem1,a(1),qsi(1),p(1),q(1),ecoh(1)
        WRITE(3,'(a12,a2,1x,5f12.7)') 'Metal B = ',elem2,a(2),qsi(2),p(2),q(2),ecoh(2)
        WRITE(3,'(a12,a2,1x,4f12.7)') 'Alloy   = ','AB',a(3),qsi(3),p(3),q(3)
     ENDIF
  !
  ELSEIF(type_potential =='lj1') THEN
     WRITE(3,*) 'LENNARD-JONES PARAMETERS for:',elem1,' and ',elem2
     WRITE(3,'(a6,1x,2f7.5)')'arete [A]=',aretebim,arete(1)
     WRITE(3,'(a2,1x,a2,1x,a4,e13.7,a1,f7.5)')&
          &  'A=',elem1,'t2m= ',t2m(1),' arete= ',arete(1)
     WRITE(3,*) 'DEPTH OF THE WELL (eV):', U0
     WRITE(3,*) 'DISTANCE BETWEEN NN (A):',R0
  !
  ELSEIF(type_potential =='gir') THEN
     WRITE(3,*) 'C ATOMS IN A BALL',Nball
     WRITE(3,*) 'DISTANCE BETWEEN NN (A):',Dnn
     WRITE(3,*) 'DIAMETR OF THE BALL (A):',dball
     WRITE(3,*) 'LJ ATTRACTIVE COSTANT (eV*A**6):',ALJ
     WRITE(3,*) 'LJ REPULSIVE COSTANT (eV*A**12)',BLJ
  !
  ENDIF

  WRITE(3,*)'IN MAIN',SIZE(x),SIZE(fx),SIZE(dfx),SIZE(u),SIZE(du),SIZE(vx)
  WRITE(3,*)'IN MAIN: number of position-output file',SIZE(filename)
  CLOSE(3)
  !
  WRITE(*,*) 'main_MD> INPUT AND POTENTIAL FILE READING IS DONE'
  WRITE(*,*) 'main_MD> THERMALIZATION IS STARTING'
  !
  !HERE I PREPARE THE INITIAL STATE AT TIME t=0
  !
  CALL preparo  !!preparo means 'preparing'
  !
  OPEN (UNIT=13,file='pr.out',status='unknown')
  WRITE(13,*) natom
  WRITE(13,'(a2,1x,a2,1x,f12.5)') elem1,elem2,etot
  DO i=1,natom
     icolor=itype(i)
     xuar=(x(i)+u(i))*aretebim
     yvar=(y(i)+v(i))*aretebim
     zwar=(z(i)+w(i))*aretebim
     WRITE(13,'(a3,1x,3f16.5,i4)') elem(i),xuar,yvar,zwar,icolor
  ENDDO
  CLOSE(13)

  !______________________________________________________________
  !                    End of initial phase
  !``````````````````````````````````````````````````````````````


  !==============================================================
  ! Initializing variables for metadynamics
  !==============================================================
  IF ((collvar_wanted).OR.(metadyn=='ya').OR.impl_env) THEN
     NG   = 0        ! Number of Gaussians added

     ! Creating the pair kind matrix
     WRITE (*,*) 'main_MD> Creating the Pair Kind Matrix (pairkindmat)...'
     DO i = 1, natom
        DO j = 1, natom
           ! Checking what kind of atoms the couple i,j is made of
           if      ((itype(i).eq.1).and.(itype(j).eq.1)) then
              pairkindmat(j,i) = 1  ! same metal A
           else if ((itype(i).eq.2).and.(itype(j).eq.2)) then
              pairkindmat(j,i) = 2  ! same metal B
           else
              pairkindmat(j,i) = 3  ! A-B interaction
           end if
        END DO
     END DO

     IF (cutoffmat_req) THEN
        ! Limits to compute the matrix cutoffmatwin for the window function in CNF
        inflim = win_in(1)-win_in(3)
        suplim = win_in(1)+win_in(3)

        ! Defining the values for CV parameters for bimetallic systems
        ! 'dist': nearest neighbours distances in arete(1) units,
        ! it is defined in bimet_rgl (read_rgl.f90), called in read_input
        ! dist(1) is relative to the case A-A
        ! dist(2) is relative to the case B-B
        ! dist(3) is relative to the case A-B

        ! CNF:
        ! parameter d0 for CNF in dependence of atom pair kind (Array Operation)
        cnf_in(6:8)   = dist * cnf_in(1) * DSQRT(2.d0)

        ! parameter r0 for CNF in dependence of atom pair kind (Array Operation)
        cnf_in(9:11)  = dist * cnf_in(2) * DSQRT(2.d0)

        ! parameter c0 for CNF in dependence of atom pair kind (Array Operation)
        cnf_in(12:14) = dist * cnf_in(3) * DSQRT(2.d0)

        ! Calculating the dcutoff for 1st coord sphere [bulk latt dist] (Arr. Op.)
        ! d0 + c0
        dcutoff = cnf_in(6:8) + cnf_in(12:14)
        WRITE(*,*) 'main_MD> cut-off distances used in the calculation of Common Neighbour Function...'
        WRITE(*,*) '- for pairs of kind A-A [bulk ref. dist.]:', dcutoff(1)
        WRITE(*,*) '- for pairs of kind B-B [bulk ref. dist.]:', dcutoff(2)
        WRITE(*,*) '- for pairs of kind A-B [bulk ref. dist.]:', dcutoff(3)
     ENDIF

     ! Calculating  pair_yesno_mat  for CN_bim
     IF (yesno_wanted) THEN
        WRITE (*,*) 'main_MD> Creating the matrix of pairs considered for CN_bim (pair_yesno_mat)'
        pair_yesno_mat = .false.
        DO i = 1, natom
           DO j = 1, natom
              IF(cn_aa .AND. (pairkindmat(j,i).EQ.1)) pair_yesno_mat(j,i) = .true.
              IF(cn_bb .AND. (pairkindmat(j,i).EQ.2)) pair_yesno_mat(j,i) = .true.
              IF(cn_ab .AND. (pairkindmat(j,i).EQ.3)) pair_yesno_mat(j,i) = .true.
           ENDDO
        ENDDO
     ENDIF

     ! Counting number atom A and B for d_com
     IF (d_com_wanted) THEN
        natoma = 0
        natomb = 0
        WRITE (*,*) 'main_MD> Counting number of atoms A and B,'
        DO i = 1, natom
           IF(itype(i).EQ.1) THEN
              natoma = natoma + 1
           ELSE
              natomb = natomb + 1
           ENDIF
        ENDDO
        IF ((natomb .EQ. 0).OR.(natoma .EQ. 0)) THEN
           WRITE(*,*) 'main_MD> Error: there are only atoms of one kind for d_com calculation!'
           STOP
        ENDIF
        WRITE (*,*) 'main_MD> Number of atoms of kind A: ', natomA
        WRITE (*,*) 'main_MD> Number of atoms of kind B: ', natomB
     ENDIF

  END IF
  !______________________________________________________________
  ! End of Initialization of variables for metadynamics
  !``````````````````````````````````````````````````````````````


  WRITE(*,*) 'main_MD> Time-loop is starting for ', TRIM(type_process)
  WRITE(*,*) 'main_MD> Now the computation of EXECUTION TIME starts'
  time_start =  OMP_GET_WTIME()
  !
  tpar = 0.d0     !averaged T
  edelta = 0.d0   !averaged energy
  nd_proc = 0     !for "caloric" is the step-number when increase/decrease T
                  !for "deposizione" is the number of deposited atom
                  !for quenching =1
  !
  nfile = 0       !for writing output-files
  !=================================================================
  !======================= START TIME-LOOP =========================
  !=================================================================
  IF((quenching=='ya').OR.(canonical=='ya').OR.(metadyn=='ya') .OR. (coalescence == 'ya') ) THEN
     ipas=0
     nd_proc=1
     OPEN(UNIT=unite,file='energy.out',status='unknown')
     !
     IF (metadyn.eq.'ya') THEN
        OPEN(UNIT=110, file='metahistory.out', status='unknown')
     ELSE
        OPEN(UNIT=unite,file='energy.out',status='unknown')
     END IF
     !
     !IF (restart_mode) THEN
     !  call read_state
     !ELSE
       IF((quenching=='ya').and.(ipas>itremp)) THEN 
       tfin=0.d0
       ENDIF
     !END IF

     !call initialize_movie

      IF ((collvar_wanted).or.(metadyn.eq.'ya').OR.mgo_substrate.OR.impl_env) THEN
             ! Calling the subroutine that calculate pair distances
             ! (It needs to be done for all the implemented CVs)
             CALL pair_distances
             IF ((collvar_wanted).or.(metadyn.eq.'ya')) THEN
             ! Calling the subroutine that calculate pair distances
             ! (It needs to be done for all the implemented CVs)
             CALL collective
             END IF
       END IF

     !=========================================
     ! Time-evolution of system:
     ! integration of Newton Equation of motion
     ! through a Velocity-Verlet algorthm
     !=========================================
     CALL ev_time
     !_________________________________________

     !!if (save_progress) call save_state(.true.)

     IF (metadyn.eq.'ya') then
        CLOSE(110)
     ELSE
        close(unite)
     END IF

     !call close_movie
  ENDIF
  !
  IF(caloric=='ya') THEN

     ! Injecting saved state
     !if (restart_mode) CALL read_state

     !call initialize_movie

     WRITE(*,*) 'main_MD> Running from loop ', start_from_nd_proc, ' to ', npassifreezing
     !
     DO nd_proc = start_from_nd_proc, npassifreezing
        !=========================================
        !Defining the new temperature
        !=========================================
        CALL thermodynamics
        !_________________________________________
        write(*,*) 'main_MD> interactions number', nd_proc, 'at Temperature=',tfin

        !=========================================
        ! Time-evolution of system:
        ! integration of Newton Equation of motion
        ! through a Velocity-Verlet algorthm
        !=========================================
        CALL ev_time
        !_________________________________________
        !if (periodic_save) call save_state(.false.)
     ENDDO

     !if (save_progress) call save_state(.true.)
     !call close_movie
  ENDIF
  !
  IF (deposizione=='ya') THEN

     OPEN(UNIT=unitd,file='depo.out',status='replace')
     initntipo2 = ntipo2
     DO nd_proc=1,ndepmax
        !=========================================
        !Subroutine to create a new atoms which is
        ! going to the cluster (deposition)
        !=========================================
        CALL growth
        !_________________________________________
        !call initialize_movie
        !=========================================
        ! Time-evolution of system:
        ! integration of Newton Equation of motion
        ! through a Velocity-Verlet algorthm
        !=========================================
        CALL ev_time
        !_________________________________________
        !if (periodic_save) call save_state(.false.)
        !call close_movie
     ENDDO
  ENDIF
  !

  !
  CLOSE(unitd)
  IF ((canonical == 'ya') .and. (vel_af) ) THEN
  call acf_routine
  call dft_routine
  end if


  !=================================================
  ! Deallocation of global variables
  !=================================================
  time_end = OMP_GET_WTIME()
  t_days    = INT((time_end-time_start)/(3600*24))
  t_hours   = INT((time_end-time_start)/3600) - t_days*24
  t_minuts  = INT((time_end-time_start)/60) - (t_days*24 + t_hours)* 60
  t_seconds = (time_end-time_start) - (t_days*24*60 + t_hours*60 + t_minuts)* 60
  WRITE(*,'(A,X,3(2X,I3,X,A),2X,F10.7,X,A)') 'main_MD> EXECUTION TIME:',&
                   & t_days, '[day],', t_hours, '[hour],', t_minuts, '[min],', t_seconds, '[s].'

  IF ((collvar_wanted).OR.(metadyn=='ya').OR.mgo_substrate.OR.impl_env) THEN
     WRITE (*,*) 'main_MD> Deallocating global variables to calculate pair distances'
     DEALLOCATE(pair_dist)
     DEALLOCATE(xyz_dist)
  ENDIF
  IF (mgo_substrate) THEN
     WRITE (*,*) 'main_MD> Deallocating global variables for MgO substrate'
     DEALLOCATE( mgo_cn, mgo_s, mgo_dS )
     DEALLOCATE( mgo_dcn_dx, mgo_dcn_dy, mgo_dcn_dz )
     DEALLOCATE( mgo_x, mgo_y, mgo_z )
     DEALLOCATE( mgo_fx, mgo_fy, mgo_fz )
     IF (metal_on_top) THEN
        DEALLOCATE( mgo_cn_mot, mgo_s_mot, mgo_dS_mot )
        DEALLOCATE( mgo_dcn_mot_dx, mgo_dcn_mot_dy, mgo_dcn_mot_dz )
     ENDIF
  ENDIF

  IF (sticky_atoms.GT.0) THEN
     WRITE (*,*) 'main_MD> Deallocating global variables for sticky_atoms'
     DEALLOCATE (sticky_fx, sticky_fy, sticky_fz)
     DEALLOCATE (sticky_labels)
  ENDIF

  IF ((collvar_wanted).OR.(metadyn=='ya').OR.impl_env) THEN
     WRITE (*,*) 'main_MD> Deallocating global variables for Metadynamics and calculation of CVs'
     DO l_cv= 1, num_cv
        IF (TRIM(collvar_name(l_cv)).EQ.'coord_number') THEN
              DEALLOCATE(s)
              DEALLOCATE(CNatom)
              DEALLOCATE(fgx, fgy, fgz)
              DEALLOCATE(dS_dx, dS_dy, dS_dz)
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'C2N') THEN
              DEALLOCATE(s2n)
              DEALLOCATE(C2Natom)
              DEALLOCATE(fg2nx, fg2ny, fg2nz)
              DEALLOCATE(dS2n_dx, dS2n_dy, dS2n_dz)
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'SFN') THEN
              DEALLOCATE(ssf)
              DEALLOCATE(SFNatom)
              DEALLOCATE(fgsfx, fgsfy, fgsfz)
              DEALLOCATE(dSsf_dx, dSsf_dy, dSsf_dz)
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'CNN') THEN
              DEALLOCATE(fgcnnx, fgcnny, fgcnnz)
              DEALLOCATE(dScnn_dx, dScnn_dy, dScnn_dz)
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'CN_bim') THEN
              DEALLOCATE(cn_s)
              DEALLOCATE(cn_CNatom)
              DEALLOCATE(cn_fgx, cn_fgy, cn_fgz)
              DEALLOCATE(cn_dS_dx, cn_dS_dy, cn_dS_dz)
              DEALLOCATE(pair_yesno_mat)
        ELSE IF (TRIM(collvar_name(l_cv)).EQ.'d_com') THEN
              DEALLOCATE(d_com_dcv_dx, d_com_dcv_dy, d_com_dcv_dz)
              DEALLOCATE(d_com_fgx, d_com_fgy, d_com_fgz)
        ENDIF
     ENDDO
     DEALLOCATE(s_of_t)
     DEALLOCATE(coll_var)
     DEALLOCATE(gausswidth)
     DEALLOCATE(deno, halfdeno, expoc, dVg_ds)
     DEALLOCATE(cutoffmat, pairkindmat)
  END IF
  WRITE (*,*) 'main_MD> Deallocating global variables'
  IF (ALLOCATED(dS)) DEALLOCATE(dS)
  DEALLOCATE(x,y,z)
  DEALLOCATE(fx,fy,fz)
  DEALLOCATE(dfx,dfy,dfz)
  DEALLOCATE(vx,vy,vz)
  DEALLOCATE(u,v,w)
  DEALLOCATE(du,dv,dw)
  DEALLOCATE(nvois,ivois)
  IF(type_potential=='rgl')  DEALLOCATE(nv4,iv4)
  IF( (caloric=='ya').OR.(quenching=='ya').OR.(deposizione=='ya').OR.(canonical=='ya')&
                                          .OR.(coalescence == 'ya') ) DEALLOCATE(filename)
IF ((canonical == 'ya') .and. (vel_af) ) THEN
    deallocate(v_acf)
    deallocate(vel_act_est)
end if

WRITE (*,*) 'main_MD> This process has finished'

END PROGRAM main
