SUBROUTINE read_input
!====================================
! The subroutine reads the input file
!====================================

USE PARACLUSTER  ! module with parameters
USE CLUSTER      ! module of global variables
USE POTENTIAL
USE ENFORCE
USE DISTANCE
USE SUBSTRATE
USE META         ! metadynamics module
USE module_sticky
USE ENVIRONMENT

IMPLICIT NONE

 ! Local variables
 REAL    :: mgo_mindist ! min distance from cluster and substrate MgO
 INTEGER :: ind,ios,l_cv
 CHARACTER(LEN=256) :: input_line, card
 LOGICAL :: tpos, end_file
 INTEGER :: nproc            ! number of chosen processes
 INTEGER :: metaperiod_def   ! default value for metaperiod (for any CVs)
 REAL :: comx0, comy0, comz0, com2x0, com2y0, com2z0, distcom
 !
 TYPE :: atoms
 CHARACTER*2 :: spx(nsiz)
 REAL        :: tau(3,nsiz)
 INTEGER     :: icolor
 ENDTYPE atoms
 TYPE(atoms) :: lb_at

NAMELIST/simul/irand,tstep,npas,scrivo,npast,tinit,&
              &type_process,mgo_substrate,metal_on_top,&
              &vnu,filepos,filepot,mgo_pot,sticky_atoms,sticky_k,cn_cal,&
              &impl_env, pot_a, pot_b, eta_a, eta_b

NAMELIST/system/type_potential,natom,fattor,elem1,elem2,clusters,wires,&
               &surface,bulk,pbcz,pbcy,pbcx,sys,lcutoff
!NAMELIST/support/sub_geom
NAMELIST/sysbim/sologeom,NA_elem
NAMELIST/canon/vel_af
NAMELIST/growth/ndepmax,lcs,at_tipo2,prob,tsorg,rad
NAMELIST/quench/itremp,tmin
NAMELIST/calor/tcaloric,deltat,tfin
NAMELIST/coal/filepos2,natom2,somedist

!----- beginning of new part (Metadynamics) -----------------------------------
NAMELIST/metalist/collvar_wanted, metaframe, gheight, metaperiod, num_cv, collvar_name, &
                & n_pwr, m_pwr, rzero, gwidth, &
                & n2n_pwr, m2n_pwr, d2n, r2n, g2nwidth, &
                & nsf_pwr, msf_pwr, dsf, rsf, gsfwidth, &
                & cnf_in, win_in, gcnnwidth, &
                & cn_aa, cn_bb, cn_ab, cn_n_pwr, cn_m_pwr, cn_rzero, cn_gwidth, &
                & d_com_gwidth
!----- end of new part --------------------------------------------------------

!----- For restarting ---
NAMELIST/restart/restart_mode, save_progress, restart_from_file, start_from, start_from_nd_proc, &
    & progress_file_label, periodic_save, periodic_save_every, periodic_save_last

!SIMUL VARIABLES
irand = 1
tstep = 0.d0
npas= 0
scrivo= 0
npast = 0
tinit =1.d0
tfin = tinit
type_process = 'no'
mgo_substrate = .false.
metal_on_top = .false.
vnu = 5.0d11
filepos = 'pos.in'
filepot = 'no-pot-name'
mgo_pot = 'no-mgo-name'
sticky_atoms = 0
sticky_k = 10
cn_cal = .false.
impl_env = .false.
pot_a = 0.d0
pot_b = 0.d0
eta_a = 0.d0
eta_b = 0.d0
output_xyz = .false.
!
!SYSTEM
natom= 0
type_potential='rgl'
elem1 = 'no'
elem2=elem1
fattor=1.d0
clusters = .TRUE.
wires =.FALSE.
surface = .FALSE.
bulk=.FALSE.
pbcx= 0.d0
pbcy= 0.d0
pbcz= 0.d0
sys='mon'
lcutoff = .true.

!
!SUPPORT
!sub_geom = 0
!
!COAL
natom2 = 0
filepos2 = 'pos2.in'
somedist = 0.0d0
!
!QUENCH
itremp = 1
tmin = 1.d-8
!
!GROWTH
ndepmax = 0
lcs = 0
at_tipo2 = 0
prob=0.5d0
tsorg=1500.d0
rad=6.d0
!
!CALOR
tcaloric= tinit
deltat=0.d0
!
!SYSBIM
sologeom = .FALSE.
NA_elem = elem1
!CANON
vel_af = .false.
!------------------
! METALIST
!------------------
metaperiod_def = 1000
metaframe = metaperiod_def
collvar_wanted = .FALSE.
gheight = 0.3d0
metaperiod = metaperiod_def
num_cv = 1
collvar_name(1) = 'coord_number'
collvar_name(2) = 'none'

! CN
n_pwr  = 6
m_pwr  = 12
rzero  = 0.147d0
gwidth  = 2.d0

! Parameters for C2N (number of 2nd neighbours)
n2n_pwr         = 6
m2n_pwr         = 12
d2n             = 1.d0          ! [bulk latt parameter]
r2n             = 0.07d0        ! [bulk latt parameter] distances distribution parameter
g2nwidth        = 1.d0          ! Gaussian width

! Parameters for SFN (Stacking fault number)
nsf_pwr         = 6
msf_pwr         = 12
dsf             = 1.354d0       ! [bulk latt parameter]
rsf             = 0.05d0        ! [bulk latt parameter] distances distribution parameter
gsfwidth        = 1.d0          ! Gaussian width

! CNN (Common Neighbour No.) - Common Neighbour Function
cnf_in(1)      = 0.70710678d0   ! d0 : 1st coord sphere [latt bulk ref]
cnf_in(2)      = 0.1d0          ! r0 : 1st coord sphere [latt bulk ref]
cnf_in(3)      = 0.4d0          ! c0 : 1st coord sphere [latt bulk ref]
cnf_in(4)      = 6              ! m  : 1st coord sphere
cnf_in(5)      = 12             ! n  : 1st coord sphere
win_in(1)      = 5.d0           ! d0 : window function  [latt bulk ref]
win_in(2)      = 0.35d0         ! r0 : window function  [latt bulk ref]
win_in(3)      = 0.70d0         ! c0 : window function  [latt bulk ref]
win_in(4)      = 5              ! m  : window function
win_in(5)      = 8              ! n  : window function
gcnnwidth      = 0.5d0          ! Gaussian width

! CN_bim (CN bimetallic)     - Sigmoid Function bimetallic
 cn_aa         = .false.        ! pairs elem1-elem1
 cn_bb         = .false.        ! pairs elem2-elem2
 cn_ab         = .true.         ! mixed pairs
 cn_n_pwr      = 6              ! n
 cn_m_pwr      = 12             ! m
 cn_rzero      = 0.147d0        ! r0 [latt bulk ref]
 cn_gwidth     = 2.d0           ! Gaussian width

! d_com (squared distance CoM atoms A and CoM atoms B)
 d_com_gwidth  = 2.d0           ! Gaussian width
!----------------

!POSITION
lb_at%spx(:)='no'
lb_at%tau(1:3,:)=0.d0

!
!read namelist
  READ(UNIT=5, nml=simul, iostat=ios)
  READ(UNIT=5, nml=system, iostat=ios)
!  IF(mgo_substrate) READ(UNIT=5, nml=support, iostat=ios)
!=======================================
!    control on system/ choice
!=======================================
IF(natom>nsiz) THEN
 WRITE(*,*) 'read_input> Error: too many atoms.'
 WRITE(*,*) 'read input> Change parameter nsiz in module_parameters.f90 and recompile the code'
 WRITE(uniterr,*) "error: too many atoms"
 STOP
ENDIF
!
IF(scrivo>npas) THEN
 WRITE(*,*) "read_input> Error: the output cannot be written"
 WRITE(uniterr,*) "read_input> Error: the output cannot be written"
 STOP
ENDIF

IF(natom==0) THEN
 WRITE(*,*) 'read_input> NO ATOMS GIVEN',natom
 WRITE(*,*) 'read_input> ######################'
 WRITE(uniterr,*) 'NO ATOMS GIVEN',natom
 WRITE(*,*) 'read_input> ######################'
 STOP
ENDIF
!
!
!----------------------------------------------------
! RANDOM NUMBER GENERATOR
WRITE(*,*) 'read_input> Initializing random number generator'
WRITE(*,*) 'read_input> Checking if irand is odd', MOD(irand,2)
IF((irand<4095).AND.(irand>1000).AND.MOD(irand,2)/=0)THEN
    irand_seed(1)=irand
    irand_seed(2)=irand_seed(1)+1
    irand_seed(3)=irand_seed(2)+1
    irand_seed(4)=2*irand_seed(3)+1
ELSE
    WRITE(uniterr,*) 'wrong irand, irand should be odd, current irand is=',irand
    STOP
ENDIF

! Checking program/type_process
type_process = TRIM(type_process)
IF (type_process == 'no' ) THEN
   WRITE(*,*) 'read_input> No type_process selected'
   WRITE(*,*) 'read_input> ######################'
   WRITE(uniterr,*) 'No type_process selected'
   WRITE(uniterr,*) '######################'
   STOP
ELSE
!default values no process selected
!---------------------------------------------
   quenching = 'no'
   deposizione = 'no'
   caloric = 'no'
   canonical = 'no'
   metadyn = 'no'
   coalescence = 'no'
!
! assigning character label for running the code
!
   WRITE(*,*) 'read_input> Process type: ', TRIM(type_process)
!
   IF(type_process == 'Quenching' .or. type_process=='quenching') quenching = 'ya'
   IF(type_process == 'microcan' .or. type_process == 'NVE' .or. type_process == 'nve' .or. type_process =='microcanonical') then
       canonical = 'ya'
       vnu = 0.d0
       tfin = tinit
       WRITE(*,*) 'read_input> Selected NVE at temperature', tinit, 'thermostat freq.',vnu
     !!! check that tinit is given
     ENDIF
   IF(type_process == 'NVT'.or. type_process=='nvt'.or. type_process=='canonical') then
         canonical = 'ya'
         tfin = tinit
         write(*,*) 'read_input> NVT at temperature',tinit,'freq thermostat',vnu
   ENDIF
   IF(type_process == 'itMD'.or. type_process=='itmd' )  then
     caloric = 'ya'
     READ(UNIT=5, nml=calor, iostat = ios) 
     Write(*,*) 'WARINING in read_input > itMD used without melting/freezing'
 ENDIF
   IF(type_process == 'melting'.or. type_process=='melt') then
        caloric='ya'
        READ(UNIT=5, nml=calor, iostat = ios) 
        write(*,*) 'tinit is=',tinit,' tcaloric is =',tcaloric,''
        IF (tcaloric < tinit) THEN
         write(*,*) 'read_input> Error: tcaloric in melting must be > tinit'
         stop
        ENDIF
   ENDIF
   IF(type_process == 'freezing'.or. type_process=='freeze') then
        caloric='ya'
        READ(UNIT=5, nml=calor, iostat = ios) 
        write(*,*) 'tinit is=',tinit,' tcaloric is =',tcaloric, ''
        IF (tcaloric > tinit) THEN
         write(*,*) 'read_input> Error: tcaloric in freezing must be < tinit'
         stop
        ENDIF
    ENDIF

   IF(type_process == 'Growth'.or. type_process=='growth' )   deposizione = 'ya'
   IF(type_process == 'Metadynamics'.or. type_process=='metadynamics' .or. type_process=='metadyn' ) metadyn = 'ya'
   IF(type_process == 'Coalescence'.or. type_process=='coalescence' .or. type_process=='coal' ) coalescence = 'ya'
ENDIF
!
!
IF(quenching=='ya') READ(UNIT=5, nml=quench, iostat = ios)
IF(canonical=='ya') THEN
   READ(UNIT=5, nml=quench, iostat = ios) 
   READ(UNIT=5, nml=canon, iostat=ios)
END IF
IF(deposizione =='ya') READ(UNIT=5, nml=growth, iostat = ios)
IF (coalescence == 'ya') READ(UNIT = 5, nml = coal, iostat = ios)

!!!!!!!!!!!!!!! ====================================
! Instruction for a substrate implemented in &system
!!!!!!!!!!!!!!! ====================================
IF (mgo_substrate) THEN
  IF (mgo_pot == 'no-mgo-name') THEN
   WRITE(*,*) 'read_input> Error: No path for substrate parameters'
   STOP
  ENDIF
!  
!   
 IF (metal_on_top)  WRITE(*,*) 'read_input> Warning: metal_on_top (MOT) NOT TESTED'

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-------------------------------------------------------------------
! Metadynamics Begins
!-------------------------------------------------------------------
READ(UNIT=5, nml=metalist, iostat = ios) ! I want to read this always because I need collvar_wanted
WRITE(*,*) 'read_input> Input for metadynamics has been read'

IF ((collvar_wanted) .OR. (metadyn == 'ya'))   THEN
   IF (metadyn == 'ya') THEN
      WRITE(*,*) 'read_input> Metadynamics has been chosen'
   ELSE
      WRITE(*,*) 'read_input> Metadynamics has not been chosen, but CV values are going to be calculated'
   END IF

   !check on number of collective variable
   IF ((num_cv.ne.1).and.(num_cv.ne.2)) THEN
      WRITE(*, *) 'read_input> Error: invalid value for number of collective variables!'
      STOP
   ENDIF

   IF (TRIM(collvar_name(1)) == TRIM(collvar_name(2)) ) THEN
      WRITE(*,*) 'read_input> Warning: collective variables are identical, setting num_cv = 1'
      num_cv = 1
      collvar_name(2) = 'none'
   ENDIF

   DO l_cv= 1, num_cv
      IF((TRIM(collvar_name(l_cv)).EQ.'coord_number').OR.(TRIM(collvar_name(l_cv)).EQ.'C2N')&
           &.OR.(TRIM(collvar_name(l_cv)).EQ.'CNN').OR.(TRIM(collvar_name(l_cv)).EQ.'SFN') ) THEN
         write(*,*) 'read_input> CV chosen: ',TRIM(collvar_name(l_cv)),'.'
      ELSE IF(TRIM(collvar_name(l_cv)).EQ.'CN_bim') THEN
         !------------------------
         ! Checks on CN_bim
         WRITE(*,*) 'read_input> CV chosen: ',TRIM(collvar_name(l_cv)),','
         IF(cn_aa) WRITE(*,*) 'read_input> considering  A-A  pairs of atoms.'
         IF(cn_bb) WRITE(*,*) 'read_input> considering  B-B  pairs of atoms.'
         IF(cn_ab) WRITE(*,*) 'read_input> considering mixed pairs of atoms.'
         IF (sys .NE. 'bim') THEN
            WRITE(*,*) 'read_input> Warning: the system is monometallic and a chemical order CV has been selected'
            IF (.NOT.(cn_aa)) THEN
               WRITE(*,*) 'read_input> Error: No pairs of NNs would be counted in CN_bim!'
               STOP
            ENDIF
         ENDIF
         IF (cn_aa .AND. cn_bb .AND. cn_ab) THEN
            WRITE(*,*) 'read_input> Warning: All the pairs of NNs will be counted in CN_bim'
         ENDIF
         IF (.NOT.(cn_aa .OR. cn_bb .OR. cn_ab)) THEN
            WRITE(*,*) 'read_input> Error: No pairs of NNs would be counted in CN_bim!'
            STOP
         ENDIF
         !------------------------
      ELSE IF(TRIM(collvar_name(l_cv)).EQ.'d_com') THEN
         !------------------------
         ! Checks on d_com
         WRITE(*,*) 'read_input> CV chosen: ',TRIM(collvar_name(l_cv)),'.'
         IF (sys .NE. 'bim') THEN
            WRITE(*,*) 'read_input> Error: the system is monometallic, d_com is for bimetallic systems!'
            STOP
         ENDIF
         !------------------------
      ELSE
         write(*,*) 'read_input> Error: a CV has been selected, that is not implemented.'
         STOP
      ENDIF
   ENDDO

   IF (metaperiod >= npas) THEN
      metaperiod = npas
      WRITE(*,*) 'read_input> Warning: chosen meta-period is too high, value resetted to: ', npas,'.'
   ENDIF

   IF ((metaperiod) .LT. 0) THEN
      WRITE(*,*) 'read_input> Error: invalid meta-period value.'
      STOP
   ENDIF
ELSE
   WRITE(*,*) 'read_input> Metadynamics has not been chosen.'
END IF
!--- End of new part (metadynamics)-----------------------------------------------------

!
!! CONTROL ON PROCEDURE CHOICE
!
IF((quenching=='ya').AND.(itremp>npas)) THEN
 WRITE(uniterr,*) "Error: quenching chosen but never runable"
 WRITE(*,*) "read_input> Error: quenching chosen but never runable because itremp>npas"
 STOP
ENDIF
!
IF (coalescence == 'ya') THEN
    IF (natom2 == 0) THEN
        WRITE(uniterr,*) 'Error: no atoms in cluster two'
        WRITE(*,*) 'read_input> Error: no atoms in cluster two'
        STOP
    END IF
    IF (somedist == 0.0d0) THEN
  WRITE(uniterr,*)'Error: possible chance of error'
  WRITE(*,*)'read_input> Error: possible fail of simulation if atoms start together'
 END IF
END IF
!
IF(deposizione=='ya') THEN
 IF(ndepmax==0) THEN
  WRITE(uniterr,*) "Error: deposition chosen but no atoms are deposited"
  WRITE(*,*) "read_input> Error: deposition chosen but no atoms are deposited"
  STOP
 ENDIF
 IF((lcs==2).AND.(prob==0.d0)) THEN
  WRITE(uniterr,*) "Error: bi-deposition chosen with non prob"
  WRITE(*,*) "read_input> Error: bi-deposition chosen with non prob"
  STOP
 ENDIF
 IF(lcs==0) THEN
  WRITE(uniterr,*) "Error: deposition chosen but no species are defined"
  WRITE(*,*) "read_input> Error: deposition chosen but no species are defined"
  STOP
 ENDIF
 IF((natom+ndepmax)>nsiz) THEN
  WRITE(uniterr,*) "Error: too many atoms"
  WRITE(*,*) 'read_input> Error: too many atoms'
  WRITE(*,*) 'read input> Change parameter nsiz in module_parameters.f90'
  STOP
 ENDIF
ENDIF
!
IF((deposizione=='ya').AND.(lcs>1).AND.(elem1==elem2)) elem2=elemd
!
IF ( deposizione == 'no' ) elemd = elem1
!
IF ( sys == 'mon' )   elem2 = elem1
IF(( sys == 'bim' ) .AND. ( elem1 == elem2 )) THEN
 WRITE(*,*) 'read_input> If the system is mono-type, sys == mon'
 WRITE(*,*) 'read_input> ###########################'
 WRITE (uniterr,*) 'If the system is mono-type, sys == mon'
 STOP
ENDIF
!
if(elem2=='NA') THEN
  READ(UNIT=5, nml=sysbim, iostat=ios)
  WRITE(*,*) 'read_input> DYNAMICS WITH ATOMS FALSE Ag ATOMS',sologeom,NA_elem
  IF((caloric == 'ya') .AND. (elem1 /= 'Ag' )) THEN
   WRITE(*,*) 'read_input> Error for False Ag',elem1,elem2
   WRITE(uniterr,*) 'Error for False Ag',elem1,elem2
   STOP
  ENDIF
ENDIF
!
! Initial position
!
WRITE(*,*) 'read_input> Reading initial positions '
OPEN(UNIT=11,file=TRIM(filepos), status = 'old')
  READ(11,*)
  READ(11,*)
IF (coalescence == 'ya') THEN
  OPEN(UNIT=12,file=TRIM(filepos2), status = 'old')
  READ(12,*)
  READ(12,*)
  DO ind = 1, (natom - natom2)
   READ(11,*) &
  &lb_at%spx(ind),lb_at%tau(1,ind),lb_at%tau(2,ind),lb_at%tau(3,ind)!!!,lb_at%icolor
  ENDDO
  DO ind = (natom - natom2 + 1), natom
   READ(12,*) &
  &lb_at%spx(ind),lb_at%tau(1,ind),lb_at%tau(2,ind),lb_at%tau(3,ind)!!!,lb_at%icolor
  END DO
 CLOSE(11)
 CLOSE(12)
ELSE
  DO ind = 1, natom
   READ(11,*) &
  &lb_at%spx(ind),lb_at%tau(1,ind),lb_at%tau(2,ind),lb_at%tau(3,ind)!!!,lb_at%icolor
  ENDDO
  CLOSE(11)
END IF
WRITE(*,*) 'read_input> Finish reading initial positions'
!!!!!!!!!!!!!!!!!!!!
  elem(1:natom)=lb_at%spx(1:natom)
  x0(1:natom)=lb_at%tau(1,1:natom)
  y0(1:natom)=lb_at%tau(2,1:natom)
  z0(1:natom)=lb_at%tau(3,1:natom)

  !Translating clusters to (0, 0, -z) and (0, 0, z)
  IF (coalescence == 'ya') THEN
   comx0 = sum(x0(1:natom-natom2)) / dble(natom-natom2)
   comy0 = sum(y0(1:natom-natom2)) / dble(natom-natom2)
   comz0 = sum(z0(1:natom-natom2)) / dble(natom-natom2)
   com2x0 = sum(x0(natom-natom2+1:natom)) / dble(natom2)
   com2y0 = sum(y0(natom-natom2+1:natom)) / dble(natom2)
   com2z0 = sum(z0(natom-natom2+1:natom)) / dble(natom2)

   x0(1:natom-natom2) = x0(1:natom-natom2) - comx0 ; y0(1:natom-natom2) = y0(1:natom-natom2) - comy0
   z0(1:natom-natom2) = z0(1:natom-natom2) - comz0 + (somedist/2.0d0)

   x0(natom-natom2+1:natom) = x0(natom-natom2+1:natom) - com2x0 ; y0(natom-natom2+1:natom) = y0(natom-natom2+1:natom) - com2y0
   z0(natom-natom2+1:natom) = z0(natom-natom2+1:natom) - com2z0 - (somedist/2.0d0)

   comx0 = sum(x0(1:natom-natom2)) / dble(natom-natom2)
   comy0 = sum(y0(1:natom-natom2)) / dble(natom-natom2)
   comz0 = sum(z0(1:natom-natom2)) / dble(natom-natom2)
   com2x0 = sum(x0(natom-natom2+1:natom)) / dble(natom2)
   com2y0 = sum(y0(natom-natom2+1:natom)) / dble(natom2)
   com2z0 = sum(z0(natom-natom2+1:natom)) / dble(natom2)

   distcom = sqrt((comx0 - com2x0)**2 + (comy0 - com2y0)**2 + (comz0 - com2z0)**2)
   WRITE(*,*)"read_input > COM distance", distcom
  END IF

!check on read position
WRITE(*,*) 'read_input> Validating positions'
DO ind=1,natom
    IF(elem(ind)=='no') THEN
        WRITE(*,*) 'read_input> ERROR: NO ELEM for atoms',ind, elem(ind)
        WRITE(uniterr,*) 'ERROR: NO ELEM for atoms',ind, elem(ind)
        STOP
    ENDIF
ENDDO
!
!--------------------------------------------------------------------
! Check if cluster is in a bad position with respect of MgO substrate
IF (mgo_substrate) THEN
   WRITE (*,*) 'read_input> MgO substrate at z=0 has been chosen'
   IF (metal_on_top) THEN
      WRITE (*,*) 'read_input> Metal-on-top effect is considered'
   ELSE
      WRITE (*,*) 'read_input> Metal-on-top effect is not considered'
   ENDIF
   mgo_mindist = MINVAL(z0(1:natom))
   IF (mgo_mindist.LE.1.d0) THEN
      WRITE(*,*) 'read_input> Error: The cluster is bad positioned w.r.t. the substrate'
      STOP
   ELSEIF (mgo_mindist.LE.2.d0) THEN
      WRITE(*,*) 'read_input> Warning: The cluster is really close to the substrate'
   ENDIF
   WRITE(*,*) 'read_input> Minimum cluster-substrate distance is [A]: ', mgo_mindist
ENDIF

!--------------------------------------------------------------
WRITE(*,*) 'read_input> Checking program'
nproc=0
IF (quenching=='ya')   nproc = nproc+1
IF (canonical=='ya')   nproc = nproc+1
IF (caloric=='ya')     nproc = nproc+1
IF (deposizione=='ya') nproc = nproc+1
IF (metadyn=='ya')     nproc = nproc+1
IF (coalescence == 'ya') nproc = nproc+1
IF (nproc /= 1) THEN
    WRITE(uniterr,*) "Error:you chosen ",nproc,"processes, please choose exacly one"
    WRITE(*,*) "read_input> Error:you chosen ",nproc,"processes, please choose exacly one, make sure you write the process as: &
quenching, microcan, melting, freezing, coalescence, metadynamics, itmd, nvt or growth"
    STOP
ENDIF

WRITE(*,*) 'read_input> Checking program (2)'

!!!!!!!control on coherence of parameter!!!!!!!!
IF(caloric == 'ya') THEN
    IF(tcaloric==0.d0) THEN
        WRITE(uniterr,*) "Error: caloric chosen but no tfin is given"
        STOP
    ENDIF
    IF (deltat==0.d0) THEN
        WRITE(uniterr,*) "Error: caloric chosen but no deltat is given"
        STOP
    ENDIF
ENDIF
!
IF( wires .AND. pbcz == 0.d0) THEN
    WRITE(*,*) 'read_input> Error: for wires give the pbc along Z'
    STOP
ENDIF
!
CALL choice_pot  !!here there is the choice for the potential
!

DO ind=1,natom
itype(ind)=1
  IF(elem(ind)==elem1) ntipo1=ntipo1+1
  IF((elem2/=elem1).AND.(elem(ind)==elem2)) THEN
     ntipo2=ntipo2+1
     itype(ind)=2
   ENDIF

!!!control on material&potential
IF((elem(ind)/=elem1).AND.(elem(ind)/=elem2))THEN
 WRITE(uniterr,*)'bad/unknown material',ind,elem(ind)
  STOP
ENDIF
ENDDO
WRITE(*,*) 'read_input> ntipo1 = ',ntipo1,' ntipo2 = ',ntipo2


IF(((imet1<9).AND.(imet1/=15)).AND.(type_potential/='rgl')) THEN
WRITE(uniterr,*)&
'ERROR:RGL-metal with a wrong potential ',elem1,' ',type_potential
  STOP
ENDIF

IF((imet1==11).AND.(type_potential/='lj1')) THEN
WRITE(uniterr,*)&
&'ERROR:Ar with a wrong potential ',elem1,' ',type_potential
  STOP

ENDIF

IF((imet1==10).AND.((type_potential=='rgl').OR.(type_potential=='lj1'))) THEN
WRITE(uniterr,*)&
&'ERROR:C60 with a wrong potential ',elem1,' ',type_potential
  STOP
ENDIF
!
IF((ntipo1+ntipo2)/=natom) THEN
 WRITE(uniterr,*) 'ERROR on counting atoms of chemical species'
 STOP
ENDIF

natinizio=natom

IF (impl_env .eqv. .TRUE.) THEN
 beta_env_angstrom = alpha_env/rc_env  ! It is used to calculate env coord number
 rc_env = rc_env/arete(1)
 beta_env = alpha_env/rc_env  ! It is used to calculate env coord number
ENDIF

start_from = 1
start_from_nd_proc = 1

restart_mode = .false.
save_progress = .false.
periodic_save = .false.
periodic_save_every = -1
periodic_save_last = 10
READ(UNIT=5, nml=restart, iostat=ios)

restart_from_file = TRIM(restart_from_file)
progress_file_label = TRIM(progress_file_label)

if (save_progress) then
    if (quenching == 'ya' .or. deposizione == 'ya') then
        write(*,*) 'read_input> The end of simulation cannot be saved at the end of the simulation for chosen process'
        stop
    end if
end if
if (periodic_save) then
   if (quenching == 'ya') then
        write(*,*) 'read_input> Periodic save not supported for quenching'
        stop
   else if (metadyn == 'ya' .or. canonical == 'ya') then
        if (periodic_save_every == -1) then
            write(*,*) 'read_input> periodic_save_every not specified'
        end if
   else
        if (periodic_save_every /= 1) then
            write(*,*) 'read_input> WARNING, periodic save for growth or caloric happens every outer cycle'
        end if
        periodic_save_every = 100000000 ! So it never happens in cycle, there must be a better way to specify it.
   end if

end if

WRITE(*,*) 'read_input> Finish reading input'

CLOSE(1)
END SUBROUTINE read_input
