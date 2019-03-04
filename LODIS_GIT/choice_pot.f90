SUBROUTINE choice_pot
  USE PARACLUSTER  !uso il modulo di definizione dei parametri
  USE CLUSTER     !uso il modulo dove definisco variabili e parametri cluster
  USE POTENTIAL
  USE ENFORCE

  IMPLICIT NONE
  INTEGER :: NA_imet
  REAL*8 :: frac
!  REAL :: rmet1,rmet2
  LOGICAL :: alcalini,legaNA
!
!Rosato-Guillope-Legrand (Gupta like) from RGL, Philos. Mag. A 1989
!
  IF(type_potential=='rgl') THEN
     !
     alcalini=.FALSE.
     legaNA = .FALSE.
     !
     IF(elem1 == 'Ni') imet1=1
     IF(elem1.EQ.'Cu') imet1=2
     IF(elem1.EQ.'Pd') imet1=3
     IF(elem1.EQ.'Ag') imet1=4
     IF(elem1.EQ.'Pt') imet1=5
     IF(elem1.EQ.'Au') imet1=6
     IF(elem1.EQ.'Co') imet1=7
     IF(elem1.EQ.'Li') THEN 
        imet1=8
        alcalini=.TRUE.
     ENDIF
     IF(elem1.EQ.'Na') THEN
        imet1=9
        alcalini=.TRUE.
     ENDIF
     !
     IF( sys == 'mon' ) THEN
      imet2=imet1
      itype(1:natom)=1
     !
     ELSE
        IF(elem2.EQ.'Ni') imet2=1
        IF(elem2.EQ.'Cu') imet2=2
        IF(elem2.EQ.'Pd') imet2=3
        IF(elem2.EQ.'Ag') imet2=4
        IF(elem2.EQ.'Pt') imet2=5
        IF(elem2.EQ.'Au') imet2=6
        IF(elem2.EQ.'Co') imet2=7
        IF(elem2.EQ.'Li') THEN
           imet2=8
           alcalini=.TRUE.
        ENDIF
        IF(elem2.EQ.'Na') THEN
           imet2=9
           alcalini=.TRUE.
        ENDIF
        !
        IF(elem2=='NA') THEN 
           legaNA=.TRUE.
           imet2=10
           IF(NA_elem == 'Ni') NA_imet = 1  
           IF(NA_elem == 'Cu') NA_imet = 2  
           IF(NA_elem == 'Pd') NA_imet = 3
           WRITE(*,*) 'WHAT NA==',NA_elem,NA_imet,imet2
        ENDIF
      ENDIF  !!on sys
!
!Reading the potential file, filepot
!
        CALL bimet_rgl
!
        dmas(imet1) = mass(1)
        dmas(imet2) = mass(2)
!
!Lennard-Jones for Argon
!
     ELSEIF(type_potential=='lj1') THEN
        Write(*,*) ' LJ is chosen'
        IF(elem1 == 'Ar' ) imet1 = 11
        leg(11)='Ar'
        ecohv(11)=0.08d0
        dmas(11)=39.948d0
        R0=3.82d0     !distance NN in A
        ratv(11)=R0/2.d0
        U0=0.01034d0 ! costant for Ar_LJ (eV)
        arete(1)=ratv(imet1)*rac8
!
!Girifalco or Pacheco for C60
!
     ELSE
        IF(elem1=='C6') imet1=12
        leg(12)='C6'
        ecohv(12)=1.743d0
        dmas(12)=720.0d0
 
       IF(type_potential=='gir') THEN
           Dnn=10.05d0        !distance NN in A
           ratv(12)=Dnn/2.d0
           Nball=60     !atom number in a ball
           dball=7.1d0  !diameter of a ball for C60 (A)
           ALJ=20.00d0  !costant as LJ  (eV*A**6)
           BLJ=34.856d3 !costant as LJ  (eV*A**12)
        ELSEIF(type_potential=='par') THEN
           dmu=10.05d0        !distance NN in A
           ratv(12)=dmu/2.d0
           epsilon=1.04d0 !!delta=1.04d0
           dM0=0.3d0
           tau=9.75d0
           DIST0=10.3d0
           C6=75600.d0
           C8=9122400.d0
           C10=2.09d8
           C12=7.78d10
        ENDIF
        arete(1)=ratv(imet1)*rac8
     ENDIF
     !
     frac=evsujoule/(angsum*angsum)       
     frac=frac*(tstep*tstep)/(2.d0*uasukg)
     !
! Calculate the tstep**2/2m as in Velocity Verlet algorithm
     !
     t2m(1)=frac/(dmas(imet1)*arete(1)*arete(1))
     if(imet2 /= imet1) t2m(2)=frac/(dmas(imet2)*arete(2)*arete(2))
     WRITE(*,*) 'choice_pot> values:'
     WRITE(*,*) 'frac = ', frac, ' arete = ', arete(1)
     WRITE(*,*) 'choice_pot> mass is==',dmas(imet1), 'for elem1= ' !,leg(imet1)
     !
   END SUBROUTINE choice_pot

