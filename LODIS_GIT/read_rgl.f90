SUBROUTINE bimet_rgl
  USE PARACLUSTER  !uso il modulo di definizione dei parametri
  USE CLUSTER     !uso il modulo dove definisco variabili e parametri cluster
  USE POTENTIAL
  USE ENFORCE

  IMPLICIT NONE

!Local variables
  Integer :: ind,i,j,k, it ,itt
  REAL*8 :: dik0,rmet1,rmet2
  REAL*8 :: u2,v2,zr,za,zap
  REAL*8 :: ar,br,cr,ab,bb,cb
  REAL*8 :: pv(n_elposs),qv(n_elposs)
  REAL*8 :: av(3),qsiv(3)
  LOGICAL :: no_alc,legaNA
  Character(2) :: metal(2)
!
  lsubstrate = .FALSE.
  no_alc=.TRUE.
  IF((imet1==7).OR.(imet1==8)) no_alc=.FALSE.
  legaNA= .FALSE.
! 
  write (*,*) 'read_rgl> START READING FILE ', filepot
  open(19,file=TRIM(filepot),status='old')
  read(19,*)
  read(19,*)
  read(19,*) metal(1), metal(2)
            !check consistency with metals declared in seed.in
            if ((metal(1).ne. elem1) .and. (metal(2) .ne. elem2))then
               write(*,*)'####################################################'
               write(*,*)'ERROR:'
               write(*,*)'Metals declared in inputfile elem1 and elem2'
               write(*,*)'MUST BE THE SAME AND IN THE SAME ORDER' 
               write(*,*)'as in the potential file.'
               write(*,*)'####################################################'
               write(uniterr,*)'##############################################'
               write(uniterr,*)'ERROR:'
               write(uniterr,*)'Metals declared in inputfile'
               write(uniterr,*)'MUST BE THE SAME AND IN THE SAME ORDER' 
               write(uniterr,*)'as in the potential file.'
               write(uniterr,*)'##############################################'
               !!close(uniterr)
               stop
            endif
            !end check consistency
  read(19,*)
  read(19,*)
  read(19,*)p(1),p(2),p(3)
  read(19,*)q(1),q(2),q(3)
  read(19,*)a(1),a(2),a(3)
  read(19,*)qsi(1),qsi(2),qsi(3)
  read(19,*)
  read(19,*)
  read(19,*)ecoh(1),ecoh(2)
  read(19,*)rmet1,rmet2
  read(19,*)mass(1),mass(2)
  read(19,*)
  read(19,*)
  read(19,*)cutoff_start,cutoff_end

!
  write (*,*) 'read_rgl> END READING FILE ', filepot
  close(19)
  !
  ratv(imet1) = rmet1
  ratv(imet2) = rmet2
  !
  !Unit conversions from AA to arete-unit
  !
  arete(1)=ratv(imet1)*dsqrt(8.d0)
  IF( sys == 'bim' ) arete(2)=ratv(imet2)*dsqrt(8.d0)
  IF( sys == 'bim' ) arete(3)=(arete(1)+arete(2))/2.0d0
  !
  !nn is the nearest neighbours distances in AA
  !
  nn(1)=arete(1)/dsqrt(2.d0)
  IF( sys == 'bim' ) nn(2)=arete(2)/dsqrt(2.d0)
  IF( sys == 'bim' ) nn(3)=arete(3)/dsqrt(2.d0)

  !dist are the nearest neighbours distances in arete(1) units
  dist(1)=1.d0/dsqrt(2.d0)
  IF( sys == 'bim' ) dist(2)=nn(2)/arete(1)
  IF( sys == 'bim' ) dist(3)=nn(3)/arete(1)

  !! IMPROVEMENT cutoff is to be logical)
  !cutoff_end and cutoff_start are converted into arete(1) units
  if( lcutoff )then
     cutoff_start=cutoff_start/arete(1)
     cutoff_end=cutoff_end/arete(1)
  else
     write(*,*) 'read_rgl> default value for cutoff is .false. Poor physical meaning'
     cutoff_start=2000.
     cutoff_end=2000.
  endif
!
! the distance to NN is stored in cutz and then used in voisin and bivoi subroutines
!
!
!IF SUBSTRATE IS ACTIVE the interaction with it is ruled by ddd
!


!  IF(legaNA) THEN
!     pv(9)=pv(4)    
!     qv(9)=qv(4)    
!     arete(1)=ratv(imet1)*rac8
!     arete(2)=ratv(imet2)*rac8
!     rmet1=1.d0/rac8
!     rmet2=rmet1*ratv(imet2)/ratv(imet1)
!     dist(1)=rmet1+rmet1
!     dist(2)=rmet1+rmet1
!     dist(3)=rmet1+rmet1 !!solo Ag      
!  ENDIF
!
!    IF( sys == 'mon' ) THEN
!     dc1 = dist(1)
!     dc11= dist(1)
!    ELSEIF( sys == 'bim' ) THEN 
!     dc1=MAX(dist(1),dist(2))
!     dc11=MIN(dist(1),dist(2))
!    ENDIF
!     dc2=rac2*dc1 
!     dc3=rac3*dc11
!     dc4=2.d0*dc3
!     dc5=rac5/rac2!!dc3+1.d0 !!rac5/rac2
!     dc33=dc3*dc3
!     dc55=dc5*dc5
!
!Cutoff parameters a5,a4,a3,x5,x4,x3
!
  do i=1,3 
        dik0=dist(i)
        ar=-a(i)*dexp(-p(i)*(cutoff_start/dik0-1.d0))/(cutoff_end-cutoff_start)**3 
        br=-(p(i)/dik0)*a(i)*dexp(-p(i)*(cutoff_start/dik0-1.d0))/(cutoff_end-cutoff_start)**2
        cr=-((p(i)/dik0)**2) &
             *a(i)*dexp(-p(i)*(cutoff_start/dik0-1.d0))/(cutoff_end-cutoff_start)
        ab=-qsi(i)*dexp(-q(i)*(cutoff_start/dik0-1.d0))/(cutoff_end-cutoff_start)**3
        bb=-(q(i)/dik0)*qsi(i)*dexp(-q(i)*(cutoff_start/dik0-1.d0))/(cutoff_end-cutoff_start)**2 
        cb=-((q(i)/dik0)**2) &
             *qsi(i)*dexp(-q(i)*(cutoff_start/dik0-1.d0))/(cutoff_end-cutoff_start)
        x5(i)=(12.d0*ab-6.d0*bb+cb)/(2.d0*(cutoff_end-cutoff_start)**2)
        x4(i)=(15.d0*ab-7.d0*bb+cb)/(cutoff_end-cutoff_start)
        x3(i)=(20.d0*ab-8.d0*bb+cb)/2.d0
        a5(i)=(12.d0*ar-6.d0*br+cr)/(2.d0*(cutoff_end-cutoff_start)**2)
        a4(i)=(15.d0*ar-7.d0*br+cr)/(cutoff_end-cutoff_start)
        a3(i)=(20.d0*ar-8.d0*br+cr)/2.d0
     enddo
!
!
  WRITE(*,*) 'read_rgl> PARAMETERS SELECTION'
  WRITE(*,*) 'FOR IMET1 :: ',elem1,imet1 !,leg(imet1)
  WRITE(*,*) 'ECOH :: ',ecoh(1),'A :: ',a(1),'QSI  :: ',qsi(1)
  WRITE(*,*) 'p :: ',p(1),'q ::',q(1)
  WRITE(*,*) 'ATOMIC RADIUS  ::',ratv(imet1)
  IF(sys == 'bim') THEN
   WRITE(*,*) 'FOR IMET2 :: ',elem2,imet2 !,leg(imet2)
   WRITE(*,*) 'ECOH :: ',ecoh(2),'A :: ',a(2),'QSI ::',qsi(2)
   WRITE(*,*) 'p :: ',p(2),'q :: ',q(2)
   WRITE(*,*) 'ATOMIC RADIUS',ratv(imet2)
  ENDIF
END SUBROUTINE BIMET_RGL

  !raccordo potenziale
  !IF(no_alc)THEN
  !   DO it=1,itt
  !      u2=EXP(-p(it)*(rac2-1.d0))
  !      v2=EXP(-2.d0*q(it)*(rac2-1.d0))
  !      zr=12.d0+6.d0*u2
  !      zrp=(12.d0+6.d0*u2*rac2)*p(it)
  !      za=SQRT(12.d0+6.d0*v2)
  !      zap=(12.d0+6.d0*v2*rac2)*q(it)/za
  !      a(it)=ecoh(it)*zap/(za*zrp-zr*zap)
  !      qsi(it)=ecoh(it)*zrp/(za*zrp-zr*zap)
  !   ENDDO
  !   IF((itt==2).AND.((imet1==6).AND.(imet2==4)))THEN
  !      a(1)=0.2096d0
  !      qsi(1)=1.8153d0
  !      a(2)=0.1031d0
  !      qsi(2)=1.1895d0
  !   ELSEIF((itt==2).AND.((imet2==6).AND.(imet1==4))) THEN
  !      a(1)=0.1031d0
  !      qsi(1)=1.1895d0
  !      a(2)=0.2096d0
  !      qsi(2)=1.8153d0
  !   ENDIF
!
  !ELSEIF((itt==1).AND.(imet1==7)) THEN
  !   a(1)=av(1)
  !   qsi(1)=qsiv(1)
  !ELSEIF ((itt==1).AND.(imet1==8))  THEN
  !   a(1)=av(2)
  !   qsi(1)=qsiv(2)
  !ELSEIF((itt==1).AND.(imet1==14)) THEN
  !  a(1)=av(3)
  !  qsi(1)=qsiv(3)
  !  p(1)=pv(10)
  !  q(1)=qv(10)
  !ENDIF

  !
  !c nanoalloy-system Cu3Au; Cu-Ag; Pd-Ag; Ni-Ag; Pt-Co
  !
  !IF(((elem1.EQ.'Cu').AND.(elem2.EQ.'Ag')).OR.&
  !     &((elem1.EQ.'Ag').AND.(elem2.EQ.'Cu'))) THEN
  !   p(3)=10.700d0
  !   q(3)=2.805d0
  !   a(3)=0.0977483d0
  !   qsi(3)=1.22746d0
  !ELSE IF(((elem1.EQ.'Pd').AND.(elem2.EQ.'Ag')).OR.&
  !     &((elem1.EQ.'Ag').AND.(elem2.EQ.'Pd'))) THEN
  !   p(3)=10.895d0
  !   q(3)= 3.492d0
  !   a(3)=0.16069328d0
  !   qsi(3)=1.559693d0
  !ELSE IF(((elem1=='Ni').AND.(elem2=='Ag')).OR.&
 !      &((elem1=='Ag').AND.(elem2=='Ni'))) THEN
 !    p(3)=11.095d0
 !    q(3)= 2.725d0
 !    a(3)=0.096d0
 !    qsi(3)=1.340d0
 ! ELSE IF(((elem1.EQ.'Cu').AND.(elem2.EQ.'Au')).OR.&
 !      &((elem1.EQ.'Au').AND.(elem2.EQ.'Cu'))) THEN
 !    p(3)=11.590d0
 !    q(3)= 3.280d0
 !    a(3)=0.1085d0
 !    qsi(3)=1.44d0
 ! ELSE IF(((elem1.EQ.'Ag').AND.(elem2.EQ.'Au')).OR.&
 !      &  ((elem1.EQ.'Au').AND.(elem2.EQ.'Ag'))) THEN!AuAg
 !    p(3)= 10.494d0
 !    q(3)=  3.607d0
 !    a(3)=  0.1488d0
 !    qsi(3)=1.4874d0
 ! ELSE IF(((elem1.EQ.'Co').AND.(elem2.EQ.'Pt')).OR.&
 !      &  ((elem1.EQ.'Pt').AND.(elem2.EQ.'Co'))) THEN
 !     p(3)=  9.97d0
 !    q(3)=  3.32d0
 !    a(3)=  0.2447d0
 !    qsi(3)= 2.3848d0
 ! ELSE  
 !    p(3)=p(1)
 !    q(3)=q(1) 
 !    a(3)=a(1) 
 !    qsi(3)=qsi(1) 
 ! ENDIF
!
!  DO i=1,3 
!     dik0=dist(i) 
!     ar=-a(i)*EXP(-p(i)*(dc2/dik0-1.))/((dc3-dc2)**3) 
!     br=-(p(i)/dik0)*a(i)*EXP(-p(i)*(dc2/dik0-1.d0))/(dc3-dc2)**2 
!     cr=-((p(i)/dik0)**2)*a(i)*EXP(-p(i)*(dc2/dik0-1.d0))/(dc3-dc2) 
!     ab=-qsi(i)*EXP(-q(i)*(dc2/dik0-1.d0))/(dc3-dc2)**3 
!     bb=-(q(i)/dik0)*qsi(i)*EXP(-q(i)*(dc2/dik0-1.d0))/(dc3-dc2)**2 
!     cb=-((q(i)/dik0)**2)*qsi(i)*EXP(-q(i)*(dc2/dik0-1.d0))/(dc3-dc2) 
!     x5(i)=(12.d0*ab-6.d0*bb+cb)/(2.d0*(dc3-dc2)**2)
!     x4(i)=(15.d0*ab-7.d0*bb+cb)/(dc3-dc2) 
!     x3(i)=(20.d0*ab-8.d0*bb+cb)/2.d0
!     a5(i)=(12.d0*ar-6.d0*br+cr)/(2.d0*(dc3-dc2)**2)
!     a4(i)=(15.d0*ar-7.d0*br+cr)/(dc3-dc2) 
!     a3(i)=(20.d0*ar-8.d0*br+cr)/2.d0 
!  ENDDO
!
