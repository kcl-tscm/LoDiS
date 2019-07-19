SUBROUTINE growth

 USE PARACLUSTER  !uso il modulo di definizione dei parametri
 USE CLUSTER     !uso il modulo dove definisco variabili e parametri cluster
 USE POTENTIAL
 USE ENFORCE
 USE DISTANCE

 IMPLICIT NONE
 INTEGER :: nat3d
 INTEGER :: is,i
 INTEGER :: ir(4)
 REAL :: perc,cb
 REAL :: partialx(nsiz),partialy(nsiz),partialz(nsiz)
 REAL :: pcmx, pcmy,pcmz, vcmx, vcmy, vcmz
 REAL :: vargausorg,szr
 REAL :: gp1,gp2,r,rvett
 REAL :: xad1,yad1,xad,yad,zad,vxad,vyad,vzad
 REAL :: dlaran

 cb = cbol

 ir(1:4)=irand_seed(1:4)

 natom=natom+1
 nat3d=natom

 ! Monometallic growth: lcs = 1
 ! Core-shell groth: lcs = 3 B is deposited over an A_core
 ! 
 ! Bimetallic versus monometallic growth:
 ! with lcs = 2 the program is allowed to deposit
 ! an atom of A or B type, with a given probability == prob
 ! element B is choosen when the random number is smaller than prob
 ! element A is choosen when the random number is bigger/= than prob 
 ! When all the max number of element B are deposited (ntipo2 - initnatom==at_tipo2)
 ! the probability to deposit B atoms is set to zero
 
 IF (lcs==1) THEN
    elem(nat3d)=elem1  
    itype(nat3d)=1
    ntipo1=ntipo1+1
    WRITE(unitd,*) 'deposition of ',elem1
 ELSEIF(lcs==2) THEN
    perc=ABS(2.d0*dlaran(ir)-1.d0)
      IF(perc>=prob) THEN
          elem(nat3d)=elem1
          itype(nat3d)=1
          ntipo1=ntipo1+1
          WRITE(unitd,*) 'deposition of ',elem1, perc 
      ELSEIF(perc<prob) THEN 
         elem(nat3d)=elem2
         itype(nat3d)=2
         ntipo2=ntipo2+1
         IF((ntipo2 -initntipo2)==at_tipo2) prob=0.d0
         WRITE(unitd,*) 'deposition of ',elem2, perc, (ntipo2-initntipo2)
      ENDIF
 ELSEIF (lcs==3) THEN
    elem(nat3d)=elem2
    itype(nat3d)=2
    ntipo2=ntipo2+1
 ENDIF

 pcmx=0.d0
 pcmy=0.d0
 pcmz=0.d0

 vcmx=0.d0
 vcmy=0.d0
 vcmz=0.d0

 DO is=1,nat3d-1
    pcmx=pcmx+x(is)+u(is)
    pcmy=pcmy+y(is)+v(is)
    pcmz=pcmz+z(is)+w(is)

    vcmx=vcmx+vx(is)
    vcmy=vcmy+vy(is)
    vcmz=vcmz+vz(is)
 ENDDO

 pcmx=pcmx/(nat3d-1)
 pcmy=pcmy/(nat3d-1)
 pcmz=pcmz/(nat3d-1)
 vcmx=vcmx/(nat3d-1)
 vcmy=vcmy/(nat3d-1)
 vcmz=vcmz/(nat3d-1)
  
 vargausorg=SQRT(2.d0*t2m(itype(nat3d))*cb*tsorg)

 DO
    gp1=2.d0*dlaran(irand_seed)-1.d0
    gp2=2.d0*dlaran(irand_seed)-1.d0
    r=gp1**2+gp2**2
    IF(r.LT.1) EXIT
 ENDDO
 xad1=rad*gp1
 yad1=rad*gp2
 xad=rad*gp1+pcmx
 yad=rad*gp2+pcmy

 szr=dlaran(irand_seed)
 IF (szr.LE.0.5d0) zad=SQRT(rad**2-xad1**2-yad1**2)+pcmz
 IF (szr.GT.0.5d0) zad=-SQRT(rad**2-xad1**2-yad1**2)+pcmz

 rvett=xad1**2+yad1**2+(zad-pcmz)**2
     
 IF(rad<rshell) THEN
    WRITE(*,*) 'growth> you need rad > rshell!'
    STOP
 ENDIF

 vxad=-vargausorg*xad1/rad
 vyad=-vargausorg*yad1/rad
 vzad=-vargausorg*(zad-pcmz)/rad

 x(nat3d)  = xad
 y(nat3d)  = yad
 z(nat3d)  = zad
 vx(nat3d) = vxad
 vy(nat3d) = vyad
 vz(nat3d) = vzad
 u(nat3d)  = 0.d0
 v(nat3d)  = 0.d0
 w(nat3d)  = 0.d0
 du(nat3d) = vx(nat3d) !==========================!
 dv(nat3d) = vy(nat3d) ! initialize! used in time !
 dw(nat3d) = vz(nat3d) !==========================!
 fx(nat3d) = 0.d0
 fy(nat3d) = 0.d0
 fz(nat3d) = 0.d0

 WRITE(unitd,*) 'atom ',nat3d,' type ',itype(nat3d),' elem ',elem(nat3d)
 WRITE(unitd,'(a10,3F10.4)') 'from ',xad, yad, zad
 WRITE(unitd,'(a10,3F10.4)') 'velocity ',vxad, vyad, vzad
      
 DO i=1,nat3d
    partialx(i)=0.d0
    partialy(i)=0.d0
    partialz(i)=0.d0
 ENDDO

 CALL bigvoi

 !WRITE(900,*) ipas, fx(natom), fy(natom), fz(natom)
 !WRITE(901,*) ipas, u(natom), v(natom), w(natom)
END SUBROUTINE growth
