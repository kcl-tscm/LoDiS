SUBROUTINE INTER
!
USE PARACLUSTER
USE CLUSTER
USE POTENTIAL
USE ENFORCE
!
IMPLICIT NONE
REAL :: pcmx,pcmy,pcmz !position of mass center
REAL :: px,py,pz !position of mass center
REAL :: Rshell1,RMAX,intr !approximated radius of each shell and 
                !largest distance beetwen atom and CM
REAL, DIMENSION(:),ALLOCATABLE ::  Rdist      !distance beetwen each atom and CM
REAL :: eps     !atoms in one shell should be in distances 
             !beetwen R1*(1+-eps)
INTEGER :: i,l,checkat !checkat is used only to check if number of atoms in shells is equal to natom
INTEGER, DIMENSION(:),ALLOCATABLE :: Nshell !number of atoms in each shell
REAL, DIMENSION(:),ALLOCATABLE ::  averr  
!Average distance between atoms in the same shell respect to the com
REAL, DIMENSION(:),ALLOCATABLE ::  intrel !Intershell relaxation
REAL, DIMENSION(:),ALLOCATABLE :: lee !local excess energy
REAL :: che
REAL, DIMENSION(:),ALLOCATABLE :: Pl !local excess energy


!! calculate the intershell relaxation  

ALLOCATE(Rdist(natom))    
     pcmx  =0.d0
     pcmy  =0.d0
     pcmz  =0.d0

       DO i=1,natom
          pcmx=pcmx+x(i)+u(i)
          pcmy=pcmy+y(i)+v(i)
          pcmz=pcmz+z(i)+w(i)
       ENDDO
         pcmx=pcmx/(natom)
         pcmy=pcmy/(natom)
         pcmz=pcmz/(natom)

        write(51,*) 'the cdm is', pcmx,pcmy,pcmz

RMAX=0.
Rshell1=1.d0
eps=2.d-2
DO i=1,natom
  
  px = (x(i)+u(i)-pcmx)*(x(i)+u(i)-pcmx)*aretebim*aretebim
  py = (y(i)+v(i)-pcmy)*(y(i)+v(i)-pcmy)*aretebim*aretebim
  pz = (z(i)+w(i)-pcmz)*(z(i)+w(i)-pcmz)*aretebim*aretebim



  Rdist(i)=sqrt(px+py+pz)

  if (Rdist(i)>RMAX) RMAX=Rdist(i)

ENDDO
   write(51,*) 'rmax values:', RMAX
!!RMAX= MAXVAL(R(:))

!nshells=int((RMAX/Rshell1)) + 1


ALLOCATE(Nshell(nshells+1))
ALLOCATE(Pl(nshells))
ALLOCATE(averr(0:nshells+1))
ALLOCATE(lee(nshells+1))
!ALLOCATE(intrel(nshells)) 

averr(:)=0.d0
lee(:)=0.d0
Nshell(:)=0
RMAX = RMAX*(1.d0+eps)
che=0
  write(51,*) '#collect in shell'
DO i=1,natom
  l=int((Rdist(i)*nshells/RMAX))+1
  Nshell(l)= Nshell(l)+1
  averr(l)=averr(l)+Rdist(i)
  lee(l)=lee(l)+potener(i)
  Pl(l)=Pl(l)+press(i)
  che=che+potener(i)
  write(51,*) i,l,Rdist(i)
ENDDO

write(51,*) "Energy sum:",che,"=",etot
write(51,*) "Ecoh", ecoh(1)


write(51,*) '#shell, #atoms #inter av(l) l-(l-1)  loc_exc_energy pressure'

checkat=0
do l = 1,nshells+1

if(Nshell(l) > 0) THEN
averr(l) =averr(l)/Nshell(l)
lee(l)=ecoh(1)+lee(l)/Nshell(l)
Pl(l) =Pl(l)/Nshell(l)
ENDIF

if(Nshell(l) == 0)  then
averr(l) = 0.d0
lee(l)= 0.d0
Pl(l)=0
endif

write(51,*) l, Nshell(l), averr(l), (averr(l) - averr(l-1)), lee(l), Pl(l)

checkat=checkat + Nshell(l)
enddo

if (checkat.ne.natom ) THEN
WRITE (*,*) "ERROR: NOT ALL ATOMS ARE IN PROPER SHELLS"
WRITE (51,*) "ERROR: NOT ALL ATOMS ARE IN PROPER SHELLS"
endif
DEALLOCATE(Rdist,Nshell,averr)
CLOSE(51)

END SUBROUTINE INTER
!   averr(1:nshells)=averr(1:nshells)/Nshell(:)

!averr(:)=0.d0
!Nshell(:)=0
!RMAX = RMAX*(1.d0+eps)

!do i= 1,natom
!do l = 1,nshells
!if((Rdist(i).ge.(l-1)*RMAX/nshells).AND.(Rdist(i).lt.l*RMAX/nshells)) THEN
! Nshell(l) = Nshell(l) + 1
! averr(l) = averr(l) + Rdist(i)
!ENDIF
!enddo
!Enddo
