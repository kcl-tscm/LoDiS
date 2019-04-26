SUBROUTINE scrittod
!
!!!FASE di SCRITTURA
!deposizione=='ya' is my choice
!
USE PARACLUSTER  
USE CLUSTER     
USE ENFORCE
USE POTENTIAL
USE DISTANCE
USE SUBSTRATE
!
IMPLICIT NONE
!
  INTEGER :: i,nt
  INTEGER :: icolor,np
  REAL :: xuar,yvar,zwar
  REAL :: emedia

IF(ipas==1) THEN
   OPEN(unite,file='energy.out',status='unknown',position='append')
   OPEN(unitm,file='movie.xyz', status='unknown', access='append')
   IF(writeheader) THEN
      WRITE(unite,*) '## Time[ps] Ntype1  Ntype2  Etot[eV]  Ekin[eV]  <Edelta>  <Etot>  <T[K]>'
                     !&tempo,ntipo1,ntipo2,etot,ecin,edelta,emedia,tpar
      writeheader = .false.
   ENDIF
ENDIF

     ecin=etot-ener
     tpar=tpar+temp
     edelta=edelta+etot      ! At this point edelta is the sum of etot over a growing number of time steps

!====================================
! Writing each scrivo-step
!

np=MOD(ipas,scrivo)

!Delta and average energy
!
 IF(np==0) THEN
     edelta=edelta/(scrivo)  ! At this point edelta is the average of etot over a growing number of time steps
     emedia=edelta           ! Now emedia is etot averaged over the number of scrivo steps
     IF(sys == 'mon') THEN
      edelta=(edelta+natom*ecoh(1))/((natom)**(2.d0/3.d0))  ! Now edelta is really edelta averaged over scrivo
     ELSE
      !edelta=(etot+ntipo1*ecoh(1)+ntipo2*ecoh(2))/((natom)**(2.d0/3.d0))   ! HOW IT WAS ??? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      edelta=(edelta+ntipo1*ecoh(1)+ntipo2*ecoh(2))/((natom)**(2.d0/3.d0))  ! Now edelta is really edelta averaged over scrivo
     ENDIF
!
! Average Temp
!
     tpar=tpar/scrivo
WRITE(unite,'(1f11.4,2i5,1x,2f14.5,1x,3f14.5)')&
&tempo,ntipo1,ntipo2,etot,ecin,edelta,emedia,tpar
WRITE(*,'(1f11.4,2i5,1x,2f14.5,1x,3f14.5)')&
&tempo,ntipo1,ntipo2,etot,ecin,edelta,emedia,tpar

WRITE(unitm,*) natom
WRITE(unitm,'(i8, 2(a3,i4), 1x, 1f14.5, 2f7.3)') &
   & ipas, elem1, ntipo1, elem2, ntipo2, etot, temp, edelta

     tpar=0.d0
     edelta=0.d0
!
     nfile=ipas/scrivo  !se npasd=10^6  1<nfile<100
                        !se npasd=30^6  1<nfile<300

      !CALL subnomi
      !OPEN(10,file=filename(nd_proc*nfile),status='unknown')
      !REWIND(10)
      nt = natom
!
      !WRITE(10,*) nt
      !WRITE(10,'(a2,1x,a2,1f14.5,i6)') elem1,elem2,etot,nd_proc
       DO i=1,nt          
         icolor=itype(i)                    
         xuar=(x(i)+u(i))*aretebim*fattor
         yvar=(y(i)+v(i))*aretebim+fattor
         zwar=(z(i)+w(i))*aretebim*fattor
       !  WRITE(10,'(a3,1x,3f16.5,i4)') elem(i),xuar,yvar,zwar,icolor
         WRITE(unitm,'(a3,1x,3f16.5,i4)') elem(i),xuar,yvar,zwar,icolor
        ENDDO
        !CLOSE(10)
ENDIF !!su np

IF(ipas==npas) THEN
   CLOSE(unite)
   close(unitm)
ENDIF

END SUBROUTINE scrittod
