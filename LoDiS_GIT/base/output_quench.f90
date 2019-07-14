SUBROUTINE scrittoq
!!!WRITING output after a quenching procedure
!
USE PARACLUSTER 
USE CLUSTER     
USE ENFORCE
USE POTENTIAL
USE DISTANCE
USE SUBSTRATE
USE ENVIRONMENT
!
IMPLICIT NONE

  INTEGER :: i,j,nt
  INTEGER :: icolor,np
  REAL :: xuar,yvar,zwar
  REAL :: emedia, save_delta
!
IF(ipas .EQ. 1) THEN
   tpar=0.d0
   IF (mgo_substrate) THEN 
      WRITE(*,*) '    Time[ps]   Etot[eV]   Epot[eV]   Ekin[eV]   <Etot>   <Edelta>   <Temp[K]>    Esub[eV]'
      WRITE(unite,*) '## Time[ps]   Etot[eV]   Epot[eV]   Ekin[eV]   <Etot>   <Edelta>   <Temp[K]>   Esub[eV]'
   ELSE
      WRITE(*,*) '    Time[ps]   Etot[eV]   Epot[eV]   Ekin[eV]   <Etot>   <Edelta>   <Temp[K]>'
      WRITE(unite,*) '## Time[ps]   Etot[eV]   Epot[eV]   Ekin[eV]   <Etot>   <Edelta>   <Temp[K]>'
   ENDIF
ENDIF
!
     tpar=tpar+temp
     edelta=edelta+etot
! 
np=MOD(ipas,scrivo)
!
!Calculating EDelta (as defined in F. Baletto, PhD thesis)
!
IF(np==0) THEN
  edelta=edelta/(scrivo)
  emedia=edelta
  IF (sys == 'mon' ) THEN
    edelta=(edelta+natom*ecoh(1))/((natom)**(2.d0/3.d0))
  ELSE
    ! edelta=(etot+ntipo1*ecoh(1)+ntipo2*ecoh(2))/((natom)**(2.d0/3.d0))  IT WAS IN THIS WAY! Like in growth
    edelta=(etot+ntipo1*ecoh(1)+ntipo2*ecoh(2))/((natom)**(2.d0/3.d0))
  ENDIF
  tpar=tpar/scrivo

!
   IF( surface .or. bulk ) THEN
      WRITE(unite,'(1f11.4,2x,4f15.8,1x,1f11.5,1x,1f9.4)')&
      & tempo, etot, ener, ecin, emedia, edelta, tpar
      WRITE(*,'(1f11.4,2x,4f15.8,1x,1f11.5,1x,1f9.4)') &
      & tempo, etot, ener, ecin, emedia, edelta, temp  ! Why not tpar?
   ELSE IF (mgo_substrate) THEN
      WRITE(unite,'(1f11.4,1x,4f16.7,2(1x,1f9.4),1f16.7)')&
      &tempo, etot, ener, ecin, emedia, edelta, tpar, ener_sub
      WRITE(*,'(1f11.4,1x,4f16.7,2(1x,1f9.4),1f16.7)')&
      &tempo, etot, ener, ecin, emedia, edelta, tpar, ener_sub
   ELSE   
      WRITE(unite,'(1f11.4,1x,4f16.7,2(1x,1f9.4))')&
      &tempo, etot, ener, ecin, emedia, edelta, tpar
      WRITE(*,'(1f11.4,1x,4f16.7,2(1x,1f9.4))')&
      &tempo, etot, ener, ecin, emedia, edelta, tpar
   ENDIF


!! For the stress
!  IF (interewr) THEN
!    filename2=filename(nd_proc*nfile)
!    write(filename2(1:1),'(A1)') 'i'
!    OPEN(unit=51,file=filename2) 
!    call inter
!  ENDIF
!
!! Intershell stress
!!  IF ( tintershell ) THEN
!!   write(*,*) 'calling inter'
!!   OPEN(unit=51,file='inter.out') 
!!   call inter 
!!  ENDIF
!
   IF(tpar<tmin) THEN
       WRITE(*,*)&
      &tempo,'QUENCHING END FOR TPAR<TMIN',etot,temp
      STOP
    ENDIF
!
     save_delta= edelta
     tpar=0.d0
     edelta=0.d0
     nfile=ipas/scrivo  !se npasd=10^6  1<nfile<100
                        !se npasd=30^6  1<nfile<300
!
! WRITE (123,*) nfile,scrivo,ipas,npas
! IF ((nfile .ge. 1).or.(nd_proc-1 .ge. 1)) THEN
    !CALL subnomi
       !OPEN(10,file=filename(nd_proc*nfile),status='unknown')
       !REWIND(10)
!
   nt = natom
   IF ( wires ) THEN 
         nt = 2*natom

!         WRITE(10,*) nt
!         WRITE(10,'(a2,1x,a2,1f15.6,1f9.4)') elem1,elem2,etot,save_delta
          DO i=1,natom          
            icolor=itype(i)                    
            xuar=(x(i)+u(i))*aretebim*fattor
            yvar=(y(i)+v(i))*aretebim*fattor
            zwar=(z(i)+w(i))*aretebim*fattor
!           WRITE(10,'(a3,1x,4f16.5,i4)') elem(i),xuar,yvar,zwar,icolor, &
!                                       & potener(i)+ener_env_atom(i) 
!            WRITE(10,'(a3,1x,3f16.5,i4)') elem(i),xuar,yvar,zwar+pbcz*aretebim,icolor+1
           ENDDO
!           CLOSE(10)
!
   ENDIF !!su np
! ENDIF

ENDIF
END SUBROUTINE scrittoq
