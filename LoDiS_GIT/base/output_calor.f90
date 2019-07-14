SUBROUTINE scrittoc
USE PARACLUSTER  
USE CLUSTER     
USE POTENTIAL
USE ENFORCE
USE DISTANCE
USE SUBSTRATE
!
IMPLICIT NONE

  INTEGER :: i,nt, unitp
  INTEGER :: icolor,np
  REAL :: xuar,yvar,zwar,xcdm,ycdm,zcdm
  REAL :: emedia


     unitp = 50
!
!!!Writing the output for a caloric curve
!!! energy.out contains the energy output
!!! out~.xyz the positions 
!! if(caloric=='ya')
!
 IF(ipas==1) THEN
    OPEN(unite,file='energy.out',status='unknown',position='append')
    OPEN(unitm,file='movie.xyz', status='unknown', access='append')
    IF(writeheader) THEN
       IF(mgo_substrate) THEN
          WRITE(unite,*) '## Time[ps]  Epot[eV]  Etot[eV]  Ekin[eV]  <Edelta>  <Etot>  <T[K]>   Esub[eV]'
       ELSE
          WRITE(unite,*) '## Time[ps]  Epot[eV]  Etot[eV]  Ekin[eV]  <Edelta>  <Etot>  <T[K]>'
       ENDIF
       writeheader = .false.
    ENDIF
 ENDIF
   
!
!Average quantities
!
     tpar  = tpar + temp
     edelta= edelta + etot
  !  emedia= emedia + etot   ! emedia is not initialized and this value is not used...


!====================================
! Writing each scrivo-step
!   
     np=MOD(ipas,scrivo)


     IF(np==0) THEN 
!
!
! Average Energy
!
     emedia=edelta/scrivo ! At this point edelta is the sum of etot over a growing number of time steps
                          ! Now emedia is etot averaged over the number of scrivo steps
!EDelta calculation
!
     edelta=edelta/(scrivo) ! Now edelta becomes etot averaged over the number of scrivo steps
     IF(sys == 'mon') THEN 
      edelta=(edelta+natom*ecoh(1))/(natom**(2.d0/3.d0))  ! Now edelta is really edelta averaged over scrivo
     ELSE
      edelta=(edelta+ntipo1*ecoh(1)+ntipo2*ecoh(2))/((natom)**(2.d0/3.d0)) ! Now edelta is really edelta averaged over scrivo
     ENDIF
!
!Average Temperature calculation
!
    tpar=tpar/(scrivo)
!  
     IF(mgo_substrate) THEN
        WRITE(unite,'(1f15.4,1x,3f14.5,1x,4f14.5)')&
        & tempo,ener ,etot ,ecin ,edelta,emedia ,tpar, ener_sub
        ! On the screen
        WRITE(*,'(1f15.4,1x,3f14.5,1x,4f14.5)') &
        & tempo,ener ,etot ,ecin ,edelta,emedia ,tpar, ener_sub
     ELSE
        WRITE(unite,'(1f15.4,1x,3f14.5,1x,3f14.5)')&
        & tempo,ener ,etot ,ecin ,edelta,emedia ,tpar
        ! On the screen
        WRITE(*,'(1f15.4,1x,3f14.5,1x,3f14.5)') &
        & tempo,ener ,etot ,ecin ,edelta,emedia ,tpar
     ENDIF
     
!
     emedia=0.d0
     edelta=0.d0
     tpar=0.d0
!
!    !!nfile=ipas/scrivo
!    !!CALL subnomi          
     !!OPEN(unitp,file=filename(nd_proc*nfile),status='unknown',access='append')
     !!OPEN(unitp,file='pos.out',status='unknown',access='append')
!
      nt = natom
!      if ( wires ) nt= 2* natom 
!
       xcdm = 0.d0  
       ycdm = 0.d0
       zcdm = 0.d0
!
      DO i =1,nt
          xcdm = xcdm + (x(i)+u(i))
          ycdm = ycdm + (y(i)+v(i))
          zcdm = zcdm + (z(i)+w(i))
      ENDDO
          xcdm = xcdm*aretebim*fattor / nt
          ycdm = ycdm*aretebim*fattor / nt
          zcdm = zcdm*aretebim*fattor / nt
!
!      open(unitp, file = filename(nd_proc*nfile), status = 'unknown')          
!      rewind(unitp)
!!      WRITE(unitp,*) nt
!      WRITE(unitp,'(a2,1x,a2,1x,1f12.6)') elem1,elem2,etot
!      close(unitp)
!
      WRITE(unitm,*) nt
      WRITE(unitm,'(i8, 2(a3,i4), 1x, 1f14.5, 1f7.3, a4, 3(f9.4,1x))') &
      & ipas, elem1, ntipo1, elem2, ntipo2, etot, temp, 'COM', xcdm, ycdm, zcdm
!
       DO i=1,nt          
         icolor=itype(i)                    
         xuar=(x(i)+u(i))*aretebim*fattor
         yvar=(y(i)+v(i))*aretebim*fattor
         zwar=(z(i)+w(i))*aretebim*fattor
!
!        WRITE(unitp,'(a3,1x,3f16.5,i4)') elem(i),xuar,yvar,zwar,icolor
         IF (mgo_substrate) THEN
            WRITE(unitm,'(a3,1x,3f16.5,i4)') elem(i),xuar,yvar,zwar,icolor
         ELSE
            WRITE(unitm,'(a3,1x,3f16.5,i4)') elem(i),xuar-xcdm,yvar-ycdm,zwar-zcdm,icolor
         ENDIF
!
!        IF( wires ) WRITE(10,'(a3,1x,3f16.5,i4)') &
!                    &elem(i),xuar,yvar,zwar+pbc*aretebim,icolor+1
!
        ENDDO

ENDIF !!su np

IF(ipas==npas) THEN
   CLOSE(unite)
   CLOSE(unitm)
ENDIF

END SUBROUTINE scrittoc
!
! For the stress
!  IF (interewr) THEN
!     call inter
!  ENDIF
