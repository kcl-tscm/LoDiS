SUBROUTINE scrittometa
  USE PARACLUSTER  !uso il modulo di definizione dei parametri
  USE CLUSTER      !uso il modulo dove definisco variabili e parametri cluster
  USE POTENTIAL
  USE ENFORCE
  USE META
  USE SUBSTRATE

  IMPLICIT NONE

  INTEGER :: i, nt
  INTEGER :: icolor, np
  REAL :: xuar, yvar, zwar, xcdm, ycdm, zcdm
  REAL :: emedia                      !LP

  !
  !!! Writing the output for a caloric curve
  !!! energy.out contains the energy output
  !!! out~.xyz the positions 
  !
      
  IF(ipas==1) THEN
     OPEN (unite,file='meta.out',status='unknown',position='append')
     IF(mgo_substrate) THEN
        WRITE(unite,'(a99)') '## STEP TIME[ps] CV(1) CV(2) Epot[eV] Etot[eV] Ekin[eV] <Edelta[eV]> <Etot> <T[K]> T[K] Esub[eV]'
     ELSE
        WRITE(unite,'(a89)') '## STEP TIME[ps] CV(1) CV(2) Epot[eV] Etot[eV] Ekin[eV] <Edelta[eV]> <Etot> <T[K]> T[K]'
     ENDIF
     OPEN (unitm,file='movie.xyz', status='unknown', access='append')
  ENDIF
 
  !Average quantities
  tpar  = tpar + temp
  edelta= edelta + etot
  ! emedia= emedia + etot !!!!!!!!!!!!!!! Not initialized!
! CONVERTING PRESSURE TO PROPER UNITS	
	pressure_i(:)=pressure_i(:)*4.d0/(3.d0*arete(1)**3)*1.602d2 !convertita in GPa)
!enddo

	pressure_total = SUM(pressure_i)
	pressure_cumulative_total =  SUM(pressure_cumulative_i)/(ipas) ! To redefine by incrementing and not overwriting pres_tot
  !------------------------------------------------------------------------------
  ! Writing each scrivo-step
  ! (This part is executed only each scrivo-step!)
  
  np=MOD(ipas,scrivo)
  IF (np==0)  THEN 
  
     ! Average Energy
     emedia=edelta/scrivo
     
     ! EDelta calculation
     edelta=edelta/(scrivo)
     
     IF(sys == 'mon') THEN 
        edelta=(edelta+natom*ecoh(1))/(natom**(2.d0/3.d0))
     ELSE
        edelta=(edelta+ntipo1*ecoh(1)+ntipo2*ecoh(2))/((natom)**(2.d0/3.d0))
     ENDIF
     
     ! Average Temperature calculation
     tpar=tpar/(scrivo)
     
     IF(mgo_substrate) THEN
        if  (num_cv.eq.1) then 
           WRITE(unite,'(1i12, 1x, 1f15.4, 1x, 2f14.5, 1x, 5f14.5, 1x, 3f14.5)')&
           & ipas, tempo, coll_var(1), 0.d0, ener, etot, ecin, edelta, emedia,tpar, temp, ener_sub
        else
           WRITE(unite,'(1i12, 1x, 1f15.4, 1x, 2f14.5, 1x, 5f14.5, 1x, 3f14.5)')&
           & ipas, tempo, coll_var(1), coll_var(2), ener, etot, ecin, edelta, emedia, tpar, temp, ener_sub
        end if
     ELSE
        if  (num_cv.eq.1) then 
           WRITE(unite,'(1i12, 1x, 1f15.4, 1x, 2f14.5, 1x, 5f14.5, 1x, 2f14.5)')&
           & ipas, tempo, coll_var(1), 0.d0, ener, etot, ecin, edelta, emedia,tpar, temp
        else
           WRITE(unite,'(1i12, 1x, 1f15.4, 1x, 2f14.5, 1x, 5f14.5, 1x, 2f14.5)')&
           & ipas, tempo, coll_var(1), coll_var(2), ener, etot, ecin, edelta, emedia, tpar, temp
        end if
     ENDIF
              
     ! Resetting variables
     emedia = 0.d0
     edelta = 0.d0
     tpar   = 0.d0
     
  ENDIF ! su np

  !------------------------------------------------------------------------------
  ! Writing the movie.xyz at each metaframe step
  if ( mod(ipas, metaframe).eq.0) then
     nt = natom
     !
     xcdm = 0.d0  
     ycdm = 0.d0
     zcdm = 0.d0
     !
     IF(mgo_substrate) THEN
        WRITE(unitm,*) nt
        WRITE(unitm,'(i12,2(a3,i4),1x,1f14.5,1x,3(f9.4,1x))') &
        & ipas,elem1,ntipo1,elem2,ntipo2,etot, xcdm, ycdm, zcdm
        DO i=1,nt          
           icolor=itype(i)                    
           xuar=(x(i)+u(i))*aretebim*fattor
           yvar=(y(i)+v(i))*aretebim*fattor
           zwar=(z(i)+w(i))*aretebim*fattor
           !         

		IF (output_pres) THEN
			write(unitm,'(a3,1x,3f16.5,i4, 1x, f7.3)') elem(i), xuar, yvar, zwar, icolor, pressure_i(i)
		ELSE		
			WRITE(unitm,'(a3,1x,3f16.5,i4)') elem(i),xuar,yvar,zwar,icolor
            	ENDIF

        ENDDO   
     ELSE
        WRITE(unitm,*) nt
        DO i =1,nt
           xcdm = xcdm + (x(i)+u(i))
           ycdm = ycdm + (y(i)+v(i))
           zcdm = zcdm + (z(i)+w(i))
        ENDDO
        xcdm = xcdm*aretebim*fattor / nt
        ycdm = ycdm*aretebim*fattor / nt
        zcdm = zcdm*aretebim*fattor / nt
        !
        WRITE(unitm,'(i12,2(a3,i4),1x,1f14.5,1x,3(f9.4,1x))') &
        & ipas,elem1,ntipo1,elem2,ntipo2,etot, xcdm, ycdm, zcdm
        DO i=1,nt          
           icolor=itype(i)                    
           xuar=(x(i)+u(i))*aretebim*fattor
           yvar=(y(i)+v(i))*aretebim*fattor
           zwar=(z(i)+w(i))*aretebim*fattor
           !
		IF (output_pres) THEN
			write(unitm,'(a3,1x,3f16.5,i4, 1x, f7.3)') elem(i), xuar-xcdm,yvar-ycdm,zwar-zcdm,icolor,&
&pressure_i(i)
		ELSE		
			WRITE(unitm,'(a3,1x,3f16.5,i4)') elem(i),xuar-xcdm,yvar-ycdm,zwar-zcdm,icolor
            	ENDIF
        ENDDO   
     ENDIF
   
  end if
  !-----------------------------------------------------------------------------
  
  if (ipas.eq.npas) then
     close(unite)
     close(unitm)
  end if

END SUBROUTINE scrittometa


