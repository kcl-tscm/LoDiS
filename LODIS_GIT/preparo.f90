SUBROUTINE preparo
USE PARACLUSTER 
USE CLUSTER
USE POTENTIAL
USE ENFORCE
USE DISTANCE
!
Implicit None
  Integer :: i
  Real :: secd,firstd,rvec!,dmax
  Real :: partialx(nsiz),partialy(nsiz),partialz(nsiz)
  Real :: delta
!
!!open(3,file='initial_phase.in',status='unknown')
!
  call init
  IF(type_potential=='rgl') then
   call bigvoi
   write(*,*) 'preparo> IN BIGVOI : LIMIT ', cutoff_end
  ENDIF
  call voisin
   write(*,*) 'preparo> IN VOISIN : LIMIT ', cutz(1,1), cutz(2,2)   !,' Neighbours',nvois(:) 
  call force_choice
   write(*,*) 'preparo> esco da force'
  call therma

  write(*,*) 'preparo> END t=0 TIME STEP'
!
! fine passo iniziale t=0
!fase di termalizzazione alla temperaturatura tinit
!
do ipas=1,npast
     tempo=ipas*tstep*1.d12        !!in picosecondi
!
if(type_potential=='rgl') then
     partialx(1:natom)=0.d0
     partialy(1:natom)=0.d0
     partialz(1:natom)=0.d0

     firstd=0.d0
     secd=0.d0

     do i=1,natom
        partialx(i)=partialx(i)+du(i)            
        partialy(i)=partialy(i)+dv(i)            
        partialz(i)=partialz(i)+dw(i) 
        delta=partialx(i)**2+partialy(i)**2+partialz(i)**2

        if ((delta.gt.secd).and.(delta.le.firstd)) then
           secd=delta
        endif
        if (delta.gt.firstd) then
           secd=firstd
           firstd=delta
        endif
     enddo
     rvec=sqrt(firstd)+sqrt(secd)
     !
     if(rvec.gt.rshell) then
!
       do i=1,natom
           partialx(i)=0.d0
           partialy(i)=0.d0
           partialz(i)=0.d0
        enddo
!
        call bigvoi     
     endif
endif
!
     call voisin !!TO CALCULATE THE NEIGHBOURS AT TIME t
!
     call force_choice !!To CALCULATE FORCE AT TIME t
!
     call vel_prep !! TO CALCULATE THE VEL. AT TIME t
!
!!!Thermalization output
!
     tpar=tpar+temp
     edelta=edelta+etot
!
     IF(sys == 'mon') THEN
      edelta=(edelta+natom*ecoh(1))/(natom**(2.d0/3.d0))
     ELSE
       edelta=(etot+ntipo1*ecoh(1)+ntipo2*ecoh(2))/((natom)**(2.d0/3.d0))
     ENDIF
      IF(ipas==npast) write(*,*)'preparo> THERMALIZATION PROCESS AT ', tpar/dfloat(npast)
      IF(ipas==npast) WRITE(*,*) 'preparo> TOTAL ENERGY, KINETIC EN., POT. ENERGY, DELTA' 
      IF(ipas==npast) &
     & write(*,*)npast,etot,ecin,ener,edelta/dfloat(npast)
!
call therma  !!TO CALCULATE THE NEW POSITIONS
!
enddo 
write(*,*) 'preparo> END THERMALISATION ON NPAST DEFINE IN PARACLUSTER MODULE'
!
END SUBROUTINE preparo
