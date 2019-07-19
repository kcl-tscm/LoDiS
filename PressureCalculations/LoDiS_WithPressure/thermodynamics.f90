subroutine thermodynamics

USE PARACLUSTER  !uso il modulo di definizione dei parametri
USE CLUSTER     !uso il modulo dove definisco variabili e parametri cluster
Implicit None

Real :: tnuova
!Real :: vargau  !Modified by LP (ifort error)
       
if((tcaloric>tinit).and.(tfin<tcaloric)) then
 !!melting process
 tnuova=tinit+nd_proc*deltat
! vargau=vargau*sqrt(tnuova/tfin)  !Modified by LP (ifort error)
 tfin=tnuova
write(*,*)'MELTING PROCESS: loop number',nd_proc,'T-loop=',tfin
elseif((tcaloric<tinit).and.(tfin>tcaloric)) then
!!freezing process
 tnuova=tinit-nd_proc*deltat
! vargau=vargau*sqrt(tnuova/tfin)  !Modified by LP (ifort error)
 tfin=tnuova
write(*,*)'FREEZING PROCESS: loop number',nd_proc,'T-loop=',tfin
endif

END subroutine thermodynamics
