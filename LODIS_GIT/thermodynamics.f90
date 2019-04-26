subroutine thermodynamics

USE PARACLUSTER 
USE CLUSTER    
Implicit None

Real :: tnuova
!Real :: vargau  
       
if((tcaloric>tinit).and.(tfin<tcaloric)) then
 !!melting process
 tnuova=tinit+nd_proc*deltat
! vargau=vargau*sqrt(tnuova/tfin)
 tfin=tnuova
write(*,*)'MELTING PROCESS: loop number',nd_proc,'T-loop=',tfin
elseif((tcaloric<tinit).and.(tfin>tcaloric)) then
!!freezing process
 tnuova=tinit-nd_proc*deltat
! vargau=vargau*sqrt(tnuova/tfin)
 tfin=tnuova
write(*,*)'FREEZING PROCESS: loop number',nd_proc,'T-loop=',tfin
endif

END subroutine thermodynamics
