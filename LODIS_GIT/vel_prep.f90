SUBROUTINE VEL_PREP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
USE CLUSTER, only : natom, itype, ipas, vnu, tstep, tfin, irand_seed, g1,g2,g3!,thermostat_off
USE POTENTIAL, only : t2m
USE ENFORCE
!
IMPLICIT NONE
REAL :: dlaran
REAL :: qx(nsiz),qy(nsiz),qz(nsiz)
REAL :: ecinet,poguy,vargau,vel2
REAL :: vdfa
REAL :: vdfae
INTEGER :: nat3
INTEGER :: ir(4)
INTEGER :: i
!
vdfa=vnu*tstep+1.
IF(ipas==1)ir(:)=irand_seed(:)
!
     ecinet=0.d0
     poguy=0.d0
     nat3=natom
!
DO i=1,nat3 

      qx(i)=(vx(i)+t2m(itype(i))*(fx(i)+dfx(i))) 
      qy(i)=(vy(i)+t2m(itype(i))*(fy(i)+dfy(i))) 
      qz(i)=(vz(i)+t2m(itype(i))*(fz(i)+dfz(i)))
!
!THERMALIZATION OF SYSTEM AT TINIT
!Through the Andersen Thermostat, 
!Honeycutt, J. D., and H. C. Andersen, 1987, J. Phys, Chem. 91, 4950.
!Velocity is velocity per time step (in arete units!)
!
!if (thermostat_off) then !skip the thermostat if its turned off
!	!write(999,*) '################SKIPPING THERMOSTAT###############'
!	vx(i)=qx(i) 
 !       vy(i)=qy(i)
  !      vz(i)=qz(i)
!else
  vdfae=dlaran(irand_seed)
     IF(vdfae > vdfa) THEN   
!
           vx(i)=qx(i) 
           vy(i)=qy(i)
           vz(i)=qz(i)
!
     ELSE 
!
!CALL TO THE ANDERSEN THERMOSTAT
!
            vargau=SQRT(2.d0*t2m(itype(i))*cbol*tfin)
            CALL gauss   

             vx(i)=g1*vargau
             vy(i)=g2*vargau
             vz(i)=g3*vargau
     ENDIF  

!endif
!
!decision of thermostat call
!
! Kinetic energy calculation
! POGUY was introduced by Guy Treglia to speed up the system

         vel2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
         poguy=vel2/(4.d0*t2m(itype(i)))
         ecinet=ecinet+poguy
!

ENDDO   !! LOOP on i

!Equipartition Theorem is applied
!
!
      ecin = ecinet
      temp=2.d0*ecin/(3.d0*(nat3)*cbol) 
      etot=ener+ecin 
!
END SUBROUTINE VEL_PREP


