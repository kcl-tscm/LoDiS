SUBROUTINE FENVIRONMENT

  USE PARACLUSTER
  USE CLUSTER 
  USE POTENTIAL
  USE ENFORCE
  USE DISTANCE
  USE ENVIRONMENT

  IMPLICIT NONE

integer :: n
 !Chooses the CN calculation method based on the input
IF (cn_cal) THEN
  CALL ENVCN
ELSE
  CALL ENVCN2
END IF

 ! Now I have cn(i) and and its derivative w.r.t. x, y, z
 ! They are called 
 ! env_cn(i)                                    adimensional
 ! env_dcn_dx(i), env_dcn_dy(i), env_dcn_dz(i)  [1/L]
 ! In this moment the unit is [1/Angstrom], and can be modified in mgo_coord_no
 ! Normalised positions  (in Angstrom and then *mgo_alat)

 ! Inizialisation of local variables
 ener_env_atom(:) = 0.d0
 !------------------------------------------------------
 ! Inizialisation of global variables

 ener_env    = 0.d0 ! This is global  
 env_fx(:)   = 0.d0
 env_fy(:)   = 0.d0
 env_fz(:)   = 0.d0

 DO n = 1, natom
 if (env_cn(n) .ge. 11.5d0) then
   env_fx(n) = 0.d0
   env_fx(n) = 0.d0
   env_fx(n) = 0.d0
 else
! if ( elem(n) .eq. elem1 ) then
    ener_env_atom(n) = - eta_a * ((12.d0 - env_cn(n))**pot_a)

!    write(30,*) n, env_cn(n), ener_env_atom
!write(130,*) n, ener_env_atom
    env_fx(n)=  eta_a * pot_a * ((12.d0 - env_cn(n))**(pot_a-1.d0)) * env_dcn_dx(n)
    env_fy(n)=  eta_a * pot_a * ((12.d0 - env_cn(n))**(pot_a-1.d0)) * env_dcn_dy(n)
    env_fz(n)=  eta_a * pot_a * ((12.d0 - env_cn(n))**(pot_a-1.d0)) * env_dcn_dz(n)
! else
!    ener_env_atom = - eta_b * ((12.d0 - env_cn(n))**pot_b)
!    env_fx(n)=  eta_b * pot_b * ((12.d0 - env_cn(n))**(pot_b-1.d0)) * env_dcn_dx(n)
!    env_fy(n)=  eta_b * pot_b * ((12.d0 - env_cn(n))**(pot_b-1.d0)) * env_dcn_dy(n)
!    env_fz(n)=  eta_b * pot_b * ((12.d0 - env_cn(n))**(pot_b-1.d0)) * env_dcn_dz(n)
 end if
! end if
    ! Converting to arete(1) units
    env_fx(n) = -env_fx(n) !*arete(1)
    env_fy(n) = -env_fy(n) !*arete(1)
    env_fz(n) = -env_fz(n) !*arete(1)
!   write(32,*) env_dcn_dx(n), env_dcn_dy(n), env_dcn_dz(n)    
!   write(31,*) env_cn(n), env_fx(n), env_fy(n), env_fz(n)
    ! Energy from the environment-metal interaction
    ener_env = sum(ener_env_atom)

 END DO

END SUBROUTINE FENVIRONMENT
