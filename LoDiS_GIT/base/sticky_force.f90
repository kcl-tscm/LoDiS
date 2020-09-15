SUBROUTINE sticky_force
  
  USE CLUSTER 
  USE ENFORCE
  USE module_sticky

  IMPLICIT NONE
  
  INTEGER :: i, k
  
  ! This is done in main_MD because only the same atoms always have these
  ! quantities different from zero 
  ! sticky_fx = 0.d0
  ! sticky_fy = 0.d0
  ! sticky_fz = 0.d0

  DO i= 1, sticky_atoms
     k = sticky_labels(i)
     ! u(:),v(:),w(:) ! displacement at time t     
     sticky_fx(k) = -sticky_k* u(k)
     sticky_fy(k) = -sticky_k* v(k)
     sticky_fz(k) = -sticky_k* w(k)        
  ENDDO

END SUBROUTINE sticky_force
