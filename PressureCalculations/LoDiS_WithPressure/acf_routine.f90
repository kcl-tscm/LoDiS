subroutine acf_routine
use CLUSTER
use ENFORCE

implicit none

integer :: i, m, n
real :: zero


do m=1,npas
	vel_act_est(m) = 0.d0
end do

do m = 1,npas-1
 do n = 1, npas-m
  do i = 1, natom
    vel_act_est(m) = vel_act_est(m) + v_acf(n+m,i,1)*v_acf(n,i,1) + v_acf(n+m,i,2)*v_acf(n,i,2) + v_acf(n+m,i,3)*v_acf(n,i,3) 
  end do
 end do
end do

do m=1,npas-1
  vel_act_est(m) = vel_act_est(m)*(npas-m)/(npas)
end do

zero = vel_act_est(1) 
do m=1, npas-1
  vel_act_est(m) = vel_act_est(m)/zero
end do

do m=1,npas-1
write(56,*) m,vel_act_est(m)
end do
close(unit=56)

end subroutine






				

!
