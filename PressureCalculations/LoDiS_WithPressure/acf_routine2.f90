subroutine acf_routine
use cluster
use enforce 
implicit none 

integer :: i,j,k,t
real :: zero

zero = 0.d0

do i=1,natom
	zero = zero + v_acf(1,i,1)*v_acf(1,i,1) + v_acf(1,i,2)*v_acf(1,i,2) + v_acf(1,i,3)*v_acf(1,i,3)
enddo

do t=1,npas
	do j=1,natom
		
