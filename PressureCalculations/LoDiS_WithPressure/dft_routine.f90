subroutine dft_routine
use enforce
use cluster


implicit none

integer :: a,c,b
real :: i
real, allocatable :: dft(:,:),vel(:)
real :: arg,freq,pie
pie = 3.14159265359
!npas = 20000
allocate(vel(npas))
allocate(dft(npas,1:2))
open(11,file = 'fort.56',status = 'old',action = 'read')
do a=1,npas-1
	read(11,*) i, vel(a)
enddo
do a = 1,npas-1
	dft(a,1) = 0.d0
	dft(a,2) = 0.d0
enddo

do b=1,npas-1
	!freq = 1.d0/(b*tstep)
	arg = -(2.d0*pie*dble(b)/npas)
	do c =1,npas-1
		dft(b,1) = dft(b,1) + cos(arg*dble(c))*vel(c)!*vel_act_est(n)
		dft(b,2) = dft(b,2) - sin(arg*dble(c))*vel(c)!*vel_act_est(n)
	enddo

	write(60,*) b,dft(b,1)*dft(b,1) + dft(b,2)*dft(b,2)
enddo

end subroutine dft_routine

