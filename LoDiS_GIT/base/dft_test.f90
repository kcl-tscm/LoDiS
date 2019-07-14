subroutine dff_test2
use enforce
use cluster
use paracluster

implicit none

integer :: i,n,k
real, allocatable :: dft(:)
real :: arg,freq

allocate(dft(npas))

do i = 1,npas
	dft(i) = 0.d0
enddo

do k=1,npas
	freq = 1.d0/(k*tstep)
	arg = -(2.d0*pi*dble(1.d0/freq)/npas)
	do n =1,npas
		dft(k) = dft(k) + cos(arg*dble(n))*vel_act_est(n)
	enddo
	write(60,*) freq,dft(k)*dft(k)
enddo

end subroutine dff_test2

!======================================================


subroutine dff_test
use ENFORCE!, only: vel_act_est
use CLUSTER
use PARACLUSTER

implicit none

!variables needed:      - N = npas + 1
!			- k = index for N
!			- n = index for chosen frequency 
!PLAN:
!	INITIALISE -- size of dft

integer :: N,N_i,freq,a,b,i
real, allocatable :: dft(:,:)
real :: arg

!N = npas + 1!SIZE(vel_acf_est)
!do i = 0,npas
!	write(*,*) vel_act_est(i)
!enddo

allocate(dft(0:npas,2))

do a = 0,npas
	dft(a,1) = 0.d0
	dft(a,2) = 0.d0
end do

!open(11,file = 'fort.59', status = 'old' , action = 'read')

!do b = 0,N
!	read(11,*) dft(b)
!end do
!open(60, file = 'fort.60', status = 'new', action = 'write')

do freq = 1,npas
	do N_i = 1,npas	
	arg = -(2*pi*dble(N_i)/npas)
	!write(13,*) N_i,arg
	dft(freq,1) = dft(freq,1) + cos(arg*dble(freq)*tstep)*vel_act_est(N_i)!/npas
	!write(14,*) N_i ,dft(N_i,1) , cos(arg) ,vel_act_est(freq)
	!dft(N_i,2) = dft(N_i,2) + sin(arg*dble(freq))*vel_act_est(freq)!/npas
	!write(15,*) N_i,dft(N_i,2) , sin(arg) ,vel_act_est(freq)
	end do
write(60,*) freq , dft(freq,1)*dft(freq,1) !+ dft(N_i,2)*dft(N_i,2)
end do

!do N_i = 0,npas
	!dft(N_i,1)=dft(N_i,1)/npas
	!dft(N_i,2)=dft(N_i,2)/npas
	!write(60,*) N_i , dft(N_i,1)*dft(N_i,1) + dft(N_i,2)*dft(N_i,2)
!end do

end subroutine dff_test
