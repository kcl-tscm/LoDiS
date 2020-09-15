subroutine suppress

use CLUSTER
use ENFORCE
use POTENTIAL

implicit none

real*8 :: group_mass
real*8, dimension(3) :: group_p
real*8, dimension(3) :: group_j
real*8, dimension(3) :: group_w
real*8, dimension(3) :: jcm
real*8, dimension(natom,3) :: j0
real*8, dimension(3,3) :: group_I , I_cm , I_inv , I_test
real*8 :: tm , tm_1
real*8, dimension(natom) :: dvx,dvy,dvz
real*8, dimension(3) :: vcm,vcm_test
real*8, dimension(natom) :: dx,dy,dz
real*8, dimension(3) :: rcm
real*8, dimension(3,3) :: tmp
real*8, dimension(natom) :: rx,ry,rz
real*8 :: xy,xz,yz
real*8 :: deter,deter_1
real*8 :: fac,rfac
real*8 :: r2,r2cm
real*8 :: sum_test

integer :: i , j , k

!initialise everything to zero 

do i = 1,3
	group_p(i) = 0.d0
	group_j(i) = 0.d0
	group_w(i) = 0.d0
	jcm(i) = 0.d0
	vcm(i) = 0.d0
	rcm(i) = 0.d0
	vcm_test(i) = 0.d0
enddo

do i = 1,natom
	dx(i) = 0.d0
	dy(i) = 0.d0
	dz(i) = 0.d0
	dvx(i) = 0.d0
	dvy(i) = 0.d0
	dvz(i) = 0.d0
	j0(i,1) = 0.d0
	j0(i,2) = 0.d0
	j0(i,3) = 0.d0
enddo

do i = 1,3
	do j = 1,3
		group_I(i,j) = 0.d0
		I_cm(i,j) = 0.d0
		tmp(i,j) = 0.d0
		I_test(i,j) = 0.d0
		I_inv(i,j) = 0.d0
	enddo
enddo

group_mass = 0.d0
tm = 0.d0
tm_1 = 0.d0
deter = 0.d0

!calculate group values 

do i = 1,natom

	!total mass
	group_mass = group_mass + mass(1)

	!total linear momentum 
	group_p(1) = group_p(1) + mass(1)*vx(i)
	group_p(2) = group_p(2) + mass(1)*vy(i) 
	group_p(3) = group_p(3) + mass(1)*vz(i)

	!mass-weighted position (to be used for calculating center-of-mass coordinates)
	rcm(1) = rcm(1) + mass(1)*x(i)
	rcm(2) = rcm(2) + mass(1)*y(i)
	rcm(3) = rcm(3) + mass(1)*z(i)

	!total angular momentum 
	j0(i,1) = (y(i)*vz(i) - z(i)*vy(i))
	j0(i,2) = (z(i)*vx(i) - x(i)*vz(i))
	j0(i,3) = (x(i)*vy(i) - y(i)*vx(i))

	group_j(1) = group_j(1) + mass(1)*j0(i,1)
	group_j(2) = group_j(2) + mass(1)*j0(i,2)
	group_j(3) = group_j(3) + mass(1)*j0(i,3)

	!update the inertia tensor 
	xy = x(i)*y(i)*mass(1)
	xz = x(i)*z(i)*mass(1)
	yz = y(i)*z(i)*mass(1)
	
	group_I(1,1) = group_I(1,1) + (y(i)*y(i) + z(i)*z(i))*mass(1)
	group_I(2,2) = group_I(2,2) + (x(i)*x(i) + z(i)*z(i))*mass(1)
	group_I(3,3) = group_I(3,3) + (x(i)*x(i) + y(i)*y(i))*mass(1)
	group_I(1,2) = group_I(1,2) - xy
	group_I(2,1) = group_I(2,1) - xy
	group_I(1,3) = group_I(1,3) - xz
	group_I(3,1) = group_I(3,1) - xz
	group_I(2,3) = group_I(2,3) - yz
	group_I(3,2) = group_I(3,2) - yz
enddo

tm = group_mass
tm_1 = 1.d0/tm

!calculate position and velocity center-of-mass coordinates
do i = 1,3
	rcm(i) = rcm(i)*tm_1
	vcm(i) = group_p(i)*tm_1
enddo

!calculate center of mass contribution to angular momentum and subtract from total angular momentum
jcm(1) = (rcm(2)*vcm(3) - rcm(3)*vcm(2))
jcm(2) = (rcm(3)*vcm(1) - rcm(1)*vcm(3))
jcm(3) = (rcm(1)*vcm(2) - rcm(2)*vcm(1))

do i = 1,3
	group_j(i) = (group_j(i) - tm*jcm(i))
enddo


!update the center-of-mass contribution to the inertia tensor 
xy = rcm(1)*rcm(2)*tm
xz = rcm(1)*rcm(3)*tm
yz = rcm(2)*rcm(3)*tm
I_cm(1,1) = I_cm(1,1) + (rcm(2)*rcm(2) + rcm(3)*rcm(3))*tm
I_cm(2,2) = I_cm(2,2) + (rcm(1)*rcm(1) + rcm(3)*rcm(3))*tm
I_cm(3,3) = I_cm(3,3) + (rcm(1)*rcm(1) + rcm(2)*rcm(2))*tm
I_cm(1,2) = I_cm(1,2) - xy
I_cm(2,1) = I_cm(2,1) - xy
I_cm(1,3) = I_cm(1,3) - xz
I_cm(3,1) = I_cm(3,1) - xz
I_cm(2,3) = I_cm(2,3) - yz
I_cm(3,2) = I_cm(3,2) - yz

!subtract center-of-mass contribution from inertia tensor 

do i = 1,3
	do j = 1,3
		group_I(i,j) = (group_I(i,j) - I_cm(i,j))
	enddo
enddo

!calculate inverse of inertia tensor 

do i = 1,3
	do j = 1,3
		tmp(i,j) = group_I(i,j)
	enddo
enddo
!

!rfac = (tmp(1,1) + tmp(2,2) + tmp(3,3))/3
!fac = 1.d0/rfac

!do i = 1,3
!	do j = 1,3
!		tmp(i,j) = tmp(i,j)*fac
!	enddo
!enddo

deter = deter + tmp(1,1)*(tmp(2,2)*tmp(3,3) - tmp(2,3)*tmp(3,2)) 
deter = deter + tmp(3,1)*(tmp(1,2)*tmp(2,3) - tmp(2,2)*tmp(1,3))
deter = deter - tmp(2,1)*(tmp(1,2)*tmp(3,3) - tmp(3,2)*tmp(1,3)) 

deter_1 = 1.d0/deter

I_inv(1,1) = deter_1*(tmp(2,2)*tmp(3,3) - tmp(3,2)*tmp(2,3))
I_inv(1,2) = -deter_1*(tmp(1,2)*tmp(3,3) - tmp(3,2)*tmp(1,3))
I_inv(1,3) = deter_1*(tmp(1,2)*tmp(2,3) - tmp(2,2)*tmp(1,3))
I_inv(2,1) = -deter_1*(tmp(2,1)*tmp(3,3) - tmp(3,1)*tmp(2,3))
I_inv(2,2) = deter_1*(tmp(1,1)*tmp(3,3) - tmp(3,1)*tmp(1,3))
I_inv(2,3) = -deter_1*(tmp(1,1)*tmp(2,3) - tmp(2,1)*tmp(1,3))
I_inv(3,1) = deter_1*(tmp(2,1)*tmp(3,2) - tmp(3,1)*tmp(2,2))
I_inv(3,2) = -deter_1*(tmp(1,1)*tmp(3,2) - tmp(3,1)*tmp(1,2))
I_inv(3,3) = deter_1*(tmp(1,1)*tmp(2,2) - tmp(2,1)*tmp(1,2))

!do i = 1,3
!	do j = 1,3
!		I_inv(i,j) = I_inv(i,j)*fac
!	enddo
!enddo


!Now inverse of I = I_inv(3,3)
!-->I_inv x J = W
!
!do i = 1,3
!	do j = 1,3
!		do k = 1,3
!			I_test(i,j) = I_test(i,j) + I_inv(i,k)*group_I(k,j)
!		enddo
!	enddo
!enddo

!do i = 1,3	
!	write(8999,*) I_test(i,1) , I_test(i,2) , I_test(i,3)
!	write(7999,*) I_inv(i,1) , I_inv(i,2) , I_inv(i,3)
!enddo
!write(8999,*) 
!write(8999,*)
!write(7999,*) 
!write(7999,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!calculate the group angular momentum 
group_w(1) = I_inv(1,1)*group_j(1) + I_inv(1,2)*group_j(2) + I_inv(1,3)*group_j(3)
group_w(2) = I_inv(2,1)*group_j(1) + I_inv(2,2)*group_j(2) + I_inv(2,3)*group_j(3)
group_w(3) = I_inv(3,1)*group_j(1) + I_inv(3,2)*group_j(2) + I_inv(3,3)*group_j(3)

!Calculate the new velocities with subtracted linear and angular velocities about the center of mass
do i = 1,natom

	vx(i) = vx(i) - vcm(1)
	vy(i) = vy(i) - vcm(2)
	vz(i) = vz(i) - vcm(3)

	dx(i) = x(i) - rcm(1)
	dy(i) = y(i) - rcm(2)
	dz(i) = z(i) - rcm(3)

	dvx(i) = (group_w(2)*z(i) - group_w(3)*y(i))
	dvy(i) = (group_w(3)*x(i) - group_w(1)*z(i))
	dvz(i) = (group_w(1)*y(i) - group_w(2)*x(i))

	vx(i) = vx(i) - dvx(i)
	vy(i) = vy(i) - dvy(i)
	vz(i) = vz(i) - dvz(i)
enddo

end subroutine
