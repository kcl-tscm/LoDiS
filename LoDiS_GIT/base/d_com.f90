SUBROUTINE d_com
!======================================
! d^2 COM
!======================================
! Finds the squared distance of the Centre of Mass for a bimetallic cluster.
! Calculates the derivatives.
! Use the d_com light version to speed up calculation for taxing simulations.
!======================================
   USE PARACLUSTER  
   USE CLUSTER   
   USE POTENTIAL
   USE ENFORCE
   USE DISTANCE
   USE META
   
   IMPLICIT NONE
   
   integer :: i
   REAL(8) :: cxa, cya, cza, cxb, cyb, czb
   REAL(8) :: d_com_x, d_com_y, d_com_z
   
   
   if (ipas.eq.1) write(*,*) 'd_com> first calculation of d_com'
   
   cxa = 0.d0
   cya = 0.d0
   cza = 0.d0
   cxb = 0.d0
   cyb = 0.d0
   czb = 0.d0
   DO i = 1, natom
      IF(itype(i).EQ.1) THEN
         cxa = cxa + u(i) +x(i)
         cya = cya + v(i) +y(i)
         cza = cza + w(i) +z(i)
      ELSE
         cxb = cxb + u(i) +x(i)
         cyb = cyb + v(i) +y(i)
         czb = czb + w(i) +z(i)
      ENDIF  
   ENDDO 
   cxa = cxa/natoma 
   cya = cya/natoma 
   cza = cza/natoma
   cxb = cxb/natomb 
   cyb = cyb/natomb
   czb = czb/natomb

   d_com_x = cxa - cxb
   d_com_y = cya - cyb
   d_com_z = cza - czb

   d_com_cv = d_com_x*d_com_x +d_com_y*d_com_y +d_com_z*d_com_z

   ! Derivatives
   d_com_dcv_dx = 2.d0 *d_com_x
   d_com_dcv_dy = 2.d0 *d_com_y
   d_com_dcv_dz = 2.d0 *d_com_z
   DO i = 1, natom
      IF(itype(i).EQ.1) THEN
         d_com_dcv_dx(i) = d_com_dcv_dx(i) /natoma
         d_com_dcv_dy(i) = d_com_dcv_dy(i) /natoma
         d_com_dcv_dz(i) = d_com_dcv_dz(i) /natoma
      ELSE
         d_com_dcv_dx(i) = -d_com_dcv_dx(i) /natomb
         d_com_dcv_dy(i) = -d_com_dcv_dy(i) /natomb
         d_com_dcv_dz(i) = -d_com_dcv_dz(i) /natomb
      ENDIF 
   ENDDO

END SUBROUTINE d_com
