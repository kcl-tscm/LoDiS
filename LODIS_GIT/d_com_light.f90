SUBROUTINE d_com_light
!======================================
! d^2 COM light version
!======================================
! Calculates the squared distance of the Centre of Mass.
! No derivatives! 
!=====================
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
   
   
   if (ipas.eq.1) write(*,*) 'd_com_light> first calculation of d_com'
   
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

END SUBROUTINE d_com_light

