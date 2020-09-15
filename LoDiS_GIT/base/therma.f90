SUBROUTINE THERMA
USE PARACLUSTER
USE CLUSTER
USE POTENTIAL
USE ENFORCE
USE DISTANCE
!
IMPLICIT NONE
INTEGER :: nat3,is
REAL :: xdmax,xdmin,xd1,xlarge
REAL :: ydmax,ydmin,yd1,ylarge
REAL :: zdmax,zdmin,zd1,zlarge

IF((deposizione=='ya').AND.(nvois(natom)==0)) THEN
        du(natom)=vx(natom)
        dv(natom)=vy(natom)
        dw(natom)=vz(natom)
        u(natom)=u(natom)+du(natom)
        v(natom)=v(natom)+dv(natom)
        w(natom)=w(natom)+dw(natom)
        nat3=natom-1
ELSE
        nat3=natom
ENDIF 
!    
DO is=1,nat3   !!EVOLUTION ACCORDING WITH VELOCITY-VERLET         
!
        du(is)=t2m(itype(is))*fx(is)+vx(is)
        dv(is)=t2m(itype(is))*fy(is)+vy(is)
        dw(is)=t2m(itype(is))*fz(is)+vz(is)    
!
        u(is)=u(is)+du(is)
        v(is)=v(is)+dv(is)
        w(is)=w(is)+dw(is)
!       
ENDDO
! 
IF (wires) THEN
!!PBC along Z
!      zdmax=(pbcz+0.1d0)  !!possible with dilation
      zdmin = MINVAL(z(:)) - 0.21213*(1.00+dilat)
      zdmax = pbcz + 0.49497*(1.00+dilat)
      zlarge=zdmax-zdmin
!      IF (ipas == 1 ) THEN
!       WRITE(*,*) 'IN THERMA compare zlarge and pbc:::'
!       WRITE(*,*) 'MAX-MIN zvalues=',zlarge,'PBC=',pbcz
!       WRITE(*,*) zdmax, zdmin
!      ENDIF
      DO is=1,natom
        zd1=w(is)+z(is)
        IF (zd1<=zdmin) zd1=zd1+pbcz
        IF (zd1>zdmax)  zd1=zd1-pbcz
        w(is)=zd1-z(is)
      ENDDO

ENDIF
!
IF( surface ) then
!      xdmax= (pbcx + 0.7)*(1.d0+dilat)/rac2 !pbcx + 0.49497d0*(1.d0+dilat)
!      ydmax= (pbcy + 0.7)*(1.d0+dilat)/rac2!pbcy + 0.49497d0*(1.d0+dilat)
!      xdmin= -0.3d0*(1.d0+dilat)/rac2 !MINVAL(u(:))-0.21213d0*(1.d0+dilat)
!      ydmin= -0.3d0*(1.d0+dilat)/rac2 !MINVAL(v(:))-0.21213d0*(1.d0+dilat)
      xdmax = (pbcx - 0.3/rac2) * (1.d0 + dilat)
      xdmin = (-0.3d0/rac2)*(1.d0 + dilat) 
      ydmax = (pbcy - 0.3/rac2) * (1.d0 + dilat)
      ydmin = (-0.3d0/rac2)*(1.d0 + dilat) 
      xlarge=xdmax-xdmin
      ylarge=ydmax-ydmin
IF(ipas==1) write(*,*) 'here Xd', xdmin,xdmax,ydmin,ydmax
IF(ipas==1) write(*,*) 'here PBC', pbcx,pbcy
IF(ipas==1) write(*,*) 'here XL', xlarge,ylarge
!
      do is=1,natom
        xd1=u(is)+x(is)
        yd1=v(is)+y(is)
        if (xd1<=xdmin) xd1=xd1+xlarge
        if (xd1> xdmax) xd1=xd1-xlarge
        if (yd1<=ydmin) yd1=yd1+ylarge
        if (yd1> ydmax) yd1=yd1-ylarge
        u(is)=xd1-x(is)
        v(is)=yd1-y(is)
      enddo
      if((ipas==1).or.(ipas==npas)) write(*,*) 'in therma', xlarge,ylarge
ENDIF

IF( bulk ) THEN
      xdmax= pbcx + 0.49497d0*(1.d0+dilat)
      ydmax= pbcy + 0.49497d0*(1.d0+dilat)
      xdmin= MINVAL(u(:))-0.21213d0*(1.d0+dilat)
      ydmin= MINVAL(v(:))-0.21213d0*(1.d0+dilat)
      zdmin = MINVAL(w(:)) - 0.21213*(1.00+dilat)
      zdmax = pbcz + 0.49497*(1.00+dilat)
      xlarge=xdmax-xdmin
      ylarge=ydmax-ydmin
      zlarge=zdmax-zdmin

      DO is=1,natom
        xd1=u(is)+x(is)
        yd1=v(is)+y(is)
        zd1=u(is)+x(is)
         yd1=v(is)+y(is)        
        if (xd1<=xdmin) xd1=xd1+xlarge
        if (xd1> xdmax) xd1=xd1-xlarge
        if (yd1<=ydmin) yd1=yd1+ylarge
        if (yd1> ydmax) yd1=yd1-ylarge
        IF (zd1<=zdmin) zd1=zd1+zlarge
        IF (zd1> zdmax) zd1=zd1-zlarge
        u(is)=xd1-x(is)
        v(is)=yd1-y(is)
        w(is)=zd1-z(is)
      ENDDO
!
ENDIF
!
END SUBROUTINE THERMA
