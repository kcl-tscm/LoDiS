SUBROUTINE force_par
USE PARACLUSTER
USE CLUSTER    !module where cluster variables and parameters are defined
USE POTENTIAL
USE ENFORCE
USE DISTANCE
!
IMPLICIT NONE
INTRINSIC sqrt
INTEGER :: nat3d,i,k
REAL :: delta,ermi,mik,dmik,fik,dfik,wik,dwik,vik
REAL :: xik,yik,zik
REAL :: dik0,dik,rik,rik2,rik23
REAL :: pmorse,eneri,wtot
REAL,ALLOCATABLE  :: fDx(:),fDy(:),fDz(:)
!
!
!      common/bimW/C6,C8,C10,C12
!      common/bimF/dmu,delta
!      common/bimM/dM0,DIST0,tau      
!
       ! parameters dmu, delta which compare in Fermi-type potential       
       !F(x)= {1+ exp(x-mu)/delta}**(-1)
!        parameter(dmu0=10.05d0)
!        parameter(delta0=1.04d0)=epsilon
        
        ! dM0, tau, DIST0  which compare in Morse potential
        ! M(x)=dM0*(exp(tau*(1-x/DIST0)){exptau*(1-x/DIST0)-2})
!        parameter(dM00=0.3d0)
!        parameter(tau0=9.75d0)
!        parameter(DIST00=10.3d0)
        
        ! parameters of van der Waal part
        !W(x)= - (C6/x**6 + C8/x**8 + C10/x**10 + C12/x**12)
!        parameter(C60=75600.d0)
!        parameter(C80=9122400.d0)
!        parameter(C100=2.09d8)
!        parameter(C120=7.78d10)

! pay attention to the units: in fact,
! we have all the routine which calculates the velocities and the positions
! on eV (for energy!) and on units of arete (for position!)
! Doye considers energy in eV BUT position in A
! the parameters are in eV and angstrom

      ! for C60 cluster => PR potential
      ener=0.d0                                                         
      eneri=0.d0
      delta=epsilon                                                        
      dik0=dist(1)  !! distanza in arete ai NN
      nat3d=natom

ALLOCATE(fDx(nat3d),fDy(nat3d),fDz(nat3d))
! conservo i vecchi valori della forza in dfx,dfy,dfz
      DO i=1,nat3d
       dfx(i)=fx(i)       
       dfy(i)=fy(i)             
       dfz(i)=fz(i)
      ENDDO
      DO i=1,nat3d                        
        fx(i)=0.d0                                                       
        fy(i)=0.d0                                                       
        fz(i)=0.d0
        fDx(i)=0.d0                                                       
        fDy(i)=0.d0                                                       
        fDz(i)=0.d0
      ENDDO  
      
!! part of PR: E=(1-F(x))*M(x)+F(x)*W(x)
        !W(x)= - (C6/x**6 + C8/x**8 + C10/x**10 + C12/x**12)
        !F(x)= {1+ exp(x-mu)/delta}**(-1)        
        ! M(x)=dM0*(exp(tau*(1-x/DIST0)){exptau*(1-x/DIST0)-2})

!       print*,tau, DIST0, dmu, dM0, delta, C6, C8, C10, C12
  DO i=1,nat3d                                                   
    DO k=1+i,nat3d                          
          
          xik=x(k)+u(k)-x(i)-u(i)                                           
          yik=y(k)+v(k)-y(i)-v(i)                                           
          zik=z(k)+w(k)-z(i)-w(i)                                           

          dik=sqrt(xik*xik+yik*yik+zik*zik)   !in unita' di arete
          rik=aretebim*dik                        !in A
          rik2=rik*rik
          rik23=rik2**3.d0
 
 !!calcolo ora le forze come Doye in eV e A
 !!poi devo passare le forze in unita' eV e arete
 !! f(i) per noi= f(i)( in eV/A)*10.05*rac2  !!
          
         ermi=EXP((rik-dmu)/delta)
         fik=1.d0/(1.d0+ermi)
         dfik=-ermi/(delta*(1.d0+ermi)*(1.d0+ermi))
              
         pmorse=EXP(tau*(1.d0-rik/DIST0))
!         write(*,*) dM0, 'who am i '
         mik=dM0*pmorse*(pmorse-2.d0)
!         write(*,*) mik,'mik for i',i,'and k',k
         dmik=(2.d0*tau*dM0*pmorse*(1.d0-pmorse))/DIST0
         wik=-(C6+(C8+(C10+C12/rik2)/rik2)/rik2)/rik23
         dwik=(6*C6+(8*C8+(10*C10+12*C12/rik2)/rik2)/rik2)/(rik*rik23)
              
    ! this is the potential
              
         vik=fik*mik+(1.d0-fik)*wik
         ener=ener+vik 
         eneri=eneri+fik*mik  !! in eV
         
!         print*,i, k, rik, vik
! this is the force  in eV/A

         wtot=mik*dfik+fik*dmik+(1.d0-fik)*dwik-dfik*wik
         wtot=-wtot
              
              fDx(i)=fDx(i)-aretebim*xik*wtot/rik
              fDy(i)=fDy(i)-aretebim*yik*wtot/rik
              fDz(i)=fDz(i)-aretebim*zik*wtot/rik

              fDx(k)=fDx(k)+aretebim*xik*wtot/rik
              fDy(k)=fDy(k)+aretebim*yik*wtot/rik
              fDz(k)=fDz(k)+aretebim*zik*wtot/rik
               
        !fx(k)=fDx(k)*arete
        !fy(k)=fDy(k)*arete
        !fz(k)=fDz(k)*arete 

ENDDO !su k

        fx(i)=fDx(i)*aretebim
        fy(i)=fDy(i)*aretebim
        fz(i)=fDz(i)*aretebim 
        ener=eneri
!!        forza=dsqrt(fx(i)**2+fy(i)**2+fz(i)**2)
ENDDO
DEALLOCATE(fDx,fDy,fDz)

END SUBROUTINE force_par
