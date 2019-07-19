SUBROUTINE force_gir
USE PARACLUSTER
USE CLUSTER     
USE POTENTIAL
USE ENFORCE
USE DISTANCE
!
!uso il modulo dove defin1:natomco variabili e parametri cluster
!
Implicit None
Intrinsic sqrt
Integer :: i,k!,j,icolor,ierr
Integer :: nat3d!,itypik
Real  :: fDx(nsiz),fDy(nsiz),fDz(nsiz)
Real :: rik,rik3,rik4,rik9,rik10,rikpiu,rikmeno
Real :: rikp3,rikp4,rikp9,rikp10
Real :: rikm3,rikm4,rikm9,rikm10
Real :: dik0,dik,xik,yik,zik,xzz,yzz,zzz!,xuar,yvar,zwar
Real :: wtot,wtotA,wtotR,wzz,conver,wconv,alpha,beta,vaik,vrik
Real :: eneri!,eri!,f,for,forsudik,denik
!!Real :: fbx,fby,fbz
!!Real :: espo,pexp,qexp
!!Real :: psudik0,dueapsudik0,dueq,qsi2,qsi2qsudik0
!!Real :: den(nsiz),frx(nsiz),fry(nsiz),frz(nsiz)
      
!aretebim=arete(1)      
! pay attention to the units: in fact,
! we have all the routine which calculates the velocities and the positions
! on eV (for energy!) and on units of arete (for position!)
! Doye considers energy in eV BUT position in A
! the parameters are in eV and angstrom
      ! for C60 cluster => girifalco potential
      ener=0.d0  
      eneri=0.d0  
!!      dik0=dc1  !! distance NN in arete 

        dik0=dist(1)

! old values for force dfx,dfy,dfz
nat3d=natom

        ! old values for force dfx,dfy,dfz

      do i=1,nat3d
       dfx(i)=fx(i)       
       dfy(i)=fy(i)             
       dfz(i)=fz(i)
      enddo
      do i=1,nat3d                                                                                                            
        fx(i)=0.d0                                                       
        fy(i)=0.d0                                                       
        fz(i)=0.d0
        fDx(i)=0.d0                                                       
        fDy(i)=0.d0                                                       
        fDz(i)=0.d0
      enddo  
!!change      Natomi=Nball/arete=aretebim
!! attractive part of Girifaclo: VA(x)        
        alpha=Nball*Nball*ALJ/(12.d0*dball**6.d0)
        
!! repulsive part of Girifaclo:  VR(x)
        beta=Nball*Nball*BLJ/(90.d0*dball**12.d0)
        conver=aretebim/dball
        wconv=aretebim/(dball*dball)
                
 do i=1,nat3d                                                   
   do  k=1+i,nat3d                          
          
          xik=x(k)+u(k)-x(i)-u(i)
          yik=y(k)+v(k)-y(i)-v(i)
          zik=z(k)+w(k)-z(i)-w(i)

          dik=sqrt(xik*xik+yik*yik+zik*zik)   !in unita' di arete
          rik=dik*conver
                    
          rik3=rik*rik*rik
          rik4=rik*rik3
          rik9=rik3*rik3*rik3
          rik10=rik*rik9
          
          rikpiu=(rik+1.d0)
          rikp3=rikpiu*rikpiu*rikpiu
          rikp4=rikp3*rikpiu
          rikp9=rikp3*rikp3*rikp3
          rikp10=rikp9*rikpiu
          
          rikmeno=(rik-1.d0)
          rikm3=rikmeno*rikmeno*rikmeno
          rikm4=rikm3*rikmeno
          rikm9=rikm3*rikm3*rikm3
          rikm10=rikm9*rikmeno
 
 !!calculus the forces in eV e A
 !!then ic onvert again in eV and arete
 !! f(i)= f(i)(in eV/A)*arete  !!
                        
        VAik=alpha*(1/(rik*rikm3)+1/(rik*rikp3)-2/rik4)      
        VRik= beta*(1/(rik*rikm9)+1/(rik*rikp9)-2/rik10)
    
    ! this is the potential
              
         eneri=eneri+VRik       !! in eV     
         ener=ener-VAik+VRik    !! in eV
         

    ! this is the force  in eV/A

         wtotA=alpha*(8/rik3-(4*rik-1)/rikm4-(4*rik+1)/rikp4)
         wtotR=beta*(8/rik9-(10*rik-1)/rikm10-(10*rik+1)/rikp10)
         wtot=-(wtotR-wtotA)

              wzz=wtot*wconv/rik3
              xzz=wzz*xik
              yzz=wzz*yik
              zzz=wzz*zik
              fDx(i)=fDx(i)-xzz 
              fDy(i)=fDy(i)-yzz
              fDz(i)=fDz(i)-zzz

              fDx(k)=fDx(k)+xzz
              fDy(k)=fDy(k)+yzz
              fDz(k)=fDz(k)+zzz
                            
!   20   continue                                                          
enddo !su 20
        fx(i)=fDx(i)*aretebim
        fy(i)=fDy(i)*aretebim
        fz(i)=fDz(i)*aretebim
!   10 continue 
enddo !su 10

END SUBROUTINE force_gir
