SUBROUTINE force_lj1
USE PARACLUSTER
USE CLUSTER
USE POTENTIAL
USE ENFORCE
USE DISTANCE
!
Implicit None
Intrinsic sqrt
Integer :: i,k!,j,icolor,ierr
Integer :: nat3d!,itypik
Real  :: fDx(nsiz),fDy(nsiz),fDz(nsiz)
Real :: rik,rik2,rik3,rik6,rik12
Real :: dik,xik,yik,zik!dik0,xuar,yvar,zwar
Real :: wtot,wtote
Real :: eneri!,eri!,f,for,forsudik,denik
!!Real :: fbx,fby,fbz

!aretebim=arete(1)      
! pay attention to the units: in fact,
! we have all the routine which calculates the velocities and the positions
! in eV (for energy!) and arete (for position!)
! Initially considers energy in eV BUT position in Angstrom

      ener=0.d0  
      eneri=0.d0  
!!      dik0=dc1  !! distance NN in arete 
!do i=1,natom-1
! do k=i+1,natom
!          if((itype(i).eq.1).and.(itype(k).eq.1)) then
!             itypik=1   ! same A
!          else if((itype(i).eq.2).and.(itype(k).eq.2)) then
!             itypik=2   ! same B
!          else
!             itypik=3  ! alloy A-B
!          endif
!        dik0=dist(itypik)
!  enddo
!enddo

! old values for force dfx,dfy,dfz
nat3d=natom
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
     
 do i=1,nat3d
   do  k=1+i,nat3d      
          
          xik=x(k)+u(k)-x(i)-u(i) 
          yik=y(k)+v(k)-y(i)-v(i) 
          zik=z(k)+w(k)-z(i)-w(i)

          dik=sqrt(xik*xik+yik*yik+zik*zik)   !Arete units
          rik=aretebim*dik         !in A
          rik=rik/R0
          
!          if((ipas.eq.1).or.(ipas.eq.3))write(*,*) rik,'rik',i,k     
          
          rik2=rik*rik
          rik3=rik2*rik
          rik6=rik3*rik3
          rik12=rik6*rik6
           
 !!calculus the forces in 1/arete
 !!then i convert again in eV and arete
 !! f(i)= f(i) in eV/arete  !!
               
               wtot=(1.d0/rik6-1.d0/rik12)  !in A
               wtote=(1.d0/rik12-2.d0/rik6) !in A
            
    ! this is the potential
         ener=ener+wtote
         
    ! this is the force  in 1/arete

                           !e' il meno della derivata 
              fDx(i)=fDx(i)+aretebim*xik*wtot/rik2  
              fDy(i)=fDy(i)+aretebim*yik*wtot/rik2
              fDz(i)=fDz(i)+aretebim*zik*wtot/rik2

              fDx(k)=fDx(k)-aretebim*xik*wtot/rik2
              fDy(k)=fDy(k)-aretebim*yik*wtot/rik2
              fDz(k)=fDz(k)-aretebim*zik*wtot/rik2
                            
        enddo

        fx(i)=fDx(i)*12.d0*U0*aretebim    !in ev/arete
        fy(i)=fDy(i)*12.d0*U0*aretebim
        fz(i)=fDz(i)*12.d0*U0*aretebim          

!        fx(k)=fDx(k)*12.d0*U0*aretebim    !in ev/arete
!        fy(k)=fDy(k)*12.d0*U0*aretebim
!        fz(k)=fDz(k)*12.d0*U0*aretebim          
       
        
enddo
   
          !ener=0.25*ener*U0   ! in eV
          ener=ener*U0   ! in eV
        

END SUBROUTINE force_lj1
