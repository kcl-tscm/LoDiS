subroutine gauss
USE PARACLUSTER
USE CLUSTER

Implicit None
Real :: dlaran
!!Real, Intrinsic :: log
Real :: gset1,gset2,gset3,gset4,r                                          
Real :: g4
Integer :: rip
Integer ::ir(4)
!!Intrinsic :: log

         rip=0
ir(1:4)=irand_seed(1:4)
do 
   gset1=2.d0*dlaran(irand_seed)-1.d0 
   gset2=2.d0*dlaran(irand_seed)-1.d0
   r=gset1*gset1+gset2*gset2 
if (r<1) then     
           g1=gset1*sqrt(-2.d0*log(r)/r)
           g2=gset2*sqrt(-2.d0*log(r)/r)        
           rip=1        
else        
           rip=0         
endif
         if (rip==1) exit        
enddo
       
do          
   gset3=2.d0*dlaran(irand_seed)-1.d0
   gset4=2.d0*dlaran(irand_seed)-1.d0
   r=gset3*gset3+gset4*gset4
if (r<1) then           
           g3=gset3*sqrt(-2.d0*log(r)/r)
           g4=gset4*sqrt(-2.d0*log(r)/r)       
           rip=2        
else        
          rip=1         
endif 
         if (rip==2) exit       
        enddo 
!open(21,file='NUME.out',status='unknown')
!   write(21,'(a5,1x,i5,4f12.5)') 'GAUSS',ipas,g1,g2,g3,g4

End Subroutine GAUSS
