SUBROUTINE bigvoi
USE PARACLUSTER !uso il modulo di definizione dei parametri
USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster
USE POTENTIAL
USE ENFORCE

IMPLICIT NONE
INTEGER :: i,j
REAL*8 :: xij,yij,zij,dij2,dc55
!
      dc55 = cutoff_end * cutoff_end
!
nv4(1:natom)=0.d0
!
      DO i=1,natom-1                                                   
        DO j=i+1,natom                                                  
           xij=x(j)+u(j)-x(i)-u(i)                                           
           yij=y(j)+v(j)-y(i)-v(i)                                           
           zij=z(j)+w(j)-z(i)-w(i)           

!------------------------------------------------------------------------------
! PBC for no_cluster system
!------------------------------------------------------------------------------
IF ( wires ) THEN
!!!  condizioni periodiche lungo la direzione del filo e' l'asse z   
          IF (ABS(zij+pbcz)<ABS(zij)) zij=zij+pbcz                         
          IF (ABS(zij-pbcz)<ABS(zij)) zij=zij-pbcz                        
ENDIF
!
IF( surface ) THEN
           if (abs(xij+pbcx)<abs(xij)) xij=xij+pbcx
           if (abs(xij-pbcx)<abs(xij)) xij=xij-pbcx
           if (abs(yij+pbcy)<abs(yij)) yij=yij+pbcy
           if (abs(yij-pbcy)<abs(yij)) yij=yij-pbcy
ENDIF
!
IF ( bulk ) THEN
           if (abs(xij+pbcx)<abs(xij)) xij=xij+pbcx
           if (abs(xij-pbcx)<abs(xij)) xij=xij-pbcx
           if (abs(yij+pbcy)<abs(yij)) yij=yij+pbcy
           if (abs(yij-pbcy)<abs(yij)) yij=yij-pbcy
           IF (ABS(zij+pbcz)<ABS(zij)) zij=zij+pbcz                       
           IF (ABS(zij-pbcz)<ABS(zij)) zij=zij-pbcz                        
ENDIF
!------------------------------------------------------------------------------
!
! Verlet cage system
! dc55 > cutz(itype(i), itype(j))
!
!nv4 == number of neighbours for atom i in a radius dc55
!iv4 == label for identify the number of the neighbour
!
!
           dij2=xij*xij+yij*yij+zij*zij
!
           IF (dij2.LT.dc55) THEN
              nv4(i)=nv4(i)+1 
              iv4(nv4(i),i)=j
           ENDIF
!
        ENDDO
!                                                                     
        IF (nv4(i).GE.nvsiz) THEN                                        
           WRITE (*,'(i4,1x,f6.3,a18)')i,nv4(i),'Too many NN in bigvoi'
           STOP                                                          
        ENDIF                    

      ENDDO                                                          
!
EndSubroutine bigvoi
