SUBROUTINE voisin
  USE PARACLUSTER  !uso il modulo di definizione dei parametri
  USE CLUSTER   !uso il modulo dove defin1:natomco variabili e parametri cluster
  USE POTENTIAL
  USE ENFORCE

  IMPLICIT NONE
  INTEGER :: i,jv,j
  REAL :: xij,yij,zij,dij2
  REAL :: dc1,dc2,dc3,dc33,dc5
  ! 
  ! distance among first/second/third neighbours in a FCC crystal

   if (sys == 'mon') dc1 = dist(1)
   if( sys == 'bim') dc1 = MAX(dist(1),dist(2)) 
   dc2 = dsqrt(2.0d0) * dc1 !MIN(dist(1),dist(2))
   dc3 = dsqrt(3.0d0) * dc1 !MIN(dist(1),dist(2)) 
   dc5 = dsqrt(2.5d0) * dc1 !MIN(dist(1),dist(2)) !!rac5/rac2 
  !
  dc33 = dc3 * dc3
!  WRITE(*,*) 'IN VOISIN :: nearest distance = ', dc1, '  Third = ', dc3
!  WRITE(*,*) 'IN VOISIN :: comparison dij with  ', dc33
  !
  DO i=1,natom
     nvois(i)=0
  ENDDO
  !
  IF(type_potential=='rgl') THEN
     DO i=1,natom-1  
        DO jv=1,nv4(i) 
!! the loop is restricted to the atoms in the Verlet cage 
           j=iv4(jv,i) 
           xij=x(j)+u(j)-x(i)-u(i) 
           yij=y(j)+v(j)-y(i)-v(i)
           zij=z(j)+w(j)-z(i)-w(i)
           !
           !FOR no_custer systems the periodic boundary condition have to be specified
           !
           IF ( wires ) THEN
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
           !
           dij2=xij*xij+yij*yij+zij*zij
           !
           !NN are defined with respect the cutoff
           !
          !!dc33 = cutz(itype(i),itype(j)) 
          IF (dij2 < dc33 ) THEN
              nvois(i)=nvois(i)+1
              nvois(j)=nvois(j)+1
              ivois(nvois(i),i)=j
              ivois(nvois(j),j)=i
           ENDIF
           !
        ENDDO      !END LOOP ON JV
        IF (nvois(i)>=nvsiz) THEN                                         
           WRITE (uniterr,'(i4,1x,i4,a18)')i,nvois(i),'too much NN in voisin'
           STOP        
        ENDIF
     ENDDO       !END LOOP ON I 
     !
  ELSE   !! for other potential types
!
     DO i=1,natom-1
        DO j=i+1,natom
           xij=x(j)+u(j)-x(i)-u(i)
           yij=y(j)+v(j)-y(i)-v(i)
           zij=z(j)+w(j)-z(i)-w(i)
           dij2=xij*xij+yij*yij+zij*zij
           !
           IF (dij2.LT.dc33) THEN
              nvois(i)=nvois(i)+1
              nvois(j)=nvois(j)+1
              ivois(nvois(i),i)=j
              ivois(nvois(j),j)=i
           ENDIF
        ENDDO      !END LOOP ON J
        IF (nvois(i).GE.nvsiz) THEN    
           WRITE (uniterr,'(i4,1x,i4,a18)')i,nvois(i),'troppi NN in voisin'
           STOP                                                          
        ENDIF

     ENDDO !su i
  ENDIF
  ! 
  !open(33,file='vic.out',status='unknown')
  !do i=1,natom                                         
  ! IF (nvois(i) > 12)  write(33,*) i,nvois(i),ivois(nvois(i),i)
  !enddo
  !
  !
EndSubroutine voisin
