SUBROUTINE force_rgl
!
  USE PARACLUSTER
  USE CLUSTER 
  USE POTENTIAL
  USE ENFORCE
  USE DISTANCE
!
  IMPLICIT NONE
  INTRINSIC sqrt
  !INTRINSIC size
  REAL(8), ALLOCATABLE :: dispx(:,:),dispy(:,:),dispz(:,:)
  INTEGER :: i,j,k,icolor,ierr
  INTEGER :: nat3,itypik
  REAL(8) :: den(nsiz),frx(nsiz),fry(nsiz),frz(nsiz)
  REAL(8) :: dikm,dikm2,dikm3,dikm4,dikm5
  REAL(8) :: dik0,dik,xik,yik,zik,xuar,yvar,zwar
  REAL(8) :: ebi,eneri,eri,f,for,forsudik,denik
  REAL(8) :: fbx,fby,fbz
  REAL(8) :: espo,pexp,qqexp
  REAL(8) :: psudik0,dueapsudik0,dueq,qsi2,qsi2qsudik0
  REAL(8) :: x1cdm, z2cdm
  
  ALLOCATE(dispx(nvsiz,natom),dispy(nvsiz,natom),dispz(nvsiz,natom))
  
  ! Initialise nsiz element local vectors (Luca)
  den(1:natom)=0.d0                                                       
  frx(1:natom)=0.d0                                                      
  fry(1:natom)=0.d0                                                       
  frz(1:natom)=0.d0  

  ener=0.d0                                                        
  
  dfx(1:natom)=fx(1:natom)
  dfy(1:natom)=fy(1:natom)
  dfz(1:natom)=fz(1:natom)
  
  IF(ipas==1) write(*,*) 'force_rgl> LIMIT ON DISTANCE ', cutoff_start
  !IF(coalescence=='ya') then
  !    x1cdm = (aretebim * fattor / float(natom - natom2)) * &
  !        sum(x(1:natom - natom2) + u(1:natom - natom2))
  !    write(*,*) 'force_rgl>',ipas,x1cdm
  !    z2cdm = (aretebim * fattor / float(natom2)) * &
  !         sum(z(natom - natom2 + 1:natom) + w(natom - natom2 + 1:natom))
  !    write(*,*) 'force_rgl>',ipas,z2cdm
  !ENDIF
  IF(ipas==1) write(*,'(a18,3F10.5)') 'force_rgl> df(1) ', dfx(1),dfy(1),dfz(1)
  IF(ipas==1) write(*,'(a21,3F10.5)') 'force_rgl> df(natom) ', dfx(natom),dfy(natom),dfz(natom)
  
  IF((deposizione =='ya').AND.(natom > natinizio).AND.(nvois(natom)==0)) THEN
     ! The added atom is moving towards the cluster, not near yet
     fx(natom)=0.d0
     fy(natom)=0.d0
     fz(natom)=0.d0       
     nat3=natom-1       
  ELSE
     ! The added atom is moving towards the cluster, it is near now
     nat3=natom
  ENDIF
  !
  DO i=1,nat3                                                     
     ebi=0.d0                                                       
     eri=0.d0                                                          
     eneri=0.d0                                                        
     IF(nvois(i)==0) THEN
        error_counter = error_counter +1
        WRITE(*,*) 'force_rgl> Detected atom with no NN at time-step', ipas
        WRITE(uniterr,*) 'AT ipas',ipas,'ATOM ', i,' TYPE ',itype(i),'HAS NO NN'
        DO ierr=1,nat3
           icolor=1       
           xuar=(x(ierr)+u(ierr))*aretebim
           yvar=(y(ierr)+v(ierr))*aretebim
           zwar=(z(ierr)+w(ierr))*aretebim         
           WRITE(uniterr,'(1x,a4,3f16.5,i4)') elem(ierr),xuar,yvar,zwar,icolor
        ENDDO !su ierr        
        WRITE(uniterr,*) 'force_rgl> atom',i,'has NO NN at time-step ',ipas
        WRITE(uniterr,*) 'FIND THE POSITION in fort.34'
        IF (error_counter .EQ. 20) THEN
           WRITE(*,*) 'force_rgl> Maximum number (20) of errors on NN count has been reached'
           WRITE(*,*) 'force_rgl> Execution stopped!' 
           STOP ! Important to avoid huge error.out
        ENDIF 
     ENDIF

     DO j=1,nvois(i) 
        k=ivois(j,i)  

        IF((itype(i).EQ.1).AND.(itype(k).EQ.1)) THEN          
           itypik=1   ! stesso metallo A          
        ELSE IF((itype(i).EQ.2).AND.(itype(k).EQ.2)) THEN             
           itypik=2   ! stesso metallo B          
        ELSE             
           itypik=3  ! interazione A-B                
        ENDIF
        !       
        dik0=dist(itypik)
        psudik0=p(itypik)/dik0
        dueapsudik0=2.d0*a(itypik)*psudik0
        dueq=2.d0*q(itypik)
        qsi2=qsi(itypik)*qsi(itypik)
        qsi2qsudik0=qsi2*q(itypik)/dik0
        !                                                              
        xik=x(k)+u(k)-x(i)-u(i)                                           
        yik=y(k)+v(k)-y(i)-v(i)                                           
        zik=z(k)+w(k)-z(i)-w(i)                                           
        !
!
!! For  no_clusters systems the PBS are imposed
!
        IF ( wires ) THEN
           IF (ABS(zik+pbcz)<ABS(zik)) zik=zik+pbcz
           IF (ABS(zik-pbcz)<ABS(zik)) zik=zik-pbcz
        ENDIF
        !
        IF ( surface ) THEN
           IF (ABS(xik+pbcx)<ABS(xik)) xik=xik+pbcx
           IF (ABS(xik-pbcx)<ABS(xik)) xik=xik-pbcx
           IF (ABS(yik+pbcy)<ABS(yik)) yik=yik+pbcy
           IF (ABS(yik-pbcy)<ABS(yik)) yik=yik-pbcy
        ENDIF
        !
        IF ( bulk ) THEN
           IF (ABS(xik+pbcx)<ABS(xik)) xik=xik+pbcx
           IF (ABS(xik-pbcx)<ABS(xik)) xik=xik-pbcx
           IF (ABS(yik+pbcy)<ABS(yik)) yik=yik+pbcy
           IF (ABS(yik-pbcy)<ABS(yik)) yik=yik-pbcy
           IF (ABS(zik+pbcz)<ABS(zik)) zik=zik+pbcz
           IF (ABS(zik-pbcz)<ABS(zik)) zik=zik-pbcz
        ENDIF
        !
        dispx(j,i)=xik
        dispy(j,i)=yik
        dispz(j,i)=zik
        !
        dik=SQRT(xik*xik+yik*yik+zik*zik)         
        !!if(pair_distance(i,k) /= dik) write(*,*) 'In force', ipas, i,k, 'distance pb', dik
        !
        !
        !IF (dik.LT.dc2) THEN
!
!!condition for distance shorter than cutoff_start limit
!                                              
        IF(dik < cutoff_start) THEN
           espo=dik/dik0-1.d0                                
           pexp=EXP(-p(itypik)*espo)
           qqexp=EXP(-dueq*espo)
           den(i)=qsi2*qqexp+den(i)
           for=pexp*dueapsudik0
           eri=eri+a(itypik)*pexp
        ELSE
!
! analitical extention of the potential up to cutoff_end
!
           dikm=dik-cutoff_end 
           !
           dikm2=dikm*dikm
           dikm3=dikm2*dikm
           dikm4=dikm3*dikm
           dikm5=dikm4*dikm  
           den(i)=den(i)+(x5(itypik)*dikm5+x4(itypik)*dikm4+ &
                & x3(itypik)*dikm3)**2
           for=-2.d0*(5.d0*a5(itypik)*dikm4+4.d0*a4(itypik)*dikm3+ &
                & 3.d0*a3(itypik)*dikm2)
           eri=eri+a5(itypik)*dikm5+a4(itypik)*dikm4+a3(itypik)*dikm3       
        ENDIF
        !
        forsudik=for/dik
        frx(i)=frx(i)-forsudik*xik                                   
        fry(i)=fry(i)-forsudik*yik                                    
        frz(i)=frz(i)-forsudik*zik  
        !                       
     ENDDO !su j
     !
     ebi=-SQRT(den(i))  
     !
     den(i)=-1.d0/ebi 
     !
     eneri=ebi+eri
     !
     ener=ener+eneri 
     !
     !MORE OUTPUT
     !
     !!potener(i)=eneri
     !!press(i)=(-q(itypik) * eri - p(itypik) * ebi)/(dik0)*dik
     !
  ENDDO !end loop on atom i
!  
  DO i=1,natom                                                  
     fbx=0.d0                                                          
     fby=0.d0                                                          
     fbz=0.d0 
     IF(nvois(i).NE.0) THEN 
        DO j=1,nvois(i) 
           k=ivois(j,i) 
           !
           IF((itype(i).EQ.1).AND.(itype(k).EQ.1)) THEN
              itypik=1   ! A-A 
           ELSE IF((itype(i).EQ.2).AND.(itype(k).EQ.2)) THEN
              itypik=2   ! B-B
           ELSE
              itypik=3   ! A-B
           ENDIF
           !                      
           dik0=dist(itypik)
           psudik0=p(itypik)/dik0
           dueapsudik0=2.d0*a(itypik)*psudik0
           dueq=2.d0*q(itypik)
           qsi2=qsi(itypik)*qsi(itypik)
           qsi2qsudik0=qsi2*q(itypik)/dik0 
           xik=dispx(j,i)
           yik=dispy(j,i)
           zik=dispz(j,i)
           dik=SQRT(xik*xik+yik*yik+zik*zik)
           !
           !          IF (dik.LT.dc2) THEN                                              
           IF(dik < cutoff_start) THEN
              espo=dik/dik0-1.d0                                
              qqexp=EXP(-dueq*espo)
              f=qsi2qsudik0*qqexp
           ELSE                                                              
              dikm=dik-cutoff_end
              dikm2=dikm*dikm
              dikm3=dikm2*dikm
              dikm4=dikm3*dikm
              dikm5=dikm4*dikm
              f=-(x5(itypik)*dikm5+x4(itypik)*dikm4+&
                   x3(itypik)*dikm3)*(5.d0*x5(itypik)*dikm4+&
                   4.d0*x4(itypik)*dikm3+3.d0*x3(itypik)*dikm2)
           ENDIF
           denik=f*(den(i)+den(k))/dik
           fbx=fbx+denik*xik                        
           fby=fby+denik*yik                          
           fbz=fbz+denik*zik                                 
        ENDDO !endloop on atom j neighbour of atom i
     ENDIF !!if atom i has neighbours
     !
     !comment on units
     !RGL force is calculated as eV/[units of dik0 which is fixed to 1/sqrt(2) thus in the so-called arete units
     !in time.f90 all forces shoudl be in the same units
     !
     fx(i)=frx(i)+fbx                                                  
     fy(i)=fry(i)+fby                                                  
     fz(i)=frz(i)+fbz 
 
  ENDDO !on atom i      

  DEALLOCATE(dispx,dispy,dispz)

END SUBROUTINE force_rgl
