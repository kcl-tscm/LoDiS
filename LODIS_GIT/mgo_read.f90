SUBROUTINE mgo_read
 !====================================
 ! The subroutine reads the *.Mgo.pot
 !====================================
 ! aretebim have to be already known

 USE PARACLUSTER
 USE CLUSTER
 USE POTENTIAL
 USE SUBSTRATE

 IMPLICIT NONE
 
 INTEGER :: i, j, l, m_i, m_index
 
 WRITE (*,*) 'mgo_read> Reading file ', TRIM(mgo_pot)
 OPEN(19, FILE=TRIM(mgo_pot), STATUS='old')

 read(19,*)
 read(19,*)
 READ(19,*) mgo_met(1), mgo_met(2)
 READ(19,*)
 READ(19,*)
 read(19,*) (((ccc(1,i,j,l),l=1,3),i=1,3),j=1,3)
 read(19,*)
 read(19,*)
 read(19,*) (((ccc(2,i,j,l),l=1,3),i=1,3),j=1,3)
 read(19,*)
 read(19,*)
 read(19,*) amgo
 read(19,*)
 read(19,*)
 read(19,*) alpha_sub, rc_sub

 CLOSE (19)
 WRITE (*,*) 'mgo_read> MgO parameters read', alpha_sub, rc_sub
 
 !---------------------------------
 ! Check
 ! Are mgo_met(i) the same as elem1, elem2 ?
 IF ((mgo_met(1).NE.elem1).OR.(mgo_met(2).NE.elem2)) THEN
    WRITE (*,*) 'mgo_read> Error: elements 1 and 2 of metal-metal and metal-MgO interactions do not match'
    STOP
 ENDIF

m_index = 1
IF (mgo_met(1) /= mgo_met(2)) m_index = 2
 DO m_i=1,m_index
    WRITE (*,*) 'mgo_read> Parameters for: ', mgo_met(m_i)
    DO j=1,3
       DO i=1,3  
          WRITE(*,'(a,i2,a,i2,a,3f11.6)') 'mgo_read> c(',i,',',j,', 1:3):  ', ccc(m_i,i,j,1), ccc(m_i,i,j,2), ccc(m_i,i,j,3)
       ENDDO
    ENDDO
 ENDDO
 !---------------------------------
 
 !WRITE (*,*) 'mgo_read> Converting parameter rc in aretebim units (what if the system is bimetallic?)'
 WRITE (*,*) 'mgo_read> aretebim is: ', aretebim, arete(1)
 beta_sub_angstrom = alpha_sub/rc_sub  ! It is used to calculate MgO coord number
 rc_sub = rc_sub/arete(1)
 beta_sub = alpha_sub/rc_sub  ! It is used to calculate MgO coord number
 mgo_alat = 2.d0*dsqrt(2.d0)*pi/amgo !! a = sqrt(2) amgo /2  then 2pi/a = 2*sqrt(2) * pi /amgo
 WRITE(*,*) 'mgo_read> a_MgO (O-O distance) [A]: ', amgo
 WRITE(*,*) 'mgo_read> Parameters to calculate CN: alpha, r_c [A], beta_sub = ', alpha_sub, rc_sub, beta_sub
 WRITE (*,*) 'mgo_read> 2*pi/a = sqrt(2)*pi/a_MgO in [1/A]: ', mgo_alat
 
END SUBROUTINE mgo_read
