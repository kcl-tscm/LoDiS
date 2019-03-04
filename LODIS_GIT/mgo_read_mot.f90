SUBROUTINE mgo_read_mot
 !====================================
 ! The subroutine reads the *.Mgo.pot
 !====================================
 ! arete(1) have to be already known

 USE paracluster
 USE cluster     !uso il modulo dove definisco variabili e parametri cluster
 USE potential
 USE substrate

 IMPLICIT NONE
 
 INTEGER :: i, j, l, m_i
 

 WRITE (*,*) 'mgo_read_mot> Reading file ', TRIM(mgo_pot)
 OPEN(19, FILE=TRIM(mgo_pot), STATUS='old')

 read(19,*)
 read(19,*)
 READ(19,*) mgo_met(1), mgo_met(2)
 READ(19,*)
 READ(19,*)
 read(19,*) (((ccc(1,i,j,l),l=1,3),i=1,3),j=1,3)
 read(19,*) ((ccc4(1,j,l),l=1,3),j=1,2)
 read(19,*)
 read(19,*)
 read(19,*) (((ccc(2,i,j,l),l=1,3),i=1,3),j=1,3)
 read(19,*) ((ccc4(2,j,l),l=1,3),j=1,2)
 read(19,*)
 read(19,*)
 read(19,*) amgo
 read(19,*)
 read(19,*)
 read(19,*) alpha_sub, rc_sub

 CLOSE (19)
 WRITE (*,*) 'mgo_read_mot> MgO parameters read'
 
 !---------------------------------
 ! Check
 ! Are mgo_met(i) the same as elem1, elem2 ?
 IF ((mgo_met(1).NE.elem1).OR.(mgo_met(2).NE.elem2)) THEN
    WRITE (*,*) 'mgo_read_mot> Error: elements 1 and 2 of metal-metal and metal-MgO interactions do not match'
    STOP
 ENDIF
 DO m_i=1,2
    WRITE (*,*) 'mgo_read_mot> Parameters for: ', mgo_met(m_i)
    DO j=1,3
       DO i=1,3  
          WRITE(*,'(a,i2,a,i2,a,3f11.6)') 'mgo_read_mot> c(',i,',',j,', 1:3):  ', ccc(m_i,i,j,1), ccc(m_i,i,j,2), ccc(m_i,i,j,3)
       ENDDO
    ENDDO
    DO j=1,2
       i=4 
       WRITE(*,'(a,i2,a,i2,a,3f11.6)') 'mgo_read_mot> c(',i,',',j,', 1:3):  ', ccc4(m_i,j,1), ccc4(m_i,j,2), ccc4(m_i,j,3)
    ENDDO
 ENDDO
 WRITE(*,*) 'mgo_read_mot> a_MgO (O-O distance) [A]: ', amgo
 WRITE(*,*) 'mgo_read_mot> Parameters to calculate CN: alpha, r_c [A] = ', alpha_sub, rc_sub
 !---------------------------------
 
 !WRITE (*,*) 'mgo_read> Converting parameter rc in arete(1) units (what if the system is bimetallic?)'
 WRITE (*,*) 'mgo_read_mot> arete(1) is: ', arete(1)
 beta_sub_angstrom = alpha_sub/rc_sub  ! It is used to calculate MgO coord number
 rc_sub = rc_sub/arete(1)
 beta_sub = alpha_sub/rc_sub  ! It is used to calculate MgO coord number
 mgo_alat = 2.d0*dsqrt(2.d0)*pi/amgo
 WRITE (*,*) 'mgo_read_mot> 2*pi/(a_MgO/sqrt(2)) is [A]: ', mgo_alat

END SUBROUTINE mgo_read_mot
