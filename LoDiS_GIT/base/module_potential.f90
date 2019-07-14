MODULE POTENTIAL
USE PARACLUSTER
USE CLUSTER
!
IMPLICIT NONE
!
! Fundamental quantity for each chemical species
!
REAL*8 :: ecoh(3), rat(3), mass(2)
REAL*8 :: t2m(3),arete(3)
REAL*8 :: dist(3)
REAL*8 :: ecohv(n_elposs),ratv(n_elposs),dmas(n_elposs)
CHARACTER*2 ::  leg(n_elposs)

!!!PARAMETERS ONLY FOR RGL
REAL*8 :: cutoff_start,cutoff_end, cutz(2,2)
REAL*8 :: znmax(2),nn(3),znmaxdist

REAL*8 :: a5(3),a4(3),a3(3),x5(3),x4(3),x3(3)
REAL*8 :: p(3),q(3),a(3),qsi(3)
CHARACTER*2 :: NA_elem
CHARACTER*3 :: sys
LOGICAL :: sologeom, lcutoff, lsubstrate

!!!PARAMETERS ONLY FOR LJ
REAL*8 :: U0,R0

!!!PARAMETERS ONLY FOR GIRIFALCO
REAL*8 :: Nball !(=60)atom number in a ball
REAL*8 :: Dnn   !distance NN in A
REAL*8 :: dball !diameter of a ball in A
REAL*8 :: ALJ   ! costante di un LJ in eV*A**6
REAL*8 :: BLJ   ! costante di un LJ in eV*A**12
!!!PARAMETERS ONLY FOR PACHECO-RAMALHO
REAL*8 :: dmu,dM0,epsilon,tau,DIST0
REAL*8 :: C6,C8,C10,C12
!
END MODULE POTENTIAL
