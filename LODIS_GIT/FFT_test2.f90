SUBROUTINE DFT (FR,FI,GR,GI,N)
!
! Subroutine to perform the discrete Fourier transform with
! FR and FI as the real and imaginary parts of the input and
! GR and GI as the corresponding  output.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,i
  REAL :: PI,X,Q
  REAL, INTENT (IN), DIMENSION (N) :: FR,FI
  REAL, INTENT (OUT), DIMENSION (N) :: GR,GI
N = 4095
open(67,file = 'fort.61' , status = 'old' , action = 'read')
do i = 1,N
	read(67,*) FR(i)
	read(67,*) FI(i)
enddo

!
  PI = 4.0*ATAN(1.0)
  X  = 2*PI/N
!
  DO I = 1, N
    GR(I) = 0.0
    GI(I) = 0.0
    DO J = 1, N
      Q = X*(J-1)*(I-1)
      GR(I) = GR(I)+FR(J)*COS(Q)+FI(J)*SIN(Q)
      GI(I) = GI(I)+FI(J)*COS(Q)-FR(J)*SIN(Q)
    END DO
  END DO
do I=1,N
write(44,*) I, GR(I)*GR(I) + GI(I)*GI(I) 
enddo
END SUBROUTINE DFT
