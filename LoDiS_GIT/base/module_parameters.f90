MODULE PARACLUSTER

IMPLICIT NONE
!
!Max values for running a MD simulation
!
!==================================================================================
INTEGER, PARAMETER :: nsiz=1200                 !Upper limit for the number of atoms
!==================================================================================
INTEGER, PARAMETER :: nvsiz=100  !!Upper limit for the number of Nearest Neighbour 
INTEGER, PARAMETER :: nomemax=10000000 !!Upper limit for writing out~~.xyz
INTEGER, PARAMETER :: n_elposs=16  !!max number of possible elements 
                                   !!(to be increased if you add a 'new' chemical species)
!
INTEGER, PARAMETER :: unite=54   !file energy
INTEGER, PARAMETER :: unitm=53   !file movie
INTEGER, PARAMETER :: unitd=52   !file depo
INTEGER, PARAMETER :: uniterr=55 !file err
!
! Conversion Factors
!
REAL, PARAMETER :: evsujoule=1.60219d-19
REAL, PARAMETER :: angsum=1.d-10
REAL, PARAMETER :: uasukg=0.16603d-26
!
!Physical Constant (Boltzman, pi)
!
REAL, PARAMETER :: cbol=8.62d-05           !in eV/K
REAL, PARAMETER :: pi=3.141592653589793d0  !pi greco
!
!Square roots
!      
REAL, PARAMETER :: rac2=1.41421356d0                                                
REAL, PARAMETER :: rac3=1.73205081d0
REAL, PARAMETER :: rac5=2.23606798d0                                                  
REAL, PARAMETER :: rac8=2.82842712d0                                                
!
!to work in double precision
!
INTEGER,PARAMETER :: kin=8
!
END MODULE PARACLUSTER
