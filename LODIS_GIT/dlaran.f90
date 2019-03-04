      FUNCTION DLARAN(ISEED)
      IMPLICIT NONE
!  DLARAN returns a random real number from a uniform (0,1)
!  distribution.

      INTEGER   ::         ISEED( 4 )
      REAL      ::      DLARAN

!      INTEGER     ::       M1, M2, M3, M4
      integer, PARAMETER   ::  M1 = 494, M2 = 322, M3 = 2508, M4 = 2549
      
 !!     real       ::   ONE = 1.0D+0
REAL :: ONE     
INTEGER, parameter         ::   IPW2 = 4096
      real ::   R
!      PARAMETER        ::  ( IPW2 = 4096, R = ONE / IPW2 )

      INTEGER          ::  IT1, IT2, IT3, IT4
      real ::  DLAR

!!    INTRINSIC        ::  DBLE, MOD

    ONE = 1.0d0  
    R=ONE / IPW2
 
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +&
           ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
!
!     return updated seed
!
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
!
!     convert 48-bit integer to a real number in the interval (0,1)
!
      DLAR = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*&
              ( DBLE( IT4 ) ) ) ) )

   !!   DLARAN=REAL(DLAR,KIND=8)
       dlaran=dlar
      RETURN
!
!     End of DLARAN
!
      END
