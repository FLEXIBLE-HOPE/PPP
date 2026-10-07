!*
SUBROUTINE oi_planet_init(PL)
!!
!*
USE const
USE orbit
IMPLICIT NONE

!*
! The arguments
!!----------------------
TYPE(PLANET_INFO) :: PL

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: i
  REAL(RL) :: au

  !*
  ! Start the exectuable code
  !!----------------------------

  PL.nplanet=12
  PL.isc=1
  PL.icb=2
  PL.isun=3
  PL.name(2)='EARTH'
  PL.name(3)='SUN'
  PL.name(4)='MOON'
  PL.name(5)='MERCURY'
  PL.name(6)='VENUS'
  PL.name(7)='MARS'
  PL.name(8)='JUPITER'
  PL.name(9)='SATURN'
  PL.name(10)='URANUS'
  PL.name(11)='NEPTUNE'
  PL.name(12)='PLUTO'

  !! km^3/s^2
  !! IERS 2010
  PL.gm(2) = GME*1.d-9
  PL.gm(3) = GMS*1.d-9
  PL.gm(4) = GME*1.d-9*MEMR
  CALL jpleph_const('AU    ', au)
  au=au*au*au/86400.d0**2
  CALL jpleph_const('GM1   ', PL.gm(5))
  CALL jpleph_const('GM2   ', PL.gm(6))
  CALL jpleph_const('GM4   ', PL.gm(7))
  CALL jpleph_const('GM5   ', PL.gm(8))
  CALL jpleph_const('GM6   ', PL.gm(9))
  CALL jpleph_const('GM7   ', PL.gm(10))
  CALL jpleph_const('GM8   ', PL.gm(11))
  CALL jpleph_const('GM9   ', PL.gm(12))
  PL.gm(5) =PL.gm(5)*au
  PL.gm(6) =PL.gm(6)*au
  PL.gm(7) =PL.gm(7)*au
  PL.gm(8) =PL.gm(8)*au
  PL.gm(9) =PL.gm(9)*au
  PL.gm(10)=PL.gm(10)*au
  PL.gm(11)=PL.gm(11)*au
  PL.gm(12)=PL.gm(12)*au

  ! radius
  PL.radius(2)=E_MAJAXIS*1.d-3
  PL.radius(3)=696000000.d-3
  PL.radius(4)=1.73814d+3
  CALL jpleph_const('RAD1  ', PL.radius(5))
  CALL jpleph_const('RAD2  ', PL.radius(6))
  CALL jpleph_const('RAD4  ', PL.radius(7))

  ! Earth is the center body
  DO i=1, 6
    PL.xj(i,2)=0.d0
    PL.xe(i,2)=0.d0
  END DO

  RETURN

END SUBROUTINE
