!*
SUBROUTINE oi_eph_planet(PL,mjd,sod,rmat,dlmat)
!!
!*
USE const
USE orbit
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: mjd
REAL(RL) :: sod,rmat(3,3),dlmat(3,3)
TYPE(PLANET_INFO) :: PL

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i,j,k
  REAL(RL) :: jd

  LOGICAL(LG) :: lfound

  !*
  ! Start the exectable code
  !!---------------------------

  jd=2400000.5d0+mjd+(sod+OFF_GPS2TT)/86400.d0

  DO i=1, PL.nplanet
    lfound=.TRUE.
    SELECT CASE(TRIM(PL.name(i)))
      CASE('SUN')
        CALL pleph(jd,11,3,PL.xj(1,i))
      CASE('MOON')
        CALL pleph(jd,10,3,PL.xj(1,i))
      CASE('MERCURY')
        CALL pleph(jd,1,3,PL.xj(1,i))
      CASE('VENUS')
        CALL pleph(jd,2,3,PL.xj(1,i))
      CASE('MARS')
        CALL pleph(jd,4,3,PL.xj(1,i))
      CASE('JUPITER')
        CALL pleph(jd,5,3,PL.xj(1,i))
      CASE('SATURN')
        CALL pleph(jd,6,3,PL.xj(1,i))
      CASE('URANUS')
        CALL pleph(jd,7,3,PL.xj(1,i))
      CASE('NEPTUNE')
        CALL pleph(jd,8,3,PL.xj(1,i))
      CASE('PLUTO')
        CALL pleph(jd,9,3,PL.xj(1,i))
    END SELECT

    PL.dist2sc(i)=0.d0
    PL.dist2cb(i)=0.d0
    IF (i .NE. PL.isc) THEN
      DO j=1, 3
        PL.dist2sc(i)=PL.dist2sc(i)+(PL.xj(j,i)-PL.xj(j,PL.isc))**2
      END DO
      PL.dist2sc(i)=DSQRT(PL.dist2sc(i))
    ENDIF
    IF (i .NE. PL.icb) THEN
      DO j=1, 3
        PL.dist2cb(i)=PL.dist2cb(i)+(PL.xj(j,i)-PL.xj(j,PL.icb))**2
      END DO
      PL.dist2cb(i)=DSQRT(PL.dist2cb(i))
    END IF
  END DO

  !! Convert J2000 Position and Velocity to earth-fixed system
  DO i=1, PL.nplanet
    IF (PL.dist2cb(i) .NE. 0.d0) THEN
      DO j=1, 3
        PL.xe(j,  i)=0.d0
        PL.xe(j+3,i)=0.d0
        DO k=1, 3
          PL.xe(j,  i)=PL.xe(j,  i)+rmat(k,j)*PL.xj(k,  i)
          PL.xe(j+3,i)=PL.xe(j+3,i)+rmat(k,j)*PL.xj(k+3,i)+dlmat(k,j)*PL.xj(k,i)
        END DO
      END DO
    END IF
  END DO

  RETURN

END SUBROUTINE
