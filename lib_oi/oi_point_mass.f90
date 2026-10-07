!*
SUBROUTINE oi_point_mass(lpart,PL,force_model,acc0,amat)
!!
!*
USE const
USE orbit
IMPLICIT NONE

!*
! The arguments
!!----------------------
LOGICAL(LG) :: lpart
CHARACTER(LEN=*) :: force_model
REAL(RL) :: acc0(1:*),amat(3,3)
TYPE(PLANET_INFO) :: PL

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) :: i,j,k
  REAL(RL) :: fac1,fac2,rat,acc(3)

  !*
  ! Start the exectable code
  !!--------------------------

  acc=0.d0

  DO i=1, PL.nplanet
    IF (i.NE.PL.icb .AND. INDEX(force_model,PL.name(i)(1:3)).EQ.0) CYCLE
    fac1=-PL.gm(i)/PL.dist2sc(i)**3
    DO j=1, 3
      acc(j)=acc(j)+(PL.xj(j,PL.isc)-PL.xj(j,i))*fac1
    END DO

    IF (i .NE. PL.icb) THEN
      fac2=-PL.gm(i)/PL.dist2cb(i)**3
      DO j=1, 3
        acc(j)=acc(j)+(PL.xj(j,i)-PL.xj(j,PL.icb))*fac2
      END DO
    END IF
    IF (.NOT. lpart) CYCLE

    ! Partial
    fac2=-3.d0/PL.dist2sc(i)**2
    DO j=1, 3
      DO k=1, j
        rat=fac2*(PL.xj(j,PL.isc)-PL.xj(j,i))*(PL.xj(k,PL.isc)-PL.xj(k,i))
        IF (j .EQ. k) rat=rat+1.d0
        rat=rat*fac1
        amat(j,k)=amat(j,k)+rat
        IF (j .NE. k) amat(k,j)=amat(k,j)+rat
      END DO
    END DO
  END DO

  DO i=1, 3
    acc0(i)=acc0(i)+acc(i)
  END DO

  RETURN

END SUBROUTINE
