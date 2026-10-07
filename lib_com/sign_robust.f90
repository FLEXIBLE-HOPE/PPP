!*
SUBROUTINE sign_robust(nxl,rxl,flg,drang,ndl)
!!
!! purpose  : 1-dimensional sign-constrained least squares robust estimation
!! parameter:
!!    input : nxl,rxl -- element array
!!    output: flg -- flag of each element
!!            ndl -- # of deleted
!! author   : Geng J
!! created  : Aug 27, 2009
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------------
INTEGER(IT) :: nxl,ndl,flg(1:*)
REAL(RL) :: drang,rxl(1:*)

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i,j,k,nmed,nmd
  REAL(RL) :: rmed(2),rmd(2),resi(nxl),pnew(nxl),p(nxl),cost_new,cost

  !*
  ! Start the exectuable code
  !!--------------------------

  !! iterate until no elements removed
  IF (COUNT(flg(1:nxl).LT.2) .GT. 2) THEN
    cost=1.d12

    !! locate median value(s) of original elements
    CALL med(nxl,rxl,flg,nmed,rmed)
    DO i=1,nmed

      !! absolute residuals
      DO j=1,nxl
        resi(j)=0.d0
        IF (flg(j) .GE. 2) CYCLE
        resi(j)=dabs(rxl(j)-rmed(i))
      END DO

      !! locate median value(s) of absolute residuals
      CALL med(nxl,resi,flg,nmd,rmd)
      DO j=1,nmd

        !! flag elements possibly contaminated & compute cost values
        cost_new=0.d0
        DO k=1,nxl
          pnew(k)=0.d0
          IF (flg(k) .GE. 2) CYCLE
          IF (resi(k).LT.drang .OR. resi(k).LE.3.d0*rmd(j)) THEN
            pnew(k)=1.d0
            cost_new=cost_new+resi(k)**2
          END IF
        END DO
        cost_new=cost_new/COUNT(pnew(1:nxl).EQ.1.d0)

        !! save flags with smallest cost value
        IF (cost_new .LT. cost) THEN
          cost=cost_new
          p(1:nxl)=pnew(1:nxl)
        END IF
      END DO
    END DO

    !! make final decision
    IF (COUNT(p(1:nxl).EQ.0.d0) .GT. COUNT(flg(1:nxl).GE.2)) THEN
      DO i=1,nxl
        IF (flg(i).LT.2 .AND. p(i).EQ.0.d0) THEN
          flg(i)=3
        END IF
      END DO
    END IF
  END IF
  ndl=COUNT(flg(1:nxl) .GE. 2)

  RETURN

END SUBROUTINE


SUBROUTINE med(nxl,rxl,flg,nd,rd)
USE par
IMPLICIT NONE

INTEGER(IT) :: i,j,nxl,nd,flg(1:*),ix
REAL(RL) :: rxl(1:*),rd(1:*),tmp(nxl),tp


  tmp(1:nxl)=rxl(1:nxl)
  DO i=1,nxl-1
    IF (flg(i).GE.2) CYCLE
    DO j=i+1,nxl
      IF (flg(j) .GE. 2) CYCLE
      IF (tmp(i) .GT. tmp(j)) THEN
        tp    =tmp(i)
        tmp(i)=tmp(j)
        tmp(j)=tp
      END IF
    END DO
  END DO

! locate median value(s)
  j=COUNT(flg(1:nxl).LT.2)
  IF (imod(j,2) .EQ. 0) THEN
    nd=2
    rd(1)=tmp(ix(j/2,nxl,flg))
    rd(2)=tmp(ix(j/2+1,nxl,flg))
  ELSE
    nd=1
    rd(1)=tmp(ix(j/2+1,nxl,flg))
  END IF

  RETURN

END SUBROUTINE


INTEGER(IT) FUNCTION ix(ind,nxl,flg)
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

INTEGER(IT) :: ind,nxl,flg(1:*),i

  i=0
  DO ix=1,nxl
    IF (flg(ix) .GE. 2) CYCLE
    i=i+1
    IF (i .EQ. ind) RETURN
  END DO
  WRITE(ERROR_UNIT,'(A)') '***ERROR(sign_robust): value not found'
  CALL exit(1)

END FUNCTION
