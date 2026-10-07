!*
LOGICAL(LG) FUNCTION chos(ndel,imax,idel)
!!
!! purpose  : candidate selection
!! parameter:
!!    input   ndel -- number of candidates to be selected
!!            imax -- maximum candidates that can be selected
!!    output: idel -- index of to-be-selected candidates
!! author   : Jianghui Geng
!! created  : Sep 6 2011
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-------------------------------
INTEGER(IT) :: ndel,imax,idel(1:*)

  !*
  ! The local variables
  !!---------------------
  INTEGER(IT) :: i,ic

  !*
  ! Start the exectuable code
  !!---------------------------

  chos=.FALSE.
  IF (imax .LT. ndel) RETURN
  IF (idel(1) .EQ. 0) THEN

    !! initialize selection array
    DO i=1,ndel
      idel(i)=i
    END DO
    chos=.TRUE.
  ELSE

    !! change selected candidates
    ic=ndel
    DO WHILE(ic .GT. 0)
      idel(ic)=idel(ic)+1
      IF (idel(ic) .GT. imax) THEN
        ic=ic-1
        CYCLE
      END IF
      i=ic+1
      DO WHILE(i .LE. ndel)
        idel(i)=idel(i-1)+1
        IF (idel(i).GT.imax) THEN
          ic=ic-1
          EXIT
        END IF
        i=i+1
      END DO
      IF (i .GT. ndel) EXIT
    END DO
    chos=.TRUE.
    IF (ic.LE.0) chos=.FALSE.
  END IF

  RETURN

END FUNCTION
