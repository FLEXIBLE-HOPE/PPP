!*
SUBROUTINE write_saved(rdone,ird)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: ird
CHARACTER(LEN=*) rdone

  !*
  ! The local variables
  !!----------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: lout

  DATA lfirst /.TRUE./
  SAVE lfirst,lout

  !*
  ! The function called
  !!----------------------
  INTEGER(IT) :: get_valid_unit


  IF (lfirst) THEN
    lfirst=.FALSE.
    lout=get_valid_unit(10)
    OPEN(lout,file='rnx_recv')
  END IF
  WRITE(lout,'(A)') rdone(1:ird)

  RETURN

END SUBROUTINE
