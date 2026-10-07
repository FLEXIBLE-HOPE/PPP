!*
SUBROUTINE rdorbh(flnorb,lfn,OH)
!!
!! purpose   : read header of PANDA orbit file
!! parameters: orbfil -- orbit file name
!!             iunit  -- unit of file
!!             OH     -- header part 1 (oi control data)
!! created   : Geng J
!!
!*
USE orbit
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: lfn
CHARACTER(LEN=*) :: flnorb
TYPE(ORBHDR) :: OH

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: ierr

  !*
  ! The function called
  !!------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!--------------------------

  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=flnorb,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=ierr)
  IF (ierr .NE. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(rdorbh): open file '//TRIM(flnorb)
    CALL exit(1)
  END IF

  !! read header of the orbit file
  READ(lfn,IOSTAT=ierr) OH
  IF (ierr .NE. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(rdorbh): read orbit '//TRIM(flnorb)
    CALL exit(1)
  END IF

  RETURN

END SUBROUTINE
