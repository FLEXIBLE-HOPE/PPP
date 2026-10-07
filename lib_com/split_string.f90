!*
SUBROUTINE split_string(lnoempty,string,c_start,c_end,seperator,nword,word)
!!
!! TO SPLIT THE STRING INTO SEPERATED WORDS
!!
!! GE MAORONG: CREATED
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! THE ARGUMENTS
!!-----------------------------
LOGICAL(LG) :: lnoempty
CHARACTER(LEN=*) :: string,word(1:*)
CHARACTER :: c_start,c_end,seperator
INTEGER(IT) :: nword

  !*
  ! THE LOCAL VARIABLES
  !!-----------------------------
  INTEGER(IT) :: i0,i1,ilast,i
  CHARACTER(LEN_STRING) :: line

  !*
  ! THE EXECTUABLE CODES
  !!-----------------------------

  line=string
  IF (LEN(string) .GT. LEN_STRING) THEN
    WRITE(ERROR_UNIT,'(A,I5)') '***ERROR(split_line): input string length > ', LEN_STRING
    CALL exit(1)
  END IF

  i0=0
  i1=LEN_TRIM(line)+1
  IF (LEN_TRIM(c_start) .NE. 0) i0=INDEX(line,c_start)
  IF (LEN_TRIM(c_end  ) .NE. 0) i1=INDEX(line,c_end  )

  nword=0
  IF (i1.EQ.0 .OR. i0.GT.i1) RETURN

  ilast=i0+1
  line(i1:i1)=seperator
  DO i=i0+1,i1
    IF (line(i:i) .EQ. seperator) THEN
      IF (LEN_TRIM(line(ilast:i-1)) .NE. 0) THEN
        nword=nword+1
        word(nword)=line(ilast:i-1)
      ELSE
        IF (.NOT.lnoempty) THEN
          nword=nword+1
          word(nword)=' '
        END IF
      END IF
      ilast=i+1
    END IF
  END DO

  RETURN

END SUBROUTINE
