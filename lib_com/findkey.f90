!*
CHARACTER(LEN=*) FUNCTION findkey(lfn,keyword,sinex_bracket)
!!
!! TO GET THE CONTENT OF THE LINE START WITH 'KEYWORD' WITHIN
!! THE SINEX_BRACKET OR IN THE WHOLE FILE IF THE SINEX_BRACKET
!! IS EMPTY. LINE START WITH '*' OR '#' IS IGNORED AS COMMENT.
!! KEYWORD AND ITS CONTENT IS SEPERATED BY '='. COMMENT ON THE
!! LINE IS SEPERATED BY '!'
!!
!*
USE par
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! THE ARGUMENTS
!!------------------------------
CHARACTER(LEN=*) :: keyword
CHARACTER(LEN=*) :: sinex_bracket
INTEGER(IT) :: lfn

  !*
  ! THE LOCAL VARIABLES
  !!------------------------------
  LOGICAL(LG) :: continous_line
  INTEGER(IT) :: i,j,k,l
  CHARACTER(LEN_STRING) :: line

  !*
  ! THE EXECTUABLE CODES
  !!------------------------------

  REWIND(lfn)
  l=LEN_TRIM(sinex_bracket)
  line=' '
  IF (l .NE. 0) THEN
    DO WHILE(INDEX(line,'+'//sinex_bracket(1:l)) .EQ. 0)
      READ(lfn,'(A)',END=100) line
    END DO
  END IF

  ! THE SYMBLE '++' AT THE END OF THE VALUE AND BEFORE '!' MEANS THAT CONTINOUS LINE FOLLOWS
  findkey='++'
  continous_Line=.FALSE.
  DO WHILE(findkey(LEN_TRIM(findkey)-1:LEN_TRIM(findkey)) .EQ. '++')
    READ(lfn,'(A)',END=100) line

    ! END OF BRACKET
    IF (l.NE.0 .AND. INDEX(line,'-'//sinex_bracket(1:l)).NE.0) GOTO 100

    ! COMMENT LINES
    IF (line(1:1).EQ.'#' .OR. line(1:1).EQ.'*') CYCLE

    ! IF KEYWORD THERE. IF NOT MUST BE CONTINOUS LINE
    j=LEN_TRIM(keyword)
    i=INDEX(line,keyword(1:j))

    ! IF SEPERATOR '=' BETWEEN KEYWORD AND ITS VALUE THERE. IF NOT FROM THE BEGINNING
    j=INDEX(line,'=')
    IF (j .EQ. 0) THEN
      CALL left_justify_string(line)
      line=' '//line
    END IF

    ! COMMENT AT THE END, STARTS WITH '!', IGNORED
    k = INDEX(line,'!')-1
    IF (i.EQ.0 .AND. .NOT.continous_line) CYCLE
    i=LEN(findkey)
    IF (k .LE. j) k=i
    findkey(LEN_TRIM(findkey)-1:)=line(j+2:k)
    continous_line=.TRUE.
  END DO

  RETURN

100 CONTINUE

  findkey='EMPTY'
  WRITE(OUTPUT_UNIT,'(A)') '###WARNING(findkey): '//TRIM(keyword)//' not found.'
  WRITE(OUTPUT_UNIT,'(A)') '            bracket: '//TRIM(sinex_bracket)

  RETURN

END FUNCTION
