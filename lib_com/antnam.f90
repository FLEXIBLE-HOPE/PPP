!*
SUBROUTINE antnam(sitnam,cin,cout,iflag)
!! GET MAIN ANTTYPE
!* 
USE par
USE tables
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The argumtnes
!!-------------------
INTEGER(IT) :: iflag
CHARACTER(LEN_SITENAME) :: sitnam
CHARACTER(LEN=*) :: cin,cout
CHARACTER(LEN_FILENAME) :: antnamfile

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: lfn,ierr
  CHARACTER(LEN_SITENAME) :: type
  CHARACTER(LEN_ANTENNA) :: cmain,c
  LOGICAL(LG) :: lexist

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!----------------------------

  iflag=0 

  antnamfile=f_tableFileName('antnam')
  lexist=.TRUE.
  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=antnamfile)

  DO WHILE(1 .EQ. 1)
    READ (lfn,'(a4,1x,a20)',end=100) type,c
    IF (type .EQ. 'MAIN') cmain=c
    IF (c(1:20) .EQ. cin(1:20)) THEN
      cout=cmain
      CLOSE(lfn)
      RETURN
    ENDIF
  ENDDO

100  WRITE (OUTPUT_UNIT,'(a4,a12,a22)') sitnam,'ANTUNKNOWN',cin(1:20)

  CLOSE (lfn)
  iflag=1

  RETURN

END SUBROUTINE
      
 
