!*
SUBROUTINE write_ics_head(flnics,mjd,sod)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!---------------------------
CHARACTER(LEN=*) :: flnics
INTEGER(IT) :: mjd
REAL(RL) :: sod

CHARACTER(LEN=*) :: pname(1:*),system,cprn
INTEGER(IT) :: npwc(1:*),npar
REAL(RL) :: val(1:*),ptime(2,1:*)

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) :: lfn
  INTEGER(IT) :: i,j,k,ipar
  SAVE lfn

  !*
  ! The function called
  !!--------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!--------------------------

  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=flnics)

  WRITE(lfn,'(A)') 'PANDA satellite ICS-file, created by write_ics'
  WRITE(lfn,'(A,I8,F13.3,2X,A)') 'Ref.Time : ',mjd,sod,'GPST'
  WRITE(lfn,'(A)') 'END of FILE'

  RETURN

ENTRY write_ics_sc(system,cprn,npar,pname,npwc,val,ptime)

  BACKSPACE(lfn)
  WRITE(lfn,'(A12,A3)') system,cprn

  ipar=0
  DO i=1, npar
    DO j=1, npwc(i)
      k=INDEX(pname(i),':')-1
      IF (k .EQ. -1) k=LEN_TRIM(pname(i))
      WRITE(lfn,'(A10,3D22.14)') pname(i)(1:k)//'          ',val(ipar+j),ptime(1,ipar+j),ptime(2,ipar+j)
    END DO
    ipar=ipar+npwc(i)
  END DO
  WRITE(lfn,'(a)') 'END of SAT'
  WRITE(lfn,'(a)') 'END of FILE'

  RETURN

ENTRY write_ics_close()
  CLOSE(lfn)

  RETURN

END SUBROUTINE
