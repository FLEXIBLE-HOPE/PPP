!*
SUBROUTINE hf_eop(mjd,utc2tt,eop_out)
!!
!*
USE tables
USE hfeop_xyu 
IMPLICIT NONE

!*
! The arguments
!!--------------------
REAL(RL) :: mjd,utc2tt
REAL(RL) :: eop_out(3)

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: n,m,i
  INTEGER(IT) :: lfn

  CHARACTER(LEN_STRING) :: line
  REAL(RL) :: arg,argd
  REAL(RL) :: fund_arg(6,2)
  ! If true, then GMST+pi is first arguement. if false, then last.
  REAL(RL) :: eop(4,2)
  LOGICAL(LG) :: knew_order
  DATA knew_order /.true./

  LOGICAL(LG) :: lexist
  LOGICAL(LG) :: lfirst
  DATA lfirst /.TRUE./
  SAVE lfirst

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!---------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.

    line=f_tableFileName('hfeop')
    INQUIRE(FILE=line,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(hf_eop): '//TRIM(line)//' is not exist.'
      CALL exit(1)
    END IF

    lfn=get_valid_unit(10)
    CALL import_tides_xyu(TRIM(line),lfn)

  END IF

  ! Get the fundamental arguments at this epoch
  CALL tide_angles(mjd,utc2tt,knew_order,fund_arg)
   
  eop=0.d0
  ! now loop over the tidal contribution.  
  DO m=1, num_tide
    ! Get the argument and the time_derivative. 
    arg=dotarg(itide_arg(1,m),fund_arg(1,1))
    ! arg = mod(arg, 2.d0*pi)
    argd=dotarg(itide_arg(1,m),fund_arg(1,2))      
    DO n=1,4
      eop(n,1)=eop(n,1)+DSIN(arg)*tide_coef(1,n,m)+DCOS(arg)*tide_coef(2,n,m)
      eop(n,2)=eop(n,2)+argd*(DCOS(arg)*tide_coef(1,n,m)-DSIN(arg)*tide_coef(2,n,m))
    END DO
  END DO

  eop_out(1)=eop(1,1)
  eop_out(2)=eop(2,1)
  eop_out(3)=eop(3,1)

  RETURN

END SUBROUTINE
