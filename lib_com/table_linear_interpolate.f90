!*
SUBROUTINE table_linear_interpolate(tname,check_only,fjd,x)
!*
USE par
USE tables
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------------------
INTEGER(IT) :: MAXVARTAB
PARAMETER(MAXVARTAB=3)
TYPE linear_table
  LOGICAL(LG) :: first
  CHARACTER(LEN_FILENAME) :: file
  CHARACTER(LEN_STRING) :: fmt
  INTEGER(IT) :: lunit
  INTEGER(IT) :: jds,jde,nvar,npnt_inline
  REAL(RL) :: dintv,unit,jdsi
  REAL(RL) :: table(MAXPNT,MAXVARTAB)
END TYPE

CHARACTER(LEN=*) :: tname
LOGICAL(LG) :: check_only
REAL(RL) :: fjd,x(1:*)
TYPE(linear_table) :: EI(5)

  !*
  ! The local variable
  !!-----------------------------
  LOGICAL(LG) :: lfirst,found
  INTEGER(IT) :: itemp,ierr,i,j,idir
  REAL(RL) :: alpha,dummy
  CHARACTER(LEN_STRING) :: line
  DATA lfirst /.TRUE./
  SAVE EI, lfirst

  !*
  ! The function called
  !!------------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start of executable code
  !!------------------------------

  !
  !! table index
  IF (tname(1:3) .EQ. 'nut') THEN
    itemp=1
  ELSE IF(tname(1:7) .EQ. 'poleut1') THEN
    itemp=2
  ELSE IF(tname(1:9) .EQ. 'geomag_kp') THEN
    itemp=3
  ELSE IF(tname(1:10) .EQ. 'solar_flux') THEN
    itemp=4
  ELSE
    WRITE(ERROR_UNIT,'(2A)') '***ERROR(table_linear_interpolate): unknwon table ',tname
    CALL exit(1)
  END IF

  !! first time call
  IF (lfirst) THEN
    DO i=1,5
      EI(i).first=.TRUE.
    END DO
    EI(1).file=f_tableFileName("nutabl")
    EI(2).file=f_tableFileName("polut1")
    EI(3).file=f_tableFileName("geomag_kp")
    EI(4).file=f_tableFileName("solar_flux")

    lfirst=.FALSE.
  END IF

  IF(EI(itemp).first .EQ. .TRUE.) THEN

    EI(itemp).lunit=get_valid_unit(10)
    OPEN(UNIT=EI(itemp).lunit,FILE=EI(itemp).file)

    READ(EI(itemp).lunit,*) EI(itemp).jds,EI(itemp).jde,EI(itemp).dintv,EI(itemp).nvar,EI(itemp).npnt_inline,EI(itemp).unit

    IF (EI(itemp).nvar .GT. MAXVARTAB) THEN
      WRITE(OUTPUT_UNIT,'(A,I5,A,I5)') '%%%MESSAGE(table_linear_interpolate): the number of variables in the table ', &
                                        EI(itemp).nvar, ' is larger than maxmium', MAXVARTAB
    END IF
    IF (EI(itemp).npnt_inline .GT. MAXPNT/2) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(table_linear_interpolate): there are too many points in each line ',&
      EI(itemp).npnt_inline, MAXPNT
    END IF
        
    READ(EI(itemp).lunit,'(A)') EI(itemp).fmt
    EI(itemp).jdsi=0.d0
    EI(itemp).first=.FALSE.
  END IF
  IF (check_only .EQ. .TRUE.) RETURN

  !! time falls in the table span
  IF (fjd.LT.EI(itemp).jds .OR. fjd.GT.EI(itemp).jde) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(table_linear_interpolate): input beyond table time'
    CALL exit(1)
  END IF

  !! check IF new DATA required
  idir=1
  ierr = 0
  DO WHILE(idir .EQ. 1)
    IF(fjd.LE.EI(itemp).jdsi .OR. EI(itemp).jdsi.EQ.0.d0) THEN
      IF(EI(itemp).jdsi .NE. 0.d0) THEN
        BACKSPACE EI(itemp).lunit
        BACKSPACE EI(itemp).lunit
        BACKSPACE EI(itemp).lunit
      END IF
      READ(EI(itemp).lunit,EI(itemp).fmt,IOSTAT=ierr) EI(itemp).jdsi,((EI(itemp).table(i,j) &
           ,j=1,EI(itemp).nvar),i=1,EI(itemp).npnt_inline)
      IF(ierr .EQ. 0) READ(EI(itemp).lunit,EI(itemp).fmt,IOSTAT=ierr) dummy,((EI(itemp).table(i,j) &
           ,j=1,EI(itemp).nvar),i=EI(itemp).npnt_inline+1,2*EI(itemp).npnt_inline)
    ELSE IF(fjd .GT. EI(itemp).jdsi+EI(itemp).dintv*(EI(itemp).npnt_inline*2-1)) THEN
      EI(itemp).jdsi=EI(itemp).jdsi+EI(itemp).npnt_inline*EI(itemp).dintv
      DO i=1,EI(itemp).npnt_inline
        DO j=1,EI(itemp).nvar
          EI(itemp).table(i,j)=EI(itemp).table(i+EI(itemp).npnt_inline,j)
        END DO
      END DO
      READ(EI(itemp).lunit,EI(itemp).fmt,IOSTAT=ierr) dummy,((EI(itemp).table(i,j), &
          j=1,EI(itemp).nvar),i=EI(itemp).npnt_inline+1,2*EI(itemp).npnt_inline)
    ELSE
      idir=0
    END IF

    IF(ierr.NE.0) THEN
      IF (EI(itemp).jdsi.LT.fjd .AND. EI(itemp).jdsi.GT.(fjd-2.5)) GOTO 100
      WRITE(ERROR_UNIT,'(3A)') '***ERROR(table_linear_interpolate): READ file ',EI(itemp).file, EI(itemp).fmt
      WRITE(ERROR_UNIT,*) fjd,EI(itemp).jdsi
      CALL exit(1)
    END IF
  END DO


  found=.FALSE.
  DO i=1,EI(itemp).npnt_inline*2-1
    IF(fjd.GE.EI(itemp).jdsi+(i-1)*EI(itemp).dintv .AND. fjd.LE.EI(itemp).jdsi+i*EI(itemp).dintv) THEN
      found=.TRUE.
      alpha=(fjd-(EI(itemp).jdsi+(i-1)*EI(itemp).dintv))/EI(itemp).dintv
      DO j=1,EI(itemp).nvar
        x(j)=EI(itemp).table(i,j)+alpha*(EI(itemp).table(i+1,j)-EI(itemp).table(i,j))
        x(j)=x(j)*EI(itemp).unit
      END DO
      EXIT
    END IF
  END DO

100 CONTINUE
  IF (ierr.NE.0 .AND. EI(itemp).jdsi.LT.fjd .AND. EI(itemp).jdsi.GT.(fjd-2.5*EI(itemp).dintv)) THEN
    found=.FALSE.
    BACKSPACE EI(itemp).lunit
    BACKSPACE EI(itemp).lunit
    BACKSPACE EI(itemp).lunit
    READ(EI(itemp).lunit,EI(itemp).fmt,iostat=ierr) EI(itemp).jdsi,((EI(itemp).table(i,j) &
        ,j=1,EI(itemp).nvar),i=1,EI(itemp).npnt_inline)
    READ(EI(itemp).lunit,EI(itemp).fmt,iostat=ierr) dummy,((EI(itemp).table(i,j)  &
        ,j=1,EI(itemp).nvar),i=EI(itemp).npnt_inline+1,2*EI(itemp).npnt_inline)
    DO i=1,EI(itemp).npnt_inline*2-1
      IF (dummy .GT. EI(itemp).jdsi+EI(itemp).dintv) EXIT
      IF(fjd.GE.EI(itemp).jdsi+(i-1)*EI(itemp).dintv .AND. fjd.LE.EI(itemp).jdsi+(i+2.5)*EI(itemp).dintv) THEN
        found=.TRUE.
        alpha=(fjd-(EI(itemp).jdsi+(i-1)*EI(itemp).dintv))/EI(itemp).dintv
        DO j=1,EI(itemp).nvar
          x(j)=EI(itemp).table(i,j)+alpha*(EI(itemp).table(i+1,j)-EI(itemp).table(i,j))
          x(j)=x(j)*EI(itemp).unit
        END DO
        EXIT
      END IF
    END DO
  END IF

  !! not found
  IF(found .EQ. .FALSE.) THEN
    WRITE(ERROR_UNIT,'(A,1X,F12.5)') '***ERROR(table_linear_interpolate): no data for '//TRIM(tname),fjd
    CALL exit(1)
  END IF

  RETURN

END SUBROUTINE

