!*
SUBROUTINE read_ion_xsy(mjd,sod,CKF,sion,sitename)
!!
!*
USE ckdctrl
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-----------------------
INTEGER(IT) :: mjd
REAL(RL) :: sod
TYPE(CKDCFG) :: CKF
REAL(RL) :: sion(MAXSAT)
CHARACTER(LEN_SITENAME) :: sitename

    !*
    ! The local variables
    !!--------------------------
    LOGICAL(LG) :: lfirst,lexist
    REAL(RL) :: sec
    REAL(RL) :: dt1,dt2,dt1_b,dt2_b,alpha

    INTEGER(IT) :: lfnion,ierr
    INTEGER(IT) :: i,j,k,iy,imon,id,ih,im,isat
    INTEGER(IT) :: mjdf(2), mjdx
    REAL(RL) :: sodf(2), sodx, Dinv
    REAL(RL) :: ion(2,MAXSAT),iontmp(MAXSAT)
    CHARACTER(LEN_STRING) :: line
    CHARACTER(LEN_PRN) :: prn

    DATA lfirst /.TRUE./
    SAVE lfirst,lexist,lfnion,mjdf,sodf,Dinv,ion

    !*
    ! The function called
    !!-----------------------------
    REAL(RL) :: timdif
    INTEGER(IT) :: modified_julday
    INTEGER(IT) :: get_valid_unit
    INTEGER(IT) :: pointer_string

    !*
    ! Start the exectuable code
    !!-----------------------------

    DO k=1, CKF.nprn
      sion(k)=0.d0
    END DO

    IF (lfirst .EQ. .TRUE.) THEN
      lfirst=.FALSE.
      ion = 0.d0
      dt1_b=0.d0
      dt2_b=0.d0
      INQUIRE(FILE=CKF.flnion,EXIST=lexist)
      IF (lexist .EQ. .TRUE.) THEN
        lfnion=get_valid_unit(10)
        OPEN(UNIT=lfnion,FILE=CKF.flnion,STATUS='OLD')
        line=' '
        DO WHILE(INDEX(line,'END OF HEADER') .EQ. 0)
            READ(lfnion,'(A)',END=200,ERR=100) line
            WRITE(*,*)line
            IF (line(61:68).EQ.'INTERVAL') THEN
              READ(line(1:9),*,ERR=200)Dinv
            END IF
        END DO
        !! read one epoch
        READ(lfnion,'(A)',END=200,ERR=100) line
        DO WHILE (INDEX(line,sitename).EQ.0)
          READ(lfnion,'(A)',END=200,ERR=100) line
        END DO
        IF (INDEX(line,'END ') .NE. 0) RETURN
        CALL yr2year(iy)
        mjdf(1)=modified_julday(id,imon,iy)
        sodf(1)=ih*3600.d0+im*60.d0+sec
        CALL timinc(mjdf(1),sodf(1),Dinv,mjdf(2),sodf(2))
        dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
        dt2 = timdif(mjd,sod,mjdf(2),sodf(2))
        
        IF (dt1 .LT. 0.d0) THEN
          WRITE(ERROR_UNIT,'(A)') '%%%WARRNING(read_ion.f90): please check the first epoch in sion file '//TRIM(CKF.flnion)
          lfirst = .TRUE.
          RETURN
        ELSE IF (dt2 .GE. 0.d0) THEN
          DO WHILE (dt2 .GE. 0.d0)
            READ(lfnion,'(A)',END=200,ERR=100) line
            IF (INDEX(line,'END ') .NE. 0) RETURN
            READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec
            CALL yr2year(iy)
            mjdf(1)=modified_julday(id,imon,iy)
            sodf(1)=ih*3600.d0+im*60.d0+sec
            CALL timinc(mjdf(1),sodf(1),Dinv,mjdf(2),sodf(2))
            dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
            dt2 = timdif(mjd,sod,mjdf(2),sodf(2))
          END DO
          DO WHILE (dt2 .LT. 0.d0 .AND. dt1 .GE. 0.d0)
            READ(line(6:8),*,IOSTAT=ierr) prn
            IF (INDEX(line,'END ') .NE. 0) RETURN
            isat = pointer_string(CKF.nprn,CKF.cprn,prn)
            READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec,j,ion(1,isat)
            READ(lfnion,'(A)',END=200,ERR=100) line
            READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec
            IF (INDEX(line,'END ') .NE. 0) RETURN
            CALL yr2year(iy)
            mjdf(1)=modified_julday(id,imon,iy)
            sodf(1)=ih*3600.d0+im*60.d0+sec
            CALL timinc(mjdf(1),sodf(1),Dinv,mjdf(2),sodf(2))
            dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
            dt2 = timdif(mjd,sod,mjdf(2),sodf(2))            
          END DO
        ELSE
          DO WHILE (dt2 .LT. 0.d0 .AND. dt1 .GE. 0.d0)
            READ(line(6:8),*,IOSTAT=ierr) prn
            IF (INDEX(line,'END ') .NE. 0) RETURN
            isat = pointer_string(CKF.nprn,CKF.cprn,prn)
            READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec,j,ion(1,isat)
            READ(lfnion,'(A)',END=200,ERR=100) line
            READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec
            IF (INDEX(line,'END ') .NE. 0) RETURN
            CALL yr2year(iy)
            mjdf(1)=modified_julday(id,imon,iy)
            sodf(1)=ih*3600.d0+im*60.d0+sec
            CALL timinc(mjdf(1),sodf(1),Dinv,mjdf(2),sodf(2))
            dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
            dt2 = timdif(mjd,sod,mjdf(2),sodf(2))            
          END DO          
        END IF
        !! read two epoch
        DO WHILE (dt1 .LT. 0.d0)
          dt1_b=dt1
          READ(line(6:8),*,IOSTAT=ierr) prn
          IF (INDEX(line,'END ') .NE. 0) RETURN
          isat = pointer_string(CKF.nprn,CKF.cprn,prn)
          READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec,j,ion(2,isat)
          READ(lfnion,'(A)',END=200,ERR=100) line
          CALL yr2year(iy)
          mjdf(1)=modified_julday(id,imon,iy)
          sodf(1)=ih*3600.d0+im*60.d0+sec
          CALL timinc(mjdf(1),sodf(1),Dinv,mjdf(2),sodf(2))
          dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
          dt2 = timdif(mjd,sod,mjdf(2),sodf(2))
          IF (dt1.NE.dt1_b) THEN
            DO isat=1,CKF.nprn
              IF (ion(1,isat).NE.0.d0 .AND.ion(2,isat).NE.0.d0) THEN
                alpha = (Dinv+dt1)/Dinv
                sion(isat) = ion(1,isat) + (ion(2,isat)-ion(1,isat))*alpha
              END IF
            END DO
            RETURN
          END IF
        END DO
      ELSE
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ion.f90): the slant ion file is not exist '//TRIM(CKF.flnion)
        CALL exit(1)
      END IF
    END IF

    IF (lexist .EQ. .FALSE.) RETURN

    dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
    dt2 = timdif(mjd,sod,mjdf(2),sodf(2))

    IF (dt1 .LT. 0.d0) THEN
      DO isat=1,CKF.nprn
        IF (ion(1,isat).NE.0.d0 .AND.ion(2,isat).NE.0.d0) THEN
          alpha = (Dinv+dt1)/Dinv
          sion(isat) = ion(1,isat) + (ion(2,isat)-ion(1,isat))*alpha
        END IF
      END DO
      RETURN
    ELSE IF (dt2 .LT. 0.d0 .AND. dt1 .GE. 0.d0) THEN
      ion(1,1:MAXSAT) = ion(2,1:MAXSAT)
      ion(2,1:MAXSAT) = 0.d0
      DO WHILE (dt2 .LT. 0.d0 .AND. dt1 .GE. 0.d0)
        READ(line(6:8),*,IOSTAT=ierr) prn
        IF (INDEX(line,'END ') .NE. 0) RETURN
        isat = pointer_string(CKF.nprn,CKF.cprn,prn)
        READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec,j,ion(2,isat)
        READ(lfnion,'(A)',END=200,ERR=100) line
        READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec
        IF (INDEX(line,'END ') .NE. 0) RETURN
        CALL yr2year(iy)
        mjdf(1)=modified_julday(id,imon,iy)
        sodf(1)=ih*3600.d0+im*60.d0+sec
        CALL timinc(mjdf(1),sodf(1),Dinv,mjdf(2),sodf(2))
        dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
        dt2 = timdif(mjd,sod,mjdf(2),sodf(2))            
      END DO
      DO isat=1,CKF.nprn
        IF (ion(1,isat).NE.0.d0 .AND.ion(2,isat).NE.0.d0) THEN
          alpha = (Dinv+dt1)/Dinv
          sion(isat) = ion(1,isat) + (ion(2,isat)-ion(1,isat))*alpha
        END IF
      END DO
      RETURN
    ELSE
      ion=0.d0
      DO WHILE (dt2 .GE. 0.d0)
        READ(lfnion,'(A)',END=200,ERR=100) line
        IF (INDEX(line,'END ') .NE. 0) RETURN
        READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec
        CALL yr2year(iy)
        mjdf(1)=modified_julday(id,imon,iy)
        sodf(1)=ih*3600.d0+im*60.d0+sec
        CALL timinc(mjdf(1),sodf(1),Dinv,mjdf(2),sodf(2))
        dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
        dt2 = timdif(mjd,sod,mjdf(2),sodf(2))
      END DO
      DO WHILE (dt2 .LT. 0.d0 .AND. dt1 .GE. 0.d0)
        READ(line(6:8),*,IOSTAT=ierr) prn
        IF (INDEX(line,'END ') .NE. 0) RETURN
        isat = pointer_string(CKF.nprn,CKF.cprn,prn)
        READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec,j,ion(1,isat)
        READ(lfnion,'(A)',END=200,ERR=100) line
        READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec
        IF (INDEX(line,'END ') .NE. 0) RETURN
        CALL yr2year(iy)
        mjdf(1)=modified_julday(id,imon,iy)
        sodf(1)=ih*3600.d0+im*60.d0+sec
        CALL timinc(mjdf(1),sodf(1),Dinv,mjdf(2),sodf(2))
        dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
        dt2 = timdif(mjd,sod,mjdf(2),sodf(2))            
      END DO

      DO WHILE (dt1 .LT. 0.d0)
        dt1_b=dt1
        READ(line(6:8),*,IOSTAT=ierr) prn
        IF (INDEX(line,'END ') .NE. 0) RETURN
        isat = pointer_string(CKF.nprn,CKF.cprn,prn)
        READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec,j,ion(2,isat)
        READ(lfnion,'(A)',END=200,ERR=100) line
        READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec
        IF (INDEX(line,'END ') .NE. 0) RETURN
        CALL yr2year(iy)
        mjdf(1)=modified_julday(id,imon,iy)
        sodf(1)=ih*3600.d0+im*60.d0+sec
        CALL timinc(mjdf(1),sodf(1),Dinv,mjdf(2),sodf(2))
        dt1 = timdif(mjd,sod,mjdf(1),sodf(1))
        dt2 = timdif(mjd,sod,mjdf(2),sodf(2))
        IF (dt1.NE.dt1_b) THEN
          DO isat=1,CKF.nprn
            IF (ion(1,isat).NE.0.d0 .AND.ion(2,isat).NE.0.d0) THEN
              alpha = (Dinv+dt1)/Dinv
              sion(isat) = ion(1,isat) + (ion(2,isat)-ion(1,isat))*alpha
            END IF
          END DO
          RETURN
        END IF    
      END DO
    END IF

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ion): read the clock file error '//TRIM(line)
  CALL exit(1)

200 WRITE(OUTPUT_UNIT,'(A)') '%%%WARNING(read_ion): end of the clock file'
  RETURN
END SUBROUTINE