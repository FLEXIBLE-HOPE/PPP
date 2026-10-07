!*
SUBROUTINE read_osb_upd(CKF,SAT)
!!
!! XU SHENGYI: CREATED [2023-05-26]  
!! PURPOSE: READ OSB FILE IN THE MODE OF UPD [EPOCH, UPD => OSB]
!*
USE const
USE ckdctrl
USE satellite
USE iso_fortran_env
IMPLICIT NONE

!*
! the argument input
!!------------------------
TYPE(CKDCFG) :: CKF
TYPE(SATE)   :: SAT(MAXSAT)

    !*
    ! the local variables
    !!------------------------
    LOGICAL(LG) :: lfirst, lexist

    INTEGER(IT) :: lfnosb,ierr
    INTEGER(IT) :: i,j,lfn,isat,isys,ind ,nosbType(MAXSYS),ifreq,indp
    INTEGER(IT) :: year,doy,sod,imon,iday,ihour,imin,cfg_year,cfg_doy,mjdf(2),mjdx,iy,im,id,ih
    REAL(RL) :: sodf(2),sodx,dt1,dt2,dt1_b,dt2_b, osbDelta,isec,osb,sec
    CHARACTER(LEN_STRING) :: line,osbTime, satSys,osbType 
    CHARACTER(LEN_PRN) :: cprn,freqtp

    DATA lfirst /.TRUE./, lexist /.TRUE./
    SAVE lfirst, lexist,lfnosb, nosbType, osbDelta,mjdf,sodf,dt1_b,dt2_b

    !*
    ! The function called with return values
    REAL(RL) :: timdif
    INTEGER(IT)  :: pointer_string
    INTEGER(IT)  :: modified_julday
    INTEGER(IT)  :: get_valid_unit

    !*
    ! Start the exectuable code
    !!-------------------------------

    ! read osb file head
    IF (lfirst .EQ. .TRUE.) THEN
        dt1_b = 0.d0
        dt2_b = 0.d0
        lfirst = .FALSE.
        CALL mjd2doy(CKF.mjd0,cfg_year,cfg_doy)

        INQUIRE(FILE = CKF.flnosbupd, EXIST = lexist)
        IF (lexist .EQ. .TRUE.) THEN
            lfnosb = get_valid_unit(10)
            OPEN(UNIT = lfnosb, FILE = CKF.flnosbupd, STATUS = 'OLD')

            line = ''
            !! 读取OSB文件第一行，判断时间是否匹配
            READ(lfnosb,'(a)',END=100) line
            READ(line,'(34X,I4,1X,I3)') year,doy
            
            IF (year .NE. cfg_year .OR. doy .NE. cfg_doy) THEN
                WRITE(ERROR_UNIT,'(A)') '***WARNING(read_osb_upd): the osb file is mismatch '//TRIM(CKF.flnosbupd)
                RETURN
            END IF

            DO WHILE (INDEX(line,'+BIAS/DESCRIPTION') .EQ. 0)
                READ(lfnosb,'(A)',END=100) line
            END DO
            READ(lfnosb,'(A)',END=100) line
            !PARAMETER_SPACING
            IF (line(1:1) .EQ. '*') READ(lfnosb,'(A)',END=100) line 
            READ(line(40:),*,IOSTAT=ierr) osbDelta

            nosbType = 0
            READ(lfnosb,'(A)',END=100) line
            DO WHILE (INDEX(line,'-BIAS/DESCRIPTION') .EQ. 0)
                READ(line(20:),*,IOSTAT=ierr) satSys,osbType

                IF (TRIM(satSys) .EQ. 'GPS') THEN
                    satSys = 'G'
                ELSEIF (TRIM(satSys) .EQ. 'GLONASS') THEN
                    satSys = 'R'
                ELSEIF (TRIM(satSys) .EQ. 'GALILEO') THEN
                    satSys = 'E'
                ELSEIF (TRIM(satSys) .EQ. 'BEIDOU') THEN
                    satSys = 'C'
                ELSEIF (TRIM(satSys) .EQ. 'QZSS') THEN
                    satSys = 'J'
                ELSEIF (TRIM(satSys) .EQ. 'LEO') THEN
                    satSys = 'L'
                ELSE
                    !WRITE(ERROR_UNIT,'(A)') '%%%WARRNING: read_osb_upd.f90, unknowning sat system '//TRIM(satSys)
                END IF

                isys = INDEX(SYS,TRIM(satSys))
                nosbType(isys) = nosbType(isys) + 1

                READ(lfnosb,'(A)',END=100) line
            END DO

            DO WHILE (INDEX(line,'+BIAS/SOLUTION') .EQ. 0)
                READ(lfnosb,'(A)',END=100) line
            END DO
            READ(lfnosb,'(A)',END=100) line

            ! read osb file contents
            READ(lfnosb,'(A)',END=100) line
            IF (INDEX(line,'-BIAS/SOLUTION') .NE. 0) RETURN
            READ(line(36:39),*,IOSTAT=ierr) year
            READ(line(41:43),*,IOSTAT=ierr) doy
            READ(line(45:49),*,IOSTAT=ierr) sodf(1)
            CALL yeardoy2monthday(year,doy,imon,iday)
            mjdf(1) = modified_julday(iday,imon,year)
            CALL timinc(mjdf(1),sodf(1),osbDelta,mjdf(2),sodf(2))

            dt1 = timdif(CKF.mjd,CKF.sod,mjdf(1),sodf(1))
            dt2 = timdif(CKF.mjd,CKF.sod,mjdf(2),sodf(2))

            IF (dt1 .LT. 0.d0) THEN
                WRITE(ERROR_UNIT,'(A)') '%%%WARRNING(read_osb_upd.f90): please check the first epoch in osb file '//TRIM(CKF.flnosbupd)
                RETURN
            END IF
        ELSE
            WRITE(ERROR_UNIT,'(A)') '***WARRING(read_osb_upd.f90): the osb file for UPD(epoch) is not exist '//TRIM(CKF.flnosbupd)
            RETURN
        END IF
    END IF

    IF (lexist .EQ. .FALSE.) RETURN

    dt1 = timdif(CKF.mjd,CKF.sod,mjdf(1),sodf(1))
    dt2 = timdif(CKF.mjd,CKF.sod,mjdf(2),sodf(2))    
    IF (dt1 .LT. 0.d0) THEN
        !WRITE(ERROR_UNIT,'(A)') '%%%WARRNING(read_osb_upd.f90): please check the first epoch in osb file '//TRIM(CKF.flnosbupd)
        RETURN
    END IF

    IF ((dt1 .GE. 0.d0 .AND. dt2 .LT. 0.d0) .AND. (dt1_b .GE. 0.d0 .AND. dt2_b .LT. 0.d0)) THEN
        RETURN
    ELSE
        !@ SMT BY XSY: 重新初始化,而不保留上历元记录,保持一致性
        ! DO isat=1,CKF.nprn
        !     SAT(isat).osbValue=0.d0
        ! END DO
    END IF

    IF (INDEX(line,'-BIAS/SOLUTION') .NE. 0) RETURN

    DO WHILE (dt2 .GE. 0.d0)
        READ(lfnosb,'(A)',END=100) line
        IF (INDEX(line,'-BIAS/SOLUTION') .NE. 0) RETURN
        READ(line(36:39),*,IOSTAT=ierr) year
        READ(Line(41:43),*,IOSTAT=ierr) doy
        READ(line(45:49),*,IOSTAT=ierr) sodf(1)
        !WRITE(*,*)year,doy,sodf(1)

        CALL yeardoy2monthday(year,doy,imon,iday)
        mjdf(1) = modified_julday(iday,imon,year)
        CALL timinc(mjdf(1),sodf(1),osbDelta,mjdf(2),sodf(2))

        dt1 = timdif(CKF.mjd, CKF.sod, mjdf(1), sodf(1))
        dt2 = timdif(CKF.mjd, CKF.sod, mjdf(2), sodf(2))
    END DO

    DO WHILE (dt1 .GE. 0.d0 .AND. dt2 .LT. 0.d0)
        READ(line,'(11X,A3,11X,A3,40X,F23.4)')cprn,freqtp,osb

        isys = index(SYS,cprn(1:1))
        IF(isys .EQ. 0) CYCLE
        
        isat = pointer_string(CKF.nprn,CKF.cprn,cprn)
        IF (isat .NE. 0) THEN
            !WRITE(*,'(1X,A3,2X,A3,F8.4)')cprn,freqtp,osb
            READ(freqtp(2:2),'(I1)')ifreq
            ind = INDEX(OBSTYPE,freqtp(3:3))

            IF (freqtp(1:1) .EQ. 'C') THEN
                SAT(isat).osbValue(1,ifreq,ind) = osb
                IF (freqtp(3:3).EQ.'W') THEN
                    indp = INDEX(OBSTYPE,'P')
                    SAT(isat).osbValue(1,ifreq,indp) = osb
                END IF
            ELSE if(freqtp(1:1) .EQ. 'L') THEN
                SAT(isat).osbValue(2,ifreq,ind) = osb
                IF (freqtp(3:3).EQ.'W') THEN
                    indp = INDEX(OBSTYPE,'P')
                    SAT(isat).osbValue(2,ifreq,indP) = osb
                END IF
            END IF 
        END IF

        READ(lfnosb,'(A)',END=100) line
        IF (INDEX(line,'-BIAS/SOLUTION') .NE. 0) RETURN
        READ(line(36:39),*,IOSTAT=ierr) year
        READ(Line(41:43),*,IOSTAT=ierr) doy
        READ(line(45:49),*,IOSTAT=ierr) sodf(1)
        CALL yeardoy2monthday(year,doy,imon,iday)
        mjdf(1) = modified_julday(iday,imon,year)
        CALL timinc(mjdf(1),sodf(1),osbDelta,mjdf(2),sodf(2))
    
        dt1 = timdif(CKF.mjd,CKF.sod,mjdf(1),sodf(1))
        dt2 = timdif(CKF.mjd,CKF.sod,mjdf(2),sodf(2))
    END DO

    dt1_b = dt1
    dt2_b = dt2

    RETURN

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_osb_upd): read the osb file error '//TRIM(line)
    CALL exit(1)

END SUBROUTINE
