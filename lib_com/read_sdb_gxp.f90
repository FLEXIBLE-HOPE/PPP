!* 
SUBROUTINE read_sdb_gxp(CKF,SDBrtype,SDB)
!!
!! Read SDB from https://www.researchgate.net/profile/Xiaopeng-Gong-2/research
!! Created by Shengyi Xu
!! 2024-01-16
!
USE ckdctrl
IMPLICIT NONE

!*
! the arguments
!!-----------------------------
TYPE(ckdcfg)  :: CKF
!!xsy: SDB correction for 50 rectype
!!    10: frequency numbers
!!    16: OBSTYPE length 
CHARACTER(LEN=50) :: SDBrtype(MAXRECTYPE)
REAL(RL) :: SDB(MAXSAT,10,16,MAXRECTYPE)

    !*
    !  The local variables
    !!-----------------------------
    LOGICAL(LG) :: lexist
    INTEGER(IT) :: i,j,ind,lfn,ierr,isat,isys,freq,ichanel,irec
    CHARACTER(LEN_STRING) :: line
    
    CHARACTER(LEN=1) :: chanel
    INTEGER(IT) :: nsignal(MAXSYS),nrectype,np
    CHARACTER(LEN_PRN) :: prn,signal(MAXSYS,50)
    REAL(RL) :: sdbValue(50)

    !*
    ! The function called
    !!-----------------------------
    INTEGER(IT) :: get_valid_unit
    INTEGER(IT) :: pointer_string
    !*
    ! Start the exectuable code
    !!-----------------------------

    INQUIRE(FILE=CKF%flnsdb, EXIST=lexist)
    IF (lexist .EQ. .TRUE.)THEN
        lfn = get_valid_unit(10)
        OPEN(UNIT=lfn,FILE=CKF%flnsdb,STATUS='OLD')
    ELSE
        WRITE(ERROR_UNIT,'(A)') '***WARNING(read_sdb_gxp.f90): the SDB file is not exist '//TRIM(CKF%flnsdb)
        RETURN
    END IF

    !! 读取SDB
    line = ' '
    READ(lfn,'(a)',END=100) line
    DO WHILE(line(1:2).EQ.'##')

        IF (line(3:5).EQ.'GPS') THEN
            isys = INDEX(SYS,'G')
            nsignal(isys) = INT((LEN(TRIM(line)) - 5)/7)
            DO ind=1,nsignal(isys)
                signal(isys,ind)=TRIM(adjustl(line(6+7*(ind-1):12+7*(ind-1))))
            END DO

        ELSE IF (line(3:5).EQ.'GAL') THEN
            isys = INDEX(SYS,'E')
            nsignal(isys) = INT((LEN(TRIM(line)) - 5)/7)
            DO ind=1,nsignal(isys)
                signal(isys,ind)=TRIM(adjustl(line(6+7*(ind-1):12+7*(ind-1))))
            END DO

        ELSE IF (line(3:5).EQ.'BDS') THEN
            isys = INDEX(SYS,'C')
            nsignal(isys) = INT((LEN(TRIM(line)) - 5)/7)
            DO ind=1,nsignal(isys)
                signal(isys,ind)=TRIM(adjustl(line(6+7*(ind-1):12+7*(ind-1))))
            END DO          

        ELSE IF (line(3:5).EQ.'QZS') THEN
            isys = INDEX(SYS,'J')
            nsignal(isys) = INT((LEN(TRIM(line)) - 5)/7)
            DO ind=1,nsignal(isys)
                signal(isys,ind)=TRIM(adjustl(line(6+7*(ind-1):12+7*(ind-1))))
            END DO            
        END IF

        READ(lfn,'(a)',END=100) line
    END DO

    np = 0
    DO WHILE(INDEX(line,'ENDSDB').EQ.0)
        IF (line(1:1) .EQ. '+') THEN
            nrectype = INT((len(TRIM(line))-4)/21)+1
            DO ind=1,nrectype
                np = np + 1
                SDBrtype(np) = TRIM(line(5+21*(ind-1):25+21*(ind-1)))
            END DO

            READ(lfn,'(a)',END=100) line
            DO WHILE(INDEX(line,'End').EQ.0)

                READ(line,'(2X,A3)')prn
                isys = INDEX(SYS,prn(1:1))
                isat = pointer_string(CKF%nprn,CKF%cprn,prn)
                IF (isat.NE.0) THEN

                    sdbValue = 0.d0
                    READ(line,'(5X,<nsignal(isys)>F7.2)')(sdbValue(i),i=1,nsignal(isys))
                    DO irec=1,nrectype
                        DO i=1,nsignal(isys)
                            READ(signal(isys,i)(2:2),'(I1)')freq
                            READ(signal(isys,i)(3:3),'(A1)')chanel
                            ! WRITE(*,'(A3,2X,A3,2X,F7.2)')prn,signal(isys,i),sdbValue(i)
                            IF (chanel.EQ.' ') THEN
                                DO ind=1,16
                                    SDB(isat,freq,ind,np-nrectype+irec) = sdbValue(i)
                                END DO
                            ELSE
                                ind = INDEX(OBSTYPE,chanel)
                                SDB(isat,freq,ind,np-nrectype+irec) = sdbValue(i)
                            END IF
                        END DO
                    END DO
                END IF
                READ(lfn,'(a)',END=100) line
            END DO
        END IF
        READ(lfn,'(a)',END=100) line
    END DO

    CLOSE(lfn)
    
    RETURN

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sdb_gxp.f90):  error '//TRIM(line)
    CALL exit(1)
END SUBROUTINE
