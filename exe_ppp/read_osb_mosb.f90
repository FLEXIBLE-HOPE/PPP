!* 
SUBROUTINE read_osb_mosb(CKF,SAT)
!!
!! Read Mult-frequency OSB, PRIDE or WUMFIN
!! Created by Shengyi Xu
!! 2023-11-29
!! 在初始历元读完所有OSB,兼容ICLK模式;目前UPD模式是逐历元读取,因此不兼容
!
USE ckdctrl
USE satellite
IMPLICIT  NONE

!*
! the arguments
!!-----------------------------
TYPE(ckdcfg)  :: CKF
TYPE(SATE)    :: SAT(MAXSAT)

    !*
    !  The local variables
    !!-----------------------------
    LOGICAL(LG) :: lfirst,lexist
    INTEGER(IT) :: i,j,k,iyear,idoy,year,doy,ifreq,isys,ind,iepo,indp
    
    INTEGER(IT) :: lfnosb,ierr
    INTEGER(IT) :: isat
    INTEGER(IT) :: timespan(2),year_0,year_1,doy_0,doy_1,sod_0,sod_1
    REAL(RL) :: osb,std,sod,slope
    CHARACTER(LEN_STRING) :: line
    
    CHARACTER(LEN_PRN) :: cprn,freqtp 

    DATA lfirst /.TRUE./, lexist /.TRUE./
    SAVE lfirst, lexist

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
  
    IF (lfirst .EQ. .TRUE.) THEN
        lfirst = .FALSE.
        CALL mjd2doy(CKF.mjd0,iyear,idoy)

        INQUIRE(FILE=CKF.flnosb, EXIST=lexist)
        IF (lexist .EQ. .TRUE.)THEN
            lfnosb = get_valid_unit(10)
            OPEN(UNIT=lfnosb,FILE=CKF.flnosb,STATUS='OLD')
            line=' '
            !! 读取OSB文件第一行，判断时间是否匹配
            READ(lfnosb,'(a)',END=100) line
            IF (line(1:1) .EQ. ' ') THEN
                !! 兼容PRIDE的错误格式
                READ(line,'(35X,I4,1X,I3)') year,doy
            ELSE
                READ(line,'(34X,I4,1X,I3)') year,doy
            END IF
            
            IF (year .NE. iyear .OR. doy .NE. idoy) THEN
                WRITE(ERROR_UNIT,'(A)') '***WARNING(read_osb_mosb.f90): the osb file is mismatch '//TRIM(CKF.flnosb)
                RETURN
            END IF

            DO WHILE(INDEX(line,'+BIAS/SOLUTION') .EQ. 0)
                READ(lfnosb,'(a)',END=100) line
            END DO
            READ(lfnosb,'(a)',END=100) line
        ELSE
            WRITE(ERROR_UNIT,'(A)') '***WARNING(read_osb_mosb.f90): the osb file is not exist '//TRIM(CKF.flnosb)
            RETURN
        END IF

        !! 读取osb
        line = ' '
        READ(lfnosb,'(a)',END=100) line
        DO WHILE(INDEX(line,'-BIAS/SOLUTION') .EQ. 0)
            IF (line(1:1) .EQ. '*' .OR. TRIM(line(15:25)).NE.'') THEN
                READ(lfnosb,'(a)',END=100) line
                CYCLE
            END IF
            READ(line,'(11X,A3,11X,A3,7X,I4,1X,I3,1X,I5,1X,I4,1X,I3,1X,I5,5X,F22.4,F12.4,F22.10)')cprn,freqtp,year_0,doy_0,sod_0,year_1,doy_1,sod_1,osb,std,slope
            ! WRITE(*,'(11X,A3,11X,A3,7X,I4,1X,I3,1X,I5,1X,I4,1X,I3,1X,I5,5X,F22.4,F12.4,F22.10)')cprn,freqtp,year_0,doy_0,sod_0,year_1,doy_1,sod_1,osb,std,slope

            isys = INDEX(SYS,cprn(1:1))
            IF(isys .EQ. 0) CYCLE

            isat = pointer_string(CKF.nprn,CKF.cprn,cprn)
            IF (isat .NE. 0) THEN

                READ(freqtp(2:2),'(I1)')ifreq
                ind = INDEX(OBSTYPE,freqtp(3:3))
                
                timespan(1) = sod_0
                timespan(2) = sod_1
                ! For WUMFIN_01D_01D.OSB
                IF (timespan(2) .EQ. 0) timespan(2) = 86400
                
                iepo = INT(sod_0/900.d0)+1
                ! WRITE(*,*)iepo,osb

                IF (freqtp(1:1) .EQ. 'C') THEN
                    SAT(isat)%osbBias(1,ifreq,ind,iepo) = osb
                    SAT(isat)%osbTimeSpan(1,ifreq,ind,iepo,1) = timespan(1)
                    SAT(isat)%osbTimeSpan(1,ifreq,ind,iepo,2) = timespan(2)
                    IF (freqtp(3:3).EQ.'W') THEN
                        indp = INDEX(OBSTYPE,'P')
                        SAT(isat)%osbBias(1,ifreq,indp,iepo) = osb
                        SAT(isat)%osbTimeSpan(1,ifreq,indp,iepo,1) = timespan(1)
                        SAT(isat)%osbTimeSpan(1,ifreq,indp,iepo,2) = timespan(2)                        
                    END IF
                    
                ELSE if(freqtp(1:1) .EQ. 'L') THEN
                    SAT(isat)%osbBias(2,ifreq,ind,iepo) = osb
                    SAT(isat)%osbTimeSpan(2,ifreq,ind,iepo,1) = timespan(1)
                    SAT(isat)%osbTimeSpan(2,ifreq,ind,iepo,2) = timespan(2)
                    SAT(isat)%osbBiasSlope(ifreq,ind,iepo) = slope               
                    IF (freqtp(3:3).EQ.'W') THEN
                        indp = INDEX(OBSTYPE,'P')
                        SAT(isat)%osbBias(2,ifreq,indp,iepo) = osb
                        SAT(isat)%osbTimeSpan(2,ifreq,indp,iepo,1) = timespan(1)
                        SAT(isat)%osbTimeSpan(2,ifreq,indp,iepo,2) = timespan(2)
                        SAT(isat)%osbBiasSlope(ifreq,indp,iepo) = slope                         
                    END IF

                END IF
            END IF

            READ(lfnosb,'(a)',END=100) line    
        END DO
        CLOSE(lfnosb)
    END IF
    
    RETURN

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_osb_mosb.f90): read the osb file error '//TRIM(line)
    CALL exit(1)
END SUBROUTINE
