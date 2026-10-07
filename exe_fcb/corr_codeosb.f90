!*
SUBROUTINE corr_codeosb(CKF,OB,SAT,SIT,isit)
!!
!! Correction for one-day Code OSB
!!
!*
USE const
USE ckdctrl
USE observation
USE satellite
USE station
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB
TYPE(SITE) :: SIT
TYPE(SATE) :: SAT(MAXSAT)
INTEGER(IT) :: isit

    !*
    ! The local variables
    !!------------------------
    INTEGER(IT) :: isat,ifreq,iy,imon,id,ih,im,ind,freq,i
    LOGICAL(LG) :: lfirst
    CHARACTER(LEN=4) :: flag
    REAL(RL) :: OsbValue(MAXSAT,2,10,16)
    REAL(RL) :: sec
    DATA lfirst /.TRUE./
    SAVE OsbValue,lfirst

    INTEGER(IT) :: pointer_string

    !*
    ! Start the exectuable code
    !!-------------------------------
    !! 读OSB文件
    IF (lfirst .EQ. .TRUE.) THEN
        OsbValue = 0.d0
        CALL read_singledayosb(CKF,OsbValue)
        lfirst = .FALSE.
    END IF

    !! time tag
    CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,sec)
    WRITE(1008,'(A3,I5,4I3,F11.7,I7,F10.2)') 'TIM',iy,imon,id,ih,im,sec,CKF.mjd,CKF.sod  
    
    DO isat=1, CKF.nprn
      DO ifreq=1, MAXFREQ
        flag = 'Y'
        IF (OB.obs(isat,ifreq+MAXFREQ) .NE. 0) THEN
          READ(OB.fob(isat,ifreq+MAXFREQ)(2:2),'(I1)')freq
          ind = INDEX(OBSTYPE,OB.fob(isat,ifreq+MAXFREQ)(3:3))

          !xsy:伪距osb没有对应通道用其他通道补全,个别频点该通道缺伪距osb,其他通道补齐会好些
          IF (OsbValue(isat,1,freq,ind) .EQ. 0.d0) THEN
            DO i=1,16
              IF ((OsbValue(isat,1,freq,i) .NE. 0.d0)) THEN
                ind = i
                flag = 'NOT'
                EXIT
              END IF
            END DO
          END IF

          IF (OsbValue(isat,1,freq,ind) .EQ. 0.d0) flag = 'NONE'

          OB.obs(isat,ifreq+MAXFREQ) = OB.obs(isat,ifreq+MAXFREQ) - (OsbValue(isat,1,freq,ind)*VEL_LIGHT*1.0d-9)
          WRITE(1008,'(2X,A5,(A),2A4,1X,F14.4,A5,F14.3,A4,F14.3)')SIT.name,'  CODE: ',CKF.cprn(isat),OB.fob(isat,ifreq+MAXFREQ),OsbValue(isat,1,freq,ind),&
                flag,OB.obs(isat,ifreq+MAXFREQ)+OsbValue(isat,1,freq,ind)*VEL_LIGHT*1.0d-9,' => ',OB.obs(isat,ifreq+MAXFREQ)
        END IF
      END DO
    END DO

  RETURN

END SUBROUTINE

!*
SUBROUTINE read_singledayosb(CKF,OsbValue)
!!
!! Read for one-day Code OSB
!!
!*
USE const
USE ckdctrl
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------------
TYPE(CKDCFG) :: CKF
!单天OSB
!   2 :1-C;2-L
!   10:freq(1,2,5,6,7,8)
!   16:channel(PWCIXSAQLDBYMZN)
REAL(RL) :: OsbValue(MAXSAT,2,10,16)

    !*
    !  The local variables
    !!-----------------------------
    LOGICAL(LG) :: lfirst,lexist
    INTEGER(IT) :: i,j,k,iyear,idoy,year,doy,ifreq,isys,ind_freq,ind_channel
    
    INTEGER(IT) :: lfnosb,ierr
    INTEGER(IT) :: mjdx,isat
    REAL(RL) :: sodx,osb,std
    CHARACTER(LEN_STRING) :: line
    
    CHARACTER(LEN_PRN) :: cprn,freqtp

    !*
    ! The function called
    !!-----------------------------
    INTEGER(IT) :: modified_julday
    INTEGER(IT) :: get_valid_unit
    INTEGER(IT) :: pointer_string
    !*
    ! Start the exectuable code
    !!-----------------------------
  
    CALL mjd2doy(CKF.mjd0,iyear,idoy)
    
    INQUIRE(FILE=CKF.flnosbcode, EXIST=lexist)
    IF (lexist .EQ. .TRUE.)THEN
      lfnosb = get_valid_unit(10)
      OPEN(UNIT=lfnosb,FILE=CKF.flnosbcode,STATUS='OLD')
      line=' '
      !! 读取OSB文件第一行，判断时间是否匹配
      READ(lfnosb,'(a)',END=100) line
      READ(line,'(34X,I4,1X,I3)') year,doy
      
      IF (year .NE. iyear .OR. doy .NE. idoy) THEN
        WRITE(ERROR_UNIT,'(A)') '***WARNING(read_osb): the osb file is mismatch '//TRIM(CKF.flnosbcode)
        RETURN
      END IF 
            
      DO WHILE(INDEX(line,'+BIAS/SOLUTION') .EQ. 0)
        READ(lfnosb,'(a)',END=100) line
      END DO
      READ(lfnosb,'(a)',END=100) line
    ELSE
      WRITE(ERROR_UNIT,'(A)') '***WARNING(read_osb): the osb file is not exist '//TRIM(CKF.flnosbcode)
      RETURN
    END IF

    !! 读取osb,获取指定通道的伪距OSB
    DO WHILE(INDEX(line,'-BIAS/SOLUTION') .EQ. 0)
      line = ' '
      READ(lfnosb,'(a)',END=100) line
      READ(line,'(11X,A3,11X,A3,40X,F23.4,1X,F11.4)')cprn,freqtp,osb,std
      isys = index(SYS,cprn(1:1))
      IF(isys .EQ. 0) CYCLE
      
      isat = pointer_string(CKF.nprn,CKF.cprn,cprn)
      IF (isat .NE. 0) THEN
        READ(freqtp(2:2),'(I1)')ind_freq
        ind_channel = INDEX(OBSTYPE,freqtp(3:3))

        IF (freqtp(1:1) .EQ. 'C') THEN
          OsbValue(isat,1,ind_freq,ind_channel) = osb
        ELSE IF (freqtp(1:1) .EQ. 'L') THEN
          OsbValue(isat,2,ind_freq,ind_channel) = osb
        END IF
      END IF      
    END DO
  
    CLOSE(lfnosb)
    RETURN

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_osb): read the osb file error '//TRIM(line)
    CALL exit(1)

END SUBROUTINE