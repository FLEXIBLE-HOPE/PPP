!* 
SUBROUTINE read_osb(CKF,SAT,osb_corr,OB,isit)
!*xsy:2022-10-27: ICLK
!
USE ckdctrl
USE observation
USE satellite
IMPLICIT  NONE
!*
! the arguments
!!-----------------------------
TYPE(ckdcfg)  :: CKF
TYPE(SATE)    :: SAT(MAXSAT)
REAL(RL) :: osb_corr(MAXSAT,MAXFREQ,4,isit)
TYPE(rnxobs)  :: OB
INTEGER(IT) :: isit

  !*
  !  The local variables
  !!-----------------------------
  LOGICAL(LG) :: lfirst,lexist
  INTEGER(IT) :: i,j,k,iyear,idoy,year,doy,ifreq,isys
  
  INTEGER(IT) :: lfnosb,ierr
  INTEGER(IT) :: mjdx,isat
  REAL(RL) :: sodx,osb,std
  CHARACTER(LEN_STRING) :: line
  
  CHARACTER(LEN_PRN) :: cprn,freqtp
  !DATA lfirst /.TRUE./
  !SAVE lfirst
  !SAVE lfnosb  

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
  
  CALL mjd2doy(CKF.mjd0,iyear,idoy)
  
  !DO K=1,CKF.nprn
  !  SAT(k).cosb = 0.d0
  !  SAT(k).posb = 0.d0
  !END DO
  
  !IF (lfirst .EQ. .TRUE.) THEN
    
  !  lfirst = .FALSE.
    INQUIRE(FILE=CKF.flnosb, EXIST=lexist)
    IF (lexist .EQ. .TRUE.)THEN
      lfnosb = get_valid_unit(10)
      OPEN(UNIT=lfnosb,FILE=CKF.flnosb,STATUS='OLD')
      line=' '
      !! 读取OSB文件第一行，判断时间是否匹配
      READ(lfnosb,'(a)',END=100) line
      READ(line,'(34X,I4,1X,I3)') year,doy
      
      IF (year .NE. iyear .OR. doy .NE. idoy) THEN
        WRITE(ERROR_UNIT,'(A)') '***WARNING(read_osb): the osb file is mismatch '//TRIM(CKF.flnosb)
        RETURN
      END IF 
            
      DO WHILE(INDEX(line,'+BIAS/SOLUTION') .EQ. 0)
        READ(lfnosb,'(a)',END=100) line
      END DO
      READ(lfnosb,'(a)',END=100) line
    ELSE
      WRITE(ERROR_UNIT,'(A)') '***WARNING(read_osb): the osb file is not exist '//TRIM(CKF.flnosb)
      RETURN
    END IF

    !WRITE(*,*)OB.code_type(:,1),' ',OB.phase_type(:,1)

    !! 读取osb   
    DO WHILE(INDEX(line,'-BIAS/SOLUTION') .EQ. 0)
      line = ' '
      READ(lfnosb,'(a)',END=100) line
      !LJQ OSB
      !READ(line,'(11X,A3,11X,A3,40X,F27.4,1X,F11.4)')cprn,freqtp,osb,std
      !PRIDE or GBM OSB
      READ(line,'(11X,A3,11X,A3,40X,F23.4,1X,F11.4)')cprn,freqtp,osb,std
    
      !WRITE(*,*)cprn,freqtp,osb,std

      isys = index(SYS,cprn(1:1))
      IF(isys .EQ. 0) CYCLE
      
      isat = pointer_string(CKF.nprn,CKF.cprn,cprn)
      
      IF (isat .NE. 0) THEN
      
        IF (freqtp(1:1) .EQ. 'C') THEN
          
          ifreq = pointer_string(CKF.nfreq(isys), OB.code_type(:,isys),freqtp)
          
          IF (ifreq .NE. 0) THEN
            SAT(isat).cosbstd(ifreq) = std
            osb_corr(isat,ifreq,1,isit) = osb
            osb_corr(isat,ifreq,2,isit) = std
          END IF
          
        ELSE if(freqtp(1:1) .EQ. 'L') THEN
        
          ifreq = pointer_string(CKF.nfreq(isys),OB.phase_type(:,isys),freqtp)
          
          IF (ifreq .NE. 0) THEN
            SAT(isat).posbstd(ifreq) = std
            osb_corr(isat,ifreq,3,isit) = osb
            osb_corr(isat,ifreq,4,isit) = std
          END IF
          
        END IF     
      END IF      
    END DO
  !END IF
  
  CLOSE(lfnosb)
  RETURN

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_osb): read the osb file error '//TRIM(line)
  CALL exit(1)

END SUBROUTINE
