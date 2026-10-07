!*
!! xsy-2022-07-26: add NLOS detection and NLOSflag
!*
SUBROUTINE read_rnxobs(leo,lfn,jd0,sod0,nprn0,cprn0,nfreq,freqused,HD,OB,ierr,SNR_limit)
!*
USE par
USE const
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: leo,lfn,jd0,nprn0,ierr,nfreq(MAXSYS),SNR_limit
REAL(RL) :: sod0
CHARACTER(LEN_PRN) :: cprn0(MAXSAT)
CHARACTER(LEN_FREQ) :: freqused(MAXFREQ,MAXSYS)
TYPE(RNXHEAD) :: HD
TYPE(RNXOBS) :: OB

  !*
  ! The local variables
  !!---------------------------------
  INTEGER(IT) :: ioerr,iy,im,id,ih,imi,nprn,    outsat,numnlos,iy0,imon0,id0,ih0,im0
  INTEGER(IT) :: iflag,i,j,l,k,nline,isat,iobs,ifreq,isys
  REAL(RL) :: sec,ds,dt,obs(MAXOBSTYP),litrg,  sec0
  CHARACTER(LEN_STRING) :: line,fmt,msg
  CHARACTER(LEN_PRN) :: cprn(MAXSAT)
  CHARACTER(LEN_OBSTYPE) :: code
  
  INTEGER(IT) :: nlosflag(MAXOBSTYP)
  INTEGER(IT) :: lli(MAXOBSTYP)

  !*
  ! The function called
  !!----------------------------------
  INTEGER(IT) :: pointer_string
  INTEGER(IT) :: modified_julday

  !*
  ! Start the exectuable code
  !!-----------------------------------
      
  ierr=0
  line=' '
  !! Maxmium differences between L1 and L2
  !! Please becareful the therhold 50.d0, it is effeced by ionosphere
  !! for LEO satellites, it is a bit too small, to be 500
  ! litrg=1000.d0 !for HISI
  litrg=50.d0
  IF (leo .NE. 0) litrg=500.d0
   
  OB.nlosflag = 0   !每个历元需要清0
  OB%lli=0

10 READ(lfn,'(a)',END=200) line
  IF (LEN_TRIM(line) .EQ. 0) GOTO 10
    cprn=''
    nprn=0
    msg=''

    !! number of satellite
    IF(HD.ver .LT. 3.d0)THEN
      READ(line(30:),'(i3)',iostat=ioerr) nprn
    ELSE
      READ(line(33:),'(i3)',iostat=ioerr) nprn
    END IF

    IF(ioerr .NE. 0)  THEN
      msg = 'READ satellite number error.'
    ELSE IF(nprn .GT. MAXSAT) THEN
      msg = 'satellite number > MAX_SAT'
    END IF
    IF(LEN_TRIM(msg) .NE. 0) GOTO 100

    IF(HD.ver .LT. 3.d0)THEN
      READ(line(27:),'(i3)',iostat=ioerr) iflag
    ELSE
      READ(line(30:),'(i3)',iostat=ioerr) iflag
    END IF
    IF(ioerr .NE. 0) THEN
      msg = 'read event flag error.'
      GOTO 100
    ELSE IF(iflag .GT. 1) THEN
      IF(HD.ver .LT. 3.d0)THEN
        msg = 'read internal antenna information error '
        DO i=1, nprn
          READ(lfn,iostat=ioerr,END=200,fmt='(a80)') line
          IF(line(61:80) .eq. 'ANTENNA: DELTA H/E/N') THEN
            READ (line,'(3f14.4)',err=100) HD.h,HD.e,HD.n
          END IF
        END DO
        GOTO 10
     END IF
   END IF

   !! initialization
   OB.obs=0.d0

   !! format of the time tag line
   IF(HD.ver .LT. 3.d0)THEN
     fmt='(5i3,f11.7,2i3,12a3)'
     msg = 'read time & svn error.'
     READ(line,fmt=fmt,err=100) iy,im,id,ih,imi,sec,iflag,nprn
   ELSE
     fmt='(2x,i4,4(1x,i2.2),f11.7,2x,i1,i3,6x,f15.12)'
     msg = 'read time & svn error.'
     READ(line,fmt=fmt,err=100) iy,im,id,ih,imi,sec,iflag,nprn,OB.dtrcv
   END IF

   IF(nprn .GT. MAXSAT) THEN
     WRITE(ERROR_UNIT,'(A)') '***ERROR(read_rnxobs): nprn > MAXSAT'
     CALL exit(1)
   END IF
   IF(HD.ver .LT. 3.d0)THEN
     READ(line(69:80),'(f12.9)',iostat=ioerr) OB.dtrcv
     IF(ioerr .NE. 0) OB.dtrcv = 0.d0
     IF (MOD(nprn,12) .NE. 0) THEN
       nline=INT(nprn/12)+1
     ELSE
       nline=INT(nprn/12)
     END IF
     IF (nprn .EQ. 0) nline=1

     BACKSPACE(lfn)
     DO i=1, nline
       READ(lfn,'(A)',ERR=100) line
       READ(line,'(32X,12A3)',err=100) (cprn(j),j=(i-1)*12+1,MIN(nprn,12*i))
     END DO

   END IF
   !! check time
   IF(im.LE.0 .OR. im.GT.12    &
       .OR. id.LE.0 .OR. id.GT.31 .OR. ih.LT.0 .OR. ih.GE.24  &
       .OR. imi.LT.0 .OR. imi.GT.60 .OR.sec .LT. 0.d0.OR.sec.GT.60.d0  &
     ) THEN
     msg = 'epoch time incorrect.'
     GOTO 100
   END IF
   CALL yr2year(iy)

   !! update start and stop time in rinex header
   IF (hd.t0(1).EQ.0) THEN
      hd.t0(1)=iy
      hd.t0(2)=im
      hd.t0(3)=id
      hd.t0(4)=ih
      hd.t0(5)=imi
      hd.t0(6)=nint(sec)
   END IF
   hd.t1(1)=iy
   hd.t1(2)=im
   hd.t1(3)=id
   hd.t1(4)=ih
   hd.t1(5)=imi
   hd.t1(6)=nint(sec)

   !! check on time tags. do not change the requested time IF there is no data
   ds=0.d0
   OB.jd=modified_julday(id,im,iy)
   OB.tsec=ih*3600.d0 + imi*60.d0 + sec
   IF(jd0.ne.0) THEN
     ds=(jd0-OB.jd)*86400.d0 + (sod0-OB.tsec)
     IF(ds .LT. -MAXWND) THEN
       OB.jd=jd0
       OB.tsec=sod0
       IF(HD.ver.LT.3.d0)THEN

         IF (MOD(nprn,12) .NE. 0) THEN
           nline=INT(nprn/12)+1
         ELSE
           nline=INT(nprn/12)
         END IF
         IF (nprn .eq. 0) nline=1
            
         DO i =1, nline
           BACKSPACE lfn
         END DO
       ELSE
         BACKSPACE lfn
       END IF
       OB.nprn=0
       RETURN
     ELSE IF(ds .GT. MAXWND) THEN
       nline = nprn
       IF(HD.ver .LT. 3.d0)THEN
         IF (MOD(HD.nobstype(1),5) .NE. 0) THEN
           nline=nprn*(INT(HD.nobstype(1)/5)+1)
         ELSE
           nline=nprn*(INT(HD.nobstype(1)/5))
         END IF
       END IF
       DO i=1,nline
         READ(lfn,'(a)',END=100) line
       END DO
       GOTO 10
     END IF
   END IF

   IF(HD.ver .LT. 3.d0)THEN
     DO i=1, nprn
       IF(cprn(i)(1:1) .EQ.' ') cprn(i)(1:1)='G'
       IF(cprn(i)(2:2) .EQ.' ') cprn(i)(2:2)='0'
     END DO

     DO i=1, nprn

       IF (MOD(HD.nobstype(1),5) .NE. 0) THEN
         l = INT(HD.nobstype(1)/5)+1
       ELSE
         l = INT(HD.nobstype(1)/5)
       END IF

       DO k=1, l
         READ(lfn,err=100,END=200,fmt='(A)') line
         READ(line,fmt='(5(f14.3,2X))',ERR=100) (obs(j),j=(k-1)*5+1,MIN(HD.nobstype(1),5*k))
       END DO

400    isat=0
       IF(nprn0 .GT. 0) THEN
         isat=pointer_string(nprn0,cprn0,cprn(i))
       ELSE
         isat=i
       END IF

      IF(isat .NE. 0) THEN

        isys=INDEX(SYS,cprn(i)(1:1))

        DO ifreq=1, nfreq(isys)

          WRITE(code,'(A1,A1)') 'P',freqused(ifreq,isys)(2:2)
          k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
          IF (k .NE. 0) THEN
            IF (obs(k) .GT. 1.d0 ) THEN
              OB.obs(isat,MAXFREQ+ifreq) = obs(k)
              OB.fob(isat,MAXFREQ+ifreq) = 'C'//freqused(ifreq,isys)(2:2)//'P'
            END IF
          END IF

          IF (OB.obs(isat,MAXFREQ+ifreq) .EQ. 0.d0) THEN
            ! Please must sure the C observations are converted to P observations
            WRITE(code,'(A1,A1)') 'C',freqused(ifreq,isys)(2:2)
            k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
            IF (k .NE. 0) THEN
              IF (obs(k) .GT. 1.d0 ) THEN
                OB.obs(isat,MAXFREQ+ifreq) = obs(k)
                OB.fob(isat,MAXFREQ+ifreq) = 'C'//freqused(ifreq,isys)(2:2)//'C'
              END IF
            END IF
          END IF

          WRITE(code,'(A1,A1)') 'L',freqused(ifreq,isys)(2:2)
          k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
          IF (k .NE. 0) THEN
            IF (DABS(obs(k)) .GT. 1.d0 ) THEN
              OB.obs(isat,ifreq) = obs(k)
              OB.fob(isat,ifreq) = 'L'//freqused(ifreq,isys)(2:2)//'P'
            END IF
          END IF

          WRITE(code,'(A1,A1)') 'D',freqused(ifreq,isys)(2:2)
          k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
          IF (k .NE. 0) THEN
            IF (DABS(obs(k)) .GT. 1.d0 ) THEN
              OB.obs(isat,2*MAXFREQ+ifreq) = obs(k)
              OB.fob(isat,2*MAXFREQ+ifreq) = 'D'//freqused(ifreq,isys)(2:2)//'P'
            END IF
          END IF

          WRITE(code,'(A1,A1)') 'S',freqused(ifreq,isys)(2:2)
          k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
          IF (k .NE. 0) THEN
            IF (DABS(obs(k)) .GT. 1.d0 ) THEN
              OB.obs(isat,3*MAXFREQ+ifreq) = obs(k)
              OB.fob(isat,3*MAXFREQ+ifreq) = 'S'//freqused(ifreq,isys)(2:2)//'P'
            END IF
            IF (obs(k) .LT. 32.d0) THEN
              OB.obs(isat,1:4*MAXFREQ) = 0.D0
              OB.fob(isat,1:4*MAXFREQ) = ' '
            END IF
          END IF
        END DO

        IF (COUNT(OB.obs(isat,1:2*MAXFREQ).NE.0.D0) .LT. 2*nfreq(isys)) THEN
          OB.obs(isat,1:4*MAXFREQ) = 0.D0
          OB.fob(isat,1:4*MAXFREQ) = ' '
        END IF

        IF (nfreq(isys).GE.2 .AND. (DABS(OB.obs(isat,MAXFREQ+1)-OB.obs(isat,MAXFREQ+2)).GT.litrg)) THEN
          OB.obs(isat,1:4*MAXFREQ) = 0.D0
          OB.fob(isat,1:4*MAXFREQ) = ' '
        END IF

      END IF
    END DO
  END IF

  IF(HD.ver .GE. 3.d0)THEN
    isat=0
    cprn(:)=''
    DO i=1,nprn
      READ (lfn,err=100,END=200,fmt='(a)') line
      !@ SMT BY XSY:  W/X/0 > L for LEOPPP
      IF (line(1:1) .EQ. 'W') line(1:1)='L'
      IF (line(1:1) .EQ. 'X') line(1:1)='L'
      IF (line(1:1) .EQ. '0') line(1:1)='L'

      isys = index(SYS,line(1:1))
      IF(isys .EQ. 0) CYCLE

      IF(nprn0 .GT. 0) THEN
        IF (line(2:2) .EQ. ' ') line(2:2)='0'
        isat=pointer_string(nprn0,cprn0,line(1:3))
        IF (isat .EQ. 0) CYCLE
      ELSE
        isat=isat+1
        cprn(isat)=line(1:3)
      END IF

      READ(line,'(3x,<MAXOBSTYP>(f14.3,2X))',ERR=100) (obs(k),k=1,HD.nobstype(isys))
      
      !@CMT BY XSY: NLOS flag [L is LOS, 0 is NLOC]
      READ(line,'(3X,<MAXOBSTYP>(14X,I2))',err=100)(nlosflag(k),k=1,HD.nobstype(isys))
      !@CMT BY XSY: LLI for cycle slip detect
      READ(line,'(3X,<MAXOBSTYP>(14X,I1,1X))',err=100)(lli(k),k=1,HD.nobstype(isys))
      
      !@CMT BY XSY: 如果所有历元不固定使用某一码通道，则取消注释，此处表示每颗卫星有什么通道用什么通道，只有对于卫星数不多的时候可以尝试
      !HD.usetype = ''

      IF(isat .NE. 0) THEN

        !! only the first type observation for the frequency has been read
        DO ifreq=1, nfreq(isys)
          !@CMT BY XSY: 首历元根据通道优先级,初始化使用的码通道,之后就会固定使用该通道类型，如果想要不限定通道类型，每个历元都将重新初始化HD.usetype(ifreq,isys)即可
          IF (LEN_TRIM(HD.usetype(ifreq,isys)) .EQ. 0) THEN
            DO iobs=1, LEN(OBSTYPE)
              WRITE(code,'(A1,A1,A1)') 'C',freqused(ifreq,isys)(2:2),OBSTYPE(iobs:iobs)
              k=0
              k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
              IF (k .NE. 0) THEN
                IF (obs(k) .GT. 1.d0 ) THEN
                  OB.obs(isat,MAXFREQ+ifreq)=obs(k)
                  OB.fob(isat,MAXFREQ+ifreq)=code
                END IF
              END IF

              WRITE(code,'(A1,A1,A1)') 'L',freqused(ifreq,isys)(2:2),OBSTYPE(iobs:iobs)
              k=0
              k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
              IF (k .NE. 0) THEN
                IF (DABS(obs(k)) .GT. 1.d0 ) THEN
                  OB.obs(isat,ifreq)=obs(k)
                  OB.fob(isat,ifreq)=code
                END IF
              END IF

              WRITE(code,'(A1,A1,A1)') 'D',freqused(ifreq,isys)(2:2),OBSTYPE(iobs:iobs)
              k=0
              k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
              IF (k .NE. 0) THEN
                IF (DABS(obs(k)) .GT. 1.d0 ) THEN
                  OB.obs(isat,2*MAXFREQ+ifreq)=obs(k)
                  OB.fob(isat,2*MAXFREQ+ifreq)=code
                END IF
              END IF

              WRITE(code,'(A1,A1,A1)') 'S',freqused(ifreq,isys)(2:2),OBSTYPE(iobs:iobs)
              k=0
              k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
              IF (k .NE. 0) THEN
                IF (DABS(obs(k)) .GT. 1.d0 ) THEN
                  OB.obs(isat,3*MAXFREQ+ifreq)=obs(k)
                  OB.fob(isat,3*MAXFREQ+ifreq)=code
                  !@ CMT BY XSY: 在SNR观测值存在的情况下去剔除低SNR观测
                  IF (DABS(obs(k)).LT.SNR_limit)THEN
                    OB.obs(isat,          ifreq) = 0.D0
                    OB.fob(isat,          ifreq) = ' ' 
                    OB.obs(isat,  MAXFREQ+ifreq) = 0.D0
                    OB.fob(isat,  MAXFREQ+ifreq) = ' ' 
                    OB.obs(isat,2*MAXFREQ+ifreq) = 0.D0
                    OB.fob(isat,2*MAXFREQ+ifreq) = ' '                                                       
                    OB.obs(isat,3*MAXFREQ+ifreq) = 0.D0
                    OB.fob(isat,3*MAXFREQ+ifreq) = ' ' 
                  END IF                  
                END IF
              END IF

              !@ CMT BY XSY: 根据通道顺序搜索OBSTYPE，频率通道的伪距和相位都存在，即认为是该频率的使用观测通道
              IF (OB.obs(isat,MAXFREQ+ifreq).NE.0.D0 .AND. OB.obs(isat,ifreq).NE.0.D0)THEN
                WRITE(code,'(A1,A1,A1)') 'L',freqused(ifreq,isys)(2:2),OBSTYPE(iobs:iobs)                
                WRITE(OB.phase_type(ifreq,isys),'(A)') code
                WRITE(code,'(A1,A1,A1)') 'C',freqused(ifreq,isys)(2:2),OBSTYPE(iobs:iobs)
                WRITE(OB.code_type(ifreq,isys),'(A)') code
                EXIT
              END IF
            END DO
            ! WRITE(*,'(A4,I3,2A4)')cprn0(isat),ifreq,OB.fob(isat,MAXFREQ+ifreq),OB.fob(isat,ifreq)
            IF (iobs .LE. LEN(OBSTYPE)) HD.usetype(ifreq,isys)=OBSTYPE(iobs:iobs)

          ELSE
            WRITE(code,'(A1,A1,A1)') 'C',freqused(ifreq,isys)(2:2),HD.usetype(ifreq,isys)
            k=0
            k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
            IF (k .NE. 0) THEN
              IF (obs(k) .GT. 1.d0 ) THEN
                OB.obs(isat,MAXFREQ+ifreq)=obs(k)
                OB.fob(isat,MAXFREQ+ifreq)=code
              END IF
            END IF

            WRITE(code,'(A1,A1,A1)') 'L',freqused(ifreq,isys)(2:2),HD.usetype(ifreq,isys)
            k=0
            k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
            IF (k .NE. 0) THEN
              IF (DABS(obs(k)) .GT. 1.d0 ) THEN
                OB.obs(isat,ifreq)=obs(k)
                OB.fob(isat,ifreq)=code
                OB%lli(isat,ifreq)=lli(k)
              END IF
            END IF

            WRITE(code,'(A1,A1,A1)') 'D',freqused(ifreq,isys)(2:2),HD.usetype(ifreq,isys)
            k=0
            k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
            IF (k .NE. 0) THEN
              IF (DABS(obs(k)) .GT. 1.d0 ) THEN
                OB.obs(isat,2*MAXFREQ+ifreq)=obs(k)
                OB.fob(isat,2*MAXFREQ+ifreq)=code
              END IF
            END IF

            WRITE(code,'(A1,A1,A1)') 'S',freqused(ifreq,isys)(2:2),HD.usetype(ifreq,isys)
            k=0
            k = pointer_string(HD.nobstype(isys), HD.obstype(:,isys), TRIM(code))
            IF (k .NE. 0) THEN
              IF (DABS(obs(k)) .GT. 1.d0 ) THEN
                OB.obs(isat,3*MAXFREQ+ifreq)=obs(k)
                OB.fob(isat,3*MAXFREQ+ifreq)=code
                !@ CMT BY XSY: 在SNR观测值存在的情况下去剔除低SNR观测
                IF (obs(k) .LT. SNR_limit) THEN
                  OB.obs(isat,          ifreq) = 0.D0
                  OB.fob(isat,          ifreq) = ' ' 
                  OB.obs(isat,  MAXFREQ+ifreq) = 0.D0
                  OB.fob(isat,  MAXFREQ+ifreq) = ' ' 
                  OB.obs(isat,2*MAXFREQ+ifreq) = 0.D0
                  OB.fob(isat,2*MAXFREQ+ifreq) = ' '                                                       
                  OB.obs(isat,3*MAXFREQ+ifreq) = 0.D0
                  OB.fob(isat,3*MAXFREQ+ifreq) = ' ' 
                END IF
              END IF
            END IF
          END IF

          !@ CMT BY XSY: 同一频率的码相观测值同时存在
          IF (OB.obs(isat,ifreq).EQ.0.D0 .OR. OB.obs(isat,MAXFREQ+ifreq).EQ.0.D0) THEN
            OB.obs(isat,          ifreq) = 0.D0
            OB.fob(isat,          ifreq) = ' ' 
            OB.obs(isat,  MAXFREQ+ifreq) = 0.D0
            OB.fob(isat,  MAXFREQ+ifreq) = ' ' 
            OB.obs(isat,2*MAXFREQ+ifreq) = 0.D0
            OB.fob(isat,2*MAXFREQ+ifreq) = ' '                                                       
            OB.obs(isat,3*MAXFREQ+ifreq) = 0.D0
            OB.fob(isat,3*MAXFREQ+ifreq) = ' '
          END IF
        END DO

        !@ CMT BY XSY: 判断nlos标识[1:LOS, 0:NLOS]
        DO K=1,HD.nobstype(isys)
          IF (nlosflag(k).EQ.1)THEN!1:LOS
            OB.nlosflag(isat)=1
            EXIT
          END IF
        END DO
        
        !@ CMT BY XSY: GPS部分卫星跟踪L5,其他系统不存在这种特性,因此要求非GPS系统指定频率的观测值都存在
        IF (line(1:1).NE.'G') THEN
          IF (COUNT(OB.obs(isat,1:2*MAXFREQ).NE.0.D0) .LT. 2*nfreq(isys)) THEN
            OB.obs(isat,1:4*MAXFREQ) = 0.D0
            OB.fob(isat,1:4*MAXFREQ) = ' '
          END IF
        END IF

        !@ CMT BY XSY: 如果大于等于双频,则前两个频率的观测值必须存在,第三频率可以不存在[GPS]
        IF (nfreq(isys).GE.2 .AND. (DABS(OB.obs(isat,MAXFREQ+1)-OB.obs(isat,MAXFREQ+2)).GT.litrg)) THEN
          OB.obs(isat,1:4*MAXFREQ) = 0.D0
          OB.fob(isat,1:4*MAXFREQ) = ' '
        END IF

      END IF
    END DO
  END IF

  IF(nprn0 .EQ. 0) THEN
    OB.nprn=nprn
    DO i=1,nprn
      OB.cprn(i)=cprn(i)
    END DO
  ELSE
    OB.nprn=nprn0
    DO i=1,nprn0
      OB.cprn(i)=cprn0(i)
    END DO
  END IF

  !@ CMT BY XSY: 判断nlos标识[1:LOS, 0:NLOS]
  outsat = 0
  numnlos = 0
  DO isat=1,nprn0
    IF (OB.obs(isat,1+MAXFREQ) .NE. 0.d0) THEN
      outsat = outsat + 1
      IF (OB.nlosflag(isat) .EQ. 0) THEN
        numnlos = numnlos + 1 !for SM, numnlos is the num of los
      END IF
    END IF
  END DO

  CALL mjd2date(jd0,sod0,iy0,imon0,id0,ih0,im0,sec0)
  WRITE(1002,'(A3,I5,4I3,F11.7,I7,F10.2,1X,I3,1X,I3,4X,A,1X,I2)') 'TIM',iy0,imon0,id0,ih0,im0,sec0,jd0,sod0,numnlos,outsat,'RAW_OBS-SNR=',SNR_limit
  DO isat=1,nprn0
    DO i=1,MAXFREQ
      IF (OB.obs(isat,MAXFREQ+i) .NE. 0) then
         WRITE(1002,'(2X,A3,2X,I3,2X,A3,4(f14.3,2X),1X,I2)')OB.cprn(isat),i,OB.fob(isat,i),OB.obs(isat,i),OB.obs(isat,i+MAXFREQ),OB.obs(isat,i+2*MAXFREQ),OB.obs(isat,i+3*MAXFREQ),OB.nlosflag(isat)
      END IF
    END DO
  END DO

  IF(ds .GT. MAXWND) GOTO 10

  RETURN

100 ierr=1
    WRITE(OUTPUT_UNIT,'(A/A/A)') '***WARNING(read_rnxobs): read file', &
                       '   line :'//TRIM(line), &
                       '   msg  :'//TRIM(msg)
  RETURN

200 ierr=1
    OB.obs=0.d0

  RETURN

END SUBROUTINE

