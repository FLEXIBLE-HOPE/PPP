!*
SUBROUTINE write_sp3_head(flnsp3,mjd0,sod0,mjd1,sod1,dintv,dused,orbt,frame,tims,nprn,cprn,acc)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-----------------------
INTEGER(IT) :: mjd0,mjd1,mjd,nprn,acc(1:*)
REAL(RL) :: sod0,sod1,sod,dintv
REAL(RL) :: x(6,1:*),clk(1:*)
CHARACTER(LEN_PRN) :: cprn(1:*)
CHARACTER(LEN=*) :: flnsp3,dused,orbt,frame,tims,type
LOGICAL(LG) :: lvel

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) :: i,isat,lfn
  INTEGER(IT) :: nepo,wk,acu(170)
  INTEGER(IT) :: iy,im,id,ih,imin

  CHARACTER(LEN_PRN) :: xprn(170)
  CHARACTER(MAXSYS) :: system

  REAL(RL) :: sec,sow
  SAVE lfn

  !*
  ! The function called
  !!-------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!-------------------------

  xprn='  0'
  xprn(1:nprn)=cprn(1:nprn)
  acu=0
  acu(1:nprn)=acc(1:nprn)

  system=''
  DO i=1, nprn
    IF (INDEX(system,cprn(i)(1:1)).EQ.0) system(LEN_TRIM(system)+1:LEN_TRIM(system)+1)=cprn(i)(1:1)
  END DO

  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=flnsp3)

  CALL mjd2date(mjd0,sod0,iy,im,id,ih,imin,sec)
  CALL mjd2wksow(mjd0,sod0,wk,sow)
  IF (tims(1:3) .EQ. 'BDT') wk=wk-1356

  nepo=NINT(((mjd1-mjd0)*86400.d0+sod1-sod0)/dintv)

  WRITE(lfn,'(A3,I4,4I3,F12.8,I8,1X,A5,1X,A5,1X,A3,1X,A4)') '#dP',iy,im,id,ih,imin,sec,nepo, &
      TRIM(dused),TRIM(frame),TRIM(orbt),' WHU'
  WRITE(lfn,'(A2,I5,F16.8,F15.8,I6,F16.9)') '##',wk,sow,dintv,mjd0,sod0
  WRITE(lfn,'("+",I5,3X,17A3,/,8("+",8X,17A3,/),"+",8X,17A3)') nprn,(xprn(i),i=1,170)
  WRITE(lfn,'(9("++",7X,17I3,/),"++",7X,17I3)') (acu(i),i=1,170)
  IF (LEN_TRIM(system) .GT. 1) THEN
    WRITE(lfn,'(A9,A3,A48)') '%c M  cc ',tims(1:3),' ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc'
  ELSE
    WRITE(lfn,'(A3,A1,A5,A3,A48)') '%c ',system(1:1),'  cc ',tims(1:3),' ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc'
  END IF
  WRITE(lfn,'(A60)') '%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc'
  WRITE(lfn,'(A60)') '%f  1.2500000  1.025000000  0.00000000000  0.000000000000000'
  WRITE(lfn,'(A60)') '%f  0.0000000  0.000000000  0.00000000000  0.000000000000000'
  WRITE(lfn,'(A60)') '%i    0    0    0    0      0      0      0      0         0'
  WRITE(lfn,'(A60)') '%i    0    0    0    0      0      0      0      0         0'
  WRITE(lfn,'(A60)') '/* GNSS RESEARCH CENTER, WUHAN UNIVERSITY (WHU), P. R. CHINA'
  WRITE(lfn,'(A60)') '/* POSITION AND NAVIGATION DATA ANALYST (PANDA) SOFTWARE    '
  WRITE(lfn,'(A60)') '/* JING GUO (EMAIL:JINGGUO@WHU.EDU.CN)                      '
  IF (orbt(1:3) .NE. 'BCT') THEN
    CALL atxwk()
    WRITE(lfn,'(A60)') '/* PCV:'//ATXWEEK//' OL/AL:FES2004  NONE     YN ORB:CoN CLK:CoN'
  ELSE
    WRITE(lfn,'(A60)') '/* BROADCAST ORBIT                                          '
  END IF

  RETURN

ENTRY write_sp3_orb(mjd,sod,lvel,nprn,cprn,x,clk,type)

  CALL mjd2date(mjd,sod,iy,im,id,ih,imin,sec)
  WRITE(lfn,'("*",I6,4I3,F12.8)') iy,im,id,ih,imin,sec

  DO isat=1, nprn
    WRITE(lfn,'("P",A3,4F14.6,15X,A1,3X,A1)') cprn(isat),(x(i,isat),i=1,3),clk(isat),type,type
    IF(lvel .EQ. .TRUE.) THEN
      WRITE(lfn,'("V",A3,4F14.6,15X,A1,3X,A1)') cprn(isat),(x(i,isat)*1.d3,i=4,6),0.d0,type,type
    END IF
  END DO

  RETURN

ENTRY write_sp3_end()

  WRITE(lfn,'(A3)') 'EOF'
  CLOSE(lfn)

  RETURN

END SUBROUTINE
