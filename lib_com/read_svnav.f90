!*
SUBROUTINE read_svnav(mjd,SAT)
!!
!! TO GET THE PHYSICAL PROPERTIES OF SATELLITES
!!
!*  
USE par
USE tables
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! THE ARGUMENTS
!!-----------------------------
REAL(RL) :: mjd
TYPE(SATE) :: SAT

  !*
  ! THE LOCAL VARIABLES
  !!-----------------------------
  INTEGER(IT) :: lfn,iq

  INTEGER(IT) :: iy,im,id,ih,imin,isec
  INTEGER(IT) :: ky,km,kd,kh,kmin,ksec

  REAL(RL) :: enu(3),bvec(3),avec(3)
  REAL(RL) :: mjds,mjde,imass

  CHARACTER(LEN_PRN) :: prn,svn
  CHARACTER(LEN_SATTYPE) :: typ
  CHARACTER(LEN=9) :: cid
  CHARACTER(LEN_STRING) :: line
  CHARACTER(LEN_STRING) :: form

  LOGICAL(LG) :: lfirst,lexist,lfind
  DATA lfirst /.TRUE./

  SAVE lfirst,lfn

  !*
  ! THE FUNCTIONS CALLED
  !!-----------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: modified_julday


  !*
  ! THE EXECTUABLE CODES
  !!-----------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    line=f_tableFileName('satnav')
    INQUIRE(FILE=line,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_svnav): no '//TRIM(line)
      CALL exit(1)
    END IF

    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=line,STATUS='OLD',ACTION='READ')
  END IF

  REWIND(lfn)
  line=''
  DO WHILE(INDEX(line,'PART 1') .EQ. 0)
    READ(lfn,'(A)') line
  END DO
  READ(lfn,'(A)') form
  READ(lfn,'(A)') line
  READ(lfn,'(A)') line
  READ(lfn,'(A)') line
  READ(lfn,'(A)') line

  line=''
  lfind=.FALSE.
  DO WHILE(.TRUE.)
    READ(lfn,'(A)',END=100,ERR=300) line
    IF (INDEX(line,'PART 2') .NE. 0) EXIT
    READ(line,form) prn,svn,typ,cid,iy,im,id,ih,imin,isec,ky,km,kd,kh,kmin,ksec,imass,iq
    IF (prn .NE. SAT.cprn) CYCLE
    mjds=modified_julday(id,im,iy)+ih/24.d0+imin/1440.d0+isec/86400.d0
    IF (ky .EQ. 0) THEN
      mjde=mjd+1.d0
    ELSE
      mjde=modified_julday(kd,km,ky)+kh/24.d0+kmin/1440.d0+ksec/86400.d0
    END IF
    IF (mjd.GE.mjds .AND. mjd.LE.mjde) THEN
      lfind=.TRUE.
      EXIT
    END IF
  END DO
100 CONTINUE
  IF (lfind .EQ. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(read_svnav): satellite lost '//SAT.cprn//' in svnav.dat file'
    CALL exit(1)
  ELSE
    SAT.csvn=svn
    SAT.type=TRIM(typ)
    SAT.ifreq=iq
    SAT.mass=imass
    READ(cid,'(I4,1X,I3)') iy,im
    CALL year2yr(iy)
    WRITE(SAT.sid,'(I2.2,I3.3,I2.2)') iy,im,ICHAR(cid(9:9))-ICHAR('A')+1
  END IF

  DO WHILE(INDEX(line,'PART 2') .EQ. 0)
    READ(lfn,'(A)') line
  END DO
  READ(lfn,'(A)') form
  READ(lfn,'(A)') line
  READ(lfn,'(A)') line
  READ(lfn,'(A)') line
  READ(lfn,'(A)') line

  SAT.offs=0.d0
  SAT.xyz0=0.d0
  line=''
  DO WHILE(.TRUE.)
    READ(lfn,form,END=200,ERR=300) prn,svn,iy,im,id,ih,imin,isec,ky,km,kd,kh,kmin,ksec,enu,bvec,avec
    IF (prn .NE. SAT.cprn) CYCLE
    mjds=modified_julday(id,im,iy)+ih/24.d0+imin/1440.d0+isec/86400.d0
    IF (ky .EQ. 0) THEN
      mjde=mjd+1.d0
    ELSE
      mjde=modified_julday(kd,km,ky)+kh/24.d0+kmin/1440.d0+ksec/86400.d0
    END IF
    IF (mjd.GE.mjds .AND. mjd.LE.mjde) THEN
      SELECT CASE(TRIM(svn))
        CASE('SLR')
          SAT.offs=enu
        CASE('ISL')
          SAT.isl=enu
          SAT.avec=avec
          SAT.bvec=bvec
        CASE('POD')
          SAT.xyz0=enu
          SAT.avec=avec
          SAT.bvec=bvec

          CALL cross(SAT.avec,SAT.bvec,SAT.rant(1:3,2))
          SAT.rant(1:3,1)=SAT.avec
          SAT.rant(1:3,3)=SAT.bvec
      END SELECT
    END IF
  END DO
200 CONTINUE
         
  RETURN

300 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_svnav): read '//TRIM(line)
  CALL exit(1)

END SUBROUTINE
