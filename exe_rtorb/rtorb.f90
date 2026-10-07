!*
PROGRAM rtorb
!!
!*
USE const
USE orbit
USE ISO_FORTRAN_ENV
IMPLICIT NONE

INTEGER(IT) :: lfnorb = 0
INTEGER(IT) :: mjd = 0
INTEGER(IT) :: mjdutc = 0
REAL(RL) :: x(6,MAXSAT) = 0.D0
REAL(RL) :: sod = 0.D0
REAL(RL) :: sodutc = 0.d0
CHARACTER(LEN_FILENAME) :: erpfile = ''

INTEGER(IT) :: mjdut1 = 0
INTEGER(IT) :: mjdtai = 0
REAL(RL) :: sodut1 = 0.d0
REAL(RL) :: sodtai = 0.d0
REAL(RL) :: utcut1r = 0.D0
REAL(RL) :: xhelp(2) = 0.D0
REAL(RL) :: mate2j(3,3) = 0.D0
REAL(RL) :: rmte2j(3,3) = 0.D0
REAL(RL) :: dxmat(3,3) = 0.D0
REAL(RL) :: dymat(3,3) = 0.D0
REAL(RL) :: gmst = 0.D0
REAL(RL) :: xpole = 0.D0
REAL(RL) :: ypole = 0.D0
REAL(RL) :: erp(3)

INTEGER(IT) :: i,j,k,iflag
TYPE(ORBHDR) :: CKF
LOGICAL(LG) :: lpost

  !*
  ! The function called
  !!------------------------
  INTEGER(IT) :: get_valid_unit
  REAL(RL) :: timdif

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! get arguements
  CALL get_rtorb_args(lpost,erpfile,CKF)
  WRITE(OUTPUT_UNIT,'(A,200(1X,A3))') 'GNSS Satellites: ',(CKF.cprn(i),i=1,CKF.nprn)

  !! open temp orbit file
  lfnorb=get_valid_unit(10)
  OPEN(UNIT=lfnorb,FILE='orb_temp',FORM='UNFORMATTED')

  !! write header of orbit file
  WRITE(lfnorb) CKF

  !! read position from sp3 file
  k  =0
  mjd=CKF.mjd0
  sod=CKF.sod0
  DO WHILE(timdif(mjd,sod,CKF.mjd1,CKF.sod1) .LE. MAXWND)
    CALL rdsp3i(mjd,sod,CKF.nprn,CKF.cprn,x,iflag)
    k=k+1
    IF (iflag .EQ. 1) THEN
      WRITE(OUTPUT_UNIT,'(A,I6,F8.2)') '###WARNING(rtorb): epoch lost ',mjd,sod
    ELSE
      IF (lpost .EQ. .FALSE.) THEN
        !! IGS ERP parameters in GPST
        CALL read_igserp(erpfile,mjd,sod,utcut1r,xhelp)
      ELSE
        !! IERS ERP parameters in UTC for EPO C04 products
        CALL timinc(mjd,sod,OFF_GPS2TAI,mjdtai,sodtai)
        !CALL iau_TAIUTC(DBLE(mjdtai+2400000.5D0),sodtai/86400.d0,erp(1),erp(2),i)
        !mjdutc=INT(erp(1)-2400000.5d0+erp(2))
        !sodutc=(erp(1)-2400000.5d0+erp(2)-mjdutc)*86400.d0
        !IF (DABS(sodutc-NINT(sodutc)) .LE. 5.d-5) sodutc=DBLE(NINT(sodutc))
        CALL taiutc(mjdtai,sodtai,mjdutc,sodutc)

        CALL table_linear_interpolate('poleut1',.FALSE.,mjdutc+sodutc/86400.d0,erp)
        xhelp=erp(1:2)
        utcut1r=erp(3)
      END IF
      CALL itrs2gcrs('IERS2010',mjd,sod,utcut1r,xhelp,mate2j,rmte2j,dxmat,dymat,mjdut1,sodut1,gmst,xpole,ypole)
      DO i=1,CKF.nprn
        !IF (CKF.cprn(i).EQ.'G04') WRITE(*,'(I8,F16.8,3(F16.6,2X))')mjd,sod,x(1:3,i)
        IF (x(1,i) .NE. 1.d15) CALL matmpy(mate2j,x(1,i),x(1,i),3,3,1) 
        !WRITE(*,*)CKF.cprn(i),x(1:3,i)
      END DO
    END IF
    WRITE(lfnorb) ((x(j,i),j=1,3),i=1,CKF.nprn)
    CALL timinc(mjd,sod,CKF.dintv,mjd,sod)
  END DO
  CALL rdsp3c()
  CLOSE(lfnorb)

  STOP

END PROGRAM
