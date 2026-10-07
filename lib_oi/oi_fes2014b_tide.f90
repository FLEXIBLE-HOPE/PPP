!*
SUBROUTINE oi_fes2014b_tide(gmst,mjd,lmax,dc,ds)
!!
!!
!*
USE tables
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------------
INTEGER(IT) :: lmax
REAL(RL) :: gmst,mjd
! corrections to spherical harmonics due to ocean tide
REAL(RL) :: dc(MAXOCNDEG,0:MAXOCNDEG)
REAL(RL) :: ds(MAXOCNDEG,0:MAXOCNDEG)

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: i,j,k,lfn,ierr
  INTEGER(IT), PARAMETER :: ntide=34
  INTEGER(IT), PARAMETER :: ndood=361
  INTEGER(IT) :: doodsonMatrix(ndood,6)

  REAL(RL) :: arg(5),beta(6),thetaf(ndood)
  REAL(RL) :: factorCos(ntide),factorSin(ntide)
  REAL(RL) :: admittance(ndood,ntide)
  REAL(RL) :: cnmCos(MAXOCNDEG,0:MAXOCNDEG,ntide)
  REAL(RL) :: cnmSin(MAXOCNDEG,0:MAXOCNDEG,ntide)
  REAL(RL) :: snmCos(MAXOCNDEG,0:MAXOCNDEG,ntide)
  REAL(RL) :: snmSin(MAXOCNDEG,0:MAXOCNDEG,ntide)

  CHARACTER(LEN=4) :: darwin(ntide)
  CHARACTER(LEN=7) :: doodson(ndood)

  DATA darwin /'om1 ','om2 ','sa  ', &
               'ssa ','mm  ','mf  ', &
               'mtm ','msq ','q1  ', &
               'o1  ','p1  ','s1  ', &
               'k1  ','j1  ','eps2', &
               '2n2 ','mu2 ','n2  ', &
               'nu2 ','m2  ','la2 ', &
               'l2  ','t2  ','s2  ', &
               'r2  ','k2  ','m3  ', &
               'n4  ','mn4 ','m4  ', &
               'ms4 ','s4  ','m6  ', &
               'm8  '/

  DATA doodson /'055.565','055.575','055.765','056.544','056.554','056.556','057.355','057.553','057.555','057.565', &
                '057.575','058.554','059.553','062.656','063.445','063.645','063.655','063.665','064.456','064.555', &
                '064.654','065.445','065.455','065.465','065.655','065.665','065.675','066.454','067.455','067.465', &
                '067.475','071.755','072.556','073.545','073.555','073.565','073.755','074.554','074.556','074.566', &
                '075.345','075.355','075.365','075.555','075.565','075.575','075.585','076.554','076.564','077.355', &
                '077.365','081.655','082.456','082.656','082.666','083.445','083.455','083.465','083.655','083.665', &
                '083.675','084.456','084.466','084.555','085.255','085.455','085.465','085.475','085.675','086.454', &
                '086.464','087.255','091.555','091.755','092.556','092.566','093.355','093.555','093.565','093.575', &
                '095.355','095.365','095.375','0a1.655','0a1.665','0a2.456','0a3.455','0a3.465','0a5.255','0a5.265', &
                '0b1.555','0b3.355','105.955','107.745','107.755','109.555','115.845','115.855','117.645','117.655', &
                '118.654','119.445','119.455','124.756','125.745','125.755','126.556','126.655','126.754','127.545', &
                '127.555','127.755','128.544','128.554','129.355','129.555','133.855','134.646','134.656','135.435', &
                '135.635','135.645','135.655','135.855','136.456','136.555','136.644','136.654','137.435','137.445', &
                '137.455','137.655','137.665','138.444','138.454','139.455','143.535','143.745','143.755','144.546', &
                '144.556','145.535','145.545','145.555','145.557','145.755','145.765','146.544','146.554','147.355', &
                '147.545','147.555','147.565','148.554','149.355','152.656','153.645','153.655','154.656','155.435', &
                '155.445','155.455','155.645','155.655','155.665','155.675','156.555','156.654','157.445','157.455', &
                '157.465','158.454','161.557','162.546','162.556','163.535','163.545','163.555','163.557','163.755', &
                '164.554','164.555','164.556','164.566','165.345','165.545','165.553','165.555','165.565','165.575', &
                '166.554','166.564','167.355','167.365','167.553','167.555','167.565','167.575','168.554','172.656', &
                '173.445','173.645','173.655','173.665','174.456','174.555','175.445','175.455','175.465','175.475', &
                '175.655','175.665','175.675','176.454','177.455','177.465','181.755','182.556','183.545','183.555', &
                '183.565','185.355','185.365','185.555','185.565','185.575','185.585','191.655','193.455','193.465', &
                '193.655','193.665','193.675','195.255','195.455','195.465','195.475','1a1.555','1a3.355','1a3.555', &
                '1a3.565','1a5.355','1a5.365','1b3.455','1b3.465','207.855','209.655','215.955','217.755','218.754', &
                '219.555','21a.554','225.845','225.855','226.656','227.645','227.655','228.654','229.455','22a.454', &
                '233.955','234.756','235.535','235.745','235.755','236.556','236.655','236.754','237.545','237.555', &
                '238.554','239.355','239.553','23a.354','243.635','243.855','244.656','245.435','245.556','245.635', &
                '245.645','245.655','245.657','246.456','246.555','246.654','247.445','247.455','247.655','247.665', &
                '248.454','252.756','253.535','253.745','253.755','254.546','254.556','254.655','255.535','255.545', &
                '255.555','255.557','255.755','255.765','256.554','257.355','257.555','257.565','257.575','262.656', &
                '263.645','263.655','264.456','264.555','265.445','265.455','265.645','265.655','265.665','265.675', &
                '267.455','267.465','271.557','272.556','273.545','273.555','273.557','274.554','274.556','275.545', &
                '275.555','275.565','275.575','276.554','277.555','283.445','283.655','283.665','285.445','285.455', &
                '285.465','285.475','293.555','293.565','295.355','295.365','295.555','295.565','295.575','2a3.455', &
                '2a5.455','2a5.465','2a5.475','355.555','435.755','445.655','455.555','473.555','491.555','655.555', &
                '855.555'/

  CHARACTER(LEN_FILENAME) :: flneot,flntid

  LOGICAL(IT) :: lfirst,lexist
  DATA lfirst /.TRUE./

  SAVE lfirst,cnmCos,cnmSin,snmCos,snmSin,admittance,doodsonMatrix

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: hex2dec
  INTEGER(IT) :: get_valid_unit


  !*
  ! Start the exectuable code
  !!---------------------------

  IF (lfirst .EQ. .TRUE.) THEN

    lfirst=.FALSE.

    flneot=f_tablefilename('fes14b')

    i=LEN_TRIM(flneot)
    IF (flneot(i:i) .NE. '/') THEN
      flneot=TRIM(flneot)//'/'
    END IF

    DO i=1, ntide
      flntid=TRIM(flneot)//'fes2014b_n180_version20170520.'//TRIM(DARWIN(i))//'.cos.gfc'
      CALL read_icgem_coef(flntid,cnmCos(1:MAXOCNDEG,0:MAXOCNDEG,i),snmCos(1:MAXOCNDEG,0:MAXOCNDEG,i))

      flntid=TRIM(flneot)//'fes2014b_n180_version20170520.'//TRIM(DARWIN(i))//'.sin.gfc'
      CALL read_icgem_coef(flntid,cnmSin(1:MAXOCNDEG,0:MAXOCNDEG,i),snmSin(1:MAXOCNDEG,0:MAXOCNDEG,i))
    END DO

    DO i=1, ndood
      doodsonMatrix(i,1)=hex2dec(doodson(i)(1:1))
      doodsonMatrix(i,2)=hex2dec(doodson(i)(2:2))-5
      doodsonMatrix(i,3)=hex2dec(doodson(i)(3:3))-5
      doodsonMatrix(i,4)=hex2dec(doodson(i)(5:5))-5
      doodsonMatrix(i,5)=hex2dec(doodson(i)(6:6))-5
      doodsonMatrix(i,6)=hex2dec(doodson(i)(7:7))-5
    END DO

    flntid=TRIM(flneot)//'fes2014b_admittance_linear_linear.txt'

    INQUIRE(FILE=flntid,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_fes2014b_tide): '//TRIM(flntid)//' is not exist'
      CALL exit(1)
    END IF

    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=flntid,STATUS='OLD',IOSTAT=ierr)
    IF (ierr .NE. 0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_fes2014b_tide): open '//TRIM(flntid)
      CALL exit(1)
    END IF

    DO i=1, ntide
      READ(lfn,*,ERR=100) (admittance(j,i),j=1,ndood)
    END DO
    CLOSE(lfn)

  END IF


  ! mjd in TT
  CALL fund_arg_nutation(mjd,arg)
  CALL doodson_arg(gmst,arg,beta)
  !CALL doodsonArguments(mjd,gmst,beta)

  ! interpolation matrix (256 tides -> 18 major tides)
  DO i=1, ndood
    thetaf(i)=0.d0
    DO j=1, 6
      thetaf(i)=thetaf(i)+doodsonMatrix(i,j)*beta(j)
    END DO
  END DO

  DO i=1, ntide
    factorCos(i)=0.d0
    factorSin(i)=0.d0
    DO j=1, ndood
      factorCos(i)=factorCos(i)+DCOS(thetaf(j))*admittance(j,i)
      factorSin(i)=factorSin(i)+DSIN(thetaf(j))*admittance(j,i)
    END DO
  END DO


  DO i=1, lmax
    DO j=0, i
      dc(i,j)=0.d0
      ds(i,j)=0.d0
      DO k=1, ntide
        dc(i,j)=dc(i,j)+factorCos(k)*cnmCos(i,j,k)+factorSin(k)*cnmSin(i,j,k)
        ds(i,j)=ds(i,j)+factorCos(k)*snmCos(i,j,k)+factorSin(k)*snmSin(i,j,k)
      END DO
    END DO
  END DO

  RETURN

100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_fes2014b_tide): read error'
  CALL exit(1)

END SUBROUTINE


INTEGER(IT) FUNCTION hex2dec(hex)
!*
!!
!*
USE par
IMPLICIT NONE

!*
! The input argumentes
!!---------------------------
CHARACTER(LEN=*) :: hex

  !*
  ! The local variables
  !!---------------------------

  INTEGER(IT) :: i,dec
  CHARACTER(LEN_STRING) :: chex

  !*
  ! Start the exectuable code
  !!---------------------------

  chex=ADJUSTL(TRIM(hex)) 
 
  hex2dec=0
  DO i=1, LEN_TRIM(chex) 
    SELECT CASE(chex(i:i)) 
      CASE('a':'z') 
        dec=ICHAR(chex(i:i))-ICHAR("a")+10 
      CASE('A':'Z') 
        dec=ICHAR(chex(i:i))-ICHAR("A")+10 
      CASE('0':'9') 
        dec=ICHAR(chex(i:i))-ICHAR("0")
    END SELECT 
    hex2dec=hex2dec*16+dec
  END DO
 
  RETURN

END FUNCTION



SUBROUTINE doodsonArguments(mjd,gmst,dood)
!*
!!
USE const
IMPLICIT NONE

!*
! The input arguments
!!---------------------------
REAL(RL) :: mjd, dood(6)

  !*
  ! The local variables
  !!---------------------------
  REAL(RL) :: Tu0,t,gmst0,gmst,r

  !*
  ! Start the exectuable code
  !!---------------------------

!  Tu0   = (NINT(mjd)-51544.5)/36525.0
!  gmst0 = (6.0/24 + 41.0/(24*60) + 50.54841/(24*60*60))
!  gmst0 = gmst0 + (8640184.812866/(24*60*60))*Tu0
!  gmst0 = gmst0 + (0.093104/(24*60*60))*Tu0*Tu0
!  gmst0 = gmst0 + (-6.2e-6/(24*60*60))*Tu0*Tu0*Tu0
!  r     = 1.002737909350795 + 5.9006e-11*Tu0 - 5.9e-15*Tu0*Tu0
!  gmst  = DMOD(2*PI*(gmst0 + r * DMOD(mjd,1.d0)), 2*PI)

  t = (mjd-51544.5)/365250.0

  dood =0.d0
  dood(2) = (218.31664562999 + (4812678.81195750 + (-0.14663889 + ( 0.00185140 +-0.00015355*t)*t)*t)*t) * pi/180
  dood(3) = (280.46645016002 + ( 360007.69748806 + ( 0.03032222 + ( 0.00002000 +-0.00006532*t)*t)*t)*t) * pi/180
  dood(4) = ( 83.35324311998 + (  40690.13635250 + (-1.03217222 + (-0.01249168 +0.00052655*t)*t)*t)*t) * pi/180
  dood(5) = (234.95544499000 + (  19341.36261972 + (-0.20756111 + (-0.00213942 +0.00016501*t)*t)*t)*t) * pi/180
  dood(6) = (282.93734098001 + (     17.19457667 + ( 0.04568889 + (-0.00001776 +-0.00003323*t)*t)*t)*t) * pi/180
  dood(1) = gmst + PI - dood(2)

  RETURN

END SUBROUTINE
