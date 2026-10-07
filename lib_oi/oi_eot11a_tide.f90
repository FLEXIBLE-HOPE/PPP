!*
SUBROUTINE oi_eot11a_tide(gmst,mjd,lmax,dc,ds)
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
  INTEGER(IT) :: doodsonMatrix(256,6)

  REAL(RL) :: arg(5),beta(6),thetaf(256)
  REAL(RL) :: factorCos(18),factorSin(18)
  REAL(RL) :: admittance(256,18)
  REAL(RL) :: cnmCos(MAXOCNDEG,0:MAXOCNDEG,18)
  REAL(RL) :: cnmSin(MAXOCNDEG,0:MAXOCNDEG,18)
  REAL(RL) :: snmCos(MAXOCNDEG,0:MAXOCNDEG,18)
  REAL(RL) :: snmSin(MAXOCNDEG,0:MAXOCNDEG,18)

  CHARACTER(LEN=3) :: darwin(18)
  CHARACTER(LEN=7) :: doodson(256)

  DATA darwin /'OM1','OM2','SA ', &
               'SSA','MM ','MF ', &
               'MTM','MSQ','Q1 ', &
               'O1 ','P1 ','K1 ', &
               '2N2','N2 ','M2 ', &
               'S2 ','K2 ','M4 '/

  DATA doodson /'055.565','055.575','056.554','056.556','057.355','057.553','057.555','057.565','057.575', &
                '058.554','059.553','062.656','063.645','063.655','063.665','064.456','064.555','065.445', &
                '065.455','065.465','065.655','065.665','065.675','066.454','067.455','067.465','071.755', &
                '072.556','073.545','073.555','073.565','074.556','075.345','075.355','075.365','075.555', &
                '075.565','075.575','076.554','077.355','077.365','081.655','082.656','083.445','083.455', &
                '083.655','083.665','083.675','084.456','085.255','085.455','085.465','085.475','086.454', &
                '091.555','092.556','093.355','093.555','093.565','093.575','095.355','095.365','107.755', &
                '109.555','115.845','115.855','117.645','117.655','118.654','119.455','125.745','125.755', &
                '126.556','126.754','127.545','127.555','128.554','129.355','133.855','134.656','135.435', &
                '135.635','135.645','135.655','135.855','136.555','136.654','137.445','137.455','137.655', &
                '137.665','138.454','139.455','143.535','143.745','143.755','144.546','144.556','145.535', &
                '145.545','145.555','145.755','145.765','146.554','147.355','147.555','147.565','148.554', &
                '153.645','153.655','154.656','155.435','155.445','155.455','155.645','155.655','155.665', &
                '155.675','156.555','156.654','157.445','157.455','157.465','158.454','161.557','162.556', &
                '163.545','163.555','163.755','164.554','164.556','165.545','165.555','165.565','165.575', &
                '166.554','167.355','167.555','167.565','168.554','172.656','173.445','173.645','173.655', &
                '173.665','174.456','174.555','175.445','175.455','175.465','175.655','175.665','175.675', &
                '176.454','182.556','183.545','183.555','183.565','185.355','185.365','185.555','185.565', &
                '185.575','191.655','193.455','193.465','193.655','193.665','195.255','195.455','195.465', &
                '195.475','207.855','209.655','215.955','217.755','219.555','225.855','227.645','227.655', &
                '228.654','229.455','234.756','235.745','235.755','236.556','236.655','236.754','237.545', &
                '237.555','238.554','239.355','243.635','243.855','244.656','245.435','245.645','245.655', &
                '246.456','246.555','246.654','247.445','247.455','247.655','248.454','253.535','253.755', &
                '254.556','255.535','255.545','255.555','255.557','255.755','255.765','256.554','257.355', &
                '257.555','257.565','257.575','262.656','263.645','263.655','264.555','265.445','265.455', &
                '265.655','265.665','265.675','267.455','267.465','271.557','272.556','273.545','273.555', &
                '273.557','274.554','274.556','275.545','275.555','275.565','275.575','276.554','277.555', &
                '283.655','283.665','285.455','285.465','285.475','293.555','293.565','295.355','295.365', &
                '295.555','295.565','295.575','455.555'/

  CHARACTER(LEN_FILENAME) :: flneot,flntid

  LOGICAL(IT) :: lfirst,lexist
  DATA lfirst /.TRUE./

  SAVE lfirst,cnmCos,cnmSin,snmCos,snmSin,admittance,doodsonMatrix

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: get_valid_unit


  !*
  ! Start the exectuable code
  !!---------------------------

  IF (lfirst .EQ. .TRUE.) THEN

    lfirst=.FALSE.

    flneot=f_tablefilename('eot11a')

    i=LEN_TRIM(flneot)
    IF (flneot(i:i) .NE. '/') THEN
      flneot=flneot//'/'
    END IF

    DO i=1, 18
      flntid=TRIM(flneot)//'eot11a.'//TRIM(DARWIN(i))//'.cos.gfc'
      CALL read_icgem_coef(flntid,cnmCos(1:MAXOCNDEG,0:MAXOCNDEG,i),snmCos(1:MAXOCNDEG,0:MAXOCNDEG,i))

      flntid=TRIM(flneot)//'eot11a.'//TRIM(DARWIN(i))//'.sin.gfc'
      CALL read_icgem_coef(flntid,cnmSin(1:MAXOCNDEG,0:MAXOCNDEG,i),snmSin(1:MAXOCNDEG,0:MAXOCNDEG,i))
    END DO

    DO i=1, 256
      READ(doodson(i)(1:1),*) doodsonMatrix(i,1)
      READ(doodson(i)(2:2),*) doodsonMatrix(i,2)
      doodsonMatrix(i,2)=doodsonMatrix(i,2)-5
      READ(doodson(i)(3:3),*) doodsonMatrix(i,3)
      doodsonMatrix(i,3)=doodsonMatrix(i,3)-5
      READ(doodson(i)(5:5),*) doodsonMatrix(i,4)
      doodsonMatrix(i,4)=doodsonMatrix(i,4)-5
      READ(doodson(i)(6:6),*) doodsonMatrix(i,5)
      doodsonMatrix(i,5)=doodsonMatrix(i,5)-5
      READ(doodson(i)(7:7),*) doodsonMatrix(i,6)
      doodsonMatrix(i,6)=doodsonMatrix(i,6)-5
    END DO

    flntid=TRIM(flneot)//'admittance.linear.txt'

    INQUIRE(FILE=flntid,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_eot11a_tide): '//TRIM(flntid)//' is not exist'
      CALL exit(1)
    END IF

    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=flntid,STATUS='OLD',IOSTAT=ierr)
    IF (ierr .NE. 0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_eot11a_tide): open '//TRIM(flntid)
      CALL exit(1)
    END IF

    DO i=1, 18
      READ(lfn,*,ERR=100) (admittance(j,i),j=1,256)
    END DO
    CLOSE(lfn)

  END IF

  ! mjd in TT
  CALL fund_arg_nutation(mjd,arg)
  CALL doodson_arg(gmst,arg,beta)

  ! interpolation matrix (256 tides -> 18 major tides)
  DO i=1, 256
    thetaf(i)=0.d0
    DO j=1, 6
      thetaf(i)=thetaf(i)+doodsonMatrix(i,j)*beta(j)
    END DO
  END DO

  DO i=1, 18
    factorCos(i)=0.d0
    factorSin(i)=0.d0
    DO j=1, 256
      factorCos(i)=factorCos(i)+DCOS(thetaf(j))*admittance(j,i)
      factorSin(i)=factorSin(i)+DSIN(thetaf(j))*admittance(j,i)
    END DO
  END DO


  DO i=1, lmax
    DO j=0, i
      dc(i,j)=0.d0
      ds(i,j)=0.d0
      DO k=1, 18
        dc(i,j)=dc(i,j)+factorCos(k)*cnmCos(i,j,k)+factorSin(k)*cnmSin(i,j,k)
        ds(i,j)=ds(i,j)+factorCos(k)*snmCos(i,j,k)+factorSin(k)*snmSin(i,j,k)
      END DO
    END DO
  END DO

  RETURN

100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_eot11a_tide): read error'
  CALL exit(1)

END SUBROUTINE


!*
SUBROUTINE read_icgem_coef(flntid,cnm,snm)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------------
CHARACTER(LEN=*) :: flntid
REAL(RL) :: cnm(MAXOCNDEG,0:MAXOCNDEG)
REAL(RL) :: snm(MAXOCNDEG,0:MAXOCNDEG)

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: i,j
  INTEGER(IT) :: lfn,ierr,max_deg
  REAL(RL) :: coef(2),EARTH_GM,EARTH_R

  CHARACTER(LEN_STRING) :: line 

  LOGICAL(LG) :: lexist

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!---------------------------

  INQUIRE(FILE=flntid,EXIST=lexist)
  IF (lexist .EQ. .FALSE.) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(read_icgem_coef): '//TRIM(flntid)//' is not exist'
    CALL exit(1)
  END IF

  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=flntid,STATUS='OLD',IOSTAT=ierr)
  IF (ierr .NE. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(read_icgem_coef): open '//TRIM(flntid)
    CALL exit(1)
  END IF

  DO WHILE(.TRUE.)
    READ(lfn,'(A)',ERR=100,END=200) line
    IF (line(1:22) .EQ. 'earth_gravity_constant') THEN
      READ(line(23:),*) EARTH_GM
    ELSE IF (line(1:6) .EQ. 'radius') THEN
      READ(line(7:),*) EARTH_R
    ELSE IF (line(1:10) .EQ. 'max_degree') THEN
      READ(line(11:),*) max_deg
      IF (max_deg .GT. MAXOCNDEG) THEN
        !WRITE(OUTPUT_UNIT,*) '###MESSAGE(read_icgem_coef): max_degree is greater than MAXOCNDEG, hence using max_degree ', max_deg
        max_deg=MAXOCNDEG
      END IF
    ELSE IF (line(1:11) .EQ. 'end_of_head') THEN
      EXIT
    END IF
  END DO

  cnm=0.d0
  snm=0.d0

  DO WHILE(.TRUE.)
    READ(lfn,'(A)',ERR=100,END=200) line
    !! note the format description to destinguish different coefficients
    READ(line(5:),*) i,j,coef(1),coef(2)
    !! do not read 0, as it is 1 or 0 for gravity and other gravitations
    IF (i .EQ. 0) CYCLE
    IF (i .GT. max_deg) CYCLE
    cnm(i,j)=coef(1)
    snm(i,j)=coef(2)
  END DO

200 CLOSE(lfn)

  RETURN
  
100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(read_icgem_coef): '//TRIM(line)
  CALL exit(1)

END SUBROUTINE
  
