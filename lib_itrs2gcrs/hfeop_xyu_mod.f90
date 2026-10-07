!*
MODULE hfeop_xyu
!*
!! Hold tidal information.
!! Maximum tides is 200.
!! So far the model with the largest is JPL which has ~160
!! History
!! 2017Nov22 JMGipson. Original version
!! 2019Oct09 Two bugs found and fixed by Michael Gerstl. Rotines calc_gmst and calc_nut_arg.
!!           Effect is ~1.0 uas in PM, ~0.01 us in UT1
!*
USE par
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  LOGICAL(LG) :: khfeop_xyu
  INTEGER(IT), PARAMETER :: max_tide=200
  ! number we find   
  INTEGER(IT) :: num_tide
  ! tidal arguments. GST+pi followed by nutation.
  INTEGER(IT) :: itide_arg(6,max_tide)
  CHARACTER(LEN=8) :: ctide(max_tide)
  CHARACTER(LEN=8) :: cdoodson(max_tide)
  ! period
  REAL(RL) :: tide_period(max_tide)
  ! Coefficients of model in order: (sin,cosine)x(X,Y,UT1,LOD)x Max_tide
  REAL(RL) :: tide_coef(2,4,max_tide)

  SAVE num_tide,itide_arg,ctide,cdoodson,tide_period,tide_coef

  CONTAINS

  SUBROUTINE import_tides_xyu(lhfeop_file,ilu)
  USE par
  USE const
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

  !*
  ! The arguments
  !!---------------------------
  ! File that contains the data 
  CHARACTER(LEN=*) :: lhfeop_file

    !*
    ! The local variables
    !!---------------------------
    INTEGER(IT) :: ilu
    ! temporary buffer to hold input line.  
    CHARACTER(LEN_STRING) :: ldum
    LOGICAL(LG) :: kexist

    !*
    ! The funcation called
    !!--------------------------
    INTEGER(IT) :: get_valid_unit

    !*
    ! Start the exectuable code
    !!--------------------------
  
    itide_arg=0
    tide_period=0.d0
    tide_coef=0.d0

    INQUIRE(EXIST=kexist,FILE=lhfeop_file)
    IF (kexist .EQ. .FALSE.) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(import_tides_xyu): did not find HF_EOP_XYU file '//TRIM(lhfeop_file)
      CALL exit(1)
    END IF
  
    num_tide=0
    ilu=get_valid_unit(10)
    OPEN(UNIT=ilu, file=lhfeop_file)

    ! Skip any lines that  begin with "#"
    ! Get the number of components in model.
    DO WHILE(num_tide .LT. max_tide)
      ! read in a line, exit when reach EOF 
      READ(ilu,'(A)',END=100) ldum
      IF (ldum(1:1) .NE. "#") THEN
        num_tide=num_tide+1  
        READ(ldum,*,ERR=110) ctide(num_tide),itide_arg(1:6,num_tide),cdoodson(num_tide), &
                tide_period(num_tide), tide_coef(1:2,1:4,num_tide)
      END IF 
    END DO 
    WRITE(OUTPUT_UNIT,'(A,I5)') '###MESSAGE(import_tides_xyu): ran out of space in import_iers_tides ', num_tide
      
100 CONTINUE    
    CLOSE(ilu)        
    
    RETURN 
      
110 CONTINUE
    WRITE(ERROR_UNIT,'(A)') "***ERROR(import_tides_xyu): error reading in line: "//TRIM(ldum)
    CALL exit(1)

  END SUBROUTINE


  SUBROUTINE tide_angles(rmjd_TT,Delta_T,knew_order,fund_arg)
  !*
  ! Calculate the tide angles
  !*
  USE par
  USE const
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE


  !*
  ! The arguments
  !!---------------------------
  ! rmjd (TT)
  REAL(RL) :: rmjd_TT
  ! TT-UT1 (seconds)
  REAL(RL) :: Delta_T
  ! True: Put the arguments in the order GMST+pi, fund_arg. (2010 IERS conventions.)
  ! False: Order is nutation, GMST+pi
  LOGICAL(LG) :: knew_order
  ! Argument and their time derivatives. 
  REAL(RL) :: fund_arg(6,2)

    !*
    ! The local variables
    !!---------------------------
    ! Fundamental arguments of nutation. Second index is time derivative.
    REAL(RL) :: nut_arg(5,2)
    ! GMST. Second index is time derivative. 
    REAL(RL) :: gmst(2)
    ! Modified julian date TT.
    REAL(RL) :: RMJD_UT

    !*
    ! Start the exectuable code
    !!---------------------------

    rmjd_UT = rmjd_TT-Delta_T/86400.d0
    CALL calc_gmst(rmjd_UT,gmst)
           
    CALL calc_nut_arg(rmjd_tt,nut_arg)
    
    ! Add pi to GMST 
    gmst(1)=DMOD(gmst(1)+PI,2.d0*PI)

    ! If knew_order=.true. then put GMST+pi as the first argument
    !              =.false. use 1996 IERS order where GMST+pi was last. 

    ! knew_order=.false. 
    IF (knew_order .EQ. .TRUE.) THEN
      fund_arg(1,1:2)=gmst
      fund_arg(2:6,1:2)=nut_arg 
    ELSE
      fund_arg(6,1:2)=gmst
      fund_arg(1:5,1:2)=nut_arg
    END IF

    RETURN

  END SUBROUTINE tide_angles
     
  SUBROUTINE calc_gmst(RMJD_UT,gmst)
  !*
  ! Compute GMST
  USE par
  USE const
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

  !*
  ! The arguments
  !!---------------------------
  REAL(RL) :: RMJD_UT
  ! GMST and it's derivative. 
  REAL(RL) :: gmst(2)
  
    !*
    ! The local variables
    !!---------------------------
    ! Do the calculation in a way that is more numerically stable 
    LOGICAL(LG) :: kuse_mod

    !UT Julian centuries since J2000. 
    REAL(RL) ::  T

    REAL(RL) :: tmp_big, tmp_small
    REAL(RL), PARAMETER :: sec_per_circ=1296000d0

    ! coefficients of expansion for GMST 
    ! An alternate form of co(1,2) is (876600d0*3600d0 + 8640184.812866d0)
    ! Below this is expanded.  
    REAL(RL) :: co(4)
    DATA co/67310.54841d0,   3164400184.812866d0, 0.093104d0, -6.2d-6/

    ! Rewrite the second coefficient in another way. 
    !    co(2)=cot(1,j)*sec_per_circ+cot(2,j)
    REAL(RL) :: cot(2)
    DATA cot/2442.d0, -431815.18734d0/ 

    !*
    ! Start the exectuable code
    !!---------------------------

    t=(RMJD_UT-51544.5d0)/36525.d0        
  
    kuse_mod=.TRUE.
    IF (kuse_mod .EQ. .TRUE.) THEN
      tmp_big=co(1)+DMOD(cot(1)*t,1.d0)*sec_per_circ+cot(2)*t
      tmp_big=DMOD(tmp_big,sec_per_circ)
    ELSE 
      tmp_big=co(1)+co(2)*T       
    END IF 

    tmp_small = co(3)*T**2 + co(4)*T**3   ! correction. Michael Gerstl

    ! Convert from time-seconds to arc-seconds
    tmp_big=tmp_big*15.d0
    tmp_small=tmp_small*15.d0

    ! If we are the boundary, may overflow
    gmst(1)=DMOD(tmp_big+tmp_small,sec_per_circ)
    IF (gmst(1) .LT. 0.d0) gmst(1)=gmst(1)+sec_per_circ
    ! convert to radians   
    gmst(1)=DMOD(gmst(1),sec_per_circ)*ARCSEC2RAD

    gmst(2)=co(2)+2.d0*co(3)*T+3.d0*co(4)*T**2
    ! convert to radians/sec
    gmst(2)=gmst(2)*15*ARCSEC2RAD/36525.d0/86400.d0
    
    RETURN

  END SUBROUTINE calc_gmst


  SUBROUTINE calc_nut_arg(RMJD_TT, arg)
  !*
  ! Compute GMST
  USE par
  USE const
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

  !*
  ! The arguments
  !!---------------------------
  ! Modified Julian Day (TT)
  REAL(RL) :: RMJD_TT 
  REAL(RL) :: arg(5,2)

    !*
    ! The local variables
    !!---------------------------
    ! use the mod calculation 
    LOGICAL(LG) :: kuse_mod 
    REAL(RL) :: t
    REAL(RL), PARAMETER :: sec_per_circ=1296000d0
    REAL(RL) :: tmp_big, tmp_small
    REAL(RL) :: tmp
    INTEGER(IT) :: i

    ! Written by JMGipson
    !   2017Oct02. 
    !
    ! Compute the fundamental arguments and their derivatives.
    ! Order is: l, lp, f,d, omega 

    ! The coefficients for the data. 
    ! The order is constant, T, T^2, T^3, T^4

    ! The values in the table below come from the routine fundarg.f 
    !  http://iers-conventions.obspm.fr/2010/2010_official/chapter8/software/FUNDARG.F

    ! the version below expands this out.

    REAL(RL) :: co(5,5), cot(2,5)
    DATA co /485868.249036d0, 1717915923.2178d0,  31.8792d0, 0.051635d0,  -0.00024470d0, &
             1287104.79305d0, 129596581.0481d0,   -0.5532d0, -0.000136d0, 0.00001149d0,  &
             335779.526232d0, 1739527262.8478d0, -12.7512d0, -0.001037d0, 0.00000417d0,  &
             1072260.70369d0, 1602961601.2090d0,  -6.3706, 0.006593d0, -0.00003169d0,    &
             450160.398036d0, -6962890.2665d0,    0.007702d0, 7.4722d0, -0.00005939d0 /

    ! The values in COT are related to co(2,*) by
    !    co(2,j)=cot(1,j)*sec_per_circ+cot(2,j)
    DATA cot / 1326.d0, -580076.722d0,  &
                100.d0, -3418.9519d0,   &
               1342.d0,  295262.8478d0, &
               1237.d0, -190398.791d0,  &
                 -5.d0, -482890.2665d0/

    !*
    ! Start the exectuable code
    !!--------------------------

    t=(RMJD_TT-51544.5)/36525.d0   
    !      write(*,*) "NEW CENT ", T
    ! The general formula is:
    !       arg(j,1)=co(1,j)+co(2,j)*T+Co(3,j)*T**2+co(4,j)*T**3+co(5,j)*T**4  
    !       arg(j,2)=co(2,j)+co(3,j)*T+Co(4,j)*T**2+Co(5,j)*T**3


    ! If this is true, do the calculation in a way that is slighltly more numerically stable. 
    kuse_mod=.TRUE.

    ! do it the following way for numerical stability.
    DO i=1,5
      IF (kuse_mod .EQ. .TRUE.) THEN
        tmp_big=co(1,i)+DMOD(cot(1,i)*T,1.d0)*sec_per_circ+cot(2,i)*T
        tmp_big=DMOD(tmp_big,sec_per_circ)
      ELSE
        tmp_big=co(1,i)+co(2,i)*T       
      END IF

      tmp_small = co(3,i)*T**2 + co(4,i)*T**3 + co(5,i)*T**4   ! Michael Gerstl

      ! Conversion from time-seconds to angle-seconds 
      arg(i,1)=DMOD(tmp_big+tmp_small,sec_per_circ)    !If we are the boundary, may overflow
      IF (arg(i,1) .LT. 0.d0) arg(i,1)=arg(i,1)+sec_per_circ 

      arg(i,1)=DMOD(arg(i,1),sec_per_circ)*ARCSEC2RAD     !convert to radians
 
      ! Now do the rates. 
      arg(i,2)=co(2,i)+2.d0*co(3,i)*T+3.d0*co(4,i)*T**2+4*co(5,i)*T**3
      !         arg(i,2)=co(2,i)    
      ! Conversion from time-seconds to angle-seconds  
      arg(i,2)=arg(i,2)*ARCSEC2RAD/36525.d0/86400.d0    !Convert to radians/per sec      
    END DO
   
    RETURN

  END SUBROUTINE calc_nut_arg

      
      
  FUNCTION dotarg(iarg,angles)
  !*
  !*
  IMPLICIT NONE

  !*
  ! The arguments
  !!---------------------------  
  ! multipliers of fundametal  arguments. 
  INTEGER(IT) :: iarg(6)
  ! Values of GMST+pi, 5 fundamental arguments
  REAL(RL) :: angles(6)
  REAL(RL) :: dotarg

    !*
    ! The local variables
    !!---------------------------
    INTEGER(IT) :: i

    !*
    ! Start the exectuable code
    !!--------------------------

    dotarg=0.d0
    DO i=1,6
      dotarg=dotarg+iarg(i)*angles(i)
    END DO
    
    RETURN

  END FUNCTION dotarg

END MODULE
