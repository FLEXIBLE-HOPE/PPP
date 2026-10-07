!*
SUBROUTINE prob_resol(est,sigma,ih,cutdev,cutsig,deci)
!
!!    Calculate the decision-function value as per Appendix A
!!    of Dong and Bock [1989].
!!       Da-nan Dong  880403
!!    Add damping factor. DND  880602
!!    Change definition of decision function based on hypothesis
!!       test theory.      DND 880929
!!    Remove arbitrary factor of 3., scale wide-lane but not narrow-lane
!!       sigmas.    King 930322
!!
!!    Input:
!!      est:     estimated (real) bias value
!!      sigma:   estimated uncertainty of bias value, scaled by nrms for
!!                 widelane, unscaled for narrow lane
!!      ih :     control for receiver ambiguity
!!                 = 1  unit=one cycle
!!                 = 2  unit=half cycle
!!      cutdev:  threshold deviation for taper function (= 0.4 in Dong and Bock;
!!                 default still 0.4 here and in FIXDRV)
!!      cutsig:  threshold sigma for taper function (=0.33 in Dong and Bock;
!!                 default still 0.4 here and in FIXDRV)
!!
!!    Output:
!!      prob : area of decision-function region, set = 1.0 for
!!             each bias on input to BDECI, reduced by the probability
!!             of an error in rounding each bias. (Since we now search
!!             only one bias at a time, the cumulative probability is
!!             not calculated by NBIASR or NBIASP.)
!!      deci : decision function d(x,sigma) = FT/Q, inverse of
!!             1. - allowable rate for a type 1 error (alpha),
!!             compared with input wlcut or nlcut in NBIASR.
!!             (But F is no longer used because we search one bias at
!!             a time.)
!
USE par
IMPLICIT NONE

REAL(RL) :: est,sigma,cutdev,cutsig,prob,dev,deci,term1,term2,c,d1,a1,bint,taper,add,cdev,csig,trun,s1,s2
REAL(RL) :: b1,b2,erfc,erfcb1,erfcb2
INTEGER(IT) :: ih,j
EXTERNAL erfc

DATA s2/1.414213562373095d0/

  !
  !! compute the deviation of the estimated value from an integer or half-integer
  add = est*dble(ih)
  bint= dint(add+0.5d0*dsign(1.d0,add))/dble(ih)
  dev = dabs(bint-est)
  !
  !! set the cutoff deviation from the input or default value
  cdev = cutdev
  !! default was 0.4 prior to 930315, now 0.15
  IF (cutdev .LT. 1.d-3) cdev = 0.15d0
  !! this was 0.6d0*cdev by mistake, prior to 930319
  IF (ih .EQ. 2) cdev = 0.5d0*cdev

  !! if the deviation is greater than the cutoff, set prob and deci and exit
  IF (dev .GT. cdev) THEN
    prob=1.d0
    deci=0.d0
    GOTO 100
  END IF

  !**old code:
  !  scale the estimated sigma by the nrms (sclerr) for both widelane and narrowlane
  !  (these should be treated differently)
  !  sigtemp = sclerr*sigma
  !  s1 = 1.0d0/(sigtemp*s2)
  !**new code
  s1 = 1.d0/(sigma*s2)
  !  numerical truncation tolerance
  trun = 1.d-9

  !  compute the taper (T)
  !  this term is (1 - dev/0.4) in Dong and Bock;
  term1 = 1.d0-dev/cdev
  !  this term is (1. - 3*scaled_sigma) in Dong and Bock; since cutsig is
  !  0.4 now by default, term2 = 1.3 - 3*scaled_sigma
  !  Dong and Bock:  term2 = (.333 - sigtemp)*3.
  !  New:
  csig = cutsig
  !  default changed from 0.4 to 0.15 930319
  IF (cutsig .LT. 1.d-3) csig = 0.15d0
  !**old code:
  !  term2 = (csig-sigtemp)*3.d0
  !**new code
  term2 = (csig-sigma)*3.d0
  !**old code      if (bcigma.lt.1.0d-3) term1 = (0.4d0-c)*3.0d0
  IF (term2 .LT. 0.d0) THEN
    prob=1.d0
    deci=0.d0
    GOTO 100
  END IF
  !  we now square the first term (linear in Dong and Bock) to
  !  achieve a greater taper
  taper = term1**2 * term2

  !  compute Q according to equation A-12 in Dong and Bock
  c=0.d0
  DO j=1,50
    a1 = dble(j)
  !  b1 = sngl((a1-dev)*s1)
  !  b2 = sngl((a1+dev)*s1)
  !  d1 = dble(erfc(b1)-erfc(b2))
    b1 = (a1-dev)*s1
    b2 = (a1+dev)*s1
    !  limit the range of erf to  avoid underflows
    IF(b1.LT.0.d0 .OR. b1.GT.15.d0) THEN
      erfcb1 = 0.d0
    ELSE
      erfcb1 = erfc(b1)
    END IF
    IF (b2.LT.0.d0 .OR. b2.GT.15.d0) THEN
      erfcb2 = 0.d0
    ELSE
      erfcb2 = erfc(b2)
    END IF
    !* d1 = erfc(b1)-erfc(b2)
    d1= erfcb1 - erfcb2
    c = c+d1
    IF (d1 .LT. trun) GOTO 440
  END DO

  ! return the decision function and reduced probability
440 CONTINUE
  prob=1.d0-c
  IF (c .LT. 1.d-9) c=1.d-9
  deci=taper/c

100 CONTINUE
  RETURN

END SUBROUTINE

!--------------------------------------------------------------
REAL(RL) FUNCTION erf(x)
USE par
IMPLICIT NONE

  REAL(RL) gammp,x,half

  half=0.5d0
  IF (x .LT. 0.d0) THEN
    erf = -gammp(half,x**2)
  ELSE
    erf = gammp(half,x**2)
  END IF

  RETURN

END FUNCTION


!--------------------------------------------------------------
REAL(RL) FUNCTION erfc(x)
USE par
IMPLICIT NONE

  REAL(RL) x,gammp,gammq,half

  half=0.5d0
  IF (x .LT. 0.d0) THEN
    erfc = 1.d0+gammp(half,x**2)
  ELSE
    erfc = gammq(half,x**2)
  END IF

  RETURN

END FUNCTION


!--------------------------------------------------------------
REAL(RL) FUNCTION gammp(a,x)
USE par
IMPLICIT NONE

  REAL(RL) :: a,x,gln,gammcf

  IF (x.LT.0.d0 .OR. a.LE.0.d0) PAUSE
  IF (x.LT.a+1.d0) THEN
    CALL gser(gammp,a,x,gln)
  ELSE
    CALL gcf(gammcf,a,x,gln)
    gammp = 1.d0-gammcf
  END IF

  RETURN

END FUNCTION


!--------------------------------------------------------------
REAL(RL) FUNCTION gammq(a,x)
USE par
IMPLICIT NONE

  REAL(RL) a,x,gamser,gln

  IF (x.LT.0.d0 .OR. a.LE.0.d0) PAUSE
  IF (x.LT.a+1.d0) THEN
    CALL gser(gamser,a,x,gln)
    gammq = 1.d0-gamser
  ELSE
    CALL gcf(gammq,a,x,gln)
  END IF

  RETURN

END FUNCTION


!--------------------------------------------------------------
REAL(RL) FUNCTION gammln(xx)
USE par
IMPLICIT NONE

  REAL(RL) :: cof(6),stp,half,one,fpf,x,xx,tmp,ser
  INTEGER(IT) :: j

  DATA cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,-1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
  DATA half,one,fpf/0.5d0,1.0d0,5.5d0/

  x   = xx-one
  tmp = x+fpf
  tmp = (x+half)*log(tmp)-tmp
  ser = one
  DO j=1,6
    x   = x+one
    ser = ser+cof(j)/x
  END DO
  gammln = tmp+LOG(stp*ser)

  RETURN

END FUNCTION


!--------------------------------------------------------------
SUBROUTINE gcf(gammcf,a,x,gln)
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  REAL(RL) :: gammcf,gammln,a,anf,x,gln,gold,a0,a1,b0,b1,fac,an,ana,g,eps
  INTEGER(IT) :: n,itmax
  PARAMETER(itmax=100,eps=3.d-7)


  gln = gammln(a)
  g = 0.d0
  gold = 0.d0
  a0 = 1.d0
  a1 = x
  b0 = 0.d0
  b1 = 1.d0
  fac = 1.d0
  DO n = 1,itmax
    an = float(n)
    ana = an-a
    a0 = (a1+a0*ana)*fac
    b0 = (b1+b0*ana)*fac
    anf = an*fac
    a1 = x*a0+anf*a1
    b1 = x*b0+anf*b1
    IF (a1 .NE. 0.D0) THEN
      fac = 1.d0/a1
      g = b1*fac
      IF (ABS((g-gold)/g) .LT. eps) THEN
        gammcf = dexp(-x+a*dlog(x)-gln)*g
        RETURN
      END IF
      gold = g
    END IF
  END DO

  WRITE(ERROR_UNIT,'(A)') '***ERROR(prob_resol/gcf): a too large, itmax too small '
  !CALL exit(1)

END SUBROUTINE


!--------------------------------------------------------------
SUBROUTINE gser(gamser,a,x,gln)
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  REAL(RL)  gamser,a,x,gln,gammln,ap,sum,del,eps
  INTEGER(IT) n,itmax
  PARAMETER(itmax=100,eps=3.d-7)

  gln = gammln(a)
  IF (x.LE.0.d0) THEN
    IF (x .LT. 0.d0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(prob_resol/gser): x < 0 '
      CALL exit(1)
    END IF
    gamser = 0.d0
    RETURN
  END IF

  ap  = a
  sum = 1.d0/a
  del = sum
  DO n = 1,itmax
    ap = ap+1.d0
    del = del*x/ap
    sum = sum+del
    IF (abs(del).LT.abs(sum)*eps) THEN
      gamser = sum*exp(-x+a*log(x)-gln)
      RETURN
    END IF
  END DO

  WRITE(ERROR_UNIT,'(A)') '***ERROR(prob_resol/gser): a too large, itmax too small'
  CALL exit(1)

END SUBROUTINE
