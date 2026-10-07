!*
SUBROUTINE fcb_wgt_mean(l_edit,x,flg,wgt,n,ndel,mean,rms,sig)
!!
!! purpose  : compute weighted mean value of an array
!!            & sign-constrained least squares
!! parameter:
!!    input : l_edit -- edit data or not
!!            x,flg  -- data array & flag array
!!            wgt    -- weight array
!!            n      -- # of data
!!    output: ndel -- # of deleted data
!!            mean -- mean value
!!            rms  -- unweighted rms of residuals
!!            sig  -- sigma of mean value
!! author   : Ge M, Geng J
!! revised  : Jan 21, 2008; Aug 26, 2009
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!----------------------
LOGICAL(LG) :: l_edit
REAL(RL) :: x(1:*),sig,rms,mean,hp,wgt(1:*)
INTEGER(IT) :: flg(1:*)
INTEGER(IT) :: n,ndel

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i,m,ndel_new,imax
  REAL(RL) :: wgt_sum,maxres

  !*
  ! Start the exectuable code
  !!--------------------------

  ndel=0
  ndel_new=1
  DO WHILE(ndel_new .NE. 0)

    !! mean
    m=0
    mean=0.d0
    wgt_sum=0.d0
    DO i=1, n
      IF (flg(i) .GE. 2) CYCLE
      mean=mean+wgt(i)*x(i)
      m=m+1
      wgt_sum=wgt_sum+wgt(i)
    END DO
    IF (m .NE. 0) THEN
      mean=mean/wgt_sum
    ELSE
      rms=999.d0
      sig=999.d0
      ndel=n
      RETURN
    END IF

    !! sigma & standard deviation
    sig=0.d0
    rms=0.d0
    DO i=1, n
      IF (flg(i) .GE. 2) CYCLE
      rms=rms+(x(i)-mean)**2*wgt(i)
    END DO
 
    IF (m .EQ. 1) THEN
      ndel=n-1
      rms=999.d0
      sig=999.d0
      RETURN
    ELSE
      rms=DSQRT(rms/(m-1))
      sig=rms/DSQRT(wgt_sum)
    END IF

    !! if achieve goal
    IF (.NOT.l_edit .OR. 3.d0*rms.LT.0.05d0) THEN
      ndel=COUNT(flg(1:n) .GE. 2)
      RETURN
    END IF

    !! edit data according to residuals
    ndel_new=0
    ndel=0

    imax=0
    maxres=0.d0
    DO i=1, n
      IF (flg(i) .GE. 2) THEN
        ndel=ndel+1
      ELSE IF (DABS(x(i)-mean) .GT. maxres) THEN
        imax=i
        maxres=DABS(x(i)-mean)
      END IF
    END DO
    IF (rms.GE.0.2d0 .OR. maxres.GE.0.3d0 .OR. maxres.GE.3.d0*rms/wgt(i)) THEN
      ndel=ndel+1
      ndel_new=ndel_new+1
      flg(imax)=2
    END IF

  END DO

  RETURN

END SUBROUTINE
