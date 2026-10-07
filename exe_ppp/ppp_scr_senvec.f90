!
!! purpose  : generate sensitivity vectors for bias or blunder
!! parameter:
!!    input : iobs -- which phase observation
!!            wgt  -- weight of this observation
!!            NM   -- information matrix
!!            infs -- one-dimensional array of infmation matrix
!!    output: v    -- sensitivity vector
!! author   : Geng J
!! created  : Nov. 24, 2007
!
SUBROUTINE ppp_scr_senvec(iobs,wgt,v,NM,infs)
!*
USE info
IMPLICIT NONE

!*
! The arguments
!!------------------
TYPE(INFM) :: NM
INTEGER(IT) :: iobs
REAL(RL) :: wgt,v(1:*),infs(1:*)

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: i,j
  REAL(RL) :: gama

  !*
  ! Start the exectuable code
  !!--------------------------

  !! initialization
  DO i=1,NM.imtx+NM.nobs
    v(i)=0.d0
  END DO
  v(iobs+NM.imtx)=wgt

  !! Householder transformation
  DO j=1,NM.imtx
    gama=0.d0
    DO i=j,NM.imtx+NM.nobs
      IF (v(i).EQ.0.d0 .OR. infs(NM.iptx(j)+i+1).EQ.0.d0) CYCLE
      gama=gama+v(i)*infs(NM.iptx(j)+i+1)
    END DO
    IF (gama .EQ. 0.d0) CYCLE
    gama=gama*infs(NM.iptx(j)+NM.imtx+NM.nobs+2)
    DO i=j,NM.imtx+NM.nobs
      IF (infs(NM.iptx(j)+i+1) .EQ. 0.d0) CYCLE
      v(i)=v(i)+gama*infs(NM.iptx(j)+i+1)
    END DO
  END DO

  RETURN

END SUBROUTINE
