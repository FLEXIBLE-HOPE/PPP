!
!! purpose  : add constraint on bias
!! parameter:
!!    input : CKF -- ckdet configuration
!!            OB -- observation struct
!!    output: NM,infs -- information struct
!! author   : Geng J
!! created  : Nov. 18, 2008
!
SUBROUTINE ppp_con_bias(CKF,PM,NM,infs)
!!
!*
USE info
USE const
USE ckdctrl
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!--------------------
TYPE(CKDCFG) :: CKF
TYPE(PRMT) :: PM(1:*)
TYPE(INFM) :: NM
REAL(RL) :: infs(1:*)

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i
  LOGICAL(LG) :: lfirst
  DATA lfirst /.TRUE./
  SAVE lfirst

  !*
  ! Start the exectuable code
  !!------------------------------

  IF (lfirst .EQ. .FALSE.) RETURN

  !! check size of information matrix
  IF (NM.imtx+NM.nobs+1 .GT. NM.nmtx) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_con_bias): S matrix too small'
    CALL exit(1)
  END IF

  !! add constraints to information matrix (all available sites)
  DO i=1,NM.imtx+1
    infs(NM.iptx(i)+NM.imtx+NM.nobs+1)=0.d0
  END DO

  DO i=1, NM.imtx
    IF (PM(i).pname(1:7).EQ.'RECCLKR' .AND. LEN_TRIM(PM(i).pname).GT.7) THEN
      infs(NM.iptx(i)+NM.imtx+NM.nobs+1)=1.d6
    END IF
  END DO

  NM.nobs=NM.nobs+1
  NM.weig(NM.nobs)=1.d6
  NM.ipob(NM.nobs,1)=0
  NM.ipob(NM.nobs,2)=0

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
  END IF

  RETURN

END SUBROUTINE
