!*
SUBROUTINE proc_xl_dirc(nxl,rxl,wgt,ifg,ndl,fxl,vxl,sxl)
!!
!! purpose  : process wide-lane directional data
!! parameter:
!!    input : nxl -- number of wide-lane fractional parts
!!            rxl -- fractional parts
!!            wgt -- weight
!!            ifg -- flag
!!    output: ndl -- number of deleted fractional parts
!!            fxl -- estimate of fractional parts
!!            vxl -- rms
!!            sxl -- sigma
!! author   : Geng J
!! created  : Apr 4, 2009
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: nxl,ndl,ifg(1:*)
REAL(RL) :: rxl(1:*),wgt(1:*),fxl,vxl,sxl

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: i,j,k,ndel,flg(nxl),loc
  REAL(RL) :: gap(3),bwl(nxl),sav(nxl),mean,rms,sig

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! compute fractional parts
  DO i=1,nxl
    bwl(i)=rxl(i)-NINT(rxl(i))
  END DO

  !! loop for each array elements
  ndl=0
  vxl=1000.d0
  sxl=0.d0
  IF (fxl .EQ. 10.d0) THEN
    DO i=1,nxl
      rxl(i)=bwl(i)
      DO j=1, nxl
        IF (j .EQ. i) CYCLE
        gap(1)=bwl(j)-rxl(i)
        gap(2)=bwl(j)+1.d0-rxl(i)
        gap(3)=bwl(j)-1.d0-rxl(i)
        ! minimum difference
        loc=MINLOC(dabs(gap),1)     !xsy: 数组中最小值的位置
        rxl(j)=gap(loc)+rxl(i)      !xsy: 和最近整数的偏差(±0.5)
      END DO
      !! get weighted mean
      flg(1:nxl)=0
      CALL fcb_wgt_mean(.TRUE.,rxl,flg,wgt,nxl,ndel,mean,rms,sig)
      !DO WHILE(nxl-ndel.GT.2 .AND. rms.GT.0.2d0)
      !  k=ndel
      !  CALL sign_robust(nxl,rxl,flg,0.3d0,k)
      !  write(*,*) rxl(1:nxl)
      !  IF (ndel .EQ. k) EXIT
      !  CALL get_wgt_mean(.FALSE.,rxl,flg,wgt,nxl,ndel,mean,rms,sig)
      !END DO
      IF (rms .LT. vxl) THEN
        sav(1:nxl)=rxl(1:nxl)
        ifg(1:nxl)=flg(1:nxl)
        ndl=ndel
        fxl=mean
        vxl=rms
        sxl=sig
      END IF
    END DO
    rxl(1:nxl)=sav(1:nxl)
  !! ambiguity fixed
  ELSE
    DO j=1,nxl
      gap(1)=bwl(j)-fxl
      gap(2)=bwl(j)+1.d0-fxl
      gap(3)=bwl(j)-1.d0-fxl
      ! minimum difference
      loc=MINLOC(dabs(gap),1)
      rxl(j)=gap(loc)+fxl
    END DO
    ifg(1:nxl)=0
    CALL get_wgt_mean(.TRUE.,rxl,ifg,wgt,nxl,ndl,fxl,vxl,sxl)
  END IF

  RETURN

END SUBROUTINE
