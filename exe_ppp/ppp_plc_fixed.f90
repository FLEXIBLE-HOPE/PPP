!
!! purpose  : sort inversed normal matrix (lower triangular part)
!! parameter:
!!    input : AB -- ambiguities
!!    output: QM,invx -- inversed normal matrix,Qxx
!! author   : Geng J
!! created  : Jan. 30, 2008
!
SUBROUTINE ppp_plc_fixed(AB,QM,invx)
!*
USE info
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!----------------------
TYPE(AMBD) :: AB(1:*)
TYPE(INVM) :: QM
REAL(RL) :: invx(1:*)

  !*
  ! The local variables
  !!-----------------------
  TYPE(AMBD) :: AX
  INTEGER(IT) :: i,j,k,ifnd
  REAL(RL) :: dump

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! place fixed ambiguities at the end of AB
  DO i=QM.nxyz+1,QM.ntot
    IF (AB(i-QM.nxyz).abst .LE. 0) CYCLE ! redundant or unfixed amb
    ifnd=1              ! find a fixed
    DO j=i+1,QM.ntot
      IF (AB(j-QM.nxyz).abst .LE. 0) THEN
        ifnd=0          ! find a redundant or unfixed
        AX=AB(i-QM.nxyz)
        AB(i-QM.nxyz)=AB(j-QM.nxyz)
        AB(j-QM.nxyz)=AX

        !! move elements in QM
        DO k=1,i-1                         ! column
          dump=invx(QM.idq(k)+i)
          invx(QM.idq(k)+i)=invx(QM.idq(k)+j)
          invx(QM.idq(k)+j)=dump
        END DO
        DO k=j+1,QM.ntot                   ! row
          dump=invx(QM.idq(i)+k)
          invx(QM.idq(i)+k)=invx(QM.idq(j)+k)
          invx(QM.idq(j)+k)=dump
        END DO
        DO k=i+1,j-1                       ! row
          dump=invx(QM.idq(i)+k)
          invx(QM.idq(i)+k)=invx(QM.idq(k)+j)
          invx(QM.idq(k)+j)=dump
        END DO
        dump=invx(QM.idq(i)+i)
        invx(QM.idq(i)+i)=invx(QM.idq(j)+j)
        invx(QM.idq(j)+j)=dump
      END IF
    END DO
    IF (ifnd .NE. 0) EXIT   !xsy: 后面找不到没固定的则停止
  END DO

  RETURN

END SUBROUTINE
