!
!! purpose  : sort inversed normal matrix (lower triangular part)
!! parameter:
!!    input : AB -- ambiguities
!!    output: QN,invx -- inversed normal matrix
!! author   : Geng J
!! created  : Jan. 30, 2008
!

SUBROUTINE ppp_plc_float(AB,QN,invx)
!*
USE info
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!-----------------------
TYPE(AMBD) :: AB(1:*)
TYPE(INVM) :: QN
REAL(RL) :: invx(1:*)

  !*
  ! The local variables
  !!-------------------------
  TYPE(AMBD) :: AX
  INTEGER(IT) :: i,j,k
  REAL(RL) :: dump

  !*
  ! Start the exectuable code
  !!---------------------------

  !! variance of ambiguities
  DO i=QN.nxyz+1,QN.ntot
    IF (AB(i-QN.nxyz).abst .EQ. 0) THEN
      AB(i-QN.nxyz).stad=invx(QN.idq(i)+i)
    ELSE IF (AB(i-QN.nxyz).abst .EQ. -1) THEN
      AB(i-QN.nxyz).stad=1.d10 ! to be placed ahead
    END IF
  END DO

  !! for unfixed amb, place amb with smaller 'stad' ahead       Qxx对角线降序排列
  DO i=QN.nxyz+1,QN.ntot
    DO j=i+1,QN.ntot
      IF (AB(i-QN.nxyz).stad .LT. AB(j-QN.nxyz).stad) THEN
        AX=AB(i-QN.nxyz)
        AB(i-QN.nxyz)=AB(j-QN.nxyz)
        AB(j-QN.nxyz)=AX

        !! move elements in QN
        DO k=1,i-1                         ! column
          dump=invx(QN.idq(k)+i)
          invx(QN.idq(k)+i)=invx(QN.idq(k)+j)
          invx(QN.idq(k)+j)=dump
        END DO
        DO k=j+1,QN.ntot                   ! row
          dump=invx(QN.idq(i)+k)
          invx(QN.idq(i)+k)=invx(QN.idq(j)+k)
          invx(QN.idq(j)+k)=dump
        END DO
        DO k=i+1,j-1                       ! row
          dump=invx(QN.idq(i)+k)
          invx(QN.idq(i)+k)=invx(QN.idq(k)+j)
          invx(QN.idq(k)+j)=dump
        END DO
        dump=invx(QN.idq(i)+i)
        invx(QN.idq(i)+i)=invx(QN.idq(j)+j)
        invx(QN.idq(j)+j)=dump
      END IF
    END DO
  END DO

  RETURN

END SUBROUTINE
