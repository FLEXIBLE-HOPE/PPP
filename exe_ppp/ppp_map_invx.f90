!
!! purpose  : map inverted normal matrix from Lc to Ln
!! parameter:
!!    output: QM,invx -- inverted normal matrix
!! author   : Jianghui Geng
!! created  : Nov 27 2011
!
SUBROUTINE ppp_map_invx(AM,SAT,AB,QM,invx)
!*
USE info
USE satellite
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!--------------------
TYPE(AMBT) :: AM(1:*)
TYPE(SATE) :: SAT(MAXSAT)
TYPE(AMBD) :: AB(1:*)
TYPE(INVM) :: QM
REAL(RL) :: invx(1:*)

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: i, j, isat, jsat

  !*
  ! Start the exectuable code
  !!--------------------------

  !! transform inverted normal matrix
  ! DO j=1,QM.ntot
  !   IF (j .GT. QM.nxyz) jsat = AM(AB(j-QM.nxyz).pab).psat
  !   DO i=j,QM.ntot
  !     IF (i .GT. QM.nxyz) isat = AM(AB(i-QM.nxyz).pab).psat
  !     IF (j.LE.QM.nxyz .OR. j.GT.QM.nxyz .AND. AB(j-QM.nxyz).abst .EQ.-1) THEN
  !       IF (i.GT.QM.nxyz .AND. AB(i-QM.nxyz).abst .NE.-1) THEN
  !         invx(QM.idq(j)+i)=invx(QM.idq(j)+i)/SAT(isat).lamdn
  !       END IF
  !     ELSE
  !       IF (AB(i-QM.nxyz).abst .EQ. -1) THEN
  !         invx(QM.idq(j)+i)=invx(QM.idq(j)+i)/SAT(jsat).lamdn
  !       ELSE
  !         invx(QM.idq(j)+i)=invx(QM.idq(j)+i)/SAT(isat).lamdn/SAT(jsat).lamdn
  !       END IF
  !     END IF
  !   END DO
  ! END DO

  DO j=1,QM.ntot
    IF (j .GT. QM.nxyz) jsat = AM(AB(j-QM.nxyz).pab).psat
    DO i=j,QM.ntot
      IF (i .GT. QM.nxyz) isat = AM(AB(i-QM.nxyz).pab).psat
      IF (j.LE.QM.nxyz) THEN
        IF (i.GT.QM.nxyz) THEN
          IF (AB(i-QM.nxyz).abst .NE.-1) THEN
            invx(QM.idq(j)+i)=invx(QM.idq(j)+i)/SAT(isat).lamdn
          END IF
        END IF
      ELSE
        IF (AB(j-QM.nxyz).abst .EQ.-1) THEN
          IF (i.GT.QM.nxyz) THEN
            IF (AB(i-QM.nxyz).abst .NE.-1) THEN
              invx(QM.idq(j)+i)=invx(QM.idq(j)+i)/SAT(isat).lamdn
            END IF
          END IF
        ELSE
          IF (AB(i-QM.nxyz).abst .EQ. -1) THEN
            invx(QM.idq(j)+i)=invx(QM.idq(j)+i)/SAT(jsat).lamdn
          ELSE
            invx(QM.idq(j)+i)=invx(QM.idq(j)+i)/SAT(isat).lamdn/SAT(jsat).lamdn
          END IF
        END IF
      END IF
    END DO
  END DO  

  RETURN

END SUBROUTINE
