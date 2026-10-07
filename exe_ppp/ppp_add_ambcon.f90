!
!! purpose  : impose already fixed ambiguities
!! parameter:
!!    input : PM -- parameter list
!!            AB -- ambiguity list (auxiliary)
!!    output: QM -- inverted normal matrix
!! author   : Jianghui Geng
!! created  : Nov 23 2011
!
SUBROUTINE ppp_add_ambcon(PM,AB,QM,invx)
!!
!*
USE info
USE ambiguity
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(PRMT) :: PM(1:*)
TYPE(AMBD) :: AB(1:*)
TYPE(INVM) :: QM
REAL(RL) :: invx(1:*)

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: i,j,nfix,nred
  REAL(RL) :: dump
  REAL(RL), POINTER :: ainc(:),pinc(:)
  REAL(RL), POINTER :: q22(:,:),q12(:,:),q21(:,:),qhelp(:,:)

  !*
  ! The function called
  !!-------------------------
  REAL(RL) :: dot

  !*
  ! Start the exectuable code
  !!--------------------------

  !! count fixed ambiguities
  nfix=0        !xsy: 固定的模糊度数量
  nred=0        !xsy: 多余的模糊度数量
  ALLOCATE(ainc(QM.ntot-QM.nxyz))
  ainc(1:QM.ntot-QM.nxyz)=0.d0
  DO i=1,QM.ntot-QM.nxyz
    IF (AB(i).abst .GT. 0) THEN             !xsy: abst=1、2,固定或者伪固定
      nfix=nfix+1
      ! ambiguity increment
      ainc(nfix)=AB(i).abfx-AB(i).abfr      !xsy: 整数模糊度和浮点模糊度偏差
    ELSE IF (AB(i).abst .EQ. -1) THEN
      ! number of redundant ambiguities
      nred=nred+1                           !xsy: 多余模糊度数量（剩下的）
    END IF
  END DO

  !! derive amb-fixed estimates
  IF (nfix .NE. 0) THEN

    !! Q22 inversed matrix，分块
    ALLOCATE(q22(nfix,nfix))
    q22(1:nfix,1:nfix)=0.d0                 !xsy: 已固定模糊度的Qxx阵
    DO j=QM.ntot-nfix+1,QM.ntot
      DO i=QM.ntot-nfix+1,QM.ntot
        IF (i .LT. j) THEN
          q22(i-QM.ntot+nfix,j-QM.ntot+nfix)=invx(QM.idq(i)+j)
        ELSE
          q22(i-QM.ntot+nfix,j-QM.ntot+nfix)=invx(QM.idq(j)+i)
        END IF
      END DO
    END DO
    CALL matinv(q22,nfix,nfix,dump)
    IF (dump .EQ. 0.d0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_add_ambcon): matrix singularity'
      CALL exit(1)
    END IF

    !! Q12 and Q21 matrix，分块
    ALLOCATE(q12(QM.ntot-nfix,nfix))            
    ALLOCATE(qhelp(nfix,QM.ntot-nfix))          
    ALLOCATE(q21(QM.ntot-nfix,QM.ntot-nfix))
    q12(1:QM.ntot-nfix,1:nfix)=0.d0
    qhelp(1:nfix,1:QM.ntot-nfix)=0.d0
    q21(1:QM.ntot-nfix,1:QM.ntot-nfix)=0.d0
    DO j=1,QM.ntot-nfix
      DO i=QM.ntot-nfix+1,QM.ntot
        q12(j,i-QM.ntot+nfix)=invx(QM.idq(j)+i)
      END DO
    END DO
    qhelp=transpose(q12)
    ! q12 becomes q12*q22^-1
    CALL matmpy(q12,q22,q12,QM.ntot-nfix,nfix,nfix)
    ! q21 becomes q12*q22^-1*q21
    CALL matmpy(q12,qhelp,q21,QM.ntot-nfix,nfix,QM.ntot-nfix)

    !! update estimates，np、nc、ns
    ALLOCATE(pinc(QM.ntot))
    pinc(1:QM.ntot)=0.d0
    CALL matmpy(q12,ainc,pinc,QM.ntot-nfix,nfix,1)
    DO i=1,QM.ntot-nfix
      IF (i .LE. QM.nxyz) THEN
        PM(i).xcor=PM(i).xcor+pinc(i)   !固定后，非模糊度参数-改正数-矫正值
        PM(i).xest=PM(i).xest+pinc(i)   !固定后，非模糊度参数-估值  -矫正值
      ELSE
        AB(i-QM.nxyz).abfr=AB(i-QM.nxyz).abfr+pinc(i)   !浮点模糊度修正，hold则作为下个历元模糊度值
      END IF
    END DO

    !! update vtpv
    !  call matmpy(ainc,q22,pinc,1,nfix,nfix)
    !  QM.vtpv=QM.vtpv+dot(nfix,pinc,ainc)

    !! update q11
    DO j=1,QM.ntot-nfix
      DO i=j,QM.ntot-nfix
        invx(QM.idq(j)+i)=invx(QM.idq(j)+i)-q21(i,j)
      END DO
    END DO

    !! clean memory
    DEALLOCATE(pinc)
    DEALLOCATE(qhelp)
    DEALLOCATE(q12)
    DEALLOCATE(q22)
    DEALLOCATE(q21)

    !! Decrease ntot, ndam after excluding nfix
    QM.ntot=QM.ntot-nfix            !xsy: 将固定的模糊度从法方程中剔除
    QM.ndam=QM.ntot-QM.nxyz-nred    !xsy: 剩下的没有固定的模糊度，剔除多余模糊度(abst=-1,没有fcb或者截至高度角或者新的模糊度)
    QM.nfix=nfix                    !xsy: 本次固定的模糊度数量
  ELSE
    QM.nfix=0
    QM.ndam=QM.ntot-QM.nxyz-nred
  END IF
  
  DEALLOCATE(ainc)

  RETURN

END SUBROUTINE
