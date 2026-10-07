!*
SUBROUTINE chos_sort(ndel,imax,iset,NM)
!!
!! purpose  :  sort and candidate selection
!! parameter:
!!    input   ndel -- number of candidates to be selected
!!            imax -- maximum candidates that can be selected(相位数)
!!    output: iset -- index of to-be-selected candidates
!! author   : shengyi xu
!! created  : 2022/3/31
!!
!*
USE info

IMPLICIT NONE

!*
! The arguments
!!-------------------------------
INTEGER(IT) :: ndel,imax,iset(1:*)
TYPE(INFM) :: NM

  !*
  ! The local variables
  !!---------------------
  INTEGER(IT) :: i,j,K,sort(1:imax)
  REAL(RL) :: resi(1:imax),res(1:imax),temp

  !*
  ! Start the exectuable code
  !!---------------------------

  IF (imax .LT. ndel) RETURN
  DO i=1,imax
    resi(i) = NM.resi(i*2-1)
    !write(*,*)NM.resi(2*i-1)
    resi(i) = ABS(resi(i))
  END DO
  res(1:imax)=resi(1:imax)           !res存储原来相位的顺序
  DO i=imax-1,1,-1
    DO j=1,i
      IF ( resi(j) < resi(j+1) ) THEN
        temp=resi(j)
        resi(j)=resi(j+1)
        resi(j+1)=temp     !resi存储排序后的顺序
      END IF
    END DO
  END DO
  k=0
  DO i=1,imax
    DO j=1,imax
       IF (resi(i) .EQ. res(j))THEN
         k=k+1
         sort(k)=j
       END IF
    END DO
  END DO
  !write(*,*)sort(1:imax)
  iset(1:ndel)=sort(1:ndel)
  iset(ndel+1:imax)=0

  RETURN

END SUBROUTINE

