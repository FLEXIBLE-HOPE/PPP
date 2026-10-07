!!
SUBROUTINE matinv(a,ndim,n,d)
!!
!! purpose    :  inversion of a matrix by means of gauss....
!!
!! parameters : a : matrix to be inverted
!!           ndim : dimension of a matrix
!!              n : dimension of a to be invert
!!              d : value of determinant a
!!
!! author     : Ge Maorong
!!              Tsinghua University
!!              Beijing 100084
!!
!! created    : 1989-12-20
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: ndim,n
REAL(RL) :: a(ndim,ndim),d

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) :: i,j,k,l(ndim),m(ndim)
  REAL(RL) :: max_a,swap

  !*
  ! Start the exectuable code
  !!--------------------------

  d=1.d0
  !
  !! loop over all row
  DO k=1,n
    l(k)=k
    m(k)=k
    max_a=a(k,k)
    !
    !! find out the diag. element with max. abs. value.
    DO i=k,n
      DO j=k,n
        IF (dabs(max_a)-dabs(a(i,j)) .LT. 0.d0) THEN
          max_a=a(i,j)
          m(k)=i
          l(k)=j
        END IF
      END DO
    END DO

    !! rank defection or bad condition
    IF (dabs(max_a) .LT. 1d-13) THEN
      d=0.d0
      RETURN
    END IF

    !! swap k and l(k) column
    IF (l(k) .GT. k) THEN
      DO i=1,n
        swap=-a(i,k)
        a(i,k)=a(i,l(k))
        a(i,l(k))=swap
      END DO
    END IF

    !! swap k and m(k) raw
    IF (m(k) .GT. k) THEN
      DO j=1,n
        swap=-a(k,j)
        a(k,j)=a(m(k),j)
        a(m(k),j)=swap
      END DO
    END IF

    !! elemination
    DO i=1,n
      IF (i .NE. k) a(k,i)=-a(k,i)/max_a
    END DO

    DO i=1,n
      DO j=1,n
        IF (i.NE.k .AND. j.NE.k) a(j,i)=a(j,i)+a(k,i)*a(j,k)
      END DO
    END DO

    !! save invert part
    DO j=1,n
      IF (j .NE. k) a(j,k)=a(j,k)/max_a
    END DO
    a(k,k)=1/max_a
  END DO

  !
  !! re-range raw and column as input
  DO k=n,1,-1
    IF (l(k) .GT. k) THEN
      DO j=1,n
        swap=a(k,j)
        a(k,j)=-a(l(k),j)
        a(l(k),j)=swap
      END DO
    END IF

    IF (m(k) .GT. k) THEN
      DO i=1,n
        swap=a(i,k)
        a(i,k)=-a(i,m(k))
        a(i,m(k))=swap
      END DO
    END IF
  END DO

  RETURN

END SUBROUTINE
