!*
SUBROUTINE MATMPY(A,B,C,NROW,NCOLA,NCOLB)
!!
!! PREMULTIPLY MATRIX B BY MATRIX A WITH RESULTS IN MATRIX C
!!
!! PARAMETERS
!!      A       I  INPUT MATRIX
!!      B       I  INPUT MATRIX
!!      C       O  OUTPUT MATRIX
!!      NROW    I  NUMBER OF ROWS OF MATRIX A
!!      NCOLA   I  NUMBER OF COLUMNS OF MATRIX A
!!      NCOLB   I  NUMBER OF COLUMNS OF MATRIX B
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-------------------------
INTEGER(IT) :: NROW,NCOLA,NCOLB
REAL(RL) :: A(NROW,NCOLA),B(NCOLA,NCOLB),C(NROW,NCOLB)

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: I,J,K
  REAL(RL) :: SUM
  REAL(RL), POINTER :: C1(:,:)

  !*
  ! Start the exectuable code
  !!--------------------------

  ALLOCATE(C1(1:NROW,1:NCOLB))

  DO I=1,NROW
    DO J=1,NCOLB
      SUM = 0.D0
      DO k=1,NCOLA
        SUM = SUM + A(I,k)*B(k,J)
      END DO
      C1(I,J) = SUM
    END DO
  END DO

  DO I=1,NROW
    DO J=1,NCOLB
      C(I,J)=C1(I,J)
    END DO
  END DO

  DEALLOCATE(C1)

  RETURN

END SUBROUTINE
