!*
SUBROUTINE corr_iondel_grid(CKF,sion)
!!
!! shenyixu@whu.edu.cn
!!      correct the ionsphere delay with grid model
!*
USE const
USE ckdctrl
USE ISO_FORTRAN_ENV

!*
! The arguments
!!------------------------
TYPE(CKDCFG) :: CKF
REAL(RL) :: sion(MAXSAT)

    !*
    ! The local variables
    !!------------------------
    LOGICAL(LG) :: lfirst,lexist
    REAL(RL) :: sec
    REAL(RL) :: coef(6),a0(2,MAXSAT)
    REAL(RL) :: a1(2,MAXSAT),a2(2,MAXSAT)
    REAL(RL) :: dt1,dt2,alpha

    INTEGER(IT) :: lfnclk,ierr
    INTEGER(IT) :: i,j,k,iy,imon,id,ih,im
    INTEGER(IT) :: mjdf(2), mjdx
    REAL(RL) :: sodf(2), sodx,dintv
    CHARACTER(LEN_STRING) :: line

    INTEGER(IT) :: nprn(2)
    CHARACTER(LEN_PRN) :: cprn(MAXSAT,2)
    DATA lfirst /.TRUE./
    SAVE lfirst,lfnclk,nprn,cprn,a0,a1,mjdf,sodf


    !*
    ! The function called
    !!-----------------------------
    REAL(RL) :: timdif
    INTEGER(IT) :: modified_julday
    INTEGER(IT) :: get_valid_unit
    INTEGER(IT) :: pointer_string

    !*
    ! Start the exectuable code
    !!-----------------------------

    DO k=1, CKF%nprn
        sion(k)=0.d0
    END DO


END SUBROUTINE