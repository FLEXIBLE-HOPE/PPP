!
!! purpose   : prediction of SRIF for p and x parameters
!! parameters: NM,PM -- Information matrix & parameters
!!             infs  -- one-dimensional information matrix
!! auther    : Geng J
!! reference : Bierman GJ  1972. pp 142 (3.2)
!!
!! OPTIMIZATION: Use LAPACK DGEMM for matrix operations
!!

SUBROUTINE fcb_tim_update(CKF,PM,NM,infs)
!!
!*
USE info
USE ckdctrl
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
LOGICAL(LG) :: lfirst
INTEGER(IT) :: i,ir,ic,j,k,iw,minut,iupc
REAL(RL) :: s,beta,gama,v(MAXPARSIT),u(MAXPARSIT)

!*
! Start the exectuable code
!!---------------------------

!! predict process parameters
DO i=1,NM.np
  PM(i).zw=0.d0
  iw=INDEX(PM(i).pname,':')
  IF (iw .NE. 0) THEN
    READ(PM(i).pname(iw+1:),*) minut
    IF ((CKF.mjd-PM(i).ptime(2))*86400.d0+CKF.sod .LT. -MAXWND) CYCLE
  END IF

  !! for receiver clock
  IF (PM(i).pname(4:6) .EQ. 'CLK') THEN
    PM(i).xini=PM(i).xini+PM(i).map*PM(i).xcor
  ELSE
    !! for ztd
    PM(i).xini=PM(i).map*PM(i).xest
  END IF
  PM(i).zw=-PM(i).rw*PM(i).map*PM(i).xcor
END DO

!! Householder transformation to pre-eliminate process parameters
DO i=1,NM.np
  iupc=0
  iw=INDEX(PM(i).pname,':')
  IF (iw .NE. 0) THEN
    READ(PM(i).pname(iw+1:),*) minut
    IF ((CKF.mjd-PM(i).ptime(2))*86400.d0+CKF.sod .LT. -MAXWND) iupc=1
  END IF
  IF (iupc .NE. 0) THEN
    !! just change rows
    DO j=1,NM.np
      v(j)=infs(NM.iptx(1)+j)
    END DO
    DO j=2,NM.np
      DO k=1,NM.np
        infs(NM.iptx(j-1)+k)=infs(NM.iptx(j)+k)
      END DO
    END DO
    DO j=1,NM.np
      infs(NM.iptx(NM.np)+j)=v(j)
    END DO
    CYCLE
  ELSE
    !! elemental transformation vector
    PM(i).iobs=0
    v(1)=-PM(i).rw*PM(i).map
    DO j=1,NM.np
      v(j+1)=infs(NM.iptx(1)+j)
    END DO
    CALL fcb_ele_hht(NM.np+1,v,s,u,beta)
    !! apply to next columns
    DO j=2,NM.imtx+1
      gama=0.d0
      DO k=1,NM.np
        IF (infs(NM.iptx(j)+k).EQ.0.d0 .OR. u(k+1).EQ.0.d0) CYCLE
        gama=gama+infs(NM.iptx(j)+k)*u(k+1)
      END DO
      IF (j .EQ. NM.imtx+1) gama=gama+PM(i).zw*u(1)
      gama=gama*beta
      ic=j
      IF (j .LE. NM.np) ic=j-1
      DO k=1,NM.np
        IF (ic.EQ.j .AND. u(k+1).EQ.0.d0) CYCLE
        infs(NM.iptx(ic)+k)=infs(NM.iptx(j)+k)+gama*u(k+1)
      END DO
    END DO
    !! move the first column to the end of np
    gama=u(1)*PM(i).rw*beta
    DO k=1,NM.np
      infs(NM.iptx(NM.np)+k)=gama*u(k+1)
    END DO
  END IF
  !! change time tags for ztd
  IF (iw .NE. 0) THEN
    PM(i).ptime(1)=PM(i).ptime(1)+minut/1440.d0
    PM(i).ptime(2)=PM(i).ptime(2)+minut/1440.d0
  END IF
END DO

RETURN

END SUBROUTINE
