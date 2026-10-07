!!
SUBROUTINE leo_interpolate_orbit(flnorb,mjd,sod,isat,cprn,npar,pname,lpos,lvel,lpart,x,v,npwc,partial,dpartial)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------------
CHARACTER(LEN=*) :: flnorb,cprn,pname(1:*)
INTEGER(IT) :: isat,mjd,npar,npwc(1:*)
REAL(RL) :: sod,x(1:*),v(1:*),partial(1:*),dpartial(1:*)
LOGICAL(LG) :: lpos,lvel,lpart

  !*
  ! The local variables
  !!-----------------------------
  INTEGER(IT) :: i,j,k,ii,jj,iics
  INTEGER(IT) :: mjd0,ipwc,mjdi
  REAL(RL) :: sod0,dt,sodi,det
  REAL(RL) :: xhelp(6),xhelp0(6),phik(6,6)

  INTEGER(IT) :: mpar(MAXSAT,MAXPWC),ipt(MAXSAT,0:MAXPWC),intv(MAXSAT,MAXPWC)
  REAL(RL) :: part(MAXEQUS),dpart(MAXEQUS)
  REAL(RL) :: mm(MAXSAT,MAXPWC,MAXPWCEQUS)


  LOGICAL(LG) :: lfirst,ldone,lsat(MAXSAT)
  DATA lfirst,lsat /.TRUE.,MAXSAT*.TRUE./
  SAVE lfirst,lsat,mjd0,sod0,mpar,ipt,mm,intv

  !*
  ! Start of executable code
  !!------------------------


  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    ipt=0
    CALL leo_everett_interpolate(flnorb,npar,pname,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., &
                       mjd0,sod0,cprn,x,v,part,dpart)
  END IF

  !! initialization for PWC parameters
  IF (lsat(isat) .EQ. .TRUE.) THEN
    lsat(isat)=.FALSE.
    ipwc=0
    DO i=1, npar
      j=INDEX(pname(i),':')
      IF (j .EQ. 0) CYCLE
      ipwc=ipwc+1
      IF (ipwc .GT. MAXPWC) WRITE(*,*) 'Too small MAXPWC'
      ipt(isat,ipwc)=i
      READ(pname(i)(j+1:),*) intv(isat,ipwc)
      mpar(isat,ipwc)=1
      DO k=1, MAXPWCEQUS
        mm(isat,ipwc,k)=0.d0
      END DO

      IF(i .LE. 6) mm(isat,ipwc,i)=1.d0
    END DO
    ipt(isat,0)=ipwc

  END IF

  !! check IF mm must be updated
  ldone=.FALSE.
  DO ipwc=1,ipt(isat,0)
    iics=ipt(isat,ipwc)
    dt=(mjd-mjd0)*1440+(sod-sod0)/60.d0
    IF (dt .LT. -MAXWND) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(leo_interpolate_orbit): PWC parameters require data should start at reference time'
      CALL exit(1)
    END IF

    !! how many parameters for this PWC
    !npwc(ipwc)=INT((dt-MAXWND)/intv(isat,ipwc))+1
    npwc(ipwc)=INT((dt-0.01d0)/intv(isat,ipwc))+1

    !! IF it is the same as what is saved
    IF(npwc(ipwc) .LE. mpar(isat,ipwc)) CYCLE

    !! other wise new parameter comes for this PWC, so we have to calculate MM
    mpar(isat,ipwc)=npwc(ipwc)
    mjdi=mjd0
    sodi=sod0+(npwc(ipwc)-1)*intv(isat,ipwc)*60.d0

    !! for this epoch we only have to DO the interpolation once  (inv(phi(ti,t0)))
    IF (.NOT. ldone) THEN
      !CALL leo_everett_interpolate(flnorb,npar,pname,.FALSE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.,mjdi,sodi,cprn,x,v,part,dpart)
      CALL leo_everett_interpolate(flnorb,npar,pname,.FALSE.,.FALSE.,.FALSE.,lpart,lpart,mjdi,sodi,cprn,x,v,part,dpart)

      ldone=.TRUE.
      DO ii=1,3
        DO jj=1,6
          phik(ii,jj)=part((jj-1)*3+ii)
          phik(ii+3,jj)=dpart((jj-1)*3+ii)
        END DO
      END DO
      CALL matinv(phik,6,6,det)
    END IF

    !!  d(x,v)/dp(ipwc) = PHIp(ti,t0), iics <= 6 for velocity or impulse
    DO ii=1,3
      IF (iics .GT. 6) THEN
        k=(iics-1)*3
        xhelp(ii)=part(k+ii)
        xhelp(ii+3)=dpart(k+ii)
      ELSE
        xhelp(ii)=0.d0
        xhelp(ii+3)=0.d0
      END IF
    END DO
    IF(iics .LE. 6) xhelp(iics)=1.d0

    !! mm(ti) = inv[PHIx(ti,t0)]*PHIp(ti,t0)
    CALL matmpy(phik,xhelp,xhelp,6,6,1)

    !! the position of MM
    DO ii=1,6
      mm(isat,ipwc,npwc(ipwc)*6-6+ii)=xhelp(ii)
    END DO
  END DO

  !! interpolating
  CALL leo_everett_interpolate(flnorb,npar,pname,.FALSE.,lpos,lvel,lpart,ipt(isat,0).NE.0,mjd,sod,cprn,x,v,part,dpart)

  IF(.NOT. lpart) RETURN

  DO ii=1,3
    DO jj=1,6
      phik(ii,jj)=part((jj-1)*3+ii)
      phik(ii+3,jj)=dpart((jj-1)*3+ii)
    END DO
  END DO

  !! partial
  ipwc=0
  jj=0
  DO i=1, npar
    j=INDEX(pname(i),':')
    IF (j .NE. 0) THEN
      ipwc=ipwc+1
      DO j=1,mpar(isat,ipwc)
        jj=jj+1
        DO k=1,6
          IF (i .GT. 6) THEN
            xhelp(k)=mm(isat,ipwc,j*6+k)-mm(isat,ipwc,(j-1)*6+k)
          ELSE
            xhelp(k)=mm(isat,ipwc,(j-1)*6+k)
          END IF
        END DO
        CALL matmpy(phik,xhelp,xhelp,6,6,1)
        DO k=1,3
          partial((jj-1)*3+k)=xhelp(k)
          dpartial((jj-1)*3+k)=xhelp(k+3)
        END DO
      END DO

      !! the last one of the PWC
      !! PHIp(t,ti) = PHIp(t,t0) - PHIx(t,t0)*(inv PHIx(ti,t0) * PHIp(ti,t0))
      IF(i .GT. 6) THEN
        DO k=1, 3
          partial((jj-1)*3+k)=part((i-1)*3+k)+xhelp(k)
          dpartial((jj-1)*3+k)=dpart((i-1)*3+k)+xhelp(k+3)
        END DO
      END IF
    ELSE
      jj=jj+1
      DO k=1,3
        partial((jj-1)*3+k)=part((i-1)*3+k)
        dpartial((jj-1)*3+k)=dpart((i-1)*3+k)
      END DO
    END IF
  END DO

  RETURN

END SUBROUTINE
