!*
SUBROUTINE leo_everett_interpolate(flnorb,npar,pname,lchk,lpos,lvel,lpart,lvpart,mjd,sod,prn,x,v,part,vpart)
!!
!*
USE orbit
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------------
CHARACTER(LEN=*) :: flnorb,prn,pname(1:*)
LOGICAL(LG) :: lchk,lpos,lvel,lpart,lvpart
INTEGER(IT) :: mjd,npar
REAL(RL) :: sod,x(1:*),v(1:*),part(1:*),vpart(1:*)

  !*
  ! The local variables
  !!----------------------
  TYPE(ORBHDR) :: OH
  TYPE(SATEPAR) :: ICS
  TYPE(ORB_INT_TAB) :: EI
  TYPE(ORB_INT_DATA) :: DT(MAXSAT)

  INTEGER(IT) :: iprn,nprn,istart,iend
  CHARACTER(LEN_PRN) :: cprn(MAXSAT)
  INTEGER(IT) :: ierr,ipos,i,j,k,l,ivar,jpos
  INTEGER(IT) :: pt(0:MAXICS,MAXSAT),icomp
  REAL(RL) :: s,u,s2,u2,sum1,sum2

  LOGICAL(LG) lfirst,update_alpha,update_memory
  DATA lfirst,EI.lunit /.TRUE.,0/

  SAVE lfirst,EI,DT,nprn,cprn,pt

  !*
  ! The function called
  !!------------------------
  REAL(RL) :: timdif
  INTEGER(IT) :: pointer_string

  !*
  ! Start of executable code
  !!------------------------

  IF (lfirst .EQ. .TRUE.) THEN

    lfirst=.FALSE.

    CALL rdorbh(flnorb,EI%lunit,OH)
    !! xsy-2023-01-08: a bug
    ! DO i=1, OH.nprn
    !   READ(EI%lunit) ICS
    ! END DO
    
    !! find out requested parameters and PWC parameters
    k=0
    pt=0
    !! there is a bug, if for different satellites with different ics, the interpolate is wrong
    DO iprn=1, OH.nprn
      DO j=1, npar
        pt(j,iprn)=j
      END DO
    END DO

    EI.mjd0=OH.mjd0
    EI.sod0=OH.sod0
    EI.dintv=OH.dintv
    EI.mjd1=OH.mjd1
    EI.sod1=OH.sod1
    EI.nvar=OH.nequ
    !EI.nrec=NINT(timdif(OH.mjd1,OH.sod1,OH.mjd0,OH.sod0)/OH.dintv)+1
    EI.nrec=((OH.mjd1-OH.mjd0)*86400.d0+OH.sod1-OH.sod0)/OH.dintv+1
    nprn=OH.nprn
    DO i=1, OH.nprn
      cprn(i)=OH.cprn(i)
    END DO
    EI.ndgr=6
    IF (EI.nvar .GT. MAXVARS) THEN
      WRITE(ERROR_UNIT,'(A,2I3)') '***ERROR(leo_everett_interpolate): EI.nvar larger than MAXVARS,',EI.nvar,MAXVARS
      CALL exit(1)
    END IF
    IF (EI.ndgr .GT. MAXDGR) THEN
      WRITE(ERROR_UNIT,'(A,2I3)') '***ERROR(leo_everett_interpolate): EI.ndgr larger than MAXDGR',EI.ndgr,MAXDGR
      CALL exit(1)
    END IF

    EI.irec_inmemory=0
    EI.nrec_inmemory=0
    EI.irec_middle_alpha =0
    CALL everett_coeff(MAXDGR,EI.ndgr,EI.ec)
  END IF

  IF(lchk .EQ. .TRUE.) THEN
    mjd=OH.rmjd
    sod=OH.rsod
    RETURN
  END IF

  !! check time
  IF (timdif(mjd,sod,EI.mjd0,EI.sod0).LT.-MAXWND .OR. timdif(mjd,sod,EI.mjd1,EI.sod1).GT.MAXWND) THEN
    WRITE(ERROR_UNIT,'(A,I5,F9.1)') '***ERROR(leo_everett_interpolate): arc not cover epoch ',mjd,sod
    CALL exit(1)
  END IF

  iprn=pointer_string(nprn,cprn,prn)
  IF(iprn .EQ. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(leo_everett_interpolate): no data for satellite '//prn
    CALL exit(1)
  END IF

  !ipos=INT((mjd-EI.mjd0)*(86400.d0/EI.dintv)+(sod-EI.sod0)/EI.dintv)+1
  ipos=((mjd-EI.mjd0)*(86400.d0/EI.dintv)+(sod-EI.sod0)/EI.dintv)+1
  IF(ipos .LE. EI.ndgr) ipos=EI.ndgr+1
  !ELSE IF (ipos+EI.ndgr+1 .GT. EI.nrec) THEN
  !  ipos=EI.irec_middle_alpha
  IF (ipos.GE.EI.nrec-EI.ndgr .AND. ipos.LE.EI.nrec) ipos=EI.nrec-EI.ndgr-1

  update_alpha=(ipos.GT.EI.irec_middle_alpha) .OR. (ipos.LT.EI.irec_middle_alpha-1)

  IF (ipos.LT.EI.irec_middle_alpha-1 .AND. ipos-EI.ndgr.LE.EI.irec_inmemory) THEN
    !! backward to the beginning
    DO i=1,EI.nrec_inmemory
      BACKSPACE(EI.lunit)
    END DO
    !! backward and replace inmemory
    jpos=EI.irec_middle_alpha-ipos+2
    IF (jpos .GE. EI.irec_inmemory) jpos=EI.irec_inmemory
    DO i=1, jpos
      !! move forward memory
      DO l=1, nprn
        DO k=1, EI.nvar
          DO j=EI.nrec_inmemory,2,-1
            DT(l).table(j,k)=DT(l).table(j-1,k)
          END DO
        END DO
      END DO
      BACKSPACE(EI.lunit)
      READ(EI.lunit,END=300) ((DT(l).table(1,j),j=1,EI.nvar),l=1,nprn)
      BACKSPACE(EI.lunit)
    END DO
    EI.irec_inmemory=EI.irec_inmemory-jpos
    EI.irec_middle_alpha=EI.irec_middle_alpha-jpos
    !! forward to the end
    DO i=1,EI.nrec_inmemory
      READ(EI.lunit,END=300)
    END DO
    !! set ipos
    IF (ipos.LT.EI.irec_middle_alpha) THEN
      ipos=EI.irec_middle_alpha
    END IF
  END IF

  update_memory=update_alpha
  DO WHILE(update_memory)
    update_memory=EI.nrec_inmemory.LT.2*EI.ndgr+2 .OR. &
          ipos+EI.ndgr+1.GT.EI.irec_inmemory+EI.nrec_inmemory
    IF (.NOT.update_memory) EXIT
    IF (EI.nrec_inmemory .GE. 2*EI.ndgr+2) THEN
      EI.nrec_inmemory=EI.nrec_inmemory-1
      EI.irec_inmemory=EI.irec_inmemory+1
      DO k=1, nprn
        DO j=1, EI.nvar
          DO i=1, EI.nrec_inmemory
            DT(k).table(i,j)=DT(k).table(i+1,j)
          END DO
        END DO
      END DO
    END IF
    READ(EI.lunit,END=300) ((DT(i).table(EI.nrec_inmemory+1,j),j=1,EI.nvar),i=1,nprn)
    EI.nrec_inmemory=EI.nrec_inmemory+1
  END DO

  !
  !! update alpha and beta
  IF (update_alpha) THEN
    EI.irec_middle_alpha=ipos
    jpos=ipos-EI.irec_inmemory
    DO i=1, nprn
      DO ivar=1, EI.nvar
        DO k=0, EI.ndgr
          DT(i).alpha(ivar,k)=0.d0
          DT(i).beta (ivar,k)=0.d0
          DO j=0, EI.ndgr
            DT(i).alpha(ivar,k)=DT(i).alpha(ivar,k)+EI.ec(k,j)*  &
                (DT(i).table(jpos-j,ivar)+DT(i).table(jpos+j,ivar))
            DT(i).beta (ivar,k)=DT(i).beta (ivar,k)+EI.ec(k,j)*  &
                (DT(i).table(jpos-j+1,ivar)+DT(i).table(jpos+j+1,ivar))
          END DO
        END DO
      END DO
    END DO
  END IF

  !! if missing ..., send a note to upper-level routine
  IF (ALL(DT(iprn).alpha(1,0:EI.ndgr) .EQ. 0.d0)) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(leo_everett_interpolate): missing data'
    CALL exit(1)
  END IF

  s=((mjd-EI.mjd0)*(86400.d0/EI.dintv)+(sod-EI.sod0)/EI.dintv)+1.d0-EI.irec_middle_alpha
  u=s-1.d0
  s2=s*s
  u2=u*u

  IF (EI.nvar.EQ.3 .AND. lpart) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(leo_everett_interpolate): orbit file without partial derivatives'
    CALL exit(1)
  END IF

  istart=1
  iend=EI.nvar
  IF(.NOT. lpos) istart=4
  IF(.NOT. lpart) iend=3
  DO ivar=istart,iend
    l=(ivar-1)/3
    IF (l.NE.0 .AND. pt(l,iprn).EQ.0) CYCLE
    icomp=ivar-3*l
    sum1=DT(iprn).beta(ivar,EI.ndgr)
    sum2=DT(iprn).alpha(ivar,EI.ndgr)
    DO k=EI.ndgr,1,-1
      sum1=sum1*s2+DT(iprn).beta(ivar,k-1)
      sum2=sum2*u2+DT(iprn).alpha(ivar,k-1)
    END DO
    sum1=sum1*s-sum2*u
    IF (l .EQ. 0) THEN
      x(ivar)=sum1
    ELSE
      IF (pt(l,iprn)*3-3+icomp .GT. MAXPWC*MAXPWCEQUS) WRITE(*,*) 'Too small MAXPWC'
      part(pt(l,iprn)*3-3+icomp)=sum1
    END IF
    IF ((l.EQ.0 .AND. lvel) .OR. lvpart) THEN
      sum1=DT(iprn).beta(ivar,EI.ndgr)*(2.d0*EI.ndgr+1.d0)
      sum2=DT(iprn).alpha(ivar,EI.ndgr)*(2.d0*EI.ndgr+1.d0)
      DO k=EI.ndgr,1,-1
        sum1=sum1*s2+DT(iprn).beta(ivar,k-1)*(2.d0*(k-1.d0)+1.d0)
        sum2=sum2*u2+DT(iprn).alpha(ivar,k-1)*(2.d0*(k-1.d0)+1.d0)
      END DO
      IF (l .EQ. 0) THEN
        v(ivar)=(sum1-sum2)/EI.dintv
      ELSE
        IF (pt(l,iprn)*3-3+icomp .GT. MAXPWC*MAXPWCEQUS) WRITE(*,*) 'Too small MAXPWC'
        vpart(pt(l,iprn)*3-3+icomp)=(sum1-sum2)/EI.dintv
      END IF
    END IF
  END DO

  RETURN

300 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(leo_everett_interpolate): end of file '
  CALL exit(1)

ENTRY leo_everett_interpolate_reset()

  lfirst=.TRUE.
  CLOSE(EI.lunit)

  RETURN

END SUBROUTINE
