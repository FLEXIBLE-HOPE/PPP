!*
SUBROUTINE get_ant_ipt(ctype,nfreq,freq,fjd_beg,fjd_end,antnam,antnum,iptatx,enu)
!!
!! purpose   : get antenna pointer in the whole table for receivers and satellites
!! parameter :
!!    input  : fjd_beg,fjd_end -- time span
!!             antnam,antnum   -- antenna name and serial number
!!    output : iptatx  -- pointer to global table
!!             enu     -- antenna offset    pco
!! author    : Geng J
!! created   : Oct. 2, 2007
!!
!*
USE atx
USE const
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: iptatx
REAL(RL) :: fjd_beg,fjd_end,enu(1:*)
CHARACTER(LEN_ANTENNA) :: antnam,antnum
CHARACTER(LEN=*) :: ctype
INTEGER(IT) :: nfreq(MAXSYS)
CHARACTER(LEN_FREQ) :: freq(MAXFREQ,MAXSYS)


  !*
  ! The local variables
  !!--------------------------
  INTEGER(IT) :: ript,sipt,isys
  INTEGER(IT) :: i,j,natx,izen,iazi
  REAL(RL) :: zeni,azim,nadir,var(MAXFREQ,MAXSYS)
  REAL(RL) :: alpha,beta,zen,azi,nad,x1,x2,dvar(0:MAXPCVDEG,MAXFREQ,MAXSYS)
  TYPE(ANTATX) :: AX(MAXSAT+MAXSIT),ATXE,ANT

  DATA natx/0/
  SAVE natx,AX

  !*
  ! Start the exectuable code
  !!---------------------------

  DO i=1, natx
    IF (antnam .NE. AX(i).antnam) CYCLE
    j=LEN_TRIM(antnum)
    IF (j.NE.0 .AND. antnum(1:j).NE.AX(i).antnum(1:j)) CYCLE
    iptatx=i
    EXIT
  END DO
  IF (iptatx .NE. 0) GOTO 5

  !! if not in memory, read from file
  ATXE.antnam=antnam
  ATXE.antnum=antnum
  CALL rdatx(nfreq,freq,fjd_beg,fjd_end,ATXE)
  natx=natx+1
  AX(natx)=ATXE
  iptatx=natx

  !! get antenna phase offset
5 CONTINUE
  IF (TRIM(ctype) .EQ. 'SITE') THEN
    DO isys=1, MAXSYS
      DO i=1, MAXFREQ
        !! Please be carefull this part
        ! if no GPS?
        IF ((isys.NE.INDEX(SYS,'G')) .AND. (LEN_TRIM(AX(iptatx).freq(i,isys)).EQ.0)) THEN
          enu(1+(i-1)*3+(isys-1)*3*MAXFREQ)=AX(iptatx).neu(2,i,INDEX(SYS,'G'))
          enu(2+(i-1)*3+(isys-1)*3*MAXFREQ)=AX(iptatx).neu(1,i,INDEX(SYS,'G'))
          enu(3+(i-1)*3+(isys-1)*3*MAXFREQ)=AX(iptatx).neu(3,i,INDEX(SYS,'G'))
        ELSE
          enu(1+(i-1)*3+(isys-1)*3*MAXFREQ)=AX(iptatx).neu(2,i,isys)
          enu(2+(i-1)*3+(isys-1)*3*MAXFREQ)=AX(iptatx).neu(1,i,isys)
          enu(3+(i-1)*3+(isys-1)*3*MAXFREQ)=AX(iptatx).neu(3,i,isys)
        END IF
      END DO
    END DO
  ELSE
    isys=INDEX(SYS,ctype(1:1))
    DO i=1, MAXFREQ
      enu(1+(i-1)*3)  =AX(iptatx).neu(1,i,isys)
      enu(2+(i-1)*3)  =AX(iptatx).neu(2,i,isys)
      enu(3+(i-1)*3)  =AX(iptatx).neu(3,i,isys)
    END DO
  END IF

  RETURN

!
!! purpose   : get antenna pcv for receivers and satellites
!! parameter :
!!    input  : ript,sipt -- pointer of receiver and satellite
!!             zeni  -- zenith angle in radian
!!             azim  -- azimuth angle in radian
!!             nadir -- nadir angle in radian
!!    output : var   -- phase center variation
!! author    : Geng J
!! created   : Oct. 2, 2007
!
ENTRY get_ant_pcv(ript,sipt,zeni,azim,nadir,var,dvar)

  var=0.d0
  dvar=0.d0

  !
  !! receiver pcv
  zen=zeni*RAD2DEG
  azi=azim*RAD2DEG
  IF (zen .GT. AX(ript).zen2) zen=AX(ript).zen2
  IF (zen .LT. AX(ript).zen1) zen=AX(ript).zen1
  IF (azi .LT.0.d0) azi=azi+360.d0

  !! azimuth dependent
  iazi=0
  IF (AX(ript).dazi.ne.0.d0) iazi=int(azi/AX(ript).dazi)+1

  !! zenith dependent
  izen=INT((zen-AX(ript).zen1)/AX(ript).dzen)+1

  ! If there is no PCVs avaliable for particular system, the GPS's will be used
  DO isys=1, MAXSYS
    DO i=1, MAXFREQ
      j=isys
      IF (LEN_TRIM(AX(ript).freq(i,isys)) .EQ. 0) j=INDEX(SYS,'G')
      x1=AX(ript).pcv(izen  ,iazi,i,j)
      x2=AX(ript).pcv(izen+1,iazi,i,j)
      IF (iazi .NE. 0) THEN
        alpha=azi/AX(ript).dazi-iazi+1
        x1=x1+(AX(ript).pcv(izen  ,iazi+1,i,j)-AX(ript).pcv(izen  ,iazi,i,j))*alpha
        x2=x2+(AX(ript).pcv(izen+1,iazi+1,i,j)-AX(ript).pcv(izen+1,iazi,i,j))*alpha
      END IF
      beta=(zen-AX(ript).zen1)/AX(ript).dzen-izen+1
      var(i,isys)=var(i,isys)+x1+(x2-x1)*beta

      !IF (LEN_TRIM(AX(ript).freq(i,isys)) .EQ. 0) CYCLE
      !x1=AX(ript).pcv(izen  ,iazi,i,isys)
      !x2=AX(ript).pcv(izen+1,iazi,i,isys)
      !IF (iazi .NE. 0) THEN
      !  alpha=azi/AX(ript).dazi-iazi+1
      !  x1=x1+(AX(ript).pcv(izen  ,iazi+1,i,isys)-AX(ript).pcv(izen  ,iazi,i,isys))*alpha
      !  x2=x2+(AX(ript).pcv(izen+1,iazi+1,i,isys)-AX(ript).pcv(izen+1,iazi,i,isys))*alpha
      !END IF
      !alpha=(zen-AX(ript).zen1)/AX(ript).dzen-izen+1
      !var(i,isys)=var(i,isys)+x1+(x2-x1)*alpha
    END DO
  END DO

  !! satellite pcv
  IF (AX(sipt).dzen .EQ. 0.d0) RETURN

  nad=nadir*RAD2DEG
  IF (nad .GT. AX(sipt).zen2) nad=AX(sipt).zen2
  IF (nad .LT. AX(sipt).zen1) nad=AX(sipt).zen1

  !! only zenith dependent
  izen=INT((nad-AX(sipt).zen1)/AX(sipt).dzen)+1
  
  !! xsy：izen会超限
  IF (izen .GT. MAXPCVDEG) RETURN

  DO isys=1, MAXSYS
    DO i=1,MAXFREQ
      IF (LEN_TRIM(AX(sipt).freq(i,isys)) .EQ. 0) CYCLE
      x1=AX(sipt).pcv(izen  ,0,i,isys)
      x2=AX(sipt).pcv(izen+1,0,i,isys)
      alpha=(nad-AX(sipt).zen1)/AX(sipt).dzen-izen+1
      var(i,isys)=var(i,isys)+x1+(x2-x1)*alpha
      dvar(izen-1,i,isys)=1.0-alpha
      dvar(izen,i,isys)=alpha
    END DO
  END DO

  RETURN

ENTRY get_atx(iptatx,ANT)

  IF (iptatx .LE. natx) THEN
    ANT=AX(iptatx)
  ELSE
    iptatx=-1
  END IF

  RETURN

END SUBROUTINE
