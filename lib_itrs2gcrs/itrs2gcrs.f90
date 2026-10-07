
SUBROUTINE itrs2gcrs(conv,mjdgps,sodgps,utcut1r,xhelp,tmat,dlmat,dxmat,dymat,mjdut1,sodut1,gmst,xp,yp)
!*
!
!  The transformation matrix to be computed to relate the International Terrestrial
!  Reference System (ITRS) to the Geocentric Celestial Reference System (GCRS) at the
!  date t of the observation (in GPST) according to IERS Conventions
!
!  This routine is part of the Position And Navigation Data Analyst (@PANDA)software.
!
!  Reference:
!
!     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
!     IERS Technical Note No. 32, BKG
!
!  This revision:  2009 April 1
!
!  PANDA release 2013-12-01
!
!  Copyright (C) 2013 Jing Guo, GNSS Research Center, Wuhan University
!
!-----------------------------------------------------------------------------------------------------
USE const
IMPLICIT NONE

!*
! The arguments
!!---------------------
CHARACTER(LEN=*) :: conv
INTEGER(IT) :: mjdgps,mjdut1
REAL(RL) :: sodgps,sodut1,tmat(3,3)
REAL(RL) :: dlmat(3,3),dxmat(3,3),dymat(3,3)
REAL(RL) :: utcut1r,xhelp(2)
REAL(RL) :: xp,yp,gmst

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: ierr
  INTEGER(IT) :: mjdtt,mjdutc,mjdtai
  REAL(RL) :: sodtt,sodutc,sodtai,ut12tt
  REAL(RL) :: erp(3),ut1r(3),ocean(3),libr(2)
  REAL(RL) :: QMAT(3,3),RMAT(3,3),WMAT(3,3),RSMA(3,3)
  REAL(RL) :: RXMA(3,3),RYMA(3,3),MATT(3,3),era,sp

  !*
  ! The function called
  !!---------------------------
  REAL(RL) :: iau_SP00,iau_ERA00,iau_GST94
  REAL(RL) :: iau_GMST06,iau_GMST00,iau_GMST82,iau_GST06A,iau_GST00A

  !*
  ! Start the exectuable code
  !!---------------------------

  ! Time system transformation
  CALL timinc(mjdgps,sodgps,OFF_GPS2TAI,mjdtai,sodtai)
  CALL timinc(mjdtai,sodtai,OFF_TAI2TT,mjdtt,sodtt)


  !! Please be careful, the iau_TAIUTC has truncation error to trans the TAI to UTC
  CALL taiutc(mjdtai,sodtai,mjdutc,sodutc)

  ! Interpolate the Earth Rotation Parameters (xpole, ypole, ut1r(s))
  !CALL tableLinearInterpolate('poleut1',.FALSE.,mjdutc+sodutc/86400.d0,erp)
  erp(1)=xhelp(1)
  erp(2)=xhelp(2)
  erp(3)=utcut1r

  ! Add the long-term of ut1 to get ut1r
  ut1r=0.d0
  CALL RG_ZONT2(conv,(mjdtt+sodtt/86400.d0-51544.50)/36525.D0, ut1r(1), ut1r(2), ut1r(3))

  ! JG: carefully, using utc2tt is better than that of ut12tt, it will caused rotation
  ut12tt=(mjdtt-mjdutc)*86400.d0+sodtt-sodutc
  CALL HF_EOP(mjdtt+sodtt/86400.d0,ut12tt,ocean)
  ocean=ocean*1.D-6

  ! libration in polar motion and ut1
  !CALL PMSDNUT2(mjdutc+sodutc/86400.d0,libr)
  CALL PMSDNUT2(mjdtt+sodtt/86400.d0,libr)
  libr = libr*1.D-6

  !CALL UTLIBR(mjdutc+sodutc/86400.d0, ut1r(2), ut1r(3))
  CALL UTLIBR(mjdtt+sodtt/86400.d0, ut1r(2), ut1r(3))
  ut1r(2) = ut1r(2)*1.D-6

  !! xp and yp in arcsec
  xp=erp(1)+ocean(1)+libr(1)
  yp=erp(2)+ocean(2)+libr(2)
  !! erp in rad
  erp(1) = xp*ARCSEC2RAD
  erp(2) = yp*ARCSEC2RAD
  erp(3) = erp(3)+ocean(3)+ut1r(1)+ut1r(2)

  ! Predict the Earth rotation angle for this UT1.
  CALL timinc(mjdutc,sodutc,erp(3),mjdut1,sodut1)


  ! Form the celestial-to-intermediate matrix for this TT.
  SELECT CASE(TRIM(conv))
    CASE('IERS2010')
      CALL iau_C2I06A(DBLE(mjdtt+2400000.5D0),sodtt/86400.d0,QMAT)
      !! CALL iau_PNM06A(DBLE(mjdtt+2400000.5D0),sodtt/86400.d0,QMAT)
    CASE('IERS2003')
      CALL iau_C2I00A(DBLE(mjdtt+2400000.5D0),sodtt/86400.d0,QMAT)
      !! CALL iau_PNM00A(DBLE(mjdtt+2400000.5D0),sodtt/86400.d0,QMAT)
    CASE('IERS1996')
      CALL iau_PNM80(DBLE(mjdtt+2400000.5D0),sodtt/86400.d0,QMAT)
  END SELECT
  ! Form the intermediate-to-celestial matrix
  QMAT=TRANSPOSE(QMAT)

  SELECT CASE(TRIM(conv))
    CASE('IERS2010')
      era = iau_ERA00(DBLE(mjdut1+2400000.5D0),sodut1/86400.d0)
      !! era = iau_GST06A(DBLE(mjdut1+2400000.5D0),sodut1/86400.d0,DBLE(mjdtt+2400000.5D0),sodtt/86400.d0)
    CASE('IERS2003')
      era = iau_ERA00(DBLE(mjdut1+2400000.5D0),sodut1/86400.d0)
      !! era = iau_GST00A(DBLE(mjdut1+2400000.5D0),sodut1/86400.d0,DBLE(mjdtt+2400000.5D0),sodtt/86400.d0)
    CASE('IERS1996')
      era = iau_GST94(DBLE(mjdut1+2400000.5D0),sodut1/86400.d0)
  END SELECT

  SELECT CASE(TRIM(conv))
    CASE('IERS2010')
      gmst = iau_GMST06(DBLE(mjdut1+2400000.5D0),sodut1/86400.d0,DBLE(mjdtt+2400000.5D0),sodtt/86400.d0)
    CASE('IERS2003')
      gmst = iau_GMST00(DBLE(mjdut1+2400000.5D0),sodut1/86400.d0,DBLE(mjdtt+2400000.5D0),sodtt/86400.d0)
  END SELECT

  ! Form the Earth rotation matrix
  CALL iau_IR(RMAT)
  CALL iau_RZ(-era, RMAT)
  dlmat=0.d0
  dlmat(1,1) = -DSIN(era)*E_ROTATE
  dlmat(2,1) =  DCOS(era)*E_ROTATE
  dlmat(1,2) = -DCOS(era)*E_ROTATE
  dlmat(2,2) = -DSIN(era)*E_ROTATE

  !*  Estimate s'
  SELECT CASE(TRIM(conv))
    CASE('IERS2003','IERS2010')
      sp = iau_SP00(DBLE(mjdtt+2400000.5D0), sodtt/86400.d0)
    CASE('IERS1996')
  END SELECT
  CALL iau_IR(RSMA)
  CALL iau_RZ(-sp, RSMA)

  CALL iau_IR(RXMA)
  CALL iau_RY(erp(1),RXMA)
  CALL iau_IR(RYMA)
  CALL iau_RX(erp(2),RYMA)

  CALL iau_IR(WMAT)
  CALL iau_RXR(RSMA,RXMA,WMAT)
  CALL iau_RXR(WMAT,RYMA,WMAT)

  CALL iau_IR(tmat)
  CALL iau_RXR(QMAT,RMAT,tmat)

  dxmat=0.d0
  dxmat(1,1) = -DSIN(erp(1))*ARCSEC2RAD
  dxmat(3,3) = -DSIN(erp(1))*ARCSEC2RAD
  dxmat(3,1) =  DCOS(erp(1))*ARCSEC2RAD
  dxmat(1,3) = -DCOS(erp(1))*ARCSEC2RAD

  CALL iau_RXR(tmat,RSMA,MATT)
  CALL iau_RXR(MATT,dxmat,dxmat)
  CALL iau_RXR(dxmat,RYMA,dxmat)

  dymat=0.d0
  dymat(2,2) = -DSIN(erp(2))*ARCSEC2RAD
  dymat(3,3) = -DSIN(erp(2))*ARCSEC2RAD
  dymat(2,3) =  DCOS(erp(2))*ARCSEC2RAD
  dymat(3,2) = -DCOS(erp(2))*ARCSEC2RAD


  CALL iau_RXR(MATT,RXMA,MATT)
  CALL iau_RXR(MATT,dymat,dymat)

  CALL iau_RXR(tmat,WMAT,tmat)

  CALL iau_RXR(QMAT,dlmat,MATT)
  CALL iau_RXR(MATT,WMAT,dlmat)

  RETURN

END SUBROUTINE
