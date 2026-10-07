!*
SUBROUTINE iono_corr_higher_order(xsit,xsat,f1,f2,stec,mjd,fmjd,xi2)
!!
!*
USE const
IMPLICIT NONE

!      
! The input arguments
! ---------------------------

! Station -> Satellite unit vecktor in TRS
REAL(RL) :: xsat2xsta(3),xsit(3),xsat(3)
! Frequency in L1 and L2 frequencies [Hz]
REAL(RL) :: f1, f2
! Slant total elecctron content (STEC) for Sta/Sat/Epoch
REAL(RL) :: stec
! Full MJD      55125
INTEGER(IT) :: mjd
! second of day   43200.0
REAL(RL) :: fmjd
! Ionosphere pierce point position in TRS in m
REAL(RL) :: xp(3)
! 2nd-order ionospheric correction term (to LC), s term
REAL(RL) :: xi2

  !
  ! The local variables
  !!---------------------------

  ! [10^16 electron per m^2], constant
  REAL(RL) :: iono_const=40.309d16

  ! IONOSPHERIC PIERCE POINT (Xp)
  ! Earth radius in km
  REAL(RL) :: re = 6378.137d0
  ! Height of Ionospere Layer in km
  REAL(RL) :: hion =   450.d0
  ! distance of the xp and its projection on xy-plane
  REAL(RL) :: norm_xipp,norm_xy_xipp
  ! Satellite -> Station unit vecktor in TRS
  REAL(RL) :: xsta2xsat(3)
  ! Colatitude (0-180) and east-longitude (0-360) of Xp 
  REAL(RL) :: colat,elong

  !
  ! DATE INTO IGRF FORMAT
  ! MJD; 4-digit-YEAR; 2-digit-YEAR; DayOfYear; #Days/year
  INTEGER(IT) :: iyyyy,idoy,ndays,ii
  ! Date in [decimal_years]
  REAL(RL) :: date

  ! 
  ! MAGNETIC FIELD AT IONOSPHERIC PIERCE POINT
  ! IGRF parameters; error code
  INTEGER(IT) :: isv, itype, ierr
  ! Altitude of Iono.Layer (re+hion)
  REAL(RL) :: alt
  ! North, East, Vertical component and total intensity of B [nT]
  REAL(RL) :: bn,be,bv,b1
  ! Magentic field at Xp in XYZ coordinates
  REAL(RL) :: bxyz(3)

  ! 
  ! MAGNETIC FIELD PROJECTION (Bproj) ALONG PROPAGETION DIRECTION
  ! Magnetic field projected toward the propagation direction [T]
  REAL(RL) :: Bproj

  !
  ! 2ND-ORDER IONOSPHERIC CORRECTION TERM (xi2)
  ! s-term
  REAL(RL) :: s2
  REAL(RL) :: rot_loc2trs(3,3), bloc(3)
  ! Ionosphere pierce point position BLH in [rad, m]
  REAL(RL) :: blh_p(3)

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: modified_julday

  !
  ! Start the exectuable code
  ! ---------------------------    

  xsat2xsta=0.d0
  DO ii=1, 3
    xsat2xsta(ii)=xsat(ii)-xsit(ii)
  END DO
  CALL unit_vector(3,xsat2xsta,xsat2xsta,s2)

  ! in km
  alt=re+hion
  ! direction
  xsta2xsat = -xsat2xsta

  s2=450000.d0/s2

  DO ii=1, 3
    xp(ii)=xsit(ii)+s2*xsat2xsta(ii)
  END DO
  
  CALL xyzblh(xp,1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,blh_p)
  CALL rot_enu2xyz(blh_p(1),blh_p(2),rot_loc2trs)

  !
  ! DATE INTO IGRF FORMAT
  ! ---------------------
  !   re-define MJD because of intend(INOUT) in SR dateconv
  
  CALL mjd2doy(mjd,iyyyy,idoy)
  ndays=-modified_julday(1,1,iyyyy)+modified_julday(31,12,iyyyy)+1
  ! Date [in decimal years]
  date = dble(iyyyy) + dble( idoy + fmjd/86400.d0) / dble(ndays)


  !
  ! MAGNETIC FIELD AT IONOSPHERIC PIERCE POINT
  ! ------------------------------------------
  !   IGRF parameters: Main field with isv=0; Geocentric distance with
  !   itype=2
  isv  =0
  itype=2

  !   
  !   Magentic field: Components of B (bn,be,bv) and Total intensity b1
  !   [all in nT]
  IF (blh_p(2) .LT. 0.d0) THEN
    CALL igrf13syn(isv,date,itype,alt,90d0-blh_p(1)*RAD2DEG,360d0+blh_p(2)*RAD2DEG,bn,be,bv,b1)
  ELSE
    CALL igrf13syn(isv,date,itype,alt,90d0-blh_p(1)*RAD2DEG,blh_p(2)*RAD2DEG,bn,be,bv,b1)
  END IF

  !
  ! Transformation of B: NEU -> XYZ [nT]
  !  call nev2xyz(xp,bn,be,bv,bxyz,ierr)
  !   write(*,'(a,3f10.3)') 'loc,be,bn,bv',be,bn,bv
  bloc(1)=be
  bloc(2)=bn
  bloc(3)=bv
  CALL matmpy(rot_loc2trs,bloc,bxyz,3,3,1)

  !   write(*,*) 'bxyz1',rot_loc2trs, bloc
  !     write(*,*) 'bxyz2', bxyz, xsta2xsat
  !
  ! MAGNETIC FIELD PROJECTION (Bproj) ALONG PROPAGETION DIRECTION
  ! -------------------------------------------------------------
  Bproj=0.d0
  DO ii=1,3
    Bproj=Bproj+bxyz(ii)*xsta2xsat(ii)
  END DO
  !/norm_xsat2xsta       ![T]
  Bproj=Bproj*1.d-9

  !
  !
  ! 2ND-ORDER IONOSPHERIC CORRECTION TERM (xi2)
  ! -------------------------------------------
  !   s term
  s2 = 7527.d0*VEL_LIGHT*Bproj*stec*1.d16*(re/(re+hion))**3

  !   
  ! I2 correction term
  xi2 = s2 / (f1*f2*(f1+f2)) ! [m]

  !     write(*,'(a,6(e20.6))')'  c,Bproj,stec,s2,f1f2f12,xi2:',vel_light,
  !     Bproj, stec*1.d16, s2, f1*f2*(f1+f2),  xi2
  !
    
  RETURN

END SUBROUTINE iono_corr_higher_order
