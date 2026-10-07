!*
SUBROUTINE oi_fes2004_tide(mjd,lmax,m_diurnal,m_long,m_radi,dc,ds)
!!
!*
USE tables
IMPLICIT NONE


!*
! The local variables
!!--------------------------------
REAL(RL) :: mjd
INTEGER(IT) :: lmax
! m_diurnal = 0: diurnal and semidiurnal constituents are considered
! m_diurnal = 1: diurnal and semidiurnal constituents are not considered
! m_long = 0: long period constituents are considered
! m_long = 1: long period constituents are not considered
! m_radi = 0: radial component constituents are considered
! m_radi = 1: radial component constituents are not considered
INTEGER(IT) :: m_diurnal, m_long, m_radi
! corrections to spherical harmonics due to ocean tide
REAL(RL) :: dC(MAXOCNDEG,0:MAXOCNDEG)
REAL(RL) :: dS(MAXOCNDEG,0:MAXOCNDEG)

  !*
  ! The local variables
  !!------------------------------
  ! total number of constituents
  INTEGER(IT) :: nr_cons
  ! file for ocean tide data
  CHARACTER(LEN_FILENAME) :: coefs_file
  ! location-dependent spherical harmonics expansions
  ! (degree,order,constituent,h1/h2)
  REAL(RL) :: Ccoefs(0:MAXOCNDEG,0:MAXOCNDEG,24,2)
  REAL(RL) :: Scoefs(0:MAXOCNDEG,0:MAXOCNDEG,24,2)

  ! the Delaunay multipliers for each constituent
  INTEGER(IT) :: mult(24,6)
  ! an additional angle for each constituent
  REAL(RL) :: chi(24)
  ! time-dependent elevation coefficients
  REAL(RL) :: C(0:lmax,0:lmax)
  REAL(RL) :: S(0:lmax,0:lmax)
  ! Love number if no radial component constituents available
  REAL(RL) :: kn(MAXOCNDEG)
  ! density of water and earth
  REAL(RL) :: rho_w,rho_e
  REAL(RL) :: GM, G
  REAL(RL) :: earth_R

  INTEGER(IT) :: i,j,l,m,err
  LOGICAL(LG) :: lfirst

  REAL(RL) :: date0,date,date_TT
  REAL(RL) :: factor, scale

  DATA lfirst /.TRUE./

  ! save data
  SAVE Ccoefs,Scoefs,mult,chi,nr_cons,kn

  !*
  ! Start the exectuable code
  !!------------------------------

  ! parameters used
  rho_w = 1025.d0             ! [kg/m^3] mean density of sea water
  rho_e = 5515.3d0            ! [kg/m^3] mean density of Earth
  earth_R = 6378136.46d0      ! Earth radius
  GM = 3.986004415e14         ! GM
  G = 6.67428e-11               ! constant of gravitation

  ! Initialization

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst = .FALSE.

    ! read location-dependent coefficieients and love numbers
    CALL get_otide_model(m_diurnal,m_long,m_radi,nr_cons,Ccoefs,Scoefs,mult,chi,kn,err)
    IF (err .NE. 0) THEN
      WRITE(*,*) '***ERROR(get_otide_model): reading the ocean tide model error'
      CALL exit(1)
    END IF
  END IF

  ! transfer time into the proper time system
  ! The input is in TT/TDT time system
  date_TT = mjd
  date0 = DBLE(int(date_TT))
  ! adding the mjd0 at [days] conventional modified Julian date
  ! offset,i.e. the number of days since midnight on November 17, 1858
  ! for example, 2002-08-01-00-00-00
  ! date0 = 2452816.5000
  ! date  = 0.0005924074 (d) = 51.184 (s)
  date = date_TT - date0
  date0 = date0 + 2400000.5d0
  CALL cal_elevation_CS(date0,date,nr_cons,lmax,Ccoefs,Scoefs,mult,chi,C,S)

  ! Notice that 'date0' and 'date' are in TT time system

  ! scaling coefficients of last step into dimensionless spherical harmonics
  dC(1:MAXOCNDEG,0:MAXOCNDEG) = 0.0d0
  dS(1:MAXOCNDEG,0:MAXOCNDEG) = 0.0d0

  factor=4.0d0*pi*G*rho_w*earth_R*earth_R/GM

  ! we only compute and use degrees larger than 2
  DO i=2, lmax
    scale = factor*(1.0d0+kn(i))/(2*i+1)
    DO j = 0,i
       dC(i,j) = scale*C(i,j)
       dS(i,j) = scale*S(i,j)
    END DO
  END DO

101 FORMAT(2i4,2e22.14)

  RETURN

END SUBROUTINE


SUBROUTINE get_otide_model(m_diurnal,m_long,m_radi, nr_cons,Ccoefs,Scoefs,mult,chi,kn,err)
!!
!*
USE par
USE tables
IMPLICIT NONE

!*
! The arguments
!!--------------------------------
! diurnal and semidiurnal ocean tides are to be modeled: 0 - no; 1 - yes
INTEGER(IT) :: m_diurnal
! long-period ocean tides are to be modeled: 0 - no; 1 - yes
INTEGER(IT)  :: m_long
! radial component ocean tides are to be modeled: 0 - no; 1 - yes
INTEGER(IT)  :: m_radi


! Number of all constituents according to model types considered
INTEGER(IT) :: nr_cons
REAL(RL) :: Ccoefs(0:MAXOCNDEG,0:MAXOCNDEG,24,2)
REAL(RL) :: Scoefs(0:MAXOCNDEG,0:MAXOCNDEG,24,2)
! the Delaunay multipliers for each constituent
INTEGER(IT) :: mult(24,6)
! an additional angle for each constituent
REAL(RL) :: chi(24)
! load love numbers if no radial spherical
! harmonics coefficients avaiable
REAL(RL) :: kn(MAXOCNDEG)
! error code: 0 = OK
INTEGER(IT) :: err

  !*
  ! The internal variables
  !!--------------------------------
  ! loop counters
  INTEGER(IT) :: i,j,k,l,m
  ! names of the files with coefficients
  ! for short, long period and radial component constituents respectively
  CHARACTER(LEN_FILENAME) :: ocean_file
  ! name of love number file if radial component
  ! coefficients are not avaiable
  CHARACTER(LEN_FILENAME) :: love_nr_file
  ! number of constituents
  INTEGER(IT) :: nr_short, nr_long, nr_radi

  !*
  ! Start the exectuable code
  !!--------------------------------

  ! Initialization:
  err = 0

  Ccoefs = 0.0d0
  Scoefs = 0.0d0
  mult = 0
  chi  = 0.0d0
  kn = 0.0d0

  nr_cons = 13
  IF (m_diurnal .EQ. 0) THEN
    nr_short = 9
    ocean_file = f_tableFileName("fessht")
    CALL read_tides(ocean_file,1,nr_short,Ccoefs,Scoefs,mult,chi,err)
  END IF

  IF (m_long .EQ. 0) THEN
    nr_long = 4
    ocean_file = f_tableFileName("feslon")
    CALL read_tides(ocean_file,10,nr_long,Ccoefs,Scoefs,mult,chi,err)
  END IF

  IF (m_radi .EQ. 0) then
    !nr_radi = 9
    !ocean_file=f_tableFileName("fesrad")
    !call read_tides(ocean_file,14,nr_radi,Ccoefs,Scoefs,mult,chi,err)
    !mult(14:22,1:5) = - mult(14:22,1:5)
    ocean_file = f_tableFileName("ldcoef")
    CALL read_love_nr(ocean_file,kn,err)
  END IF

  RETURN

END SUBROUTINE get_otide_model


SUBROUTINE read_tides(file_in,cons_beg,cons_add,Ccoefs,Scoefs,mult,chi,err)
!
! Be aware that the input files must be the one prepared by Xianglin in PANDA format
! the order of constituents is extremely important, otherwise, it may cause problems.
!
! the orders for short, long, and radial component constituents are described
! in the header files provided by Xianglin Liu
!
USE par
IMPLICIT NONE

!*
! The arguments
!!--------------------------------
! names of the files with coefficients
! for short or long period or radial component constituents
CHARACTER(LEN_FILENAME) :: file_in
! the beginning constituent in Ccoefs and Scoefs
INTEGER(IT) :: cons_beg
! number of constituents considered
INTEGER(IT) :: cons_add

REAL(RL) :: Ccoefs(0:MAXOCNDEG,0:MAXOCNDEG,24,2)
REAL(RL) :: Scoefs(0:MAXOCNDEG,0:MAXOCNDEG,24,2)
! the Delaunay multipliers for each constituent
INTEGER(IT) :: mult(24,6)
! an additional angle for each constituent
REAL(RL) :: chi(24)
! error code: 0 = OK
INTEGER(IT) :: err


  !*
  ! The internal variables
  !!--------------------------------
  INTEGER(IT) :: i,j,k,l,m
  INTEGER(IT) :: lun,l_in,m_in
  CHARACTER(LEN_STRING) :: name_ot

  !*
  ! The function called
  !!--------------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!--------------------------------

  ! Initialization
      err = 0

! An attempt to open the file:

      lun=get_valid_unit(10)
      open(lun,file=file_in,status='old',action='read',err=1001)

! Skipping header:

      do i=1,4
          read(lun,*)
      enddo

! Reading coefficients from file for per constituent:

      ! loop over constituents
      do j = cons_beg, cons_beg+cons_add-1

        ! read h1 coefficients
        read(lun,99) name_ot, mult(j,1:6), chi(j)
        do l=0,80
          do m=0,l
            read(lun,101,err=1002) l_in,m_in,Ccoefs(l,m,j,1),Scoefs(l,m,j,1)
            if (l_in.ne.l.or.m_in.ne.m) then
              write(*,*) '***ERROR(read_tides): wrong order or degree:',l_in,m_in
              call exit(1)
            endif
            enddo
        enddo

        ! read h2 coefficients they are the same with that of h1
        read(lun,99) name_ot, mult(j,1:6), chi(j)
        do l=0,80
          do m=0,l
            read(lun,101,err=1002) l_in,m_in,Ccoefs(l,m,j,2),Scoefs(l,m,j,2)
            if (l_in.ne.l.or.m_in.ne.m) then
              write(*,*) '***ERROR(read_tides): wrong order or degree:',l_in,m_in
              call exit(1)
            endif
            enddo
        enddo

      enddo
      close(lun)

99    format(a15,6(1x,i3),1x,f20.15)
101   format(2i4,2e22.14)
      return

1001    err = 1
        write(*,*) '***ERROR(read_tides): open file ',trim(file_in)
      call exit(1)
        return

1002    err = 2
      write(*,*) '***ERROR(read_tides): read file ',trim(file_in)
        close(lun)
      call exit(1)
      return

      end subroutine read_tides


      subroutine read_love_nr(file_in,kn,ierr)
      use par
      implicit none

! Input:
      character*(*) file_in         ! names of the files with coefficients
                                    ! for short or long period or radial component

! Output

      real(rl) :: kn(MAXOCNDEG)     ! load love numbers if no radial spherical
                                         ! harmonics coefficients avaiable

      integer(it) :: ierr       ! error code: 0 = OK

! Internal:

      integer(it) :: i,j        ! loop counters
      integer(it) :: lun
      character(len_string) :: text             ! to store some text

! Function called
      integer(it) :: get_valid_unit

! Initialization:

      ierr = 0

! An attempt to open the file:

      lun=get_valid_unit(10)
      open(lun,file=file_in,status='old',action='read',err=1001)

      text = '  '
      do while(index(text,'+load coefficient').eq.0)
        read(lun, '(a)') text
      enddo
      do while(index(text,'-load coefficient').eq.0)
         if(text(1:1).eq.' ') then
           read(text,*,iostat=ierr) i, kn(i)
           if(ierr.ne.0) then
              write(*,*) '***ERROR(read_love_nr): read the following line error'
              write(*,*) text
              call exit(1)
           endif
         endif
         read(lun,'(a)') text
      enddo
      close(lun)

      return

1001    ierr = 1
        write(*,*) '***ERROR(read_love_nr): cannot find or open file ',trim(file_in)
      call exit(1)

      end subroutine read_love_nr

!-----------------------------------------------------------------------
! cal_elevation_CS: this routine performs the following summation over
! the in-phase (h1) and quadrature (h2) parts of all constituents
! for every degree/order. If certain a constituent is not considered,
! the corresponding coefficients are set as zeros in the previous step.
!
!       C = sum_j C_h1_j*cos(argument_j) + C_h2_j*sin(argument_j),
!       S = sum_j S_h1_j*cos(argument_j) + S_h2_j*sin(argument_j),
!
! where:
!       j = current constituent number (see below)
!
!       argument_j = sum_i N_ji*F_i(t)
!       N_ji   = six Delaunay multipliers for constituent j
!       F_i(t) = six Delaunay variables for time t
!
! Note: The sixth argument is GMST+pi of which the multiplier is
!       denoted CHI in the IERS conventions 2003, chapter 8.
!
! The output 'C/S' is the resulting elevation coefficients that is
! induced by the constituents. This effect is returned in the same
! units as the in-phase and quadrature parts (m).
!
! The constituents should be given in order of frequency, for example
! in FES2004 model, the order of constituents is arranged as follows
!    j         j          j           j
!    1 = Q1    5 = N2     9 = 2N2    13 = MSqm
!    2 = O1    6 = M2    10 = Mm
!    3 = P1    7 = S2    11 = Mf
!    4 = K1    8 = K2    12 = Mtm
!
! The radial component constiuents follow as
!    j         j          j
!   14 = Q1   18 = N2    22 = 2N2
!   15 = O1   19 = M2
!   16 = P1   20 = S2
!   17 = K1   21 = K2
!
! Modified by Xianglin Liu (x.l.liu@tudelft.nl) from
! an original subroutine called 'tidal_effect' written by Sander
! to compute the elevation (time-dependent) coefficients from
! location-dependent spherical harmonics instead of compute the
! tide effect from input potentials/accelerations/tensors.
! The original subroutine is as
! subroutine tidal_effect(date0,date,nr,h1,h2,tide)
!
!
      subroutine cal_elevation_CS(date0,date,nr_cons,lmax,Ccoefs,Scoefs,mult,chi,C,S)
      use par
      implicit none


! Input:

      ! The Julian date: JD = date0 + date [TAI]
      ! It was TAI time system in Sander's subroutine
      ! Here, the input is in TT time system
      real(rl) :: date0,date
      ! number of constituents considered
      integer(it) :: nr_cons
      ! maximum degree defined by user
      integer(it) :: lmax
        ! Location-dependent spherical harmonics expansion computed
      ! by SHA for each constituent
      real(rl) :: Ccoefs(0:MAXOCNDEG,0:MAXOCNDEG,24,2)
      real(rl) :: Scoefs(0:MAXOCNDEG,0:MAXOCNDEG,24,2)
      ! the Delaunay multipliers for each constituent
      integer(it) :: mult(24,6)
        ! an additional angle for each constituent
      real(rl) :: chi(24)

! Output:

      ! Time-dependent elevation coefficients
      ! looped over all constituents
      real(rl) :: C(0:lmax,0:lmax)
      real(rl) :: S(0:lmax,0:lmax)

! Internal:

      ! loop variable
      integer(it) :: i,j,k
        ! the Delaunay variables: l, l', F, D, Omega, GMST + p
      real(rl) :: F(6)
        ! arguments computed by multiplying the Delaunay
      ! variables with the multipliers for each constituent
      real(rl) :: args(lmax)

      ! Compute Delaunay variables:
      call Delaunay(date0,date,F)

      ! Compute the resulting arguments:
      do j = 1, nr_cons
         args(j) = mult(j,1)*F(1) + mult(j,2)*F(2) + &
                   mult(j,3)*F(3) + mult(j,4)*F(4) + &
                   mult(j,5)*F(5) + mult(j,6)*F(6)

         ! set arguments to 0 =< args < 2*pi
         !args(j) = mod(args(j),2*pi)
         !if(args(j).lt.0.d0) args(j) = args(j) + 2*pi
      enddo

! Compute elevation coefficeints:
      C = 0.d0
      S = 0.d0

      do i = 0,lmax
         do j = 0,i
            do k = 1,nr_cons
               C(i,j) = C(i,j) + Ccoefs(i,j,k,1)*dcos(args(k)+chi(k)) &
                               + Ccoefs(i,j,k,2)*dsin(args(k)+chi(k))
               S(i,j) = S(i,j) + Scoefs(i,j,k,1)*dcos(args(k)+chi(k)) &
                               + Scoefs(i,j,k,2)*dsin(args(k)+chi(k))
            enddo
         enddo
      enddo

      end subroutine cal_elevation_CS

!-----------------------------------------------------------------------
! Delaunay: This subroutine returns the 5 Delaunay variables in radians:
! F1 = l; F2 = l'; F3 = F; F4 = D; F5 = Omega; and the additional
! rotation angle: F6 = GMST + p.
! The expressions for these arguments are taken from the IERS 2003
! Conventions, chapter 5, which can be downloaded from:
!
!         http://maia.usno.navy.mil/conv2000.html
!
! Note: The following approximation is used: TDB = TT
!       This introduces an error of < 10^-5 mas, which
!       is negligible for all but the most high precision
!       computations.
!
! Programmer: Sander van Eck van der Sluijs (s.vaneck@citg.tudelft.nl)
! Date      : September 2003
!

      subroutine Delaunay(date0,date,F)
      use par
      use const
      implicit none

! Input:

      ! The Julian date: JD = date0 + date [TAI]
      real(rl) :: date0,date

! Output:

      ! the Delaunay variables in [rad]
      real(rl) :: F(6)

! Internal:

      real(rl), parameter :: mjd2000 = 2451545.0d0
    ! input date converted to TT time
      real(rl) :: date_TT
    ! time in Julian centuries of TDB plus its powers
      real(rl) :: t,t2,t3,t4

! Functions:

!     real*8 iau_GMST82

! Compute t in Julian centuries:

      ! if input is in TAI; convert to TT
      ! (see IERS 2003 conventions, chapter 10):
!       date_TT = date + 32.184d0/86400.d0

      ! if input is in TT
        date_TT = date

      ! compute time in Julian centuries w.r.t. J2000.0
        t = ((date0 - mjd2000) + date_TT)/36525.d0

        t2 = t*t
        t3 = t*t2
        t4 = t*t3

      ! Compute Delaunay variables (IERS 2003 conventions, chapter 5, equ. 40, p48):

      ! l : Mean anomaly of the Moon
      F(1)= 134.96340251d0*deg2rad + 1717915923.2178d0*ARCSEC2RAD*t &
            + 31.8792d0*ARCSEC2RAD*t2 + 0.051635d0*ARCSEC2RAD*t3 &
            - 0.00024470d0*ARCSEC2RAD*t4

      ! l': Mean anomaly of the Sun
      F(2)= 357.52910918d0*deg2rad + 129596581.0481d0*ARCSEC2RAD*t &
            - 0.5532d0*ARCSEC2RAD*t2 + 0.000136d0*ARCSEC2RAD*t3 &
            - 0.00001149d0*ARCSEC2RAD*t4

      ! L - Omega:  L = Mean longitude of the Moon
      F(3)= 93.27209062d0*deg2rad + 1739527262.8478d0*ARCSEC2RAD*t &
            - 12.7512d0*ARCSEC2RAD*t2 - 0.001037d0*ARCSEC2RAD*t3 &
            + 0.00000417d0*ARCSEC2RAD*t4

      ! D : Mean elongation of the Moon from the Sun
      F(4)= 297.85019547d0*deg2rad + 1602961601.2090d0*ARCSEC2RAD*t &
            - 6.3706d0*ARCSEC2RAD*t2 + 0.006593d0*ARCSEC2RAD*t3 &
            - 0.00003169d0*ARCSEC2RAD*t4

      ! Omega : Mean longitude of the ascending node of the Moon
      F(5)= 125.04455501d0*deg2rad - 6962890.5431d0*ARCSEC2RAD*t &
            + 7.4722d0*ARCSEC2RAD*t2 + 0.007702d0*ARCSEC2RAD*t3 &
            - 0.00005939d0*ARCSEC2RAD*t4

      ! Rotation angle  GMST + p (p = pi)
      F(6) = ((67310.54841d0 + (876600.d0*3600.d0 + 8640184.812866d0)*t &
              + 0.093104d0*t2 - 6.2d0-6*t3 )*15.d0 + 648000.d0)*ARCSEC2RAD

      ! alternative
      !F(6) = iau_GMST82(date0,date) + pi


      end subroutine Delaunay
