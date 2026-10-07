!*
SUBROUTINE gpt3(mjd,lat,lon,h_ell,n_stat,itr,p,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w )
!!
!
! (c) Department of Geodesy and Geoinformation, Vienna University of
! Technology, 2018
!
!
! This subroutine determines pressure, temperature, temperature lapse rate, 
! mean temperature of the water vapor, water vapour pressure, hydrostatic 
! and wet mapping function coefficients ah and aw, water vapour decrease
! factor, geoid undulation and empirical tropospheric gradients for 
! specific sites near the earth's surface. All output values are valid for
! the specified ellipsoidal height h_ell.
! GPT3_1 is based on a 1°x1° external grid file ('gpt3_1.grd') with mean 
! values as well as sine and cosine amplitudes for the annual and
! semiannual variation of the coefficients.
! Because the .grd file has to be opened anew every time this function is 
! called, the process is fairly time-consuming for a longer set of 
! mjd's. 
!
! In order to get the respective mapping function, the ah and aw coefficients 
! must be input to vmf3_ht.f90.
! In order to get the respective zenith hydrostatic delay, pressure p must be
! input to saasthyd.f.
! In order to get an approximate value for the respective zenith wet delay,
! water vapor pressure e, mean temperature Tm and water vapor decrease factor
! la must be input to asknewet.f.
!
!
! Reference:
! D. Landskron, J. Böhm (2018), VMF3/GPT3: Refined Discrete and Empirical Troposphere Mapping Functions, 
! J Geod (2018) 92: 349., doi: 10.1007/s00190-017-1066-2. 
! Download at: https://link.springer.com/content/pdf/10.1007%2Fs00190-017-1066-2.pdf
!
!
! Input parameters:
!
! mjd:   modified Julian date (scalar, only one epoch per call is possible)
! lat:   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
! lon:   longitude in radians [-pi:pi] or [0:2pi] (vector)
! h_ell: ellipsoidal height in m (vector)
! it:    case 1: no time variation but static quantities
!        case 0: with time variation (annual and semiannual terms)
! 
! Output parameters:
!
! p:    pressure in hPa (vector) 
! T:    temperature in degrees Celsius (vector)
! dT:   temperature lapse rate in degrees per km (vector)
! Tm:   mean temperature weighted with the water vapor in degrees Kelvin (vector) 
! e:    water vapour pressure in hPa (vector)
! ah:   hydrostatic mapping function coefficient (VMF3) (vector)
! aw:   wet mapping function coefficient (VMF3) (vector)
! la:   water vapour decrease factor (vector)
! undu: geoid undulation in m (vector)
! Gn_h: hydrostatic north gradient in m (vector)
! Ge_h: hydrostatic east gradient in m (vector)
! Gn_w: wet north gradient in m (vector)
! Ge_w: wet east gradient in m (vector)
!
!
! File created by Daniel Landskron, 2018/02/22
!
!*
USE const
USE tables
IMPLICIT NONE
 
!*
! The input arguments
!!---------------------------
INTEGER(IT), INTENT(IN) :: n_stat
INTEGER(IT), INTENT(IN) :: itr
REAL(RL), INTENT(IN) :: mjd
REAL(RL), DIMENSION(n_stat), INTENT(IN) :: lat
REAL(RL), DIMENSION(n_stat), INTENT(IN) :: lon
REAL(RL), DIMENSION(n_stat), INTENT(IN) :: h_ell

!*
! Output arguments
!!---------------------------
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: p
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: T
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: dT
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: Tm
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: e
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: ah
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: aw
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: la
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: undu
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: Gn_h
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: Ge_h
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: Gn_w
REAL(RL), DIMENSION(n_stat), INTENT(OUT) :: Ge_w

  !*
  ! The local variables
  !!---------------------------
  ! file open
  INTEGER(IT) :: lfn
  INTEGER(IT) :: i_grid
  CHARACTER(len_string) :: line
  CHARACTER(len_filename) :: grdfile

  ! variables for the gpt3_5.grd
  REAL(RL), DIMENSION(64800) :: lat_grid
  REAL(RL), DIMENSION(64800) :: lon_grid
  REAL(RL), DIMENSION(64800,5) :: p_grid
  REAL(RL), DIMENSION(64800,5) :: T_grid
  REAL(RL), DIMENSION(64800,5) :: Q_grid
  REAL(RL), DIMENSION(64800,5) :: dT_grid
  REAL(RL), DIMENSION(64800) :: u_grid
  REAL(RL), DIMENSION(64800) :: Hs_grid
  REAL(RL), DIMENSION(64800,5) :: ah_grid
  REAL(RL), DIMENSION(64800,5) :: aw_grid
  REAL(RL), DIMENSION(64800,5) :: la_grid
  REAL(RL), DIMENSION(64800,5) :: Tm_grid
  REAL(RL), DIMENSION(64800,5) :: Gn_h_grid
  REAL(RL), DIMENSION(64800,5) :: Gn_w_grid
  REAL(RL), DIMENSION(64800,5) :: Ge_h_grid
  REAL(RL), DIMENSION(64800,5) :: Ge_w_grid

  LOGICAL(LG) :: lfirst
  DATA lfirst /.TRUE./

  SAVE lfirst,lat_grid,lon_grid,p_grid,T_grid,Q_grid,dT_grid,u_grid,Hs_grid,&
       ah_grid,aw_grid,la_grid,Tm_grid,Gn_h_grid,Gn_w_grid,Ge_h_grid,Ge_w_grid


  ! variable for the conversion from mjd to doy
  INTEGER(IT) :: year,doy,i_stat

  ! constants
  REAL(RL) :: gm
  REAL(RL) :: dMtr
  REAL(RL) :: Rg

  ! coefficients for amplitudes
  REAL(RL) :: cosfy
  REAL(RL) :: coshy
  REAL(RL) :: sinfy
  REAL(RL) :: sinhy

  ! variables for coordinate conversions
  REAL(RL) :: plon, ppod, difflon, diffpod
  INTEGER(IT) :: ilon, ipod, ilon1, ipod1
  REAL(RL) :: dnpod1, dnpod2, dnlon1, dnlon2

  ! variable to decide whether bilinear or nearest neighbor shall be applied
  INTEGER(IT) :: bilinear

  ! variable for the indices in the .grd file
  INTEGER(IT), DIMENSION(4) :: indx
  INTEGER(IT) :: ix,i_ind

  ! further variables
  REAL(RL) :: hgt, redh, c, Hs1
  REAL(RL) :: Tv, T0, p0, Q, e0
  REAL(RL), DIMENSION(4) :: undul, Ql, dTl, pl, Tl, lal, ahl, awl, Tml, el, Gn_hl, Ge_hl, Gn_wl, Ge_wl
  REAL(RL) :: R1, R2

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!---------------------------

  if (lfirst .eq. .true.) then
    lfirst=.false.
    
    ! read gridfile
    grdfile = f_tablefilename('gptgrd')
    lfn = get_valid_unit(10)
    open(lfn,file=grdfile)

    ! read first comment line
    read (lfn,'(a80)') line

    ! loop over grid points
    do i_grid = 1,64800
        
      read(lfn,*) lat_grid(i_grid), lon_grid(i_grid), &
                  p_grid(i_grid,1:5),  & ! pressure in Pascal
                  T_grid(i_grid,1:5), & ! temperature in Kelvin
                  Q_grid(i_grid,1:5), & ! specific humidity in kg/kg
                  dT_grid(i_grid,1:5), & ! temperature lapse rate in Kelvin/m
                  u_grid(i_grid), & ! geoid undulation in m
                  Hs_grid(i_grid), & ! orthometric grid height in m
                  ah_grid(i_grid,1:5), & ! hydrostatic mapping function coefficient, dimensionless
                  aw_grid(i_grid,1:5), & ! wet mapping function coefficient, dimensionless
                  la_grid(i_grid,1:5), & ! water vapor decrease factor, dimensionless
                  Tm_grid(i_grid,1:5), & ! mean temperature in Kelvin
                  Gn_h_grid(i_grid,1:5), & ! hydrostatic north gradient in m
                  Ge_h_grid(i_grid,1:5), & ! hydrostatic east gradient in m
                  Gn_w_grid(i_grid,1:5), & ! wet north gradient in m
                  Ge_w_grid(i_grid,1:5) ! wet east gradient in m
                
    end do
    close (lfn)

    ! divide by constants in order to get the real values
    Q_grid = Q_grid/1000
    dT_grid = dt_grid/1000
    ah_grid = ah_grid/1000
    aw_grid = aw_grid/1000
    Gn_h_grid = Gn_h_grid/100000
    Ge_h_grid = Ge_h_grid/100000
    Gn_w_grid = Gn_w_grid/100000
    Ge_w_grid = Ge_w_grid/100000

  end if

  ! convert mjd to doy
  call mjd2doy(INT(mjd),year,doy)


  ! mean gravity in m/s**2
  gm = 9.80665d0
  ! molar mass of dry air in kg/mol
  dMtr = 28.965d-3
  ! universal gas constant in J/K/mol
  Rg = 8.3143d0

  ! factors for amplitudes
  ! then  constant parameters
  if (itr==1) then
    cosfy = 0
    coshy = 0
    sinfy = 0
    sinhy = 0
  else 
    cosfy = dcos(dble(doy)/365.25d0*2*pi)   ! coefficient for A1
    coshy = dcos(dble(doy)/365.25d0*4*pi)   ! coefficient for B1
    sinfy = dsin(dble(doy)/365.25d0*2*pi)   ! coefficient for A2
    sinhy = dsin(dble(doy)/365.25d0*4*pi)   ! coefficient for B2
  end if


  ! loop over all stations
  do i_stat = 1,n_stat

    ! only positive longitude in degrees
    if (lon(i_stat) .lt. 0.d0) then
      plon = (lon(i_stat) + 2*pi)*180/pi
    else
      plon = lon(i_stat)*180/pi
    end if
    ! transform to polar distance in degrees
    ppod = (-lat(i_stat) + pi/2)*180/pi
    
    ! find the index (line in the grid file) of the nearest point
    ipod = floor(ppod+1.d0)
    ilon = floor(plon+1.d0)
    
    ! normalized (to one) differences, can be positive or negative
     diffpod = ppod - (ipod - 0.5d0)
     difflon = plon - (ilon - 0.5d0)
     if (ipod == 181) then
       ipod = 180
     end if
     if (ilon == 361) then
       ilon = 1
     end if
     if (ilon == 0) then
       ilon = 360
     end if

     ! get the number of the corresponding line
     indx(1) = (ipod - 1)*360 + ilon
    
     ! near the poles: nearest neighbour interpolation, otherwise: bilinear
     bilinear = 0
     if (ppod.gt.0.5d0 .AND. ppod.lt.179.5d0 ) then
       bilinear = 1          
     end if
    
    
    ! case of nearest neighbourhood
    if (bilinear .eq. 0) then

      ix = indx(1)
        
      ! transforming ellipsoidial height to orthometric height
      undu(i_stat) = u_grid(ix)
      hgt = h_ell(i_stat)-undu(i_stat)
            
      ! pressure, temperature at the heigtht of the grid
      T0 = T_grid(ix,1) + &
           T_grid(ix,2)*cosfy + T_grid(ix,3)*sinfy + &
           T_grid(ix,4)*coshy + T_grid(ix,5)*sinhy
      p0 = p_grid(ix,1) + &
           p_grid(ix,2)*cosfy + p_grid(ix,3)*sinfy + &
           p_grid(ix,4)*coshy + p_grid(ix,5)*sinhy
         
      ! specific humidity
      Q = Q_grid(ix,1) + &
          Q_grid(ix,2)*cosfy + Q_grid(ix,3)*sinfy + &
          Q_grid(ix,4)*coshy + Q_grid(ix,5)*sinhy
            
      ! lapse rate of the temperature
      dT(i_stat) = dT_grid(ix,1) + &
                   dT_grid(ix,2)*cosfy + dT_grid(ix,3)*sinfy + &
                   dT_grid(ix,4)*coshy + dT_grid(ix,5)*sinhy 

      ! station height - grid height
      redh = hgt - Hs_grid(ix)

      ! temperature at station height in Celsius
      T(i_stat) = T0 + dT(i_stat)*redh - 273.15d0
        
      ! temperature lapse rate in degrees / km
      dT(i_stat) = dT(i_stat)*1000

      ! virtual temperature in Kelvin
      Tv = T0*(1+0.6077d0*Q)
        
      c = gm*dMtr/(Rg*Tv)
        
      ! pressure in hPa
      p(i_stat) = (p0*exp(-c*redh))/100.d0
            
      ! hydrostatic and wet coefficients ah and aw 
      ah(i_stat) = ah_grid(ix,1) + &
                   ah_grid(ix,2)*cosfy + ah_grid(ix,3)*sinfy + &
                   ah_grid(ix,4)*coshy + ah_grid(ix,5)*sinhy
      aw(i_stat) = aw_grid(ix,1) + &
                   aw_grid(ix,2)*cosfy + aw_grid(ix,3)*sinfy + &
                   aw_grid(ix,4)*coshy + aw_grid(ix,5)*sinhy

      ! water vapour decrease factor la
      la(i_stat) = la_grid(ix,1) + &
                   la_grid(ix,2)*cosfy + la_grid(ix,3)*sinfy + &
                   la_grid(ix,4)*coshy + la_grid(ix,5)*sinhy

      ! mean temperature of the water vapor Tm
      Tm(i_stat) = Tm_grid(ix,1) + &
                   Tm_grid(ix,2)*cosfy + Tm_grid(ix,3)*sinfy + &
                   Tm_grid(ix,4)*coshy + Tm_grid(ix,5)*sinhy
            
      ! north and east gradients (total, hydrostatic and wet)
      Gn_h(i_stat) = Gn_h_grid(ix,1) + &
                   Gn_h_grid(ix,2)*cosfy + Gn_h_grid(ix,3)*sinfy + &
                   Gn_h_grid(ix,4)*coshy + Gn_h_grid(ix,5)*sinhy
      Ge_h(i_stat) = Ge_h_grid(ix,1) + &
                   Ge_h_grid(ix,2)*cosfy + Ge_h_grid(ix,3)*sinfy + &
                   Ge_h_grid(ix,4)*coshy + Ge_h_grid(ix,5)*sinhy
      Gn_w(i_stat) = Gn_w_grid(ix,1) + &
                   Gn_w_grid(ix,2)*cosfy + Gn_w_grid(ix,3)*sinfy + &
                   Gn_w_grid(ix,4)*coshy + Gn_w_grid(ix,5)*sinhy
      Ge_w(i_stat) = Ge_w_grid(ix,1) + &
                   Ge_w_grid(ix,2)*cosfy + Ge_w_grid(ix,3)*sinfy + &
                   Ge_w_grid(ix,4)*coshy + Ge_w_grid(ix,5)*sinhy

      ! water vapor pressure in hPa
      ! on the grid
      e0 = Q*p0/(0.622d0+0.378d0*Q)/100.d0
      ! on the station height - (14) Askne and Nordius, 1987
      e(i_stat) = e0*(100.d0*p(i_stat)/p0)**(la(i_stat)+1)
    ! bilinear interpolation
    else
   
      ipod1 = ipod + int(sign(1.d0,diffpod))
      ilon1 = ilon + int(sign(1.d0,difflon))
      if (ilon1 == 361) then
        ilon1 = 1
      end if
      if (ilon1 == 0) then
        ilon1 = 360
      end if
        
      ! get the number of the line
      indx(2) = (ipod1-1)*360 + ilon    ! along same longitude
      indx(3) = (ipod-1)*360 + ilon1   ! along same polar distance
      indx(4) = (ipod1-1)*360 + ilon1   ! diagonal
                
      do i_ind = 1,4
              
        ! transforming ellipsoidal height to orthometric height:
        ! Hortho = -N + h_ell
        undul(i_ind) = u_grid(indx(i_ind))
        hgt = h_ell(i_stat)-undul(i_ind)
        
        ! pressure, temperature at the height of the grid
        T0 = T_grid(indx(i_ind),1) + T_grid(indx(i_ind),2)*cosfy + T_grid(indx(i_ind),3)*sinfy + T_grid(indx(i_ind),4)*coshy + T_grid(indx(i_ind),5)*sinhy
        p0 = p_grid(indx(i_ind),1) + p_grid(indx(i_ind),2)*cosfy + p_grid(indx(i_ind),3)*sinfy + p_grid(indx(i_ind),4)*coshy + p_grid(indx(i_ind),5)*sinhy

        ! humidity 
        Ql(i_ind) = Q_grid(indx(i_ind),1) + Q_grid(indx(i_ind),2)*cosfy + Q_grid(indx(i_ind),3)*sinfy + Q_grid(indx(i_ind),4)*coshy + Q_grid(indx(i_ind),5)*sinhy
 
        ! reduction = stationheight - gridheight
        Hs1 = Hs_grid(indx(i_ind))
        redh = hgt - Hs1

        ! lapse rate of the temperature in degree / m
        dTl(i_ind) = dT_grid(indx(i_ind),1) + dT_grid(indx(i_ind),2)*cosfy + dT_grid(indx(i_ind),3)*sinfy + dT_grid(indx(i_ind),4)*coshy + dT_grid(indx(i_ind),5)*sinhy 

        ! temperature reduction to station height
        Tl(i_ind) = T0 + dTl(i_ind)*redh - 273.15d0

        ! virtual temperature
        Tv = T0*(1+0.6077d0*Ql(i_ind))  
        c = gm*dMtr/(Rg*Tv)
            
        ! pressure in hPa
        pl(i_ind) = (p0*exp(-c*redh))/100.d0
            
        ! hydrostatic coefficient ah
        ahl(i_ind) = ah_grid(indx(i_ind),1) + &
                     ah_grid(indx(i_ind),2)*cosfy + ah_grid(indx(i_ind),3)*sinfy + &
                     ah_grid(indx(i_ind),4)*coshy + ah_grid(indx(i_ind),5)*sinhy
            
        ! wet coefficient aw
        awl(i_ind) = aw_grid(indx(i_ind),1) + &
                     aw_grid(indx(i_ind),2)*cosfy + aw_grid(indx(i_ind),3)*sinfy + &
                     aw_grid(indx(i_ind),4)*coshy + aw_grid(indx(i_ind),5)*sinhy

        ! water vapor decrease factor la - added by GP
        lal(i_ind) = la_grid(indx(i_ind),1) + &
                     la_grid(indx(i_ind),2)*cosfy + la_grid(indx(i_ind),3)*sinfy + &
                     la_grid(indx(i_ind),4)*coshy + la_grid(indx(i_ind),5)*sinhy

        ! mean temperature of the water vapor Tm - added by GP
        Tml(i_ind) = Tm_grid(indx(i_ind),1) + &
                     Tm_grid(indx(i_ind),2)*cosfy + Tm_grid(indx(i_ind),3)*sinfy + &
                     Tm_grid(indx(i_ind),4)*coshy + Tm_grid(indx(i_ind),5)*sinhy

        ! water vapor pressure in hPa - changed by GP
        e0 = Ql(i_ind)*p0/(0.622d0+0.378d0*Ql(i_ind))/100.d0      ! on the grid
        el(i_ind) = e0*(100.d0*pl(i_ind)/p0)**(lal(i_ind)+1.d0)   ! on the station height  (14) Askne and Nordius, 1987
            
        ! gradients
        Gn_hl(i_ind) = Gn_h_grid(indx(i_ind),1) + &
                       Gn_h_grid(indx(i_ind),2)*cosfy + Gn_h_grid(indx(i_ind),3)*sinfy + &
                       Gn_h_grid(indx(i_ind),4)*coshy + Gn_h_grid(indx(i_ind),5)*sinhy
        Ge_hl(i_ind) = Ge_h_grid(indx(i_ind),1) + &
                       Ge_h_grid(indx(i_ind),2)*cosfy + Ge_h_grid(indx(i_ind),3)*sinfy + &
                       Ge_h_grid(indx(i_ind),4)*coshy + Ge_h_grid(indx(i_ind),5)*sinhy
        Gn_wl(i_ind) = Gn_w_grid(indx(i_ind),1) + &
                       Gn_w_grid(indx(i_ind),2)*cosfy + Gn_w_grid(indx(i_ind),3)*sinfy + &
                       Gn_w_grid(indx(i_ind),4)*coshy + Gn_w_grid(indx(i_ind),5)*sinhy
        Ge_wl(i_ind) = Ge_w_grid(indx(i_ind),1) + &
                       Ge_w_grid(indx(i_ind),2)*cosfy + Ge_w_grid(indx(i_ind),3)*sinfy + &
                       Ge_w_grid(indx(i_ind),4)*coshy + Ge_w_grid(indx(i_ind),5)*sinhy
            
      end do

            
      dnpod1 = abs(diffpod)   ! distance nearer point
      dnpod2 = 1 - dnpod1     ! distance to distant point
      dnlon1 = abs(difflon)
      dnlon2 = 1 - dnlon1
        
      ! pressure
      R1 = dnpod2*pl(1)+dnpod1*pl(2)
      R2 = dnpod2*pl(3)+dnpod1*pl(4)
      p(i_stat) = dnlon2*R1+dnlon1*R2
            
      ! temperature
      R1 = dnpod2*Tl(1)+dnpod1*Tl(2)
      R2 = dnpod2*Tl(3)+dnpod1*Tl(4)
      T(i_stat) = dnlon2*R1+dnlon1*R2
        
      ! temperature in degree per km
      R1 = dnpod2*dTl(1)+dnpod1*dTl(2)
      R2 = dnpod2*dTl(3)+dnpod1*dTl(4)
      dT(i_stat) = (dnlon2*R1+dnlon1*R2)*1000.d0
            
      ! water vapor pressure in hPa
      R1 = dnpod2*el(1)+dnpod1*el(2)
      R2 = dnpod2*el(3)+dnpod1*el(4)
      e(i_stat) = dnlon2*R1+dnlon1*R2
            
      ! ah and aw
      R1 = dnpod2*ahl(1)+dnpod1*ahl(2)
      R2 = dnpod2*ahl(3)+dnpod1*ahl(4)
      ah(i_stat) = dnlon2*R1+dnlon1*R2
      R1 = dnpod2*awl(1)+dnpod1*awl(2)
      R2 = dnpod2*awl(3)+dnpod1*awl(4)
      aw(i_stat) = dnlon2*R1+dnlon1*R2
        
      ! undulation
      R1 = dnpod2*undul(1)+dnpod1*undul(2)
      R2 = dnpod2*undul(3)+dnpod1*undul(4)
      undu(i_stat) = dnlon2*R1+dnlon1*R2
      
      ! water vapour decrease factor
      R1 = dnpod2*lal(1)+dnpod1*lal(2)
      R2 = dnpod2*lal(3)+dnpod1*lal(4)
      la(i_stat) = dnlon2*R1+dnlon1*R2
        
      ! gradients
      R1 = dnpod2*Gn_hl(1)+dnpod1*Gn_hl(2)
      R2 = dnpod2*Gn_hl(3)+dnpod1*Gn_hl(4)
      Gn_h(i_stat) = (dnlon2*R1 + dnlon1*R2)
      R1 = dnpod2*Ge_hl(1)+dnpod1*Ge_hl(2)
      R2 = dnpod2*Ge_hl(3)+dnpod1*Ge_hl(4)
      Ge_h(i_stat) = (dnlon2*R1 + dnlon1*R2)
      R1 = dnpod2*Gn_wl(1)+dnpod1*Gn_wl(2)
      R2 = dnpod2*Gn_wl(3)+dnpod1*Gn_wl(4)
      Gn_w(i_stat) = (dnlon2*R1 + dnlon1*R2)
      R1 = dnpod2*Ge_wl(1)+dnpod1*Ge_wl(2)
      R2 = dnpod2*Ge_wl(3)+dnpod1*Ge_wl(4)
      Ge_w(i_stat) = (dnlon2*R1 + dnlon1*R2)

      ! mean temperature of the water vapor Tm
      R1 = dnpod2*Tml(1)+dnpod1*Tml(2)
      R2 = dnpod2*Tml(3)+dnpod1*Tml(4)
      Tm(i_stat) = dnlon2*R1+dnlon1*R2

    end if

  end do

end subroutine gpt3
