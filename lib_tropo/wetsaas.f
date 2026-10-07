cs:: @(#)wetsaas.f	1.1  05/09/29 ::NGS


cb::wetchao
c
      subroutine wetsaas( rh, tc, stlat, hgt, wetzen )
c
c********1*********2*********3*********4*********5*********6*********7**
c
c name:         wetchao
c version:      9101.09
c written by:   t. l. davis
c purpose:      computes wet zenith delay using saastamoinen formula
c
c
c input parameters:
c -----------------
c rh                relative humidity
c tc                temperature [c]
c stlat             site latitude [rad]
c hgt               ellipsoidal height [m]
c
c output parameters:
c ------------------
c wetzen()          zenith delay [m]
c
c
c local variables and constants:
c ------------------------------
c
c global variables and constants:
c ------------------------------
c
c
c calls:
c ------------------------------
c phifun            computes the decrease in gravity due to height/latitude
c wetchao           Chao wet mapping function
c wetpp             the wet partial pressure from temperature and humidity
c
c include files:
c ------------------------------
c
c common blocks:
c ------------------------------
c
c references:
c ------------------------------
c
c comments:
c ------------------------------
c in computing the rate, we currently
c neglect the temperature derivative
c
c********1*********2*********3*********4*********5*********6*********7**
c::modification history
c::8702.15, tld, creation
c::9101.09, MSS, standard header.  remove vlbi specific code
c********1*********2*********3*********4*********5*********6*********7**
ce::wetchao
c
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
      double precision wetzen(2)
c
c........  0.0  explicit initialization
c
      parameter( t0= 273.15d0 )
      parameter( rhrate= 0.d0 )
      parameter( tcrate= 0.d0 )
c
      character*64 sccsid
      sccsid='@(#)wetsaas.f	1.1  05/09/29 ::NGS'//char(0) 
c
c........  1.0  gravity function
c
      fphif= phifun( stlat, hgt )
c
c........  2.0  water vapor partial pressure
c
      pp= wetpp( tc, 1.0d0 )
c
c........  3.0  calculate parameters
c               note: 0.11862 = 6.11*5.3/273
c
      dwetdt= ( -0.11862*( ( tc + t0 )/273 )**(-6.3)
     $    *dexp( 25.2d0*tc/( tc + t0 ) )
     $    + 25.2d0*pp/( tc + t0 )
     $    *( 1.d0 - tc/( tc + t0 ) ) )
     $    *tcrate
c
c........  4.0  zenith delay
c
      wetzen(1)= 0.0022768d0*( 1255.d0/( tc + t0 ) + 0.05d0 )*
     $    rh*pp/fphif
c
      wetzen(2)= 0.0022768d0*( 1255.d0/( tc + t0 ) + 0.05d0 )*
     $    ( rhrate*pp + rh*dwetdt )/fphif
c
      return
      end
