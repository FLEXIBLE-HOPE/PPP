cs:: @(#)wetpp.f	1.1  05/09/29 ::NGS


cb::wetpp
c
      real(rl) function wetpp( tc, rh )
c
c********1*********2*********3*********4*********5*********6*********7**
c
c name:         wetpp
c version:      9101.09
c written by:   M. Schenewerk
c purpose:      calculates the wet partial pressure in mbars given the
c               temperature and relative humidity
c
c input parameters:
c -----------------
c tc                temperature [c]
c rh                relative humidity
c
c output parameters:
c ------------------
c wetpp            wet parial pressure [millibars]
c
c
c local variables and constants:
c ------------------------------
c t0               centigrade - kelvins offset [degress temperature]
c tk               temperature [k]
c
c global variables and constants:
c ------------------------------
c
c
c calls:
c ------------------------------
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
c
c********1*********2*********3*********4*********5*********6*********7**
c::modification history
c::8707.13, jld, creation
c::8801.00, whd, correct wetpp calc as instructed by dsr after he
c::              brought the problem up with  tah
c::9101.09, MSS, standard header
c********1*********2*********3*********4*********5*********6*********7**
ce::wetpp
c
      use par
      use tables
      implicit real(rl) (a-h, o-z)
      implicit integer(it) (i-n)
c
c........  0.0  explicit initialization
c
      parameter ( t0= 273.15d0 )
c
      character sccsid*64
      sccsid='@(#)wetpp.f	1.1  05/09/29 ::NGS'//char(0) 
c
c........  1.0  absolute temperature
c
      tk = tc + t0
c
c........  2.0  saturation pressure in mbars
c
      wetpp = 6.11d0*( tk/t0 )**(-5.3)
     $    *dexp( 25.2d0*tc/tk )
c
c........  3.0  water vapor partial pressure in mbars
c
      wetpp = wetpp * rh
c
      return
      end
