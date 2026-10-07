      subroutine drysaas( bp, stlat, hgt, dryzen )
c
c********1*********2*********3*********4*********5*********6*********7**
c
c name:         drysaas
c version:      9101.09
c written by:   J. L. Davis
c purpose:      computes the dry zenith delay using the saastamoinen
c               formula
c
c input parameters:
c -----------------
c bp                barometric pressure [millibars]
c hgt               elliptic height [m]
c stlat             site latitude [rad]
c
c output parameters:
c ------------------
c dryzen()          dry atmosphere zenith path delay [m]
c
c
c local variables and constants:
c ------------------------------
c fphih             gravity correction to the delay
c
c global variables and constants:
c ------------------------------
c
c
c calls:
c ------------------------------
c phifun            computes the decrease in gravity due to
c                   height/latitude
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
c::8702.15, jld, creation
c::9101.09, MSS, standard header.  remove vlbi specific code
c********1*********2*********3*********4*********5*********6*********7**
ce::drysaas
c
      use par
      implicit real(rl) (a-h, o-z)
      implicit integer(it) (i-n)
c
      real(rl) dryzen(2)
c
c........  0.0  explicit initialization
c
      parameter( bprate= 0.d0 )
c
      character*64 sccsid
      sccsid='@(#)drysaas.f	1.1  05/09/29 ::NGS'//char(0) 
c
c........  1.0  gravity function
c
      fphih= phifun( stlat, hgt )
c
c........  2.0  path delay
c
      dryzen(1)= ( 0.0022768d0*bp )/fphih
      dryzen(2)= ( 0.0022768d0*bprate )/fphih
c
      return
      end
c
c
c
      real(rl) function phifun( stlat, hgt )
c
c********1*********2*********3*********4*********5*********6*********7**
c
c name:         phifun
c version:      9101.09
c written by:   j. l. davis
c purpose:      computes saastamoinen's function for the decrease
c               in gravity due to height/latitude
c
c input parameters:
c -----------------
c hgt               ellipsoidal height [m]
c stlat             site latitude [rad]
c
c output parameters:
c ------------------
c phifun            gravity correciton
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
c::9101.09, MSS, standard header.  remove vlbi specific code
c********1*********2*********3*********4*********5*********6*********7**
ce::phifun
c
      use par
      implicit real(rl) (a-h, o-z)
      implicit integer(it) (i-n)
c
      character*64 sccsid
      sccsid='@(#)phifun.f	1.1  05/09/29 ::NGS'//char(0) 
c
      phifun = 1.d0 - 0.00266d0*dcos(2.d0*stlat) - 0.00028d-03*hgt
c
      return
      end
