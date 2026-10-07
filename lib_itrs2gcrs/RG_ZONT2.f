      SUBROUTINE RG_ZONT2 (CONV, T, DUT, DLOD, DOMEGA )
*+
*  - - - - - - - - - - -
*   R G _ Z O N T 2
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine evaluates the effects of zonal Earth tides on the
*  rotation of the Earth.  The model used is a combination of Yoder
*  et al. (1981) elastic body tide, Wahr and Bergen (1986) inelastic
*  body tide, and Kantha et al. (1998) ocean tide models 
*  as recommended by the IERS Conventions (2010).  Refer to
*  Chapter 8 pp. xx - xx.  The latest version of the model is located
*  at http://tai.bipm.org/iers/convupdt/convupdt_c8.html.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status:  Class 3 model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as
*     a Class 1, 2, or 3 model.
*
*  Given:
*     T           d      TT, Julian centuries since J2000 (Note 1)
*
*  Returned:
*     DUT         d      Effect on UT1 (Note 2)
*     DLOD        d      Effect on excess length of day (LOD) (Note 3)
*     DOMEGA      d      Effect on rotational speed (Note 4)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use
*     TT, which makes no significant difference.  Julian centuries since
*     J2000 is (JD - 2451545.0)/36525.
*
*  2) The expression used is as adopted in IERS Conventions (2010).
*     DUT is expressed in seconds and is double precision.
*
*  3) The expression used is as adopted in IERS Conventions (2010).
*     DLOD is the excess in LOD and is expressed in seconds per day
*     and is double precision.  The phrase 'per day' is generally
*     understood, so it has been omitted commonly in speech and
*     literature.  
*     See: Stephenson, F. R., Morrison, L. V., Whitrow, G. J., 1984,
*     "Long-Term Changes in the Rotation of the Earth: 700 B. C. to
*     A. D. 1980 [and Discussion]", Phil. Trans. Roy. Soc. of London.
*     Series A, 313, pp. 47 - 70.
* 
*  4) The expression used is as adopted in IERS Conventions (2010).
*     Rotational speed is expressed in radians per second and is
*     double precision.
*  
*  Called:
*     FUNDARG      Computation of the fundamental lunisolar arguments
*
*  Test case:
*     given input: T = .07995893223819302 Julian centuries since J2000
*                  (MJD = 54465)
*     expected output: DUT    =  7.983287678576557467E-002 seconds
*                      DLOD   =  5.035303035410713729E-005 seconds / day
*                      DOMEGA = -4.249711616463017E-014 radians / second
*
*  References:
*
*     Yoder, C. F., Williams, J. G., and Parke, M. E., (1981),
*     "Tidal Variations of Earth Rotation," J. Geophys. Res., 86,
*     pp. 881 - 891.
*
*     Wahr, J. and Bergen, Z., (1986), "The effects of mantle 
*     anelasticity on nutations, Earth tides, and tidal variations
*     in rotation rate," Geophys. J. Roy. astr. Soc., 87, pp. 633 - 668.
*
*     Kantha, L. H., Stewart, J. S., and Desai, S. D., (1998), "Long-
*     period lunar fortnightly and monthly ocean tides," J. Geophys.
*     Res., 103, pp. 12639 - 12647.
*
*     Gross, R. S., (2009), "Ocean tidal effects on Earth rotation,"
*     J. Geodyn., 48(3-5), pp. 219 - 225.
* 
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. xx, BKG (to be issued 2010)
*
*  Revisions:  
*  2008 January 18 B.E. Stetzler  Initial changes to header
*               and used 2PI instead of PI as parameter
*  2008 January 25 B.E. Stetzler Additional changes to header
*  2008 February 21 B.E. Stetzler Definition of (excess) LOD clarified
*  2008 March   12 B.E. Stetzler Applied changes to wording of notes.
*  2008 March   14 B.E. Stetzler Further changes applied to code.
*  2008 April   03 B.E. Stetzler Provided example test case
*  2009 February 11 B.E. Stetzler Updated test case due to changes made
*                                 to FUNDARG.F subroutine
*  2009 April   10 B.E. Stetzler DLOD corrected to say it is expressed
*                                in seconds per day
*  2009 May     04 B.E. Stetzler Code formatting changes based on 
*                                client recommendations
*  2009 May     07 B.E. Stetzler Updated test case due to above changes
*  2010 February 19 B.E. Stetzler Replaced Conventions 2003 recommended
*                                 model with Conventions 2010 model
*  2010 February 22 B.E. Stetzler Provided example test case
*  2010 February 23 B.E. Stetzler Updated values to two decimal places
*  2010 February 23 B.E. Stetzler Split fundamental arguments and
*                                 coefficients for four decimal place
*                                 precision
*  2010 February 25 B.E. Stetzler Recalculation of fundamental arguments
*  2010 March    01 B.E. Stetzler Updated table values to four decimal
*                                 places and double precision
*  2010 March    12 B.E. Stetzler Applied changes to wording of notes.
*  2010 March    22 B.E. Stetzler Corrected DOMEGA output for test case
*-----------------------------------------------------------------------
      USE tables
      IMPLICIT NONE

      INTEGER I, J
      DOUBLE PRECISION T, DUT, DLOD, DOMEGA, L, LP, F, D, OM, ARG, D2PI

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  2Pi
      PARAMETER (D2PI= 6.283185307179586476925287D0)

*  ----------------------
*  Zonal Earth tide model
*  ----------------------

*  Number of terms in the zonal Earth tide model  
      INTEGER NZONT
      PARAMETER ( NZONT = 62 )

*  Coefficients for the fundamental arguments
      INTEGER NFUND(5,NZONT)

*  Zonal tide term coefficients
      DOUBLE PRECISION TIDE(6,NZONT)

      CHARACTER(256) LINE
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
      SAVE LFIRST, NFUND, TIDE

      CHARACTER(LEN=*) CONV
      TYPE(T_FILETABLE) FT
      CHARACTER(LEN_FILENAME) UT1FILE
      INTEGER FLN, GET_VALID_UNIT

*  -----------------------------------------
*   Read the tabluar coefficients form file
*  -----------------------------------------
      IF (LFIRST .EQ. .TRUE.) THEN
        LFIRST = .FALSE.

        SELECT CASE(TRIM(CONV))
          CASE('IERS2010')
            UT1FILE = f_tablefilename('ut1t10')
          CASE('IERS1996','IERS2003')
            UT1FILE = f_tablefilename('ut1t03')
          CASE DEFAULT
        END SELECT

        FLN=GET_VALID_UNIT(10)
        OPEN(UNIT=FLN,FILE=UT1FILE,ACTION='READ')
        DO WHILE(INDEX(LINE,'-----------------').EQ.0)
          READ(FLN,'(A)',ERR=30) LINE
        END DO

        I=0
        DO WHILE(.TRUE.)
          READ(FLN,'(A)',ERR=30,END=30) LINE
          I=I+1
          READ(LINE(4:),*) (NFUND(J,I),J=1,5),(TIDE(J,I),J=1,6)
        END DO
30      CLOSE(FLN)
      END IF

*  -------------------------------------
*   Computation of fundamental arguments
*  -------------------------------------
      CALL FUNDARG(T,L,LP,F,D,OM)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*  Set initial values to zero.
      DUT    = 0.0D0
      DLOD   = 0.0D0
      DOMEGA = 0.0D0

*  Sum zonal tide terms.
      DO 10 I = 1, NZONT 
*     Formation of multiples of arguments.
         ARG =      MOD ( DBLE ( NFUND(1,I) ) * L 
     .       +            DBLE ( NFUND(2,I) ) * LP 
     .       +            DBLE ( NFUND(3,I) ) * F 
     .       +            DBLE ( NFUND(4,I) ) * D 
     .       +            DBLE ( NFUND(5,I) ) * OM, D2PI )

         IF (ARG.LT.0D0) ARG = ARG + D2PI

*     Evaluate zonal tidal terms.
         DUT    = DUT    + TIDE(1,I) *DSIN(ARG) + TIDE(2,I) *DCOS(ARG)
         DLOD   = DLOD   + TIDE(3,I) *DCOS(ARG) + TIDE(4,I) *DSIN(ARG)
         DOMEGA = DOMEGA + TIDE(5,I) *DCOS(ARG) + TIDE(6,I) *DSIN(ARG)
10    CONTINUE

*  Rescale corrections so that they are in units involving seconds.

      DUT    = DUT    * 1.0D-4
      DLOD   = DLOD   * 1.0D-5
      DOMEGA = DOMEGA * 1.0D-14

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
