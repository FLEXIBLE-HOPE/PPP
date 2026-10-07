C*
      SUBROUTINE SRPFBOXW(BLKTYP,CPRN,MASS,MJD,YSAT,SUN,ACCEL)
CC
CC NAME       :  SRPFBOXW
CC
CC PURPOSE    :  COMPUTATION OF SOLAR RADIATION PRESSURE ACTING ON A
CC               BOW-WING SATELLITE
CC
CC PARAMETERS :
CC         IN :  BLKTYP    : SATELLITE TYPE
CC         IN :  BLKNUM     : BLOCK NUMBER
CC                             1 = GPS-I
CC                             2 = GPS-II
CC                             3 = GPS-IIA
CC                             4 = GPS-IIR
CC                             5 = GPS-IIR-A
CC                             6 = GPS-IIR-B
CC                             7 = GPS-IIR-M
CC                             8 = GPS-IIF
CC                             9 = GPS-IIIA
CC                            10 = GLONASS
CC                            11 = GLONASS-M
CC                            12 = GLONASS-K1 (not yet avaliable)
CC                            13 = GALILEO-0A (not yet avaliable)
CC                            14 = GALILEO-0B (not yet avaliable)
CC                            15 = GALILEO-1 (not yet avaliable)
CC                            16 = GALILEO-2 (not yet avaliable)
CC                            17 = BEIDOU-2G
CC                            18 = BEIDOU-2I
CC                            19 = BEIDOU-2M
CC                            20 = BEIDOU-3IS-SECM
CC                            21 = BEIDOU-3IS-CAST
CC                            22 = BEIDOU-3MS-SECM
CC                            23 = BEIDOU-3MS-CAST
CC                            24 = BEIDOU-3G-SECM
CC                            25 = BEIDOU-3G-CAST
CC                            26 = BEIDOU-3I-SECM
CC                            27 = BEIDOU-3I-CAST
CC                            28 = BEIDOU-3M-SECM
CC                            29 = BEIDOU-3M-CAST
CC                            30 = QZSS (not yet avaliable)
CC                            31 = QZSS-2I (not yet avaliable)
CC                            32 = QZSS-2G (not yet avaliable)
CC                            33 = QZSS-2A (not yet avaliable)
CC                            34 = IRNSS-1IGSO (not yet avaliable)
CC                            35 = IRNSS-1GEO (not yet avaliable)
CC               CPRN      : SPACE VEHICLE NUMBER
CC               MASS      : MASS OF SATELLITE [kg]
CC               MJD       : MODIFIED JULIAN DAY
CC               YSAT      : SATELLITE POSITION [km] AND VELOCITY [km/s] (INERTIAL),
CC                           (POSITION,VELOCITY) = (RX,RY,RZ,VX,VY,VZ)
CC               SUN       : SUN POSITION VECTOR [km] (INERTIAL)
CC
CC        OUT : ACCEL      : ACCELERATION VECTOR [km/s^2]
CC
CC AUTHOR     : C.J. RODRIGUEZ-SOLANO
CC              rodriguez@bv.tum.de
CC
CC VERSION    : 1.0 (OCT 2010)
CC
CC CREATED    : 2010/10/18
CC
CC MODIFICATION : JING GUO
CC                jingguo@whu.edu.cn
CC
C*
C
      USE CONST
      IMPLICIT NONE
C
C*
C THE ARGUMENTS
C -------------------------------------
      CHARACTER(LEN=*) :: BLKTYP,CPRN
      REAL(RL) :: MASS,MJD,YSAT(1:*),SUN(1:*),ACCEL(1:*)

C*
C THE LOCAL VARIABLES
C --------------------------------------
      INTEGER(IT) :: IBLK,II,JJ,K,INDB
      INTEGER(IT), PARAMETER :: MAXBLK=36

C
C CONST
      REAL(RL) :: AU,S0,TOA,ALB
C
C PROPRETIES OF SATELLITES
      REAL(RL) :: AREA(4,2,MAXBLK),REFL(4,2,MAXBLK)
      REAL(RL) :: DIFU(4,2,MAXBLK),ABSP(4,2,MAXBLK)
      REAL(RL) :: AREA2(4,2,MAXBLK),REFL2(4,2,MAXBLK)
      REAL(RL) :: DIFU2(4,2,MAXBLK),ABSP2(4,2,MAXBLK)
      REAL(RL) :: REFLIR(4,2,MAXBLK),DIFUIR(4,2,MAXBLK)
      REAL(RL) :: ABSPIR(4,2,MAXBLK)
      REAL(RL) :: AREAS(4,2),REFLS(4,2),DIFUS(4,2),ABSPS(4,2)
      REAL(RL) :: AREA2S(4,2),REFL2S(4,2),DIFU2S(4,2),ABSP2S(4,2)
      REAL(RL) :: REFLIRS(4,2),DIFUIRS(4,2),ABSPIRS(4,2)
C
C     ATTITUDE
      REAL(RL) :: RADVEC(3),ALGVEC(3),CRSVEC(3),ESUN(3)
      REAL(RL) :: RSUN,ABSPOS,ABSCRS,ABSY0V,Y0SAT(3)
      REAL(RL) :: Z_SAT(3),D_SUN(3),Y_SAT(3),B_SUN(3),X_SAT(3)
      REAL(RL) :: ATTSURF(3,4),FORCE(3)

C
      REAL(RL) :: ABSNCFVI,ABSNCFIR,ALBFAC,PHASEVI,PHASEIR
      REAL(RL) :: NCFVEC(3),PSI,PSIDOT,ABSSUN,ANTFORCE,ANTPOW

      LOGICAL(LG) :: LFIRST

      DATA LFIRST /.TRUE./
      SAVE LFIRST,AREA,REFL,DIFU,ABSP,AREA2,REFL2,DIFU2,ABSP2,
     1     REFLIR,DIFUIR,ABSPIR

C*
C START THE EXECTUABLE CODE
C -------------------------------------

C Constants needed
      AU = 149597870.691D0
C Solar Constant(W/m2)(W=kg*m^2/s^3)
      S0 = 1367D0

C Initialization of force vector
      FORCE(1) = 0D0
      FORCE(2) = 0D0
      FORCE(3) = 0D0

C ----------------------------------------
C LOAD SATELLITE PROPERTIES
C ----------------------------------------
      IF (LFIRST .EQ. .TRUE.) THEN

        LFIRST=.FALSE.

C     PROPERTIES FOR ALL SATELLITES BLOCKS
        DO IBLK=1, MAXBLK

          CALL PROPBOXW(IBLK,AREAS,REFLS,DIFUS,ABSPS,AREA2S,REFL2S,
     1                  DIFU2S,ABSP2S,REFLIRS,DIFUIRS,ABSPIRS)

          DO II = 1,4
            DO JJ = 1,2
              AREA(II,JJ,IBLK) = AREAS(II,JJ)
              REFL(II,JJ,IBLK) = REFLS(II,JJ)
              DIFU(II,JJ,IBLK) = DIFUS(II,JJ)
              ABSP(II,JJ,IBLK) = ABSPS(II,JJ)

              AREA2(II,JJ,IBLK) = AREA2S(II,JJ)
              REFL2(II,JJ,IBLK) = REFL2S(II,JJ)
              DIFU2(II,JJ,IBLK) = DIFU2S(II,JJ)
              ABSP2(II,JJ,IBLK) = ABSP2S(II,JJ)

              REFLIR(II,JJ,IBLK) = REFLIRS(II,JJ)
              DIFUIR(II,JJ,IBLK) = DIFUIRS(II,JJ)
              ABSPIR(II,JJ,IBLK) = ABSPIRS(II,JJ)
            ENDDO
          ENDDO
        ENDDO
      ENDIF


C --------------------------
C NOMINAL SATELLITE ATTITUDE
C --------------------------
C JG: THIS SHOULD BE CHANGED FOR USING ATTITUDE MODE

      ABSPOS = DSQRT(YSAT(1)**2+YSAT(2)**2+YSAT(3)**2)
      DO K=1,3
         RADVEC(K) = YSAT(K)/ABSPOS
      ENDDO

      CRSVEC(1) = YSAT(2)*YSAT(6)-YSAT(3)*YSAT(5)
      CRSVEC(2) = YSAT(3)*YSAT(4)-YSAT(1)*YSAT(6)
      CRSVEC(3) = YSAT(1)*YSAT(5)-YSAT(2)*YSAT(4)
      ABSCRS = DSQRT(CRSVEC(1)**2+CRSVEC(2)**2+CRSVEC(3)**2)
      DO K=1,3
         CRSVEC(K) = CRSVEC(K)/ABSCRS
      ENDDO

      ALGVEC(1) = CRSVEC(2)*RADVEC(3)-CRSVEC(3)*RADVEC(2)
      ALGVEC(2) = CRSVEC(3)*RADVEC(1)-CRSVEC(1)*RADVEC(3)
      ALGVEC(3) = CRSVEC(1)*RADVEC(2)-CRSVEC(2)*RADVEC(1)

C     DISTANCE FROM SATELLITE TO SUN
      RSUN = DSQRT((YSAT(1)-SUN(1))**2+(YSAT(2)-SUN(2))**2+
     1       (YSAT(3)-SUN(3))**2)

C     D VECTOR AND Z VECTOR
      DO K=1,3
        D_SUN(K) = (SUN(K)-YSAT(K))/RSUN
        Z_SAT(K) = -YSAT(K)/ABSPOS
      ENDDO

C     Y VECTOR
      Y0SAT(1) = Z_SAT(2)*D_SUN(3)-Z_SAT(3)*D_SUN(2)
      Y0SAT(2) = Z_SAT(3)*D_SUN(1)-Z_SAT(1)*D_SUN(3)
      Y0SAT(3) = Z_SAT(1)*D_SUN(2)-Z_SAT(2)*D_SUN(1)
      ABSY0V = DSQRT(Y0SAT(1)**2 + Y0SAT(2)**2 + Y0SAT(3)**2)
      DO K=1,3
        Y_SAT(K) = Y0SAT(K)/ABSY0V
      ENDDO

C     B VECTOR, -B WITH BERN
      B_SUN(1) = Y_SAT(2)*D_SUN(3) - Y_SAT(3)*D_SUN(2)
      B_SUN(2) = Y_SAT(3)*D_SUN(1) - Y_SAT(1)*D_SUN(3)
      B_SUN(3) = Y_SAT(1)*D_SUN(2) - Y_SAT(2)*D_SUN(1)

C     X VECTOR
      X_SAT(1) = Y_SAT(2)*Z_SAT(3) - Y_SAT(3)*Z_SAT(2)
      X_SAT(2) = Y_SAT(3)*Z_SAT(1) - Y_SAT(1)*Z_SAT(3)
      X_SAT(3) = Y_SAT(1)*Z_SAT(2) - Y_SAT(2)*Z_SAT(1)

      DO K=1,3
        ATTSURF(K,1) = Z_SAT(K)
        ATTSURF(K,2) = Y_SAT(K)
        ATTSURF(K,3) = X_SAT(K)
        ATTSURF(K,4) = D_SUN(K)
      ENDDO

C ----------------------------
C OPTICAL PROPERTIES PER BLOCK
C ----------------------------
      SELECT CASE(TRIM(BLKTYP))
        CASE('BLOCK I')
          INDB=1
        CASE('BLOCK II')
          INDB=2
        CASE('BLOCK IIA')
          INDB=3
        CASE('BLOCK IIR')
          INDB=4
        CASE('BLOCK IIR-A')
          INDB=5
        CASE('BLOCK IIR-B')
          INDB=6
        CASE('BLOCK IIR-M')
          INDB=7
        CASE('BLOCK IIF')
          INDB=8
        CASE('BLOCK IIIA')
          INDB=9
        CASE('GLONASS')
          INDB=10
        CASE('GLONASS-M')
          INDB=11
        CASE('GLONASS-K1')
          INDB=12
        CASE('GALILEO-0A')
          INDB=13
        CASE('GALILEO-0B')
          INDB=14
        CASE('GALILEO-1')
          INDB=15
        CASE('GALILEO-2')
          INDB=16
        CASE('BEIDOU-2G')
          INDB=17
        CASE('BEIDOU-2I')
          INDB=18
        CASE('BEIDOU-2M')
          INDB=19
        CASE('BEIDOU-3IS-SECM')
          INDB=20
        CASE('BEIDOU-3IS-CAST')
          INDB=21
        CASE('BEIDOU-3MS-SECM')
          INDB=22
        CASE('BEIDOU-3MS-CAST')
          INDB=23
        CASE('BEIDOU-3G-SECM')
          INDB=24
        CASE('BEIDOU-3G-CAST')
          INDB=25
        CASE('BEIDOU-3I-SECM')
          INDB=26
        CASE('BEIDOU-3I-CAST')
          INDB=27
        CASE('BEIDOU-3M-SECM')
          INDB=28
        CASE('BEIDOU-3M-CAST')
          INDB=29
        CASE('QZSS')
          INDB=30
        CASE('QZSS-2I')
          INDB=31
        CASE('QZSS-2G')
          INDB=32
        CASE('QZSS-2A')
          INDB=33
        CASE('IRNSS-1IGSO')
          INDB=34
        CASE('IRNSS-1GEO')
          INDB=35
      END SELECT

      DO II = 1,4
        DO JJ = 1,2
          AREAS(II,JJ) = AREA(II,JJ,INDB)
          REFLS(II,JJ) = REFL(II,JJ,INDB)
          DIFUS(II,JJ) = DIFU(II,JJ,INDB)
          ABSPS(II,JJ) = ABSP(II,JJ,INDB)

          AREA2S(II,JJ) = AREA2(II,JJ,INDB)
          REFL2S(II,JJ) = REFL2(II,JJ,INDB)
          DIFU2S(II,JJ) = DIFU2(II,JJ,INDB)
          ABSP2S(II,JJ) = ABSP2(II,JJ,INDB)

          REFLIRS(II,JJ) = REFLIR(II,JJ,INDB)
          DIFUIRS(II,JJ) = DIFUIR(II,JJ,INDB)
          ABSPIRS(II,JJ) = ABSPIR(II,JJ,INDB)
        ENDDO
      ENDDO

C ----------------------
C EARTH RADIATION MODELS
C ----------------------

      ABSSUN = DSQRT((SUN(1)-YSAT(1))**2 + (SUN(2)-YSAT(2))**2 +
     1               (SUN(3)-YSAT(3))**2)

      S0 = -S0*(AU/ABSSUN)**2

C     ANALYTICAL MODEL

      NCFVEC(1) = D_SUN(1)
      NCFVEC(2) = D_SUN(2)
      NCFVEC(3) = D_SUN(3)

      ABSNCFVI = S0/VEL_LIGHT
      ABSNCFIR = 0D0

      CALL SURFBOXW(AREAS,REFLS,DIFUS,ABSPS,
     1              AREA2S,REFL2S,DIFU2S,ABSP2S,
     2              REFLIRS,DIFUIRS,ABSPIRS,
     3              ABSNCFVI,ABSNCFIR,NCFVEC,ATTSURF,FORCE)

C     CONVERSION TO ACCELERATION
      DO K=1,3
        ACCEL(K) = ACCEL(K)+FORCE(K)/MASS*1.D-3
      ENDDO
c      write(3000,"(F23.12,1X,A3,1X,E23.12)") mjd, cprn,
c     1   dsqrt(force(1)**2+force(2)**2+force(3)**2)/mass*1.d-3

      RETURN

      END SUBROUTINE
