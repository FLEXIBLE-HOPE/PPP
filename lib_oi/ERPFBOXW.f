C*
      SUBROUTINE ERPFBOXW(BLKTYP,CPRN,CSVN,MASS,MJD,ANT,YSAT,SUN,ACCEL)
CC
CC NAME       :  ERPFBOXW
CC
CC PURPOSE    :  COMPUTATION OF EARTH RADIATION PRESSURE ACTING ON A
CC               BOW-WING SATELLITE
CC
CC PARAMETERS :
CC         IN :  BLKTYP    : SATELLITE TYPE
CC                             1 = GPS-I
CC                             2 = GPS-II
CC                             3 = GPS-IIA
CC                             4 = GPS-IIR
CC                             5 = GPS-IIR-A
CC                             6 = GPS-IIR-B
CC                             7 = GPS-IIR-M
CC                             8 = GPS-IIF
CC                             9 = GPS-III (not yet avaliable)
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
CC               CSVN      : SPACE VEHICLE NUMBER
CC               MASS      : MASS OF SATELLITE [kg]
CC               MJD       : MODIFIED JULIAN DAY
CC               ANT       : 0 = NO ANTENNA THRUST
CC                         : 1 = WITH ANTENNA THRUST
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
      CHARACTER(LEN=*) :: BLKTYP,CSVN,CPRN
      INTEGER(IT) :: ANT
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
C Top of Atmosphere for CERES
      TOA = 6371.000D0 + 30.000D0
C Albedo of the Earth
      ALB = 0.3D0

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

          CALL PROPBOXW(IBLK,CSVN,AREAS,REFLS,DIFUS,ABSPS,AREA2S,REFL2S,
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

C     B VECTOR
      B_SUN(1) = Y_SAT(2)*D_SUN(3) - Y_SAT(3)*D_SUN(2)
      B_SUN(2) = Y_SAT(3)*D_SUN(1) - Y_SAT(1)*D_SUN(3)
      B_SUN(3) = Y_SAT(1)*D_SUN(2) - Y_SAT(2)*D_SUN(1)

C     X VECTOR
      X_SAT(1) = Y_SAT(2)*Z_SAT(3) - Y_SAT(3)*Z_SAT(2)
      X_SAT(2) = Y_SAT(3)*Z_SAT(1) - Y_SAT(1)*Z_SAT(3)
      X_SAT(3) = Y_SAT(1)*Z_SAT(2) - Y_SAT(2)*Z_SAT(1)

      CALL rot_scfix2j2000(INT(MJD),(MJD-INT(MJD))*86400.d0,CPRN,CSVN,
     1                  BLKTYP,YSAT,SUN,X_SAT,Y_SAT,Z_SAT)

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

      ABSSUN = DSQRT(SUN(1)**2 + SUN(2)**2 + SUN(3)**2)
      DO K=1,3
        ESUN(K) = SUN(K)/ABSSUN
      ENDDO

      PSIDOT = ESUN(1)*RADVEC(1)+ESUN(2)*RADVEC(2)+ESUN(3)*RADVEC(3)
      IF(DABS(PSIDOT).GT.(1D0-1D-6))THEN
        PSI = 0D0
      ELSE
        PSI = DACOS(PSIDOT)
      ENDIF
      S0 = S0*(AU/ABSSUN)**2

C     ANALYTICAL MODEL

      NCFVEC(1) = RADVEC(1)
      NCFVEC(2) = RADVEC(2)
      NCFVEC(3) = RADVEC(3)
      ALBFAC = (PI*TOA**2)*(S0/VEL_LIGHT)/(ABSPOS**2)
      PHASEVI = (2*ALB/(3*PI**2))*((PI-PSI)*DCOS(PSI)+DSIN(PSI))
      PHASEIR = (1-ALB)/(4*PI)
      ABSNCFVI = ALBFAC*PHASEVI
      ABSNCFIR = ALBFAC*PHASEIR

      CALL SURFBOXW(AREAS,REFLS,DIFUS,ABSPS,
     1              AREA2S,REFL2S,DIFU2S,ABSP2S,
     2              REFLIRS,DIFUIRS,ABSPIRS,
     3              ABSNCFVI,ABSNCFIR,NCFVEC,ATTSURF,FORCE)


C     ANTENNA POWER OF GPS SATELLITES (IN WATTS)
C     IGS MODEL (JIM RAY, 2011)
C     FOR IGS REPROCESSION 3, THE DIFFERENT VALUES ARE USED

      SELECT CASE(TRIM(BLKTYP))
C     GPS BLOCK IIA (ASSUMED THE SAME FOR BLOCK I AND II)
        CASE('BLOCK I','BLOCK II','BLOCK IIA')
          ANTPOW = 50D0
C     GPS BLOCK IIR
        CASE('BLOCK IIR','BLOCK IIR-A','BLOCK IIR-B')
          ANTPOW = 60D0
C     GPS BLOCK IIR-M
        CASE('BLOCK IIR-M')
          ANTPOW = 145D0
C     GPS BLOCK IIF
        CASE('BLOCK IIF')
          ANTPOW = 240D0
C         NO M-CODE FOR SVN62/PRN25 STARTING 05APR2011; NANU 2011026
C          IF(CSVN(1:3).EQ.'062' .AND. (MJD.GE.55656D0))THEN
C            ANTPOW=154D0
C          ENDIF
C     GPS BLOCK IIIA
        CASE('BLOCK IIIA')
          ANTPOW = 300D0
C     NO ANTENNA POWER INFORMATION FOR GLONASS SATELLITES
C     THE ASSUMED VALUES USED AS CODE
        CASE('GLONASS')
          ANTPOW=100D0
        CASE('GLONASS-M')
          IF (CSVN(1:3).EQ.'735') THEN
            IF (MJD .LT. 57420.d0) THEN
              ANTPOW=20D0
            ELSE
              ANTPOW=25D0
            END IF
          ELSE IF (CSVN(1:3).EQ.'715' .OR. CSVN(1:3).EQ.'721' .OR.
     1             CSVN(1:3).EQ.'733' .OR. CSVN(1:3).EQ.'734' .OR.
     2             CSVN(1:3).EQ.'736') THEN
            ANTPOW=25D0
          ELSE IF (CSVN(1:3).EQ.'719') THEN
            ANTPOW=40D0
          ELSE IF (CSVN(1:3).EQ.'716' .OR. CSVN(1:3).EQ.'720') THEN
            ANTPOW=60D0
          ELSE IF (CSVN(1:3).EQ.'717' .OR. CSVN(1:3).EQ.'730' .OR.
     1             CSVN(1:3).EQ.'732') THEN
            ANTPOW=65D0
          ELSE IF (CSVN(1:3).EQ.'731' .OR. CSVN(1:3).EQ.'742' .OR.
     1             CSVN(1:3).EQ.'743' .OR. CSVN(1:3).EQ.'744' .OR.
     2             CSVN(1:3).EQ.'745' .OR. CSVN(1:3).EQ.'747') THEN
            ANTPOW=85D0
          ELSE IF (CSVN(1:3).EQ.'851' .OR. CSVN(1:3).EQ.'852' .OR.
     1             CSVN(1:3).EQ.'853' .OR. CSVN(1:3).EQ.'854' .OR.
     2             CSVN(1:3).EQ.'857') THEN
            ANTPOW=70D0
          ELSE IF (CSVN(1:3).EQ.'855') THEN
            ANTPOW=100D0
          ELSE IF (CSVN(1:3).EQ.'856') THEN
            ANTPOW=120D0
          ELSE IF (CSVN(1:3).EQ.'858') THEN
            ANTPOW=110D0
          ELSE
            ANTPOW=50D0
          END IF
        CASE('GLONASS-K1')
          IF (CSVN(1:3).EQ.'801') THEN
            ANTPOW=135D0
          ELSE IF (CSVN(1:3).EQ.'802') THEN
            ANTPOW=105D0
          ELSE
            ANTPOW=105D0
          END IF
        CASE('GALILEO-1')
          ANTPOW=155D0
          IF (CSVN(1:3).EQ.'101' .AND. (MJD.GE.56839D0)) THEN
            ANTPOW=135D0
          END IF
          IF (CSVN(1:3).EQ.'102' .AND. (MJD.GE.56839D0)) THEN
            ANTPOW=135D0
          END IF
          IF (CSVN(1:3).EQ.'103' .AND. (MJD.GE.56839D0)) THEN
            ANTPOW=95D0
          END IF
        CASE('GALILEO-2')
          ANTPOW=260D0
        CASE('QZSS')
          ANTPOW=250D0
        CASE('QZSS-2I')
          IF (CPRN(1:3) .EQ. 'J02') THEN
            ANTPOW=500D0
          ELSE
            ANTPOW=550D0
          END IF
        CASE('QZSS-2G')
          ANTPOW=550D0
        CASE('QZSS-2A')
          ANTPOW=460D0
        CASE('BEIDOU-2I')
          ANTPOW=185D0
        CASE('BEIDOU-2M')
          ANTPOW=130D0
        CASE('BEIDOU-3MS-CAST')
          ANTPOW=250D0
        CASE('BEIDOU-3M-CAST')
          ANTPOW=310D0
        !!!!!! JG: we notice a bias in SLR validation, so 
        CASE('BEIDOU-3M-SECM')
          ANTPOW=280D0
          ! ANTPOW=0.D0
        CASE DEFAULT
          ANTPOW=0.D0
      END SELECT


C     NAVIGATION ANTENNA THRUST (SIMPLE MODEL)
      IF (ANT .EQ. 1)THEN
        ANTFORCE = ANTPOW/VEL_LIGHT
        FORCE(1) = FORCE(1) + ANTFORCE*RADVEC(1)
        FORCE(2) = FORCE(2) + ANTFORCE*RADVEC(2)
        FORCE(3) = FORCE(3) + ANTFORCE*RADVEC(3)
      ENDIF

C     CONVERSION TO ACCELERATION
      DO K=1,3
        ACCEL(K) = ACCEL(K)+FORCE(K)/MASS*1.D-3
      ENDDO
c      write(3000,"(F23.12,1X,A3,1X,E23.12)") mjd, cprn, 
c     1   dsqrt(force(1)**2+force(2)**2+force(3)**2)/mass*1.d-3

      RETURN

      END SUBROUTINE
