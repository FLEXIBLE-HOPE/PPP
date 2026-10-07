C*
      SUBROUTINE PROPBOXW(BLKNUM,CSVN,AREA,REFL,DIFU,ABSP,AREA2,REFL2,DIFU2,
     1                    ABSP2,REFLIR,DIFUIR,ABSPIR)
CC
CC NAME       :  PROPBOXW
CC
CC PURPOSE    :  SATELLITE DIMENSIONS AND OPTICAL PROPERTIES FROM ROCK MODELS
CC               BOX-WING MODELS FOR GPS AND GLONASS SATELLITES
CC
CC REFERENCES :  Fliegel H, Gallini T, Swift E (1992) Global Positioning System Radiation
CC                  Force Model for Geodetic Applications. Journal of Geophysical Research
CC                  97(B1): 559-568
CC               Fliegel H, Gallini T (1996) Solar Force Modelling of Block IIR Global
CC                  Positioning System satellites. Journal of Spacecraft and Rockets
CC                  33(6): 863-866
CC               Ziebart M (2001) High Precision Analytical Solar Radiation Pressure
CC                  Modelling for GNSS Spacecraft. PhD Thesis, University of East London
CC               [IGEXMAIL-0086] GLONASS S/C mass and dimension
CC               [IGSMAIL-5104] GLONASS-M dimensions and center-of-mass correction
CC               http://acc.igs.org/orbits/IIF_SV_DimensionsConfiguration.ppt
CC
CC PARAMETERS : 
CC         IN :  BLKNUM     : BLOCK NUMBER
CC                             1 = GPS-I
CC                             2 = GPS-II
CC                             3 = GPS-IIA
CC                             4 = GPS-IIR
CC                             5 = GPS-IIR-A
CC                             6 = GPS-IIR-B
CC                             7 = GPS-IIR-M
CC                             8 = GPS-IIF
CC                             9 = GPS-IIIA (not yet avaliable)
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
CC                            33 = QZSS-2G (not yet avaliable)
CC                            34 = IRNSS-1IGSO (not yet avaliable)
CC                            35 = IRNSS-1GEO (not yet avaliable)
CC        OUT :  AREA(I,J)  : AREAS OF FLAT SURFACES [m^2]
CC               REFL(I,J)  : REFLEXION COEFFICIENT
CC               DIFU(I,J)  : DIFFUSION COEFFICIENT
CC               ABSP(I,J)  : ABSORPTION COEFFICIENT
CC                            I = 1 +Z (TOWARDS THE EARTH)
CC                            I = 2 +Y (ALONG SOLAR PANELS BEAMS)
CC                            I = 3 +X (ALONG BUS DIRECTION ALWAYS ILLUMINATED BY THE SUN)
CC                            I = 4 SOLAR PANELS
CC                            J = 1 POSITIVE DIRECTION
CC                            J = 2 NEGATIVE DIRECTION
CC               ????2(I,J) : INDEX 2 INDICATES AREAS AND OPTICAL PROPERTIES OF CYLINDRICAL
CC                            SURFACES, SAME MEANING OF (I,J) AS BEFORE
CC               ????IR(I,J): OPTICAL PROPERTIES IN THE INFRARED (ASSUMED), NO SEPARATION
CC                            BETWEEN FLAT AND CYLINDRICAL SURFACES
CC
CC AUTHOR     :  C.J. RODRIGUEZ-SOLANO
CC               rodriguez@bv.tum.de
CC
CC VERSION    :  1.0 (OCT 2010)
CC
CC CREATED    :  2010/10/18             LAST MODIFIED :  17-MAR-11
CC
CC CHANGES    :  17-MAR-11 : CR: ADD BOX-WING MODELS FOR GLONASS, GLONASS-M, GPS II-F
C*
      USE PAR
      IMPLICIT NONE
C
      INTEGER(IT) BLKNUM,II,JJ,KK,SS
      CHARACTER(LEN=*) :: CSVN
C
      REAL(RL) AREA(4,2),REFL(4,2),DIFU(4,2),ABSP(4,2)
      REAL(RL) AREA2(4,2),REFL2(4,2),DIFU2(4,2),ABSP2(4,2)
      REAL(RL) REFLIR(4,2),DIFUIR(4,2),ABSPIR(4,2)
      REAL(RL) X_SIDE(5,4),Z_SIDE(4,4),S_SIDE(2,4),Y_SIDE(2,4)
      REAL(RL) SURFALL(14,4)
      REAL(RL) ALLREFL(14),ALLDIFU(14),REFL_AREA(14),DIFU_AREA(14)
      REAL(RL) G_AREA1(5),G_REFL1(5),G_DIFU1(5),G_ABSP1(5)
      REAL(RL) G_AREA2(5),G_REFL2(5),G_DIFU2(5),G_ABSP2(5)
      REAL(RL) BUSFAC

C SATELLITE PROPERTIES FROM ROCK MODEL
C VALUES GIVEN FOR +X, -Z, +Z AND SOLAR PANELS
C X_SIDE(I,J) :   I     DIFFERENT SURFACE COMPONENTS
C                 J = 1 AREA
C                 J = 2 SPECULARITY
C                 J = 3 REFLECTIVITY
C                 J = 4 SHAPE: 1 = flat, 2 = cylindrical

      DO JJ = 1,4
         DO II = 1,5
            X_SIDE(II,JJ) = 0D0
         ENDDO
         DO II = 1,4
            Z_SIDE(II,JJ) = 0D0
         ENDDO
         DO II = 1,2
            S_SIDE(II,JJ) = 0D0
         ENDDO
         DO II = 1,2
            Y_SIDE(II,JJ) = 0D0
         ENDDO
      ENDDO


C     -----------
C     GPS BLOCK I
C     -----------
C     SEE FLIEGEL ET AL (1992)
      IF(BLKNUM.EQ.1)THEN
          
C     +X SIDE
          X_SIDE(1,1) = 1.055D0
          X_SIDE(1,2) = 0.80D0
          X_SIDE(1,3) = 0.50D0
          X_SIDE(1,4) = 1D0

C     ENGINE SIDE
          X_SIDE(2,1) = 0.570D0
          X_SIDE(2,2) = 0.75D0
          X_SIDE(2,3) = 0.86D0
          X_SIDE(2,4) = 2D0

C     TT&C ANTENNA SIDE
          X_SIDE(3,1) = 0.055D0
          X_SIDE(3,2) = 0.05D0
          X_SIDE(3,3) = 0.28D0
          X_SIDE(3,4) = 2D0

C     TT&C ANTENNA TIP
          X_SIDE(4,1) = 0.019D0
          X_SIDE(4,2) = 0.85D0
          X_SIDE(4,3) = 0.28D0
          X_SIDE(4,4) = 2D0

C     EACH NAVIGATIONAL ANTENNA ADAPTER
          X_SIDE(5,1) = 0.029D0
          X_SIDE(5,2) = 0.75D0
          X_SIDE(5,3) = 0.36D0
          X_SIDE(5,4) = 2D0

C     BODY AFT END (-Z)
          Z_SIDE(1,1) = 0.816D0
          Z_SIDE(1,2) = 0.80D0
          Z_SIDE(1,3) = 0.86D0
          Z_SIDE(1,4) = 1D0

C     ENGINE AFT END (-Z)                                                         
          Z_SIDE(2,1) = 0.694D0
          Z_SIDE(2,2) = 0D0
          Z_SIDE(2,3) = 0D0
          Z_SIDE(2,4) = 1D0

C     FORWARD END (+Z)                                                                     
          Z_SIDE(3,1) = 1.510D0
          Z_SIDE(3,2) = 0.75D0
          Z_SIDE(3,3) = 0.86D0
          Z_SIDE(3,4) = 1D0

C     ALL SOLAR PANELS
          S_SIDE(1,1) = 5.583D0
          S_SIDE(1,2) = 0.85D0
          S_SIDE(1,3) = 0.23D0
          S_SIDE(1,4) = 1D0

C     SOLAR PANELS MASTS                                                                                     
          S_SIDE(2,1) = 0.470D0
          S_SIDE(2,2) = 0.85D0
          S_SIDE(2,3) = 0.85D0
          S_SIDE(2,4) = 2D0

C     -----------------
C     GPS BLOCK II, IIA
C     -----------------
C     SEE FLIEGEL ET AL (1992)
      ELSEIF((BLKNUM.EQ.2).OR.(BLKNUM.EQ.3))THEN

C     +X SIDE                                                                                           
          X_SIDE(1,1) = 1.553D0
          X_SIDE(1,2) = 0.20D0
          X_SIDE(1,3) = 0.56D0
          X_SIDE(1,4) = 1D0

C     ENGINE SIDE (INCLUDES PLUME SHILED 0.22*1.84 m^2)                                            
          X_SIDE(2,1) = 1.054D0
          X_SIDE(2,2) = 0.20D0
          X_SIDE(2,3) = 0.56D0
          X_SIDE(2,4) = 2D0

C     TT&C ANTENNA SIDE                                                                       
          X_SIDE(3,1) = 0.105D0
          X_SIDE(3,2) = 0.20D0
          X_SIDE(3,3) = 0.28D0
          X_SIDE(3,4) = 2D0

C     EACH NAVIGATIONAL ANTENNA ADAPTER                                                               
          X_SIDE(5,1) = 0.181D0
          X_SIDE(5,2) = 0.20D0
          X_SIDE(5,3) = 0.36D0
          X_SIDE(5,4) = 2D0

C     BODY AFT END (-Z)                                                                               
          Z_SIDE(1,1) = 2.152D0
          Z_SIDE(1,2) = 0.20D0
          Z_SIDE(1,3) = 0.56D0
          Z_SIDE(1,4) = 1D0

C     ENGINE AFT END (-Z)                                                                          
          Z_SIDE(2,1) = 0.729D0
          Z_SIDE(2,2) = 0D0
          Z_SIDE(2,3) = 0D0
          Z_SIDE(2,4) = 1D0

C     FORWARD END (+Z)                                                                                    
          Z_SIDE(3,1) = 2.881D0
          Z_SIDE(3,2) = 0.20D0
          Z_SIDE(3,3) = 0.56D0
          Z_SIDE(3,4) = 1D0

C     ALL SOLAR PANELS                                                                                     
          S_SIDE(1,1) = 10.886D0
          S_SIDE(1,2) = 0.85D0
          S_SIDE(1,3) = 0.23D0
          S_SIDE(1,4) = 1D0

C     SOLAR PANELS MASTS                                                                    
          S_SIDE(2,1) = 0.985D0
          S_SIDE(2,2) = 0.41D0
          S_SIDE(2,3) = 0.52D0
          S_SIDE(2,4) = 2D0

C     ----------------------------------
C     GPS BLOCK IIR, IIR-A, IIR-B, IIR-M
C     ----------------------------------
C     SEE FLIEGEL AND GALINI (1996)
      ELSEIF((BLKNUM.GE.4).AND.(BLKNUM.LE.7))THEN

C     + AND -X FACES                                                                                
         X_SIDE(1,1) = 3.05D0
         X_SIDE(1,2) = 0D0
         X_SIDE(1,3) = 0.06D0
         X_SIDE(1,4) = 1D0

C     PLUME SHILED                                                                                  
         X_SIDE(2,1) = 0.17D0
         X_SIDE(2,2) = 0D0
         X_SIDE(2,3) = 0.06D0
         X_SIDE(2,4) = 2D0

C     ANTENNA SHROUD                                                                              
         X_SIDE(3,1) = 0.89D0
         X_SIDE(3,2) = 0D0
         X_SIDE(3,3) = 0.06D0
         X_SIDE(3,4) = 2D0

C     -Z FACE                                                                                        
         Z_SIDE(1,1) = 3.75D0
         Z_SIDE(1,2) = 0D0
         Z_SIDE(1,3) = 0.06D0
         Z_SIDE(1,4) = 1D0

C     -Z W-SENSOR
         Z_SIDE(2,1) = 0.50D0
         Z_SIDE(2,2) = 0D0
         Z_SIDE(2,3) = 0.06D0
         Z_SIDE(2,4) = 1D0

C     +Z FACE                                                                                              
         Z_SIDE(3,1) = 3.75D0
         Z_SIDE(3,2) = 0D0
         Z_SIDE(3,3) = 0.06D0
         Z_SIDE(3,4) = 1D0

C     +Z W-SENSOR                                                                                             
         Z_SIDE(4,1) = 0.50D0
         Z_SIDE(4,2) = 0D0
         Z_SIDE(4,3) = 0.06D0
         Z_SIDE(4,4) = 1D0

C     ALL SOLAR PANELS                                                                                
         S_SIDE(1,1) = 13.60D0
         S_SIDE(1,2) = 0.85D0
         S_SIDE(1,3) = 0.28D0
         S_SIDE(1,4) = 1D0

C     SOLAR PANELS BEAMS                                                                     
         S_SIDE(2,1) = 0.32D0
         S_SIDE(2,2) = 0.85D0
         S_SIDE(2,3) = 0.85D0
         S_SIDE(2,4) = 1D0

C     -------------
C     GPS BLOCK IIF
C     -------------
C     SEE IIF PRESENTATION
C     OPTICAL PROPERTIES NOT KNOWN, SAME ASSUMPTIONS AS ZIEBART (2001)
      ELSEIF(BLKNUM.EQ.8)THEN
C     +/- X SIDE
         X_SIDE(1,1) = 5.72D0
         X_SIDE(1,2) = 0.20D0
         X_SIDE(1,3) = 0.56D0
         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 7.01D0
         Y_SIDE(1,2) = 0.20D0
         Y_SIDE(1,3) = 0.56D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 5.40D0
         Z_SIDE(1,2) = 0D0
         Z_SIDE(1,3) = 0D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 5.40D0
         Z_SIDE(3,2) = 0D0
         Z_SIDE(3,3) = 0D0
         Z_SIDE(3,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 22.25D0
         S_SIDE(1,2) = 0.85D0
         S_SIDE(1,3) = 0.23D0
         S_SIDE(1,4) = 1D0

C Added from acc_albedo_propboxw.f
C     -------------
C     GPS BLOCK IIIA
C     -------------
C     DIMENSIONS ACCORDING TO PETER STEIGENBERGER (IGS-ACS-1219)
C     (mass = 2161 kg / wide = 2.46 m deep = 1.78 m high = 3.40 m)
C     SOLAR PANEL CAN BE TAKEN FROM (span = 5.33 m)
C     http://www.deagel.com/Space-Systems/GPS-Block-III_a000238005.aspx
C     BUT SO FAR,SOLAR PANEL BLOCK IIF VALUES ARE USED
C     OPTICAL PROPERTIES NOT KNOWN, COPIED FROM IIF
      ELSEIF(BLKNUM.EQ.9)THEN

  
C      +/- X SIDE
          X_SIDE(1,1) = 6.05D0
          X_SIDE(1,2) = 0.20D0
          X_SIDE(1,3) = 0.56D0
          X_SIDE(1,4) = 1D0

C      +/- Y SIDE
          Y_SIDE(1,1) = 8.36D0
          Y_SIDE(1,2) = 0.20D0
          Y_SIDE(1,3) = 0.56D0
          Y_SIDE(1,4) = 1D0

C      -Z SIDE
          Z_SIDE(1,1) = 4.38D0
          Z_SIDE(1,2) = 0D0
          Z_SIDE(1,3) = 0D0
          Z_SIDE(1,4) = 1D0

C      +Z SIDE
          Z_SIDE(3,1) = 4.38D0
          Z_SIDE(3,2) = 0D0
          Z_SIDE(3,3) = 0D0
          Z_SIDE(3,4) = 1D0

C      SOLAR PANELS
          S_SIDE(1,1) = 22.25D0
          S_SIDE(1,2) = 0.85D0
          S_SIDE(1,3) = 0.23D0
          S_SIDE(1,4) = 1D0

C     -------
C     GLONASS
C     -------
C     SEE ZIEBART (2001)
C     OPTICAL PROPERTIES NOT KNOWN, SAME ASSUMPTIONS AS ZIEBART (2001)
C     FROM PHD DISSERTATION OF RODRIGUEZ-SOLANO
      ELSEIF(BLKNUM .EQ. 10)THEN

C      +/- X SIDE (FLAT) TAH VALUES: +X   1.265   0.200   0.560  1
          X_SIDE(1,1) = 1.258D0
          X_SIDE(1,2) = 0.20D0
          X_SIDE(1,3) = 0.56D0
          X_SIDE(1,4) = 1D0

C      +/- X SIDE (CYLINDRICAL) TAH VALUES: +X   2.065   0.200   0.560  2
          X_SIDE(2,1) = 2.052D0
          X_SIDE(2,2) = 0.20D0
          X_SIDE(2,3) = 0.56D0
          X_SIDE(2,4) = 2D0

C      +/- Y SIDE (FLAT)  TAH VALUES:  +Y   2.592   0.200   0.560  1
          Y_SIDE(1,1) = 2.591D0
          Y_SIDE(1,2) = 0.20D0
          Y_SIDE(1,3) = 0.56D0
          Y_SIDE(1,4) = 1D0

C      +/- Y SIDE (CYLINDRICAL)  TAH VALUES:  +Y   2.531   0.200   0.560  2 
          Y_SIDE(2,1) = 2.532D0
          Y_SIDE(2,2) = 0.20D0
          Y_SIDE(2,3) = 0.56D0
          Y_SIDE(2,4) = 2D0

C      -Z SIDE (BUS)   TAH VALUES:  -Z   1.662   0.200   0.295  1
          Z_SIDE(1,1) = 0.877D0
          Z_SIDE(1,2) = 0.20D0
          Z_SIDE(1,3) = 0.56D0
          Z_SIDE(1,4) = 1D0

C      -Z SIDE (APOGEE ENGINE)
          Z_SIDE(2,1) = 0.785D0
          Z_SIDE(2,2) = 0D0
          Z_SIDE(2,3) = 0D0
          Z_SIDE(2,4) = 1D0

C      +Z SIDE (BUS)   TAH VALUES: +Z   1.662   0.391   0.626  1
          Z_SIDE(3,1) = 1.412D0
          Z_SIDE(3,2) = 0.20D0
          Z_SIDE(3,3) = 0.56D0
          Z_SIDE(3,4) = 1D0

C      +Z SIDE (RETRO REFLECTOR ARRAY)
          Z_SIDE(4,1) = 0.250D0
          Z_SIDE(4,2) = 1D0
          Z_SIDE(4,3) = 1D0
          Z_SIDE(4,4) = 1D0

C      SOLAR PANELS  TAH VALUES: SP  23.616   0.848   0.230  1
          S_SIDE(1,1) = 23.616D0
          S_SIDE(1,2) = 0.85D0
          S_SIDE(1,3) = 0.23D0
          S_SIDE(1,4) = 1D0

C     -------------------------------------
C     GLONASS-M
C     -------------------------------------
C     FROM PHD DISSERTATION OF RODRIGUEZ-SOLANO
      ELSEIF(BLKNUM .EQ. 11)THEN

C      +/- X SIDE (FLAT)
          X_SIDE(1,1) = 1.232D0
          X_SIDE(1,2) = 0.20D0
          X_SIDE(1,3) = 0.56D0
          X_SIDE(1,4) = 1D0

C      +/- X SIDE (CYLINDRICAL)
          X_SIDE(2,1) = 3.298D0
          X_SIDE(2,2) = 0.20D0
          X_SIDE(2,3) = 0.56D0
          X_SIDE(2,4) = 2D0

C      +/- Y SIDE (FLAT)
          Y_SIDE(1,1) = 2.700D0
          Y_SIDE(1,2) = 0.20D0
          Y_SIDE(1,3) = 0.56D0
          Y_SIDE(1,4) = 1D0

C      +/- Y SIDE (CYLINDRICAL)
          Y_SIDE(2,1) = 3.330D0
          Y_SIDE(2,2) = 0.20D0
          Y_SIDE(2,3) = 0.56D0
          Y_SIDE(2,4) = 2D0

C      -Z SIDE (BUS)
          Z_SIDE(1,1) = 2.120D0
          Z_SIDE(1,2) = 0.20D0
          Z_SIDE(1,3) = 0.21D0
          Z_SIDE(1,4) = 1D0

C      -Z SIDE (APOGEE ENGINE)
          Z_SIDE(2,1) = 0.000D0
          Z_SIDE(2,2) = 0D0
          Z_SIDE(2,3) = 0D0
          Z_SIDE(2,4) = 1D0

C      +Z SIDE (BUS)
          Z_SIDE(3,1) = 2.12D0
          Z_SIDE(3,2) = 0.30D0
          Z_SIDE(3,3) = 0.59D0
          Z_SIDE(3,4) = 1D0

C      +Z SIDE (RETRO REFLECTOR ARRAY)
          Z_SIDE(4,1) = 0.00D0
          Z_SIDE(4,2) = 1D0
          Z_SIDE(4,3) = 1D0
          Z_SIDE(4,4) = 1D0

C      SOLAR PANELS
          S_SIDE(1,1) = 30.850d0 
          S_SIDE(1,2) = 0.85D0
          S_SIDE(1,3) = 0.23D0
          S_SIDE(1,4) = 1D0

C     -------------------------------------
C     GLONASS-K1
C     -------------------------------------
C     FROM PHD DISSERTATION OF RODRIGUEZ-SOLANO
      ELSEIF (BLKNUM .EQ. 12) THEN

C      +/- X SIDE
          X_SIDE(1,1) = 2.12D0
          X_SIDE(1,2) = 0.20D0
          X_SIDE(1,3) = 0.56D0
          X_SIDE(1,4) = 1D0

C      +/- Y SIDE
          Y_SIDE(1,1) = 4.35D0
          Y_SIDE(1,2) = 0.20D0
          Y_SIDE(1,3) = 0.56D0
          Y_SIDE(1,4) = 1D0

C      -Z SIDE
          Z_SIDE(1,1) = 1.73D0
          Z_SIDE(1,2) = 0.20D0
          Z_SIDE(1,3) = 0.46D0
          Z_SIDE(1,4) = 1D0

C      +Z SIDE
          Z_SIDE(3,1) = 1.73D0
          Z_SIDE(3,2) = 0.32D0
          Z_SIDE(3,3) = 0.60D0
          Z_SIDE(3,4) = 1D0

C      SOLAR PANELS
          S_SIDE(1,1) = 16.960D0
          S_SIDE(1,2) = 0.85D0
          S_SIDE(1,3) = 0.23D0
          S_SIDE(1,4) = 1D0

C     -------------------------------------
C     GALILEO-0A, DIFFERENT ORIENTATION FOR X/Y AS IGS
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 13) THEN

C     -------------------------------------
C     GALILEO-0B, DIFFERENT ORIENTATION FOR X/Y AS IGS
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 14) THEN

C     -------------------------------------
C     GALILEO-1, DIFFERENT ORIENTATION FOR X/Y AS IGS
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 15) THEN

C     - X SIDE (IGS +X)
         X_SIDE(1,1) = 1.320D0
         X_SIDE(1,2) = 0.00D0
         X_SIDE(1,3) = 0.06D0
         X_SIDE(1,4) = 1D0

C     OMITTED +X SIDE (IGS -X)
C     IMPACT ON POD IN YAW MANEUVERS

C     - Y SIDE (IGS +Y) (Material 1)
         Y_SIDE(1,1) = 1.03D0
         Y_SIDE(1,2) = 0.00D0
         Y_SIDE(1,3) = 0.06D0
         Y_SIDE(1,4) = 1D0

C     - Y SIDE (IGS +Y) (Material 2)
         Y_SIDE(2,1) = 1.97D0
         Y_SIDE(2,2) = 0.80D0
         Y_SIDE(2,3) = 0.90D0
         Y_SIDE(2,4) = 1D0

C     -Z SIDE (Material 1)
         Z_SIDE(1,1) = 3.00D0
         Z_SIDE(1,2) = 0.00D0
         Z_SIDE(1,3) = 0.06D0
         Z_SIDE(1,4) = 1D0


C     +Z SIDE (Material 1)
         Z_SIDE(3,1) = 1.72D0
         Z_SIDE(3,2) = 0.00D0
         Z_SIDE(3,3) = 0.06D0
         Z_SIDE(3,4) = 1D0

C     +Z SIDE (Material 2)
         Z_SIDE(4,1) = 1.28D0
         Z_SIDE(4,2) = 0.51D0
         Z_SIDE(4,3) = 0.43D0
         Z_SIDE(4,4) = 1D0

C     SOLAR PANELS (Sollar Cells)
         S_SIDE(1,1) = 7.76D0
         S_SIDE(1,2) = 1.00D0
         S_SIDE(1,3) = 0.08D0
         S_SIDE(1,4) = 1D0

C     SOLAR PANELS (Kapton HN - Insulation layer)
         S_SIDE(2,1) = 3.06D0
         S_SIDE(2,2) = 1.00D0
         S_SIDE(2,3) = 0.10D0
         S_SIDE(2,4) = 1D0

C     -------------------------------------
C     GALILEO-2
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 16) THEN

C     - X SIDE (IGS +X)
         X_SIDE(1,1) = 1.320D0
         X_SIDE(1,2) = 0.00D0
         X_SIDE(1,3) = 0.07D0
         X_SIDE(1,4) = 1D0

C     - Y SIDE (IGS +Y) (Material 1)
         Y_SIDE(1,1) = 1.244D0
         Y_SIDE(1,2) = 0.00D0
         Y_SIDE(1,3) = 0.07D0
         Y_SIDE(1,4) = 1D0

C     - Y SIDE (IGS +Y) (Material 2)
         Y_SIDE(2,1) = 1.539D0
         Y_SIDE(2,2) = 0.793D0
         Y_SIDE(2,3) = 0.92D0
         Y_SIDE(2,4) = 1D0

C     -Z SIDE (Material 1)
         Z_SIDE(1,1) = 2.077D0
         Z_SIDE(1,2) = 0.00D0
         Z_SIDE(1,3) = 0.07D0
         Z_SIDE(1,4) = 1D0

C     -Z SIDE (Material 2)
         Z_SIDE(2,1) = 0.959D0
         Z_SIDE(2,2) = 0.793D0
         Z_SIDE(2,3) = 0.92D0
         Z_SIDE(2,4) = 1D0

C     +Z SIDE (Material 1)
         Z_SIDE(3,1) = 1.053D0
         Z_SIDE(3,2) = 0.00D0
         Z_SIDE(3,3) = 0.07D0
         Z_SIDE(3,4) = 1D0

C     +Z SIDE (Material 2)
         Z_SIDE(4,1) = 1.969D0
         Z_SIDE(4,2) = 0.51D0
         Z_SIDE(4,3) = 0.43D0
         Z_SIDE(4,4) = 1D0

C     SOLAR PANELS (Sollar Cells)
         S_SIDE(1,1) = 7.76D0
         S_SIDE(1,2) = 1.00D0
         S_SIDE(1,3) = 0.08D0
         S_SIDE(1,4) = 1D0

C     SOLAR PANELS (Kapton HN - Insulation layer)
         S_SIDE(2,1) = 3.06D0
         S_SIDE(2,2) = 1.00D0
         S_SIDE(2,3) = 0.10D0
         S_SIDE(2,4) = 1D0

C     -------------------------------------
C     BEIDOU-2G
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 17) THEN

C     +/- X SIDE
         X_SIDE(1,1) = 3.784D0
         X_SIDE(1,2) = 1.00D0
         X_SIDE(1,3) = 0.65D0
         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 4.4D0
         Y_SIDE(1,2) = 1.00D0
         Y_SIDE(1,3) = 0.856D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 3.44D0
         Z_SIDE(1,2) = 1.0D0
         Z_SIDE(1,3) = 0.65D0
         Z_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(2,1) = 3.00D0
         Z_SIDE(2,2) = 1.0D0
         Z_SIDE(2,3) = 0.65D0
         Z_SIDE(2,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 3.44D0
         Z_SIDE(3,2) = 1.0D0
         Z_SIDE(3,3) = 0.65D0
         Z_SIDE(3,4) = 1D0

C     +Z SIDE
         Z_SIDE(4,1) = 3.00D0
         Z_SIDE(4,2) = 1.0D0
         Z_SIDE(4,3) = 0.65D0
         Z_SIDE(4,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 22.704D0
         S_SIDE(1,2) = 1.0D0
         S_SIDE(1,3) = 0.28D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     BEIDOU-2I, BEIDOU-2M
C     -------------------------------------
      ELSEIF (BLKNUM.EQ.18 .OR. BLKNUM.EQ.19) THEN

C     +/- X SIDE
         X_SIDE(1,1) = 3.784D0
         X_SIDE(1,2) = 1.00D0
         X_SIDE(1,3) = 0.65D0
         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 4.4D0
         Y_SIDE(1,2) = 1.00D0
         Y_SIDE(1,3) = 0.856D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 3.44D0
         Z_SIDE(1,2) = 1.0D0
         Z_SIDE(1,3) = 0.65D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 3.44D0
         Z_SIDE(3,2) = 1.0D0
         Z_SIDE(3,3) = 0.65D0
         Z_SIDE(3,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 22.704D0
         S_SIDE(1,2) = 1.0D0
         S_SIDE(1,3) = 0.28D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     BEIDOU-3IS-SECM
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 20) THEN

C     -------------------------------------
C     BEIDOU-3IS-CAST
C     -------------------------------------	  
      ELSEIF (BLKNUM .EQ. 21) THEN

C     +/- X SIDE
         X_SIDE(1,1) = 7.56D0
         X_SIDE(1,2) = 1.00D0
         X_SIDE(1,3) = 0.65D0
         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 9.00D0
         Y_SIDE(1,2) = 1.00D0
         Y_SIDE(1,3) = 0.856D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 5.25D0
         Z_SIDE(1,2) = 1.00D0
         Z_SIDE(1,3) = 0.08D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 5.25D0
         Z_SIDE(3,2) = 1.0D0
         Z_SIDE(3,3) = 0.65D0
         Z_SIDE(3,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 40.56D0
         S_SIDE(1,2) = 1.0D0
         S_SIDE(1,3) = 0.08D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     BEIDOU-3MS-SECM
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 22) THEN


C     -------------------------------------
C     BEIDOU-3MS-CAST
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 23) THEN

C     +/- X SIDE
         X_SIDE(1,1) = 2.88D0
         X_SIDE(1,2) = 1.00D0
         X_SIDE(1,3) = 0.65D0
         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 4.32D0
         Y_SIDE(1,2) = 1.00D0
         Y_SIDE(1,3) = 0.856D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 2.30D0
         Z_SIDE(1,2) = 1.00D0
         Z_SIDE(1,3) = 0.65D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 2.30D0
         Z_SIDE(3,2) = 1.0D0
         Z_SIDE(3,3) = 0.65D0
         Z_SIDE(3,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 20.43D0
         S_SIDE(1,2) = 1.0D0
         S_SIDE(1,3) = 0.08D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     BEIDOU-3G-SECM
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 24) THEN

C     -------------------------------------
C     BEIDOU-3G-CAST
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 25) THEN

C     +/- X SIDE
         X_SIDE(1,1) = 8.4960D0
         X_SIDE(1,2) = 1.00D0
         X_SIDE(1,3) = 0.65D0
         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 7.56D0
         Y_SIDE(1,2) = 1.00D0
         Y_SIDE(1,3) = 0.865D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 4.956D0
         Z_SIDE(1,2) = 1.0D0
         Z_SIDE(1,3) = 0.13D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 4.956D0
         Z_SIDE(3,2) = 1.0D0
         Z_SIDE(3,3) = 0.13D0
         Z_SIDE(3,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 35.4D0
         S_SIDE(1,2) = 1.00D0
         S_SIDE(1,3) = 0.08D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     BEIDOU-3I-SECM
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 26) THEN

C     -------------------------------------
C     BEIDOU-3I-CAST
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 27) THEN

C     +/- X SIDE
         X_SIDE(1,1) = 8.4960D0
         X_SIDE(1,2) = 1.00D0
         X_SIDE(1,3) = 0.65D0
         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 7.56D0
         Y_SIDE(1,2) = 1.00D0
         Y_SIDE(1,3) = 0.865D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 4.956D0
         Z_SIDE(1,2) = 1.0D0
         Z_SIDE(1,3) = 0.13D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 4.956D0
         Z_SIDE(3,2) = 1.0D0
         Z_SIDE(3,3) = 0.13D0
         Z_SIDE(3,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 35.4D0
         S_SIDE(1,2) = 1.00D0
         S_SIDE(1,3) = 0.08D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     BEIDOU-3M-SECM
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 28) THEN

C     +/- X SIDE
         X_SIDE(1,1) = 1.250D0
         X_SIDE(1,2) = 1.00D0
         X_SIDE(1,3) = 0.80D0
         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 3.13D0
         Y_SIDE(1,2) = 1.00D0
         Y_SIDE(1,3) = 0.80D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 2.59D0
         Z_SIDE(1,2) = 1.0D0
         Z_SIDE(1,3) = 0.80D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 2.59D0
         Z_SIDE(3,2) = 1.0D0
         Z_SIDE(3,3) = 0.80D0
         Z_SIDE(3,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 10.8D0
         S_SIDE(1,2) = 1.00D0
         S_SIDE(1,3) = 0.08D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     BEIDOU-3M-CAST
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 29) THEN

C     + X SIDE
         X_SIDE(1,1) = 2.860D0
         X_SIDE(1,2) = 1.00D0
         X_SIDE(1,3) = 0.65D0
         X_SIDE(1,4) = 1D0

C     - X SIDE
C         X_SIDE(1,1) = 1.750D0
C         X_SIDE(1,2) = 1.00D0
C         X_SIDE(1,3) = 0.08D0
C         X_SIDE(1,4) = 1D0

C     - X SIDE
C         X_SIDE(1,1) = 1.110D0
C         X_SIDE(1,2) = 1.00D0
C         X_SIDE(1,3) = 0.865D0
C         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 3.63D0
         Y_SIDE(1,2) = 1.00D0
         Y_SIDE(1,3) = 0.865D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 2.18D0
         Z_SIDE(1,2) = 1.0D0
         Z_SIDE(1,3) = 0.65D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 2.18D0
         Z_SIDE(3,2) = 1.0D0
         Z_SIDE(3,3) = 0.08D0
         Z_SIDE(3,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 20.44D0
         S_SIDE(1,2) = 1.0D0
         S_SIDE(1,3) = 0.08D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     QZSS Michibiki (Montenbruck, QZSS SRP)
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 30) THEN

C     +/- X SIDE (black)
         X_SIDE(1,1) = 9.90D0
         X_SIDE(1,2) = 0.00D0
         X_SIDE(1,3) = 0.06D0
         X_SIDE(1,4) = 1D0

C     +/- X SIDE (sliver)
         X_SIDE(2,1) = 2.30D0
         X_SIDE(2,2) = 0.178D0
         X_SIDE(2,3) = 0.56D0
         X_SIDE(2,4) = 1D0

C     +Y SIDE (black)
         Y_SIDE(1,1) = 4.6D0
         Y_SIDE(1,2) = 0.00D0
         Y_SIDE(1,3) = 0.06D0
         Y_SIDE(1,4) = 1D0

C     +Y SIDE (radiator)
         Y_SIDE(2,1) = 5.3D0
         Y_SIDE(2,2) = 1.00D0
         Y_SIDE(2,3) = 0.94D0
         Y_SIDE(2,4) = 1D0

C THIS SHOULD BE ADDED TO OTHER CONTRIBUTION
C OMITTED THIS WILL CAUSED ERRORS !!!!!!!!!!
C     +Y SIDE (sliver)
C         Y_SIDE(3,1) = 2.7D0
C         Y_SIDE(3,2) = 0.178D0
C         Y_SIDE(3,3) = 0.56D0
C         Y_SIDE(3,4) = 1D0

C     -Z SIDE (black)
         Z_SIDE(1,1) = 6.00D0
         Z_SIDE(1,2) = 0.00D0
         Z_SIDE(1,3) = 0.06D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE (black)
         Z_SIDE(3,1) = 2.00D0
         Z_SIDE(3,2) = 0.00D0
         Z_SIDE(3,3) = 0.06D0
         Z_SIDE(3,4) = 1D0

C     +Z SIDE (sliver)
         Z_SIDE(4,1) = 4.00D0
         Z_SIDE(4,2) = 0.178D0
         Z_SIDE(4,3) = 0.56D0
         Z_SIDE(4,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 45.408D0
         S_SIDE(1,2) = 1.0D0
         S_SIDE(1,3) = 0.28D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     QZSS-2I
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 31) THEN

C     +/- X SIDE (MLI)
         X_SIDE(1,1) = 10.10D0
         X_SIDE(1,2) =  0.473D0
         X_SIDE(1,3) =  0.074D0
         X_SIDE(1,4) =  1D0

C     +Y SIDE (MLI)
         Y_SIDE(1,1) = 5.7D0
         Y_SIDE(1,2) = 0.473D0
         Y_SIDE(1,3) = 0.074D0
         Y_SIDE(1,4) = 1D0

C     +Y SIDE (radiator)
         Y_SIDE(2,1) = 4.4D0
         Y_SIDE(2,2) = 0.979D0
         Y_SIDE(2,3) = 0.974D0
         Y_SIDE(2,4) = 1D0

C     -Z SIDE (MLI)
         Z_SIDE(1,1) = 5.60D0
         Z_SIDE(1,2) = 0.473D0
         Z_SIDE(1,3) = 0.074D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE (MLI)
         Z_SIDE(3,1) = 2.70D0
         Z_SIDE(3,2) = 0.473D0
         Z_SIDE(3,3) = 0.074D0
         Z_SIDE(3,4) = 1D0

C     +Z SIDE (L-ANT Cover+L1S/L5S-ANT Covers)
         Z_SIDE(4,1) = 5.83D0
         Z_SIDE(4,2) = 0.224D0
         Z_SIDE(4,3) = 0.447D0
         Z_SIDE(4,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 29.80D0
         S_SIDE(1,2) = 0.883D0
         S_SIDE(1,3) = 0.077D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     QZSS-2G
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 32) THEN

C     +/- X SIDE (MLI)
         X_SIDE(1,1) = 10.10D0
         X_SIDE(1,2) =  0.473D0
         X_SIDE(1,3) =  0.074D0
         X_SIDE(1,4) =  1D0

C     +Y SIDE (MLI)
         Y_SIDE(1,1) = 4.8D0
         Y_SIDE(1,2) = 0.473D0
         Y_SIDE(1,3) = 0.074D0
         Y_SIDE(1,4) = 1D0

C     +Y SIDE (radiator)
         Y_SIDE(2,1) = 5.3D0
         Y_SIDE(2,2) = 0.979D0
         Y_SIDE(2,3) = 0.974D0
         Y_SIDE(2,4) = 1D0

C     -Z SIDE (MLI)
         Z_SIDE(1,1) = 5.60D0
         Z_SIDE(1,2) = 0.473D0
         Z_SIDE(1,3) = 0.074D0
         Z_SIDE(1,4) = 1D0

C     -Z SIDE (Reflector)
         Z_SIDE(2,1) = 9.10D0
         Z_SIDE(2,2) = 0.013D0
         Z_SIDE(2,3) = 0.381D0
         Z_SIDE(2,4) = 1D0

C +Z SHOULD CONTAIN OTHERS MLI, L-ANT Cover, AND L1S/Sb/L5S-ANT Covers
C THE FOLLOWING IS AVERAGE VALUE

C     +Z SIDE (MLI)
         Z_SIDE(3,1) = 5.50D0
         Z_SIDE(3,2) = 0.747D0 ! 0.2016
         Z_SIDE(3,3) = 0.270D0 ! 0.7304
         Z_SIDE(3,4) = 1D0

C     +Z SIDE (Reflector)
         Z_SIDE(4,1) = 9.1D0
         Z_SIDE(4,2) = 0.013D0
         Z_SIDE(4,3) = 0.381D0
         Z_SIDE(4,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 29.80D0
         S_SIDE(1,2) = 0.883D0
         S_SIDE(1,3) = 0.077D0
         S_SIDE(1,4) = 1D0

C     -------------------------------------
C     QZSS-2A
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 33) THEN

C     +/- X SIDE (MLI)
         X_SIDE(1,1) = 10.80D0
         X_SIDE(1,2) =  0.473D0
         X_SIDE(1,3) =  0.074D0
         X_SIDE(1,4) =  1D0

C     +Y SIDE (MLI)
         Y_SIDE(1,1) = 5.1D0
         Y_SIDE(1,2) = 0.473D0
         Y_SIDE(1,3) = 0.074D0
         Y_SIDE(1,4) = 1D0

C     +Y SIDE (radiator)
         Y_SIDE(2,1) = 4.9D0
         Y_SIDE(2,2) = 0.979D0
         Y_SIDE(2,3) = 0.974D0
         Y_SIDE(2,4) = 1D0

C     -Z SIDE (MLI)
         Z_SIDE(1,1) = 5.60D0
         Z_SIDE(1,2) = 0.473D0
         Z_SIDE(1,3) = 0.074D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE (MLI)
         Z_SIDE(3,1) = 3.40D0
         Z_SIDE(3,2) = 0.473D0
         Z_SIDE(3,3) = 0.074D0
         Z_SIDE(3,4) = 1D0

C     +Z SIDE (L-ANT Cover+L1S/L5S-ANT Covers)
         Z_SIDE(4,1) = 2.60D0
         Z_SIDE(4,2) = 0.823D0
         Z_SIDE(4,3) = 0.507D0
         Z_SIDE(4,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 29.80D0
         S_SIDE(1,2) = 0.883D0
         S_SIDE(1,3) = 0.077D0
         S_SIDE(1,4) = 1D0


C     -------------------------------------
C     IRNSS-1IGSO
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 34) THEN

C     -------------------------------------
C     IRNSS-1IGSO
C     -------------------------------------
      ELSEIF (BLKNUM .EQ. 35) THEN

      ENDIF


C ---------------------------
C BOX-WING MODEL APROXIMATION
C ---------------------------

C NAVIGATION ANTENNAS
C FOR BLOCK I/II/IIA
      X_SIDE(5,1) = X_SIDE(5,1)*12

C MATRIX WITH ALL SURFACES
      DO II = 1,14
         DO JJ = 1,4
            IF(II.LE.5)THEN
               SURFALL(II,JJ) = X_SIDE(II,JJ)
            ELSEIF((II.GE.6).AND.(II.LE.9))THEN
               SURFALL(II,JJ) = Z_SIDE(II-5,JJ)
            ELSEIF((II.GE.10).AND.(II.LE.11))THEN
               SURFALL(II,JJ) = S_SIDE(II-9,JJ)
            ELSEIF(II.GE.12)THEN
               SURFALL(II,JJ) = Y_SIDE(II-11,JJ)
            ENDIF
         ENDDO
      ENDDO
      

C COMPUTATION OF FRACTION OF REFLECTED, DIFFUSSED AND ABSORVED PHOTONS
C FROM SPECULARITY AND REFLECTIVITY
      DO II = 1,14
         ALLREFL(II) = SURFALL(II,2)*SURFALL(II,3)
         ALLDIFU(II) = SURFALL(II,3)*(1-SURFALL(II,2))
      ENDDO

C MULTIPLICATION OF OPTICAL PROPERTIES WITH SURFACE AREA
      DO II = 1,14
         REFL_AREA(II) = ALLREFL(II)*SURFALL(II,1)
         DIFU_AREA(II) = ALLDIFU(II)*SURFALL(II,1)
      ENDDO
      
C AVERAGE PER SURFACE (SEPARATE FOR FLAT AND CYLINDRICAL SURFACES)
      DO SS = 1,5
         G_AREA1(SS) = 0D0
         G_REFL1(SS) = 0D0
         G_DIFU1(SS) = 0D0
         G_ABSP1(SS) = 0D0
         G_AREA2(SS) = 0D0
         G_REFL2(SS) = 0D0
         G_DIFU2(SS) = 0D0
         G_ABSP2(SS) = 0D0
      ENDDO


      DO II = 1,14
C        +X SIDE 
         IF(II.LE.5)THEN
            SS = 1
C        -Z SIDE
         ELSEIF((II.GE.6).AND.(II.LE.7))THEN
            SS = 2
C        +Z SIDE
         ELSEIF((II.GE.8).AND.(II.LE.9))THEN
            SS = 3
C        SOLAR PANELS
         ELSEIF((II.GE.10).AND.(II.LE.11))THEN
            SS = 4
C        +/- Y SIDE
         ELSEIF(II.GE.12)THEN
            SS = 5
         ENDIF

         IF(SURFALL(II,4).EQ.1D0)THEN
            G_AREA1(SS) = G_AREA1(SS) + SURFALL(II,1)
            G_REFL1(SS) = G_REFL1(SS) + REFL_AREA(II)
            G_DIFU1(SS) = G_DIFU1(SS) + DIFU_AREA(II)
         ELSEIF(SURFALL(II,4).EQ.2D0)THEN
            G_AREA2(SS) = G_AREA2(SS) + SURFALL(II,1)
            G_REFL2(SS) = G_REFL2(SS) + REFL_AREA(II)
            G_DIFU2(SS) = G_DIFU2(SS) + DIFU_AREA(II)
         ENDIF
      ENDDO

C AVERAGE OF OPTICAL PROPERTIES ACCORDING TO SURFACES AREA
      DO SS = 1,5
         IF(G_AREA1(SS).GT.0D0)THEN
            G_REFL1(SS) = G_REFL1(SS)/G_AREA1(SS)
            G_DIFU1(SS) = G_DIFU1(SS)/G_AREA1(SS)
            G_ABSP1(SS) = 1D0 - (G_REFL1(SS) + G_DIFU1(SS))
         ENDIF
         IF(G_AREA2(SS).GT.0D0)THEN
            G_REFL2(SS) = G_REFL2(SS)/G_AREA2(SS)
            G_DIFU2(SS) = G_DIFU2(SS)/G_AREA2(SS)
            G_ABSP2(SS) = 1D0 -(G_REFL2(SS) + G_DIFU2(SS))
         ENDIF
      ENDDO


C --------------------------------------------
C ARRAYS OF OPTICAL PROPERTIES AND DIMENSIONS
C --------------------------------------------

      DO SS = 1,4
         DO KK = 1,2
            AREA(SS,KK) = 0D0
            REFL(SS,KK) = 0D0
            DIFU(SS,KK) = 0D0
            ABSP(SS,KK) = 0D0
            AREA2(SS,KK) = 0D0
            REFL2(SS,KK) = 0D0
            DIFU2(SS,KK) = 0D0
            ABSP2(SS,KK) = 0D0
         ENDDO
      ENDDO

C     +Z FACE      
      AREA(1,1) = G_AREA1(3)
      REFL(1,1) = G_REFL1(3)
      DIFU(1,1) = G_DIFU1(3)
      ABSP(1,1) = G_ABSP1(3)

C     -Z FACE
      AREA(1,2) = G_AREA1(2)
      REFL(1,2) = G_REFL1(2)
      DIFU(1,2) = G_DIFU1(2)
      ABSP(1,2) = G_ABSP1(2)

C     +X FACE
      AREA(3,1) = G_AREA1(1)
      REFL(3,1) = G_REFL1(1)
      DIFU(3,1) = G_DIFU1(1)
      ABSP(3,1) = G_ABSP1(1)

C     +X CYLINDRICAL
      AREA2(3,1) = G_AREA2(1)
      REFL2(3,1) = G_REFL2(1)
      DIFU2(3,1) = G_DIFU2(1)
      ABSP2(3,1) = G_ABSP2(1)

C     -X FACE (ASSUMED FOR ALL BLOCKS)
      AREA(3,2) = G_AREA1(1)
      REFL(3,2) = G_REFL1(1)
      DIFU(3,2) = G_DIFU1(1)
      ABSP(3,2) = G_ABSP1(1)

C     -X CYLINDRICAL (ASSUMED FOR ALL BLOCKS)                                               
      AREA2(3,2) = G_AREA2(1)
      REFL2(3,2) = G_REFL2(1)
      DIFU2(3,2) = G_DIFU2(1)
      ABSP2(3,2) = G_ABSP2(1)

C     +Y FACE (ASSUMED FOR BLKNUM = 1...7)
      IF((BLKNUM.GE.1).AND.(BLKNUM.LE.7))THEN
C        DUE TO BLOCK IIR W-SENSOR
         IF((BLKNUM.GE.4).AND.(BLKNUM.LE.7))THEN      
            AREA(2,1) = ((AREA(1,1)-0.5D0) + AREA(3,1))/(2D0)
         ELSE
            AREA(2,1) = (AREA(1,1) + AREA(3,1))/(2D0)
         ENDIF
         REFL(2,1) = (REFL(1,1) + REFL(1,2) + REFL(3,1) + REFL(3,2))/4
         DIFU(2,1) = (DIFU(1,1) + DIFU(1,2) + DIFU(3,1) + DIFU(3,2))/4
         ABSP(2,1) = (ABSP(1,1) + ABSP(1,2) + ABSP(3,1) + ABSP(3,2))/4

C        CYLINDRICAL PART EQUAL TO +/- X CYLINDRICAL
         AREA2(2,1) = AREA2(3,1)
         REFL2(2,1) = REFL2(3,1)
         DIFU2(2,1) = DIFU2(3,1)
         ABSP2(2,1) = ABSP2(3,1)

C     FOR OTHER SATELLITE BLOCKS
      ELSE
         AREA(2,1) = G_AREA1(5)
         REFL(2,1) = G_REFL1(5)
         DIFU(2,1) = G_DIFU1(5)
         ABSP(2,1) = G_ABSP1(5)

         AREA2(2,1) = G_AREA2(5)
         REFL2(2,1) = G_REFL2(5)
         DIFU2(2,1) = G_DIFU2(5)
         ABSP2(2,1) = G_ABSP2(5)
      ENDIF

C     -Y FACE (ASSUMED FOR ALL BLOCKS)
      AREA(2,2) = AREA(2,1)
      REFL(2,2) = REFL(2,1)
      DIFU(2,2) = DIFU(2,1)
      ABSP(2,2) = ABSP(2,1)

C     -Y CYLINDRICAL (ASSUMED FOR ALL BLOCKS)
      AREA2(2,2) = AREA2(2,1)
      REFL2(2,2) = REFL2(2,1)
      DIFU2(2,2) = DIFU2(2,1)
      ABSP2(2,2) = ABSP2(2,1)

C     FRONT OF SOLAR PANELS
      AREA(4,1) = G_AREA1(4)
      REFL(4,1) = G_REFL1(4)
      DIFU(4,1) = G_DIFU1(4)
      ABSP(4,1) = G_ABSP1(4)

C     FRONT OF SOLAR PANELS (CYLINDER)
      AREA2(4,1) = G_AREA2(4)
      REFL2(4,1) = G_REFL2(4)
      DIFU2(4,1) = G_DIFU2(4)
      ABSP2(4,1) = G_ABSP2(4)

C     BACK OF SOLAR PANELS (ASSUMED FOR ALL BLOCKS)
      AREA(4,2) = AREA(4,1)
      REFL(4,2) = 0.055D0
      DIFU(4,2) = 0.055D0
      ABSP(4,2) = 0.890D0


C     BACK OF SOLAR PANELS (CYLINDER)
      AREA2(4,2) = G_AREA2(4)
      REFL2(4,2) = G_REFL2(4)
      DIFU2(4,2) = G_DIFU2(4)
      ABSP2(4,2) = G_ABSP2(4)


C     OPTICAL PROPERTIES IN THE INFRARED
C     ASSUMED FOR ALL BLOCKS AND ALL SURFACES
      DO SS=1,4
         DO II=1,2
            REFLIR(SS,II) = 0.10D0
            DIFUIR(SS,II) = 0.10D0
            ABSPIR(SS,II) = 0.80D0
         ENDDO
      ENDDO

C     NOT ENOUGH INFORMATION AVAILABLE FOR GLONASS-M
C     SEE IGSMAIL-5104
      IF(BLKNUM.EQ.11)THEN
         BUSFAC = 4.2D0/3.31D0
         DO SS=1,3
            DO II=1,2
               AREA(SS,II) = BUSFAC*AREA(SS,II)
               AREA2(SS,II) = BUSFAC*AREA2(SS,II)
            ENDDO
         ENDDO
         AREA(4,1) = 30.85D0
         AREA(4,2) = 30.85D0
      ENDIF

C     NO (YET) INFORMATION AVAILABLE FOR BEIDOU,IRNSS
C      IF((BLKNUM.GE.20 .AND. BLKNUM.LE.29) .OR. BLKNUM.EQ.33 .OR. BLKNUM.EQ.34)THEN
      IF((BLKNUM.GE.20 .AND. BLKNUM.LE.27) .OR. BLKNUM.EQ.17 .OR.
     1    BLKNUM.EQ.18 .OR. (BLKNUM.GE.30 .AND. BLKNUM.LE.35)) THEN
C      IF(BLKNUM.EQ.33 .OR. BLKNUM.EQ.34)THEN
         DO SS=1,4
            DO II=1,2
               AREA(SS,II) = 0D0
               REFL(SS,II) = 0D0
               DIFU(SS,II) = 0D0
               ABSP(SS,II) = 0D0

               AREA2(SS,II) = 0D0
               REFL2(SS,II) = 0D0
               DIFU2(SS,II) = 0D0
               ABSP2(SS,II) = 0D0

               REFLIR(SS,II) = 0D0
               DIFUIR(SS,II) = 0D0
               ABSPIR(SS,II) = 0D0
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE
