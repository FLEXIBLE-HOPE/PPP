C*
      FUNCTION DJUL(J1,M1,T)
CC
CC NAME       :  DJUL
CC
CC    XX = DJUL(J1,M1,T)
CC
CC PURPOSE    :  COMPUTES THE MODIFIED JULIAN DATE (MJD)
CC               FROM YEAR,MONTH AND DAY
CC               MJD = JULIAN DATE - 2400000.5
CC
CC PARAMETERS :
CC         IN :  J1     : YEAR (E.G. 1984)                    I*4
CC               M1     : MONTH(E.G. 2 FOR FEBRUARY)          I*4
CC               T      : DAY OF MONTH                        R*8
CC        OUT :  DJUL   : MODIFIED JULIAN DATE                R*8
CC
CC SR CALLED  :  ---
CC
CC REMARKS    :  ---
CC
CC AUTHOR     :  G.BEUTLER
CC
CC VERSION    :  3.4  (JAN 93)
CC
CC CREATED    :  87/10/30 17:23        LAST MODIFIED :  88/11/21 16:42
CC
CC COPYRIGHT  :  ASTRONOMICAL INSTITUTE
CC      1987      UNIVERSITY OF BERNE
CC                    SWITZERLAND
CC
C*
      REAL*8 DJUL,T
      J=J1
      M=M1
      IF(M.GT.2)GO TO 1
      J=J-1
      M=M+12
1     I=J/100
      K=2-I+I/4
      DJUL=(365.25D0*J-DMOD(365.25D0*J,1.D0))-679006.D0
      DJUL=DJUL+AINT(30.6001*(M+1))+T+K
      RETURN
      END
