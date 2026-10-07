!*
MODULE qcflag
!*
USE par
IMPLICIT NONE

  !*
  ! Length parameters
  !!-----------------

  ! For qc
  ! Bit flag. (no flag means good data)
  ! 0 - 15 (ok), 16 -31 (del)
  INTEGER(IT), PARAMETER :: FLAG_NODATA   = 31
  INTEGER(IT), PARAMETER :: FLAG_NO4      = 30
  INTEGER(IT), PARAMETER :: FLAG_LOWELE   = 29
  INTEGER(IT), PARAMETER :: FLAG_SHORT    = 28
  INTEGER(IT), PARAMETER :: FLAG_LWBAD    = 25
  INTEGER(IT), PARAMETER :: FLAG_LGBAD    = 24
  INTEGER(IT), PARAMETER :: FLAG_LCCHECK  = 23
  INTEGER(IT), PARAMETER :: FLAG_PC1MS    = 18
  INTEGER(IT), PARAMETER :: FLAG_PCBAD    = 16
  INTEGER(IT), PARAMETER :: FLAG_LWCONN   = 6
  INTEGER(IT), PARAMETER :: FLAG_LGJUMP   = 2
  INTEGER(IT), PARAMETER :: FLAG_LWJUMP   = 3
  INTEGER(IT), PARAMETER :: FLAG_GAP      = 4
  INTEGER(IT), PARAMETER :: FLAG_LLI      = 5
  INTEGER(IT), PARAMETER :: FLAG_BIGSD    = 7

  ! For edtres
  ! Not AMB BAD DEL NODATA
  INTEGER(IT), PARAMETER :: GOOD          = 0
  ! ambiguity flag from log file
  INTEGER(IT), PARAMETER :: OLDAMB        = 1
  ! newly found ambiguity
  INTEGER(IT), PARAMETER :: NEWAMB        = 2
  ! lower elevation
  INTEGER(IT), PARAMETER :: LOWELE        = 3
  ! insufficent observations
  INTEGER(IT), PARAMETER :: NOTALL        = 4
  ! removed as short piece
  INTEGER(IT), PARAMETER :: DELSHRT       = 5
  ! removed because of large residual
  INTEGER(IT), PARAMETER :: DELBAD        = 6
  ! removed by rmstest
  INTEGER(IT), PARAMETER :: DELRMS        = 7
  ! temporerily used in check_jump
  INTEGER(IT), PARAMETER :: NEWBAD        = 8
  ! no data
  INTEGER(IT), PARAMETER :: NODATA        = 9

END MODULE qcflag

