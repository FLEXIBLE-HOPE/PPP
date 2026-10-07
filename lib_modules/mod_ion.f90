!*
MODULE ion
!!
USE const
USE par
IMPLICIT NONE

  TYPE ionhead
    ! 文件第一行 版本
    REAL(RL) :: ver=0.d0
    CHARACTER :: typ=''
    CHARACTER(LEN=3) :: sys

    ! 文件第二行 数据来源
    CHARACTER(LEN=20) :: pgm=''
    CHARACTER(LEN=20) :: run=''
    CHARACTER(LEN=20) :: dat=''

    ! 第一个和最后一个电离层图的历元
    INTEGER(IT) :: mjd0,mjd1
    REAL(RL) :: sod0,sod1

    ! 历元间隔（秒）
    INTEGER(IT) :: intv=0

    ! 电离层图的个数
    INTEGER(IT) :: nmaps=0

    ! 电离层图的投影函数
    CHARACTER(LEN_STRING) :: mapfunc=''

    ! 高度角
    REAL(RL) :: elev=0.d0

    ! 站点数和卫星个数
    INTEGER(IT) :: nsat=0
    INTEGER(IT) :: nsit=0

    ! 地球半径
    REAL(RL) :: radius=0.d0

    ! 电离层图的维数
    INTEGER(IT) :: ndim=0

    ! 始高程、终高程、高程增量
    REAL(RL) :: hgt1=0
    REAL(RL) :: hgt2=0
    REAL(RL) :: dhgt=0

    ! 始纬度、终纬度、纬度增量
    REAL(RL) :: lat1=0
    REAL(RL) :: lat2=0
    REAL(RL) :: dlat=0

    ! 始经度、终经度、经度增量
    REAL(RL) :: lon1=0
    REAL(RL) :: lon2=0
    REAL(RL) :: dlon=0

    ! 指数
    INTEGER(IT) :: unt=0

    INTEGER(IT) :: nprn=0
    REAL(RL) :: dcb(MAXSAT)=0.d0
    REAL(RL) :: rms(MAXSAT)=0.d0
    CHARACTER(LEN_PRN) :: cprn(MAXSAT)=''

  END TYPE


  TYPE iondatalat
    REAL(RL) :: lat
    REAL(RL) :: hgt

    ! data and rms
    ! 一般经度增量为5°,360/5=72，所以有73个点，其中第一个点和最后一个重合*/
    INTEGER(IT) :: tec(MAXIONLON)=0
    INTEGER(IT) :: rms(MAXIONLON)=9999  
  END TYPE

  TYPE iondata
    ! ID for GIM maps
    INTEGER(IT) :: ind=0
    ! Epoch of current GIM maps
    INTEGER(IT) :: mjd
    REAL(RL) :: sod

    ! 一般纬度增量为2.5°,87.5*2/2.5=70，所以有71个纬度数据带
    TYPE(IONDATALAT) :: ionlat(MAXIONLAT,MAXIONHGT)
  END TYPE

  !! IONEX file
  TYPE ionex
    ! IONEX head
    TYPE(IONHEAD) :: hd

    INTEGER(IT) :: nlat=0
    INTEGER(IT) :: nlon=0
    INTEGER(IT) :: nhgt=0

    ! Data record of IONDEX file, normally 13 records for daily file
    TYPE(IONDATA) :: iondat(MAXIONMAPS)
  END TYPE


END MODULE
