#=============================================================================
# PPP 的 CMake 构建 (供 Trae / CMake Tools 直接使用)
#
# 2026-10-07 by xsy: 参考 /home/xfd/ckdrt_irc_cr_orb_sdb/CMakeLists.txt 编写。
#   编译选项与该参考程序对齐; 数值与环境事实全部取自本目录原有 Makefile,
#   逐项对应关系见各段注释。原 Makefile 保留, 两条构建路线并存。
#=============================================================================

cmake_minimum_required(VERSION 3.16)

#-----------------------------------------------------------------------------
# 编译器: 与 Makefile 构建保持一致
#   依据 OBJ_INTEL/dtm_model.mod 头部 "2021.5.0" = ifort 2022.0.1 (classic)。
#   注意: 只 pin 编译器, 不改任何源码。
#-----------------------------------------------------------------------------
set(CMAKE_Fortran_COMPILER "/opt/intel/oneapi/compiler/2022.0.1/linux/bin/intel64/ifort")
set(CMAKE_C_COMPILER       "/opt/intel/oneapi/compiler/2022.0.1/linux/bin/intel64/icc")

#-----------------------------------------------------------------------------
# Intel 运行库 —— 必须在 project() 之前
#   project() 要对编译器做"能否编译简单程序"的 ABI 自检, 那一步就会调 ld;
#   而 ld 加载 icx-lto.so 插件时需要 libimf.so。Trae 继承的环境未必 source 过
#   oneAPI, 不先设这两行会在 project() 阶段直接失败:
#     ld: .../icx-lto.so: error loading plugin: libimf.so: 无法打开共享目标文件
#-----------------------------------------------------------------------------
set(INTEL_LIB_DIR "/opt/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin")
set(_ORIG_LIBRARY_PATH    "$ENV{LIBRARY_PATH}")
set(_ORIG_LD_LIBRARY_PATH "$ENV{LD_LIBRARY_PATH}")
set(ENV{LIBRARY_PATH}    "${INTEL_LIB_DIR}:${_ORIG_LIBRARY_PATH}")
set(ENV{LD_LIBRARY_PATH} "${INTEL_LIB_DIR}:${_ORIG_LD_LIBRARY_PATH}")

project(PPP LANGUAGES Fortran C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Release/Debug" FORCE)
endif()

#-----------------------------------------------------------------------------
# 产物输出目录 —— 对齐参考程序 /home/xfd/ckdrt_irc_cr_orb_sdb/CMakeLists.txt
#   目的: 让目录名由 build type 决定, 而不是由 CMake 的 binary dir 决定。
#   Trae/CMake Tools 默认 binary dir 是 ${workspaceFolder}/build, 若不管的话
#   产物落在 build/ppp, 与 .vscode/launch.json 里的 build_debug/ppp、
#   build_release/ppp 对不上。这里显式改写后, 只要在 Trae 里切换
#   Debug/Release 变体, ppp 就会自动落到对应目录, launch.json 直接可用。
#   注意: 二者是源码树下的目录, 已在 .gitignore 中排除。
#-----------------------------------------------------------------------------
if(CMAKE_BUILD_TYPE MATCHES "Release")
    set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/build_release")
    set(LIBRARY_OUTPUT_PATH    "${CMAKE_SOURCE_DIR}/build_release/Lib")
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/build_debug")
    set(LIBRARY_OUTPUT_PATH    "${CMAKE_SOURCE_DIR}/build_debug/Lib")
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/build_relwithdebinfo")
    set(LIBRARY_OUTPUT_PATH    "${CMAKE_SOURCE_DIR}/build_relwithdebinfo/Lib")
endif()

#-----------------------------------------------------------------------------
# MKL: 取自 exe_ppp/Makefile 的 MKLROOT (2024.2, 非参考程序的 2025.0)
#-----------------------------------------------------------------------------
set(MKLROOT     "/opt/intel/oneapi/mkl/2024.2")
set(MKL_lib     "${MKLROOT}/lib/intel64")
set(MKL_include "${MKLROOT}/include/intel64/lp64")

include_directories("${MKL_include}" "${MKLROOT}/include")
link_directories("${MKL_lib}" "${INTEL_LIB_DIR}")

#-----------------------------------------------------------------------------
# 链接时的 Intel 运行库环境 (关键)
#   编译 (-c) 和 ar 归档都不调 ld, 所以只有 final link 会踩这个坑: ld 要 dlopen
#   icx-lto.so 插件, 插件 NEEDED libimf.so。
#   实测结论: 上面的 set(ENV{...}) 只影响 configure 进程, 对 cmake --build 派生
#   出的链接子进程无效; 而 -L 与 -Wl,-rpath-link 都无法让 ld 找到插件的依赖,
#   只有给链接命令本身带上 LD_LIBRARY_PATH 才行。
#   因此把 Fortran 链接规则改成 cmake -E env 包装 —— 这样在 Trae 里直接 Build
#   即可, 不需要先 source oneAPI 的 setvars.sh。
#-----------------------------------------------------------------------------
set(_INTEL_LINK_ENV "LD_LIBRARY_PATH=${INTEL_LIB_DIR}")
if(NOT "${_ORIG_LD_LIBRARY_PATH}" STREQUAL "")
    set(_INTEL_LINK_ENV "${_INTEL_LINK_ENV}:${_ORIG_LD_LIBRARY_PATH}")
endif()
set(CMAKE_Fortran_LINK_EXECUTABLE
    "<CMAKE_COMMAND> -E env ${_INTEL_LINK_ENV} <CMAKE_Fortran_COMPILER> <FLAGS> <CMAKE_Fortran_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")

#-----------------------------------------------------------------------------
# 编译标志
#   Release 对齐参考程序的: -w -132 -traceback -O3 -xHost -qopt-zmm-usage=high
#                           -nocpp
#   未启用参考程序的三项 (与本次 Makefile 改动保持同一判断):
#     -qopenmp -save   : 本程序无 !$OMP 并行区, 无收益; 参考程序自身记录过
#                        -qopenmp 隐含 -auto 会压低 WL 固定率, 需配 -save 才
#                        恢复, 属被动修补, 不是收益。
#     -heap-arrays     : 不给阈值时把全部自动数组/临时量挪到堆, 会拖慢热循环,
#                        本程序无栈溢出证据。
#     -fimplicit-none  : 20 年老代码, 逐文件校验成本高, 未加。
#   注意: -xHost 会改变向量化路径与 FMA 收缩, 浮点结果末位会变, 需重跑验证。
#         本机 (i9-13900H) 无 AVX-512, -qopt-zmm-usage=high 实际不生效,
#         保留只为与参考程序一致。
#-----------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -w -132 -traceback -O3 -xHost -qopt-zmm-usage=high -nocpp")
    set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -w -132 -traceback -g -O0 -check bounds -nocpp")
endif()

#-----------------------------------------------------------------------------
# .mod 统一输出目录 (对应 Makefile 里各库共用的 -module ../OBJ_INTEL)
#-----------------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")
file(MAKE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}")
include_directories("${CMAKE_Fortran_MODULE_DIRECTORY}")

# lib_de405/pleph.f 需要 INCLUDE 'jpleph.h'; lib_comm 的 C 头文件
include_directories("${CMAKE_SOURCE_DIR}/lib_de405" "${CMAKE_SOURCE_DIR}/lib_comm")

#-----------------------------------------------------------------------------
# SOFA: 沿用预编译库
#   参考程序同样不在 CMake 里编 sofa —— 其注释写明: SOFA 不是标准语法, CMake 的
#   Release/Debug 结果有差异, 因此直接沿用 makefile 编出的 libsofa.a。
#   本程序同构处理, 依赖 OBJ_INTEL/libsofa.a (= 顶层 makefile 的 sofalib)。
#-----------------------------------------------------------------------------
set(SOFA_LIB "${CMAKE_SOURCE_DIR}/OBJ_INTEL/libsofa.a")
if(NOT EXISTS "${SOFA_LIB}")
    message(FATAL_ERROR "缺少 ${SOFA_LIB}; 请先执行: cd lib_sofa && make && make test")
endif()

#-----------------------------------------------------------------------------
# MKL 链接库 (照抄 exe_ppp/Makefile 的 Link_mul)
#   lp64 + intel_thread + core + iomp5, 与编译层的 -mkl=parallel 语义等价。
#-----------------------------------------------------------------------------
set(MKL_LIBS
    "${MKL_lib}/libmkl_blas95_lp64.a" "${MKL_lib}/libmkl_lapack95_lp64.a"
    -Wl,--start-group
    "${MKL_lib}/libmkl_intel_lp64.a" "${MKL_lib}/libmkl_intel_thread.a" "${MKL_lib}/libmkl_core.a"
    -Wl,--end-group
    iomp5 pthread m dl)

#=============================================================================
# 静态库 (对应顶层 makefile 的 moduleslib/de405lib/itrs2gcrslib/commlib/
#          comlib/lambdalib/tropolib/tidelib/ionlib/oilib)
#=============================================================================
file(GLOB MOD_SRC    CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_modules/*.f90")
file(GLOB DE405_SRC  CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_de405/*.f" "${CMAKE_SOURCE_DIR}/lib_de405/*.f90")
file(GLOB ITRS_SRC   CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_itrs2gcrs/*.f90" "${CMAKE_SOURCE_DIR}/lib_itrs2gcrs/*.f")
file(GLOB COMM_SRC   CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_comm/*.c")
file(GLOB COM_SRC    CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_com/*.f90")
file(GLOB LAMBDA_SRC CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_lambda/*.f")
file(GLOB TROPO_SRC  CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_tropo/*.f90" "${CMAKE_SOURCE_DIR}/lib_tropo/*.f")
file(GLOB TIDE_SRC   CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_tide/*.f90"  "${CMAKE_SOURCE_DIR}/lib_tide/*.f")
file(GLOB ION_SRC    CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_ion/*.f90"   "${CMAKE_SOURCE_DIR}/lib_ion/*.for")
file(GLOB OI_SRC     CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib_oi/*.f90"    "${CMAKE_SOURCE_DIR}/lib_oi/*.f")
# 说明: lib_modules / lib_com / exe_* 下的 *.i90 是老的固定格式死副本,
#       Makefile 不编译 (无 .i90 规则), 这里靠扩展名天然排除。

add_library(modules   STATIC ${MOD_SRC})
add_library(de405     STATIC ${DE405_SRC})
add_library(itrs2gcrs STATIC ${ITRS_SRC})
add_library(comm      STATIC ${COMM_SRC})
add_library(com       STATIC ${COM_SRC})
add_library(lambda    STATIC ${LAMBDA_SRC})
add_library(tropo     STATIC ${TROPO_SRC})
add_library(tide      STATIC ${TIDE_SRC})
add_library(ion       STATIC ${ION_SRC})
add_library(oi        STATIC ${OI_SRC})

# lib_com 全局 -fp-model precise (沿用 lib_com/Makefile 的 FFLAGS)
target_compile_options(com PRIVATE -fp-model precise)

# 库间依赖: 照抄顶层 makefile 的依赖关系 (静态库链接顺序 + 模块生成顺序)
# SOFA 是预编译静态库 (不是 CMake target), 单遍扫描下必须排在它的消费者之后。
# 直接写进 exe 的链接列表会被 CMake 排到消费者之前 (先 issue 直接项, 再 issue
# target 的传递依赖), 导致 itrs2gcrs 的 iau_* 解析不到。挂到消费者身上才能
# 借依赖图拿到正确顺序。引用 iau_* 的是 itrs2gcrs / oi / tide 三个库。
target_link_libraries(de405     modules)
target_link_libraries(itrs2gcrs modules ${SOFA_LIB})
target_link_libraries(tropo     modules)
target_link_libraries(tide      modules ${SOFA_LIB})
target_link_libraries(lambda    modules)
target_link_libraries(com       modules)
target_link_libraries(ion  de405 itrs2gcrs tropo tide com modules)
target_link_libraries(oi   de405 itrs2gcrs com modules ${SOFA_LIB})

#=============================================================================
# 可执行程序 (ppp / fcb / rtorb)
#=============================================================================
# --disable-new-dtags: 强制产出 DT_RPATH 而不是 DT_RUNPATH。
#   链接期只靠它不够 (见上), 但运行期必须要: libiomp5.so 自身 NEEDED
#   libintlc.so.5, 而现代 ld 默认发的 DT_RUNPATH 不会传给非直接依赖, 于是
#     ./ppp: error while loading shared libraries: libintlc.so.5
#   DT_RPATH 会沿整条依赖链生效, 这样产物不设 LD_LIBRARY_PATH 也能直接跑。
set(RPATH_FLAGS
    "-Wl,--disable-new-dtags"
    "-Wl,-rpath,${MKL_lib}"
    "-Wl,-rpath,${INTEL_LIB_DIR}")

#--------------------------------------------------------------- ppp ---------
# ⚠ exe_ppp/Makefile 的 ComFLAG = -i8 只作用于"库对象"; 主程序 ppp.f90 在链接
#   命令里是单独编的, 不带 -i8。这个整型 kind 分裂必须保留, 否则 ABI 会变。
#   因此把 exe_ppp 拆成两个 target: 库对象(带 -i8) + 主程序(不带 -i8)。
#   配套约束: 库对象里调 LAPACK 时整型出参必须写成 -1_IT 这类显式 4 字节 kind
#   (MKL lp64 是 4 字节整型)。
file(GLOB PPP_LIB_SRC CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/exe_ppp/*.f90")
list(REMOVE_ITEM PPP_LIB_SRC "${CMAKE_SOURCE_DIR}/exe_ppp/ppp.f90")

add_library(ppp_lib STATIC ${PPP_LIB_SRC})
target_compile_options(ppp_lib PRIVATE -i8)
target_link_libraries(ppp_lib oi modules de405 itrs2gcrs tropo ion tide comm com lambda)

add_executable(ppp "${CMAKE_SOURCE_DIR}/exe_ppp/ppp.f90")
target_link_libraries(ppp ppp_lib ${SOFA_LIB} ${MKL_LIBS})
target_link_options(ppp PRIVATE ${RPATH_FLAGS})

#--------------------------------------------------------------- fcb ---------
# exe_fcb/Makefile 的 FFLAGS 含全局 -fp-model precise, 主程序与库对象都带。
file(GLOB FCB_SRC CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/exe_fcb/*.f90")
list(REMOVE_ITEM FCB_SRC "${CMAKE_SOURCE_DIR}/exe_fcb/fcb.f90")

add_library(fcb_lib STATIC ${FCB_SRC})
target_compile_options(fcb_lib PRIVATE -fp-model precise)
# 顶层 makefile 的 fcblib 依赖 (无 lambda/oi)。必须显式写, 否则 CMake 不知道
# fcb_lib 要等 modules 先产出 ckdctrl.mod 等模块文件 —— 并行构建时 fcb_lib 会
# 抢在 modules 之前编译, 读到尚未写全的 .mod, 报
#   error #6460: This is not a component name ... [CPRN]
target_link_libraries(fcb_lib modules de405 itrs2gcrs tropo ion tide comm com)

add_executable(fcb "${CMAKE_SOURCE_DIR}/exe_fcb/fcb.f90")
target_compile_options(fcb PRIVATE -fp-model precise)
target_link_libraries(fcb fcb_lib modules de405 itrs2gcrs tropo ion tide comm com ${SOFA_LIB} ${MKL_LIBS})
target_link_options(fcb PRIVATE ${RPATH_FLAGS})

#------------------------------------------------------------- rtorb ---------
# exe_rtorb/Makefile 把 rtorb.f90 与 get_rtorb_args.f90 一起编, 都带 precise。
add_executable(rtorb
    "${CMAKE_SOURCE_DIR}/exe_rtorb/rtorb.f90"
    "${CMAKE_SOURCE_DIR}/exe_rtorb/get_rtorb_args.f90")
target_compile_options(rtorb PRIVATE -fp-model precise)
target_link_libraries(rtorb modules de405 itrs2gcrs com ${SOFA_LIB} ${MKL_LIBS})
target_link_options(rtorb PRIVATE ${RPATH_FLAGS})