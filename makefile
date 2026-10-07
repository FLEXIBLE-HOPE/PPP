#
#  Top Level Makefile for PANDA
#  Version 4.0.1
#  April 2025
#

OBJ_DIR = OBJ_INTEL

# 伪目标声明
.PHONY: all clean
.PHONY: moduleslib de405lib itrs2gcrslib sofalib tropolib ionlib tidelib commlib comlib lambdalib oilib
.PHONY: ppplib rtorblib fcblib

all: moduleslib de405lib itrs2gcrslib sofalib tropolib ionlib tidelib commlib comlib lambdalib oilib \
	ppplib rtorblib fcblib

# 确保 OBJ_INTEL 目录存在
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# ======== 底层库 (无依赖) ========

moduleslib: | $(OBJ_DIR)
	( cd lib_modules; $(MAKE) )

de405lib: | $(OBJ_DIR)
	( cd lib_de405; $(MAKE) )

sofalib: | $(OBJ_DIR)
	( cd lib_sofa; $(MAKE); $(MAKE) test )

tropolib: | $(OBJ_DIR)
	( cd lib_tropo; $(MAKE) )

tidelib: | $(OBJ_DIR)
	( cd lib_tide; $(MAKE) )

commlib: | $(OBJ_DIR)
	( cd lib_comm; $(MAKE) )

comlib: | $(OBJ_DIR)
	( cd lib_com; $(MAKE) )

lambdalib: | $(OBJ_DIR)
	( cd lib_lambda; $(MAKE) )

# ======== 上层库 (依赖底层库) ========

ionlib: moduleslib de405lib itrs2gcrslib sofalib tropolib tidelib comlib | $(OBJ_DIR)
	( cd lib_ion; $(MAKE) )

itrs2gcrslib: moduleslib | $(OBJ_DIR)
	( cd lib_itrs2gcrs; $(MAKE) )

oilib: moduleslib de405lib itrs2gcrslib sofalib comlib | $(OBJ_DIR)
	( cd lib_oi; $(MAKE) )

# ======== 可执行程序 ========

ppplib: oilib moduleslib de405lib itrs2gcrslib sofalib tropolib ionlib tidelib commlib comlib lambdalib
	( cd exe_ppp; $(MAKE) )

rtorblib: moduleslib de405lib itrs2gcrslib sofalib comlib
	( cd exe_rtorb; $(MAKE) )

fcblib: moduleslib de405lib itrs2gcrslib sofalib tropolib ionlib tidelib commlib comlib
	( cd exe_fcb; $(MAKE) )

# ======== 清理 ========

clean:
	-rm -rf $(OBJ_DIR)
	- ( cd lib_modules; $(MAKE) clean )
	- ( cd lib_de405; $(MAKE) clean )
	- ( cd lib_itrs2gcrs; $(MAKE) clean )
	- ( cd lib_sofa; $(MAKE) clean )
	- ( cd lib_tropo; $(MAKE) clean )
	- ( cd lib_ion; $(MAKE) clean )
	- ( cd lib_tide; $(MAKE) clean )
	- ( cd lib_comm; $(MAKE) clean )
	- ( cd lib_com; $(MAKE) clean )
	- ( cd lib_lambda; $(MAKE) clean )
	- ( cd lib_oi; $(MAKE) clean )
	- ( cd exe_ppp; $(MAKE) clean )
	- ( cd exe_rtorb; $(MAKE) clean )
	- ( cd exe_fcb; $(MAKE) clean )
