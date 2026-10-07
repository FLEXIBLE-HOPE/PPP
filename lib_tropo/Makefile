# Makefile for modulesLib Sep 01 13:09:26 2009

all: tropoLib


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FFLAGS7  = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB = $(PATH_OBJ)/libtropo.a

tropoLib : $(LIB)

$(LIB) : \
	$(LIB)(drysaas.o) \
        $(LIB)(gpt2.o) \
        $(LIB)(gpt3.o) \
        $(LIB)(gpt2_1w.o) \
        $(LIB)(gpt.o) \
        $(LIB)(gmf.o) \
	$(LIB)(VMF1_HT.o) \
	$(LIB)(vmf1.o) \
	$(LIB)(wetpp.o) \
	$(LIB)(wetsaas.o) \
	$(LIB)(wpress.o) \
	$(LIB)(FCUL_A.o) \
	$(LIB)(FCUL_B.o) \
	$(LIB)(FCUL_ZD_HPA.o) \

%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o
 
%.o : %.f
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB)
