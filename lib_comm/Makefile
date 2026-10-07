# Makefile for bnc2panda

all: commLib

PATH_OBJ = ../OBJ_INTEL
CFLAGS   = -O3 -g3 -Wall
CC       = gcc

LIB = $(PATH_OBJ)/libcomm.a

commLib : $(LIB)

$(LIB) : \
	$(LIB)(bridge.o) \
	$(LIB)(read_rtcm3.o) \
        $(LIB)(read_saved.o) \
        $(LIB)(shm_ckdrt.o) \
        $(LIB)(shm_writer.o) \
        $(LIB)(shm_reader.o) \
        $(LIB)(sem_init.o) \
        $(LIB)(sem_lock.o) \
        $(LIB)(sem_free.o) \
        $(LIB)(pth_brdck.o) \
        $(LIB)(pth_orbit.o) \
        $(LIB)(pth_cc2nc.o) \
        $(LIB)(pth_acqcf.o)

$(LIB)(shm_ckdrt.o): const.h clkobs.h sem.h
$(LIB)(shm_writer.o): const.h clkobs.h clkfcb.h rwsyn.h sem.h
$(LIB)(shm_reader.o): const.h clkfcb.h rwsyn.h sem.h
$(LIB)(pth_brdck.o): const.h
$(LIB)(pth_orbit.o): const.h
$(LIB)(pth_cc2nc.o): const.h
$(LIB)(pth_acqcf.o): const.h clkfcb.h

%.o : %.c
	$(CC) -c $(CFLAGS) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB)
