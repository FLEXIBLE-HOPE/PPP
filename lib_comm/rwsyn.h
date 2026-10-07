#ifndef RWSYN_H
#define RWSYN_H

#include "clkfcb.h"

//
// struct for writer-reader synchronization
typedef struct {
	int     semid;      // semaphore for writer and readers
	int     mutid;      // mutex for readers
	int     semts;      // semaphore in case of writer starvation
	int     nreader;    // number of readers in the shared memory
	CLK_FCB CF;         // satellite clocks and FCBs
} RWSYN;

#endif
