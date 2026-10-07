#ifndef CLKOBS_H
#define CLKOBS_H

#include "const.h"

typedef struct {
	int    mjd;
	double sod;
	int    nsit;
	int    semid;
	char   snam[MAXSIT][4];
	char   cprn[MAXSAT][3];
	int    ifreq[MAXSAT];
	double est[MAXSAT];
	double obs[MAXSIT][2*MAXFREQ][MAXSAT];
} CLK_OBS;

#endif
