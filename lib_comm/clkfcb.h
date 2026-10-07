#ifndef CLKFCB_H
#define CLKFCB_H

#include "const.h"

typedef struct {
	int    mjd;
	double sod;
	char   cprn[MAXSAT][3];
	int    ifreq[MAXSAT];
	double est[MAXSAT];
	double fwl[MAXSAT];
	double swl[MAXSAT];
	double fnl[MAXSAT];
	double snl[MAXSAT];
} CLK_FCB;

#endif
