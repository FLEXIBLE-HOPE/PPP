#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"
#include "clkfcb.h"
//
// functions to be used
double timdif_(int*,double*,int*,double*);
_Bool shm_reader(CLK_FCB* CF);
//
// Global variables
static CLK_FCB CF[MAXREC];
static pthread_mutex_t pthmt = PTHREAD_MUTEX_INITIALIZER;
//
// Function to be called
void* pthFuncAcq(void* arg)
{
  int     i;
  CLK_FCB tmCF;
  while( shm_reader(&tmCF) ) {
    pthread_mutex_lock(&pthmt);
    if( fabs(timdif_(&tmCF.mjd,&tmCF.sod,&CF[0].mjd,&CF[0].sod)) > 1.0e-6 ) {
      for( i=MAXREC-1; i>0; i-- ) {
        memcpy(&CF[i],&CF[i-1],sizeof(CLK_FCB));
      }
      memcpy(&CF[0],&tmCF,sizeof(CLK_FCB));
    }
    pthread_mutex_unlock(&pthmt);
  }

  return NULL;
}
//
// Main
void pth_acqcf_()
{
  pthread_t pth;
  if( pthread_create(&pth,NULL,pthFuncAcq,NULL) != 0 ) {
    printf("***ERROR(pth_acqcf): thread fails");
    exit(1);
  } else {
    printf("Thread for clock&FCB acquisition is created successfully!\n");
    if( pthread_detach(pth) == 0 ) printf("Thread detached successfully!\n");
  }
}
//
// locate clocks and FCBs
_Bool pth_acqcf_read_(int* mjd, double* sod, char cprn[MAXSAT][3], int ifreq[MAXSAT], double est[MAXSAT],
                      double fwl[MAXSAT], double swl[MAXSAT], double fnl[MAXSAT], double snl[MAXSAT])
{
  int i;
  pthread_mutex_lock(&pthmt);
  for( i=0; i<MAXREC; i++ ) {
    //printf("%d,%d,%f,%d,%f\n",i,*mjd,*sod,CF[i].mjd,CF[i].sod);
    if( fabs(timdif_(mjd,sod,&CF[i].mjd,&CF[i].sod)) < 1.0e-6 ) {
      memcpy(cprn,CF[i].cprn,sizeof(char)*3*MAXSAT);
      memcpy(ifreq,CF[i].ifreq,sizeof(int)*MAXSAT);
      memcpy(est,CF[i].est,sizeof(double)*MAXSAT);
      memcpy(fwl,CF[i].fwl,sizeof(double)*MAXSAT);
      memcpy(swl,CF[i].swl,sizeof(double)*MAXSAT);
      memcpy(fnl,CF[i].fnl,sizeof(double)*MAXSAT);
      memcpy(snl,CF[i].snl,sizeof(double)*MAXSAT);
      break;
    }
  }
  pthread_mutex_unlock(&pthmt);
  if( i == MAXREC ) return 0;

  return 1;
}
