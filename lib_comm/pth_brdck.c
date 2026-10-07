#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"
//
// Fortran functions to be used
void   mjd2doy_(int*, int*, int*);
int    modified_julday_(int*, int*, int*);
double timdif_(int*, double*, int*, double*);
void   brdtime_(char*, int*, double*);
//
// broadcast clocks
typedef struct {
  int    mjd;
  double sod;
  double a[3];
} BRDCK;
// argument struct
typedef struct {
  int    mjd;
  double sod;
  int    nprn;
  char   cprn[MAXSAT][3];
  BRDCK*  bk;
} BKARG;
// struct of first line for RINEX 3.X
typedef struct {
  char   cprn[3];
  char   cyr[5],cmon[3],cdat[3],chr[3],cmin[3],csec[3];
  char   ca0[19],ca1[19],ca2[19];
} LBRD;
//
// mutex
static BKARG bkarg;
static pthread_mutex_t pthmt = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t pthnb = PTHREAD_MUTEX_INITIALIZER;
//
// function to be called
void* pthFuncBrd(void* arg)
{
  if( pthread_mutex_lock(&pthnb) != 0 ) {
    perror("Thread broadcast function lock: ");
    exit(1);
  }
  int    i,mjd,isat,iy,imon,id,ih,im,istat;
  double sec,sod,dt,tmp[3];
  char   cprn[3];
  char   ccmd[80];
  char   line[100];
  LBRD*  plin;
  BKARG* parg = (BKARG*) arg;
//
// download broadcast ephemeris
  mjd=parg->mjd;
  do {
    mjd2doy_(&mjd,&iy,&id);
    iy=iy-2000;
    sprintf(ccmd,"rtgbrd %3.3d %2.2d",id,iy);
    istat=system(ccmd);
    //istat=0;
    if( istat == 0 ) {
     	break;
    }
  } while( mjd-- == parg->mjd );
//
// open file and skip header part
  if( istat == 0 ) {
    char* pch;
    FILE* pfil;
    if( (pfil = fopen("brdc.eph","r")) != NULL ) {
      memset(line,' ',100);
      while( (pch=strstr(line, "END OF HEADER")) == NULL ) fgets(line,100,pfil);
//
// read the whole file, but record the latest ones
      if( pthread_mutex_lock(&pthmt) != 0 ) {perror("Mutex lock in pthFuncBrd: ");exit(1);}
      while( (pch=fgets(line,100,pfil)) != NULL ) {
        if( strncmp(line,"  ",2) != 0 ) {
          do {
            if( *pch == 'D' || *pch == 'd' ) *pch = 'e';
          } while( *pch++ != '\0' );
          plin = (LBRD*) line;
          strncpy(ccmd,plin->cprn,3); ccmd[3] ='\0'; memcpy(cprn,ccmd,sizeof(char)*3);
          strncpy(ccmd,plin->cyr, 5); ccmd[5] ='\0'; iy     = atoi(ccmd);
          strncpy(ccmd,plin->cmon,3); ccmd[3] ='\0'; imon   = atoi(ccmd);
          strncpy(ccmd,plin->cdat,3); ccmd[3] ='\0'; id     = atoi(ccmd);
          strncpy(ccmd,plin->chr, 3); ccmd[3] ='\0'; ih     = atoi(ccmd);
          strncpy(ccmd,plin->cmin,3); ccmd[3] ='\0'; im     = atoi(ccmd);
          strncpy(ccmd,plin->csec,3); ccmd[3] ='\0'; sec    = atof(ccmd);
          strncpy(ccmd,plin->ca0,19); ccmd[19]='\0'; tmp[0] = atof(ccmd);
          strncpy(ccmd,plin->ca1,19); ccmd[19]='\0'; tmp[1] = atof(ccmd);
          strncpy(ccmd,plin->ca2,19); ccmd[19]='\0'; tmp[2] = atof(ccmd);
          //printf("%s %f",cprn,tmp[0]);
          //if( fabs(tmp[1]) < 1.0e-15 ) continue;
          for (i=0;i<MAXSAT;i++){
        	isat = strncmp(&parg->cprn[i][0],cprn,3);
        	if (isat == 0) {
        	  isat=i;
        	  break;
        	}
          }
          if (i == MAXSAT) isat = -1;
          if( isat >= 0 ) {
            mjd = modified_julday_(&id,&imon,&iy);
            sod = ih*3600.0+im*60.0+sec;
            brdtime_(cprn,&mjd,&sod);
            dt  = timdif_(&parg->mjd,&parg->sod,&mjd,&sod);
            if (strncmp(&cprn[0],"R",1) == 0) {
            	tmp[1]=0.0;
            	tmp[2]=0.0;
            }
// JG: 6 hours, read the latest epoch, but the time differences should be less then 6 hours
// maybe this hould be changed
            if( dt > 1.0e-15 && dt < 21600.0 ) {
              parg->bk[isat].mjd = mjd;
              parg->bk[isat].sod = sod;
              memcpy(parg->bk[isat].a,tmp,sizeof(double)*3);
            }
          }
        }
      }
      pthread_mutex_unlock(&pthmt);
      fclose(pfil);
    }
  }
  pthread_mutex_unlock(&pthnb);

  return NULL;
}
//
// Main
void pth_brdck_(int* mjd, double* sod, int* nprn, char cprn[MAXSAT][3], BRDCK bk[MAXSAT])
{
  if( pthread_mutex_trylock(&pthnb) == 0 ) {
    bkarg.mjd  = *mjd;
    bkarg.sod  = *sod;
    bkarg.nprn = *nprn;
    bkarg.bk   = bk;
    memcpy(bkarg.cprn,cprn,sizeof(char)*MAXSAT*3);
//
// create thread
    pthread_t pth;
    if( pthread_create(&pth,NULL,pthFuncBrd,(void*)&bkarg) != 0 ) {
      printf("***ERROR(pth_brdck): thread fails");
      exit(1);
    } else {
      printf("Thread for broadcast ephemeris is created successfully!\n");
      if( pthread_detach(pth) == 0 ) printf("Thread detached successfully!\n");
    }
    pthread_mutex_unlock(&pthnb);
  }
}
//
// check whether broadcast clocks are zero
void pth_brdck_chk_(BRDCK bk[MAXSAT], int* nprn, double obs[2*MAXFREQ][MAXSAT])
{
  int i,isat;
  if( pthread_mutex_lock(&pthmt) != 0 ) {perror("Mutex lock in pth_brdck_chk: ");exit(1);}
  for( i=0; i<MAXFREQ; i++ ) {
    for( isat=0; isat<*nprn; isat++ ) {
      if( fabs(bk[isat].a[0]) < 1.0e-15 ) {
    	  obs[i][isat]=0.0;
    	  obs[MAXFREQ+i][isat]=0.0;
      }
    }
  }
  pthread_mutex_unlock(&pthmt);
}
//
// calculate satellite clocks
void pth_brdck_clk_(BRDCK* bk, int* mjd, double* sod, double* dsatclk)
{
  double dt;
  if( pthread_mutex_lock(&pthmt) != 0 ) {perror("Mutex lock in pth_brdck_clk: ");exit(1);}
  dt = (*mjd - bk->mjd) * 86400.0 + *sod - bk->sod;
  *dsatclk = bk->a[0] + (bk->a[1] + bk->a[2] * dt) * dt;
  pthread_mutex_unlock(&pthmt);
}
