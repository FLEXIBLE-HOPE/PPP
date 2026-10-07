#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"

//
// function called
int modified_julday_(int*, int*, int*);
double timdif_(int*, double*, int*, double*);
void timinc_(int*, double*, double*, int*, double*);
//
// argument struct
typedef struct {
  int     mjd;
  int*    mjd_sav;
  double  sod;
  double* sod_sav;
  double* bias;
} CCARG;
//
// mutex
static CCARG ccarg;
static pthread_mutex_t pthmt = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t pthnb = PTHREAD_MUTEX_INITIALIZER;
//
// function to be called
void* pthFunc2nc(void* arg)
{
  if( pthread_mutex_lock(&pthnb) != 0 ) {
    perror("Thread cc2nc function lock: ");
    exit(1);
  }
  CCARG* parg = (CCARG*) arg;
// download predicted p1c1bias file
  int  i,j,istat,iy,imon,id;
  char line[100];
  istat=system("rm -f p1c1bias.2000p");
  istat=system("wget -nv ftp://ftp.aiub.unibe.ch/bcwg/cc2noncc/p1c1bias.2000p");
  istat=0;
  if(istat == 0) {
    FILE* pfil;
    char* pch;
    pfil=fopen("p1c1bias.2000p","r");
    if( pthread_mutex_lock(&pthmt) != 0 ) {perror("Mutex lock in pthFunc2nc: ");exit(1);}
    while( (pch=fgets(line,100,pfil)) != NULL ) {
      if( (i=strncmp(pch," h ",3)) == 0 ) {
        sscanf(pch+2,"%d%d%d",&iy,&imon,&id);
        *parg->mjd_sav=modified_julday_(&id,&imon,&iy);
        for( i=0; i<8; i++ ) {
          pch=fgets(line,100,pfil);
          do {
            if( *pch == 'D' || *pch == 'd' ) *pch = 'e';
          } while( *pch++ != '\0' );
          for( j=0; j<5; j++){
        	  sscanf(line+17+10*j,"%lf",&parg->bias[i*5+j]);
          }
          /*sscanf(line+17,"%lf,%lf,%lf,%lf,%lf",&parg->bias[i*5],&parg->bias[i*5+1],
                          &parg->bias[i*5+2],&parg->bias[i*5+3],&parg->bias[i*5+4]);*/
        }
      }
    }
    fclose(pfil);
    while( timdif_(parg->mjd_sav,parg->sod_sav,&parg->mjd,&parg->sod) < -1.0e-15 ) {
      double tmp = 86400.0;
      timinc_(parg->mjd_sav,parg->sod_sav,&tmp,parg->mjd_sav,parg->sod_sav);
    }
    pthread_mutex_unlock(&pthmt);
  } else {
    printf("***ERROR(pth_cc2nc): p1c1bias.2000p missing\n");
    exit(1);
  }
  pthread_mutex_unlock(&pthnb);

  return NULL;
}
//
// Main
void pth_cc2nc_(int* mjd_sav, double* sod_sav, int* mjd, double* sod, double bias[MAXSAT])
{
  if( pthread_mutex_trylock(&pthnb) == 0 ) {
// whether update p1c1bias.2000p
    double dt;
    dt = timdif_(mjd_sav,sod_sav,mjd,sod);
    if( dt < -1.0e-15 ) {
      ccarg.mjd     = *mjd;
      ccarg.sod     = *sod;
      ccarg.mjd_sav = mjd_sav;
      ccarg.sod_sav = sod_sav;
      ccarg.bias    = bias;
//
// create thread
      pthread_t pth;
      if( pthread_create(&pth,NULL,pthFunc2nc,(void*)&ccarg) != 0 ) {
        printf("***ERROR(pth_cc2nc): thread fails");
        exit(1);
      } else {
        printf("Thread for p1c1bias.2000p is created successfully!\n");
        if( pthread_detach(pth) == 0 ) printf("Thread detached successfully!\n");
      }
    }
    pthread_mutex_unlock(&pthnb);
  }
}
//
// correct pseudorange, or clear observations
void pth_cc2nc_corr_(char cprn[MAXSAT][3], double bias[MAXSAT], double obs[2*MAXFREQ][MAXSAT], char fob[2*MAXFREQ][MAXSAT][3])
{
  int i,k;
  if( pthread_mutex_lock(&pthmt) != 0 ) {perror("Mutex lock in pth_cc2nc_corr: ");exit(1);}
  for( i=0; i<MAXSAT; i++ ) {
// JG: -9.9d9 for none
    if( bias[i] > -9.0e8 ) break;
  }
  if( i == MAXSAT ) {
    memset(obs,0,sizeof(double)*2*MAXFREQ*MAXSAT);
  } else {
    for( i=0; i<MAXSAT; i++ ) {
    	if ( strncmp(cprn[i],"G",1) == 0) {
    		if( fabs(obs[MAXFREQ][i]) > 1.0e-15 && fabs(obs[MAXFREQ+1][i]) > 1.0e-15 ) {
    			k=atol(&cprn[i][1]);
    			if( bias[k-1] > -9.0e8 ) {
    				if( strncmp(fob[MAXFREQ][i],"C1C",3) == 0 ) obs[MAXFREQ][i]+=bias[k-1]/3.335641;
    				if( strncmp(fob[MAXFREQ+1][i],"C2D",3) == 0 ) obs[MAXFREQ+1][i]+=bias[k-1]/3.335641;
    			}
    		}
    	}
    }
  }
  pthread_mutex_unlock(&pthmt);
}
