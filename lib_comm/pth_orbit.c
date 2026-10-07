#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"
//
// binary orbit header
typedef struct {
  int    nprn;
  char   cprn[MAXSAT][3];
  int    mjd0,mjd1,rmjd,nequ;
  double sod0,sod1,rsod;
  double dintv;
} ORBHD;
//
// Fortran functions to be used
void read_igserp_(char*,int*,double*,double*,double*);
void rdorbh_(char*,int*,ORBHD*);
void everett_reset_();
void igserp_reset_();
//
// argument struct
typedef struct {
  char   cweek[5];
  char   cwkd[2];
  char   chour[3];
  char*  flnorb;
  char*  flnerp;
} OBARG;
//
// mutex
static OBARG obarg;
static pthread_mutex_t pthmt = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t pthnb = PTHREAD_MUTEX_INITIALIZER;
//
// function to be called
void* pthFuncOrb(void* arg)
{
  if( pthread_mutex_lock(&pthnb) != 0 ) {
    perror("Thread orbit function lock: ");
    exit(1);
  }
  int    istat;
  char*  pch;
  char   str[20];
  char   ccmd[80] = "rtgorb ";
  OBARG* parg = (OBARG*) arg;
// download requested orbit and erp files
  strncpy(str,parg->cweek,4); str[4]='\0';
  strncat(str,parg->cwkd,1);  str[5]='\0';
  strncat(str,"_",1);
  strncat(str,parg->chour,2); str[8]='\0';
  strncat(ccmd,parg->cweek,4);
  strncat(ccmd," ",1);
  strncat(ccmd,parg->cwkd,1);
  strncat(ccmd," ",1);
  strncat(ccmd,parg->chour,2);
  ccmd[16]='\0';
  if( strlen(parg->flnorb) == 4 ) {
    if( (istat=system(ccmd)) == 0 ) {
	  //istat=0;
      //if (istat == 0) {
      if( pthread_mutex_lock(&pthmt) != 0 ) {
        perror("Mutex lock in pthFuncOrb: ");
        exit(1);
      }
      strcat(parg->flnorb,str);
      strcat(parg->flnerp,str);
      pthread_mutex_unlock(&pthmt);
    } else {
      printf("***ERROR(pth_orbit): no ultra-rapid orbit or ERP files");
      exit(1);
    }
  } else if( strstr(parg->flnorb,str) == NULL ) {
    if( (istat=system(ccmd)) == 0 ) {
      if( pthread_mutex_lock(&pthmt) != 0 ) {
        perror("Mutex lock in pthFuncOrb: ");
        exit(1);
      }
      pch = parg->flnorb;
      while(*pch++ != '_') ;
      *pch='\0';
      strcat(parg->flnorb,str);
      pch=parg->flnerp;
      while(*pch++ != '_') ;
      *pch='\0';
      strcat(parg->flnerp,str);
      everett_reset_();
      igserp_reset_();
      pthread_mutex_unlock(&pthmt);
    }
  }
  pthread_mutex_unlock(&pthnb);

  return NULL;
}
//
// Main
void pth_orbit_(char cweek[4], char cwkd[1], char chour[2], char flnorb[1024], char flnerp[1024])
{
  if( pthread_mutex_trylock(&pthnb) == 0 ) {
//
// copy arguments into struct
    obarg.flnorb = flnorb;
    obarg.flnerp = flnerp;
    strncpy(obarg.cweek,cweek,4); obarg.cweek[4]='\0';
    strncpy(obarg.cwkd, cwkd, 1); obarg.cwkd[1] ='\0';
    strncpy(obarg.chour,chour,2); obarg.chour[2]='\0';
//
// create thread
    char  str[20];
    strcpy(str,obarg.cweek);
    strcat(str,obarg.cwkd);
    strcat(str,"_");
    strcat(str,obarg.chour);
    if( strlen(flnorb) == 4 || strlen(flnerp) == 4 || strstr(flnorb,str) == NULL ) {
      pthread_t pth;
      if( pthread_create(&pth,NULL,pthFuncOrb,(void*)&obarg) != 0 ) {
        printf("***ERROR(pth_orbit): thread fails");
        exit(1);
      } else {
        printf("Thread for satellite orbit is created successfully!\n");
        if( pthread_detach(pth) == 0 ) printf("Thread detached successfully!\n");
      }
    }
    pthread_mutex_unlock(&pthnb);
  }
}
//
// add '\0' to flnorb and  flnerp
// 1024 should be consistent with LEN_FILENAME in mod_par.f90
void pth_orbit_null_(char flnorb[1024], char flnerp[1024])
{
  if( pthread_mutex_lock(&pthmt) != 0 ) {
    perror("Mutex lock in pth_orbit_null: ");
    exit(1);
  }
  flnorb[4]='\0';
  flnerp[4]='\0';
  pthread_mutex_unlock(&pthmt);
}
//
// check whether orbit and erp files are present
// 1024 should be consistent with LEN_FILENAME in mod_par.f90
void pth_orbit_chk_(char flnorb[1024], char flnerp[1024], int* nprn, double obs[2*MAXFREQ][MAXSAT])
{
  if( pthread_mutex_lock(&pthmt) != 0 ) {
    perror("Mutex lock in pth_orbit_chk: ");
    exit(1);
  }
  if( strlen(flnorb) == 4 || strlen(flnerp) == 4 ) memset(obs,0,sizeof(double)*2*MAXFREQ*MAXSAT);
  pthread_mutex_unlock(&pthmt);
}
//
// orbit
_Bool pth_orbit_first_(_Bool* first)
{
  if( pthread_mutex_lock(&pthmt) != 0 ) {
    perror("Mutex lock in pth_orbit_rdorbh: ");
    exit(1);
  }
  _Bool iret = *first;
  if( *first != '\000' ) *first = '\000';
  pthread_mutex_unlock(&pthmt);
  return iret;
}
//
// read orbit header
void pth_orbit_rdorbh_(char* orbfil, int* iunit, ORBHD* OH)
{
  if( pthread_mutex_lock(&pthmt) != 0 ) {
    perror("Mutex lock in pth_orbit_rdorbh: ");
    exit(1);
  }
  printf("%s",orbfil);
  rdorbh_(orbfil,iunit,OH);
  pthread_mutex_unlock(&pthmt);
}
//
// read erp
void pth_orbit_igserp_(char* erpfil, int* jd, double* sod, double* taiut1r, double xhelp[2])
{
  if( pthread_mutex_lock(&pthmt) != 0 ) {
    perror("Mutex lock in pth_orbit_igserp: ");
    exit(1);
  }
  read_igserp_(erpfil,jd,sod,taiut1r,xhelp);
  pthread_mutex_unlock(&pthmt);
}
