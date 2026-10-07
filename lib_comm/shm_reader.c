#include <stdio.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <errno.h>
#include <signal.h>
#include <string.h>
#include "const.h"
#include "sem.h"
#include "clkfcb.h"
#include "rwsyn.h"

//
// Global variables
static int    iexit = 0; // flag for SIGINT
//
// SIGINT signal handler
void sighandler(int signum)
{
  iexit = 1;
  return;
}
//
// Main
_Bool shm_reader(CLK_FCB* CF)
{
  static RWSYN* pshm = NULL;  // pointer to shared memory
  static int    iattach = 0;
  int           shmid;
  key_t         mem_key;
// get shared memory
  if( iattach == 0 ) {
    printf("\nGetting the shared memory ...\n");
    mem_key = ftok("/dev/null", 'c');
    if( (shmid =shmget(mem_key, sizeof(RWSYN), IPC_CREAT | 0666)) == -1) {
      perror("shmget error: ");
      printf("shmid : %d errno %d %d\n",shmid,errno,EACCES);
    }
    if ( (pshm = (RWSYN*)shmat(shmid, NULL, 0)) == (void *)(-1) ) {
      perror("shmat error: ");
      exit(1);
    } else {
      printf("Shared memory is successfully attached %d\n", shmid);
    }
    iattach=1;
    signal(SIGINT, sighandler);
  }
//
// read clocks and FCBs
  sem_lock(pshm->semts,0);
  sem_free(pshm->semts,0);
  sem_lock(pshm->mutid,0);
  pshm->nreader += 1;
  if( pshm->nreader == 1 ) sem_lock(pshm->semid,0);
  sem_free(pshm->mutid,0);
  CF->mjd = pshm->CF.mjd;
  CF->sod = pshm->CF.sod;
  memcpy(CF->cprn,pshm->CF.cprn,sizeof(char)*3*MAXSAT);
  memcpy(CF->ifreq,pshm->CF.ifreq,sizeof(int)*MAXSAT);
  memcpy(CF->est,pshm->CF.est,sizeof(double)*MAXSAT);
  memcpy(CF->fwl,pshm->CF.fwl,sizeof(double)*MAXSAT);
  memcpy(CF->swl,pshm->CF.swl,sizeof(double)*MAXSAT);
  memcpy(CF->fnl,pshm->CF.fnl,sizeof(double)*MAXSAT);
  memcpy(CF->snl,pshm->CF.snl,sizeof(double)*MAXSAT);
  sem_lock(pshm->mutid,0);
  pshm->nreader -= 1;
  if( pshm->nreader == 0 ) sem_free(pshm->semid,0);
  sem_free(pshm->mutid,0);
// Detach Memory address and exit
  if( iexit == 1 ) {
    if( shmdt(pshm) == -1 ) {
      perror("shmdt failed in shm_reader.c for detaching memory\n");
      exit(1);
    }
    exit(0);
  }
  return 1;
}
