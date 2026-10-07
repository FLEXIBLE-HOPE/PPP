#include <stdio.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <errno.h>
#include <signal.h>
#include <string.h>
#include "const.h"
#include "clkobs.h"
#include "sem.h"
//
// Global variables
static int      iexit = 0; // flag for SIGINT exit
static CLK_OBS* pshm  = NULL;
//
// Signal handler
void sighandler(int signum)
{
  iexit = 1;
  if( semctl(pshm->semid, 0, GETVAL) == 0 ) sem_free(pshm->semid,0);
  return;
}
//
// Main
void shm_writ_clkobs_(int* mjd, double* sod, int* nsit, char snam[MAXSIT][4],
                      char cprn[MAXSAT][3], int ifreq[MAXSAT], double est[MAXSAT], double obs[MAXSIT][2*MAXFREQ][MAXSAT])
{
  static int      shmid = 0;
  static int      icreate = 0;
  key_t           mem_key;
// create shared memory
  if( icreate == 0 ) {
    printf("\nSetting up shared memory for FCB generation ...\n");
    mem_key = ftok("/dev/null", 'a');
    if( (shmid=shmget(mem_key, sizeof(CLK_OBS), IPC_CREAT | 0666)) == -1) {
      perror("shmget");
      printf("shmid : %d errno %d %d\n",shmid,errno,EACCES);
    }
    if( (pshm=(CLK_OBS*)shmat(shmid, NULL, 0)) == (void *)(-1) ) {
      perror("shmat");
      exit(1);
    } else {
      printf("Shared memory is successfully attached %d\n", shmid);
    }
    memset(pshm, 0, sizeof(CLK_OBS));
    icreate=1;
    pshm->semid=sem_init_l(2);
    signal(SIGINT, sighandler);
  }
// write GPS satellite clocks
  sem_lock(pshm->semid,0);
  pshm->mjd  = *mjd;
  pshm->sod  = *sod;
  pshm->nsit = *nsit;
  memcpy(pshm->snam,snam,sizeof(char)*4*MAXSIT);
  memcpy(pshm->cprn,  cprn,sizeof(char)*3*MAXSAT);
  memcpy(pshm->ifreq,ifreq,sizeof(int)*MAXSAT);
  memcpy(pshm->est,  est,sizeof(double)*MAXSAT);
  memcpy(pshm->obs,  obs,sizeof(double)*MAXSIT*2*MAXFREQ*MAXSAT);
  sem_free(pshm->semid,1);
//
// whether exit
  if( iexit == 1 ) {
    if( semctl(pshm->semid, 0, IPC_RMID) < 0 ) {
      perror("semctl(IPC_RMID) failed in shm_writ_clkobs.c\n");
      exit(1);
    }
    if( shmdt(pshm) == -1 ) {
      perror("shmdt failed in shm_writ_clkobs.c\n");
      exit(1);
    }
    if( shmctl(shmid, IPC_RMID, 0) == -1 ) {
      perror("shmctl(IPC_RMID) failed in shm_writ_clkobs.c\n");
      exit(1);
    }
    exit(0);
  }
}
