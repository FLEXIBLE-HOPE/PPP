#include <stdio.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"
#include "clkobs.h"
#include "clkfcb.h"
#include "rwsyn.h"
#include "sem.h"

//
// Global variables
static int      iexit = 0;    // flag for SIGINT exit
static CLK_OBS* pcob  = NULL; // pointer to shared memory for clock and observations
static RWSYN*   pshm  = NULL; // pointer to shared memory for clock and FCBs
static int      shmid = 0;    // ID for shared memory
//
// SIGINT signal handler
void sighandler(int signum)
{
  iexit = 1;
  return;
}
//
// destruction
void shm_clean()
{
  if( semctl(pshm->semid, 0, IPC_RMID) < 0 ) {
    perror("semctl(IPC_RMID) failed in shm_writer.c for semid\n");
    exit(1);
  }
  if( semctl(pshm->mutid, 0, IPC_RMID) < 0 ) {
    perror("semctl(IPC_RMID) failed in shm_writer.c for mutid\n");
    exit(1);
  }
  if( semctl(pshm->semts, 0, IPC_RMID) < 0 ) {
    perror("semctl(IPC_RMID) failed in shm_writer.c for semts\n");
    exit(1);
  }
  if( shmdt(pcob) == -1 ) {
    perror("shmdt failed in shm_writer.c for detaching memory of pcob\n");
    exit(1);
  }
  if( shmdt(pshm) == -1 ) {
    perror("shmdt failed in shm_writer.c for detaching memory of pshm\n");
    exit(1);
  }
  if( shmctl(shmid, IPC_RMID, 0) == -1 ) {
    perror("shmctl(IPC_RMID) failed in shm_writer.c for deleting memory\n");
    exit(1);
  }
  exit(0);
}
//
// Read shared memory from ckdrt
void shm_read_clkobs_(int* mjd, double* sod, int* nsit, char snam[MAXSIT][4],
                      char cprn[MAXSAT][3], int ifreq[MAXSAT], double est[MAXSAT], double obs[MAXSIT][2*MAXFREQ][MAXSAT])
{ 
  static int      iattach=0;
  int             shmck;
  key_t           mem_key;
// get shared memory
  if( iattach == 0 ) {
    printf("\nGetting the shared memory ...\n"); 
    mem_key = ftok("/dev/null", 'a');
    if( (shmck=shmget(mem_key, sizeof(CLK_OBS), IPC_CREAT | 0666)) == -1 ) {
      perror("shmget");
      printf("shmid : %d errno %d %d\n",shmck,errno,EACCES);
    }
    if ( (pcob=(CLK_OBS*)shmat(shmck, NULL, 0)) == (void *)(-1) ) {
      perror("shmat error: ");
      exit(1);
    } else {
      printf("Shared memory is successfully attached %d\n", shmck);
    }
    iattach=1;
  }
//
// read clocks
  sem_lock(pcob->semid,1);
  *mjd  = pcob->mjd;
  *sod  = pcob->sod;
  *nsit = pcob->nsit;
  memcpy(snam,pcob->snam,sizeof(char)*4*MAXSIT);
  memcpy(cprn,pcob->cprn, sizeof(char)*3*MAXSAT);
  memcpy(ifreq,pcob->ifreq,sizeof(int)*MAXSAT);
  memcpy(est, pcob->est, sizeof(double)*MAXSAT);
  memcpy(obs, pcob->obs, sizeof(double)*MAXSIT*2*MAXFREQ*MAXSAT);
  sem_free(pcob->semid,0);
//
// whether exit
  if( iexit == 1 ) shm_clean();
}
//
// Main writer
void shm_writer_(int* mjd, double* sod, char cprn[MAXSAT][3], int ifreq[MAXSAT], double est[MAXSAT],
                double fwl[MAXSAT], double swl[MAXSAT], double fnl[MAXSAT], double snl[MAXSAT])
{
  static int      icreate = 0;
  key_t           mem_key;
// create shared memory
  //printf("%d %f\n",*mjd,*sod);
  if( icreate == 0 ) {
    printf("\nCreating shared memory for satellite clocks and FCBs ...\n");
    mem_key = ftok("/dev/null", 'c');
    if( (shmid = shmget(mem_key, sizeof(RWSYN), IPC_CREAT | 0666)) == -1) {
      perror("shmget error: ");
      printf("shmid : %d errno %d %d\n",shmid,errno,EACCES);
    }
    if( (pshm = (RWSYN*)shmat(shmid, NULL, 0)) == (void *)(-1) ) {
      perror("shmat error: ");
      exit(1);
    } else {
      printf("Shared memory is successfully attached %d\n", shmid);
    }
    memset(pshm, 0, sizeof(RWSYN));
    icreate=1;           // created already
    pshm->semid=sem_init_l(1); // for writing
    pshm->mutid=sem_init_l(1); // for readers
    pshm->semts=sem_init_l(1); // for reading
    signal(SIGINT, sighandler);
  }
//
// write GPS satellite clocks and FCBs
  sem_lock(pshm->semts,0);
  sem_lock(pshm->semid,0);
  pshm->CF.mjd = *mjd;
  pshm->CF.sod = *sod;
  memcpy(pshm->CF.cprn,cprn,sizeof(char)*3*MAXSAT);
  memcpy(pshm->CF.ifreq,ifreq,sizeof(int)*MAXSAT);
  memcpy(pshm->CF.est,est,sizeof(double)*MAXSAT);
  memcpy(pshm->CF.fwl,fwl,sizeof(double)*MAXSAT);
  memcpy(pshm->CF.swl,swl,sizeof(double)*MAXSAT);
  memcpy(pshm->CF.fnl,fnl,sizeof(double)*MAXSAT);
  memcpy(pshm->CF.snl,snl,sizeof(double)*MAXSAT);
  sem_free(pshm->semts,0);
  sem_free(pshm->semid,0);
  //printf("%d %f\n",pshm->CF.mjd,pshm->CF.sod);
//
// whether exit
  //printf("%d\n",iexit);
  if( iexit == 1 ) shm_clean();
}
