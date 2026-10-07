#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <errno.h>

int sem_init_l(int nsems)
{
  int semid;
  if( nsems > 2 ) {
    printf("At most 2 semaphores could be created\n");
    exit(1);
  }
// create a new semaphore set of `nsems' semaphores
  if( (semid = semget(IPC_PRIVATE, nsems, IPC_CREAT | 0666)) < 0 ) {
    perror( "error in creating semaphore" );
    printf( "semid : %d errno %d %d\n",semid,errno,EACCES );
    exit(1);
  } else {
    printf("Semaphore created successfully: %d %d\n",semid,nsems);
  }
//
// initialization of semaphores
// free buffer
  if( semctl(semid, 0, SETVAL, 1) < 0 ) {
    perror( "error in initializing first semaphore" );
    exit(1);
  }
// full buffer
  if( nsems > 1 ) {
    if( semctl(semid, 1, SETVAL, 0) < 0 ) {
      perror( "error in initializing second semaphore" );
      exit(1);
    }
  }
  return semid;
}
