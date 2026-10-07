#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <errno.h>

void sem_free(int semid, int semnum)
{
  struct sembuf sops;  /* define operation on semaphore with given index */

  sops.sem_num = semnum;/* define operation on semaphore with given index */
  sops.sem_op  = 1;    /* add 1 to value for V operation */
  sops.sem_flg = 0;    /* type "man semop" in shell window for details */

  if( semop(semid, &sops, 1) == -1 ) {
    perror( "error in semaphore free operation" );
    exit(1);
  }
}
