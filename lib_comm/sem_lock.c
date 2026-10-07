#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <errno.h>

void sem_lock(int semid, int semnum)
{
  struct sembuf sops;  /* only one semaphore operation to be executed */

  sops.sem_num = semnum; /* define operation on semaphore with given index */
  sops.sem_op  = -1;   /* subtract 1 to value for P operation */
  sops.sem_flg = 0;    /* type "man semop" in shell window for details */

  int istat;
  while( (istat=semop(semid, &sops, 1)) == -1 && errno == EINTR ) {
    printf("EINTR %d %d\n",semid,semnum);
  }
  if( istat == -1 ) {
    perror( "error in semaphore lock operation" );
    exit(1);
  }
}
