#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>

void read_saved_(char* one_obs, int* bytesRecv)
{
  static FILE* fin = NULL;
  char* ch;

  if( fin == NULL ) {
	  if ( access("rnx_recv",F_OK) != -1 ){
		  fin = fopen("rnx_recv","r");
	  } else {
		  printf("rnx_recv is not exist\n");
		  exit(EXIT_FAILURE);
	  }
  }
  fgets(one_obs,1024,fin);
  ch = one_obs;
  while( *ch != '\0' ) {
    if( *ch == '\n' ) {
      *ch = '\0';
      break;
    }
    ch++;
  }
  *bytesRecv = strlen(one_obs);
}
