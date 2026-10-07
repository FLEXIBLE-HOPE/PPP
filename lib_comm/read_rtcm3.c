#include <sys/types.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

void bridge_(char* portno, int* sock);

void read_rtcm3_(int* sock, char* one_obs, char* portno, int* bytesRecv)
{
  char* ch;
  static char buff[1024];


  while(1) {
    memset(buff,'\0',sizeof(buff));
    *bytesRecv = recv(*sock, &buff, sizeof(buff)-1, 0);
    if(*bytesRecv < 0) {
      perror("Reading rtcm3 observations problems! ");
      close(*sock);
      exit(1);
    } else if(*bytesRecv == 0) {
      close(*sock);
      sleep(0.1);
      bridge_(portno, sock);
      continue;
    } else {
      buff[*bytesRecv]='\0';
      ch = buff;
      while(ch != NULL) {
        ch = strstr(ch,"\n");
        if(ch != NULL) *ch++ = '%';
      }
      strncpy(one_obs,buff,*bytesRecv);
      break;
    }
  }
}
