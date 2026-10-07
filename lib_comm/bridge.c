#include <sys/types.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

void bridge_(char* portno, int* sock)
{
  int                   port;
  struct sockaddr_in    server;
  struct hostent       *hp;

   *sock = socket(AF_INET, SOCK_STREAM, 0);
  if (*sock < 0) {
    perror("error opening stream socket!");
    exit(1);
  }
  server.sin_family = AF_INET;
  hp = gethostbyname("localhost");
  printf("%s is a known host! address = %d.%d.%d.%d = \n", hp->h_name,hp->h_addr_list[0][0],
                         hp->h_addr_list[0][1],hp->h_addr_list[0][2],hp->h_addr_list[0][3]);
  port = atoi(portno);

  bcopy(hp->h_addr, &server.sin_addr, hp->h_length);

  server.sin_port = htons(port);

  printf("Trying to connect...\n");
  if(connect(*sock, (struct sockaddr*)&server, sizeof(struct sockaddr_in)) < 0) {
     perror("connecting stream socket");
     close(*sock);
     exit(1);
  } else {
    printf("Connection succeeds!\n");
  }
}
