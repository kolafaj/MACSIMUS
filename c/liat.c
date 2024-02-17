/* cc -O2 -o liat liat.c
 */
#include "../gen/include.h"

#define MAXLEN 64000
char line[MAXLEN];

int main(int narg,char **arg)
{
  int N,ibuf;
  char **buf;

  if (narg<2) {
    fprintf(stderr,"\
Print a stream but the last N lines. Call by:\n\
  liat N < INPUT > OUTPUT\n\
Output is empty if there are N or less lines on input (1 is returned)\n\
See also:\n\
  tail head line itail\n");
    exit(0); }

  getsbufsize=MAXLEN;

  N=atoi(arg[1]);
  
  if (N<=0) {
    while (gets(line)) puts(line); }
  else {
    allocarray(buf,N);
  
    loop (ibuf,0,N) {
      if (!gets(line)) return 1;
      buf[ibuf]=strdup(line); }
    
    ibuf=0;
    while (gets(line)) {
      puts(buf[ibuf%N]);
      free(buf[ibuf%N]);
      buf[ibuf%N]=strdup(line); 
      ibuf++; } }

  return 0;
}
