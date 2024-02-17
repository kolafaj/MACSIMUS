/* cc -O2 -o sum sum.c
 */
#include "../gen/include.h"

#define M 202

double sum[M],x[M];
int maxn=0;

int main(int narg, char **arg)
{
  char line[M*16];
  char *fmt="%.8g",*c;
  int n,i;

  if (narg>1 && !strcmp(arg[1],"-h")) {
    fprintf(stderr,"\
Sum a table of data by columns. Call by:\n\
  sum [FORMAT] < DATA > LINE_OF_RESULTS\n\
FORMAT = C-style format, default=\"%s\"\n\
See also:\n\
  sumetc\n",fmt);
    exit(0); }

  if (narg>1) fmt=arg[1];

  c=strchr(fmt,'%');
  if (!c) Error("no %% in format, sum -h for help");
  else if (strchr(c+1,'%')) Error("multiple %% in format, sum -h for help");

  getsbufsize=M*16;
  
  while (gets(line))
    if (!strchr("!#",line[0])) {
      char *tok=strtok(line,"\t \r,");
     
      for (n=0;;n++) {
        if (n>=M) Error("sum: too many columns, increase #define M and recompile");
        if (!tok || sscanf(tok,"%lf",&x[n])!=1) break;
        tok=strtok(NULL,"\t \r,"); }
      loop (i,0,n) sum[i]+=x[i];
      Max(maxn,n) }

  loop (i,0,maxn) {
    printf(fmt,sum[i]);
    putchar(" \n"[i==maxn-1]); } 

  return 0;
}
