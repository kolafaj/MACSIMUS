/* cc -O2 -o runsum runsum.c
 */
#include "../gen/include.h"

#define M 64

double sum[M],x[M];
int maxn=0;

int main(int narg, char **arg)
{
  char line[M*16];
  char *c;
  char *FMT=getenv("FMT");
  int n,i;
  int nlines=0;
  int av=0;
  int from=0;

  if (!FMT || !strchr(FMT,'%')) FMT="%.8g";

  if (narg<2) fprintf(stderr,"[FMT=%%C-format] [FROM=column] runsum [-] < DATA > RUNNING SUMS[-AVERAGES]\n");

  if (getenv("FROM")) from=atoi(getenv("FROM"))-1;
  if (from<0) from=0;

  if (narg>1) av++;

  getsbufsize=M*32;

  while (gets(line)) {
    if (!strchr("!#",line[0])) {
      char *tok=strtok(line,"\t \r,");
     
      for (n=0;;n++) {
	if (n>=M) Error("too many columns");
	if (!tok || sscanf(tok,"%lf",&x[n])!=1) break;
	tok=strtok(NULL,"\t \r,"); }
      loop (i,0,from) sum[i]=x[i];
      loop (i,from,n) sum[i]+=x[i];
      Max(maxn,n)
      nlines++; }

    loop (i,0,maxn) {
      if (av) printf(FMT,sum[i]/nlines);
      else printf(FMT,sum[i]);
      putchar(" \n"[i==maxn-1]); } }
      
  return 0;    
}
