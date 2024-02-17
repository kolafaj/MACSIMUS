/* make lr
 */
#include "ground.h"
#include "linregr.h"

int main(int narg,char **arg)
{
  double x,y,s=1,DX=0;
  char line[1024];
  unsigned n=0;

  initscroll(0);
  if (narg>1) {
  fprintf(stderr,"\
Linear regression (filter). Call by:\n\
  [LRDX=] lr < DATA\n\
Format of DATA:\n\
  x  y  [err] # LRDX not defined or not number or 0\n\
  y  [err]    # x=n*LRDX\n\
If err is given, weight of 1/err^2 is used (but err is NOT interpreted as\n\
the standard error, i.e., only weights are used).\n\
# and ! in 1st column is comment\n\
See also:\n\
  lineareq lincomb linfit\n");
  exit(0); }

  if (getenv("LRDX")) DX=atof(getenv("LRDX")); 
  
  while (gets(line))
    if (!strchr("!#",line[0])) {
      int err;
      
      if (DX) {
        err=sscanf(line,"%lf%lf",&y,&s)+1;
        x=n*DX; }
      else
        err=sscanf(line,"%lf%lf%lf",&x,&y,&s);
      
      if (s==0 || err==2) s=1;
      if (err>=2) LRAdd("stdin",1/Sqr(s),x,y); 
      n++; }

  LRPrint(' ');

  return 0;
}
