/* cc -o tab tab.c -lm
*/
#include "../gen/include.h"

int main(int narg,char **arg)
{
  double x,from,to,by=1,sg;
  char *fmt="%8g",*p;
  int intf=0,nfmt=0,n=0;

  p=getenv("FMT");
  if (p) fmt=p;

  if (narg<3) {
    fprintf(stderr,"\
Make a column of numbers FROM,FROM+BY,FROM*2*BY,...,TO (incl.):\n\
  tab FROM TO [BY [FORMAT]]\n\
Make a column of N equidistant real numbers from FROM to TO-BLOCK (incl.):\n\
  tab FROM TO [N:BLOCK [FORMAT]]\n\
Make a column of N (approx.) equidistant integers from FROM to TO-BLOCK+1\n\
(the block ends at TO incl.). Optional FORMAT must be integer (e.g., %%d):\n\
  tab FROM TO [-N:BLOCK [FORMAT]]\n\
Environment: FMT (overrides the default, will be chaged by FORMAT)\n\
Defaults: BY=1, FORMAT=%s.\n\
Max. 2 formats (%%) accepted, int-type (cxXodiu) or double-type (gGeEf)\n\
Example: Calculate sum 0.1+0.2+ ... + 99.9:\n\
  tab .1 99.9 .1 | sum\n\
Example: Call script test.sh on 10 files testA.dat .. testJ.dat, w.arg.:\n\
  tab 65 74 1 \"test.sh test%%c.dat %%.2f\" | sh\n\
See also:\n\
  mergetab tabproc plot\n",fmt);
    exit(0); }

  from=atof(arg[1]);
  to=atof(arg[2]);
  if (narg>3) {
    char *sem=strchr(arg[3],':');
    by=atof(arg[3]);

    if (sem) n=atoi(arg[3]),by=atof(sem+1);
    else by=atof(arg[3]);

    if (abs(n)==1) by=1; }

  if (narg>4) fmt=arg[4];
  if (by==0) return -1;
  if (by<0) sg=-1; else sg=1;

  if (n && by<0) Error("cannot use negative block with the N:BLOCK form");

  for (p=fmt;;) {
    p=strchr(p,'%');
    if (!p) break;

    if (p[1]=='%')
      p+=2;
    else {
      do {
	p++;
	if (!*p) Error("bad format");
        } while (!isalpha(*p));
      if (*p=='l') Error("l (long) specifier in format not supported");
      if (*p=='h') Error("h (short) specifier in format not supported");
      if (strchr("cxXodiu",*p)) intf|=1<<(2*nfmt);      // 1=int
      else if (strchr("gGeEf",*p)) intf|=1<<(2*nfmt+1); // 2=double
      else Error("unknown or unsupported format");
      nfmt++; } }
  if (!nfmt) Error("no valid format detected");
  if (nfmt>2) Error("more than 2 formats detected");

  if (n==0)
    to=(to+by*1e-9)*sg;
  else if (n>0) {
    to=to-by;
    by=(to-from)/(n-1);
    if (by<=0) Error("identical numbers"); }
  else {
    n=-n;
    to=to-by+1+1e-9;
    by=(to-from)/(n-1);
    if (by<=0) Error("identical numbers"); }

  if (n==1) by=(to-from)*2;

  for (x=from; x*sg<=to; x+=by) {
    switch (intf) {
      case 1: printf(fmt,(int)x); break;
      case 2: printf(fmt,x); break;
      case 4|1: printf(fmt,(int)x,(int)x); break;
      case 4|2: printf(fmt,x,(int)x); break;
      case 8|1: printf(fmt,(int)x,x); break;
      case 8|2: printf(fmt,x,x); break;
      default: Error("more than 2 formats or wrong format"); }
    printf("\n"); }

  return 0;
}
