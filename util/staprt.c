/* \make staprt

06/2009 -N changed
10/2008 bug fixed: .stx not recognized
*/

#include "ground.h"
#include "alloc.h"
#include "statics.h"

static char how[2048];
int isname=0,iarg,files=0,tcf=0,iiarg,prec=0;
int CHAR; /* 0: change / -> _ only ; 1: remove all problematic chars */

void howcpy(char *fn) /********************************************** howcpy */
{
  char *p=NULL,*c;

  if (strlen(fn)>1024) ERROR(("too long filename"))
  if (strlen(fn)<5) ERROR(("%s: too short filename",fn))

  strcpy(how,fn);
  for (c=how; *c; c++) if (*c=='.') p=c;
  if (!p) ERROR(("%s: no extension",fn))
  if (strcmp(p,".sta")) WARNING(("%s: not standard .sta extension",fn))
  *p=0;
}

void prefix(char *p) /*********************************************** prefix */
{
  memmove(how+strlen(p),how,strlen(how)+1);
  memcpy(how,p,strlen(p));
}

char *addext(char *fn) /********************************************* addext */
{
  char *c,*dot=NULL;

  for (c=fn; *c; c++) if (*c=='.') dot=c;

  if (dot && (!strcmp(dot,".sta") || !strcmp(dot,".stx") || !strcmp(dot,".sta~")))
    c=fn;
  else {
    alloc(c,strlen(fn)+5);
    sprintf(c,"%s.sta",fn); }
  return c;
}

void process(char *fn) /******************************************** process */
{
  char *c=addext(fn);

  StaLoad(c);

  printf("====== %s ======\n",fn);
  how[0]=0;
  if (files) howcpy(c);
  if (tcf && !files) prefix("-"); /* SHOULD NEVER HAPPEN */
  //  if (tcf==2) prefix('*');
  if (tcf==1) prefix(CHAR?"-T":"-t");
  if (tcf==2) prefix(CHAR?"-C":"-c");

  if (prec && !files && !tcf) how[0]='+';
}

char *key=NULL;
char *fmt=NULL;
char *ifmt=NULL;
double Qexp=-3e33;
double Q=1;

void one(char *fn,char *name) /***************************************** one */
{
  char *k,*c=addext(fn);
  double Q0;

  if (!fmt) fmt = prec ? " %18.12g" : " %.8g";
  if (!ifmt) ifmt=" %7d";

  StaLoad(c);

  Q0=0;
  for (c=fn; *c; c++) if (isdigit(*c) || *c=='.' && isdigit(c[1])) {
    Q0=atof(c);
    break; }
  if (Qexp>-2e33) Q=pow(Q0,Qexp);

  for (k=key; *k; k++)
    switch (*k) {
      case 'a': printf(fmt,Q*StaMean(name)); break;
      case 's': printf(fmt,Q*StaMean(name)*StaN(name)); break;
      case 'e': printf(fmt,Q*StaStdErr(name)); break;
      case 'r': printf(fmt,StaStdErr(name)/StaMean(name)); break;
      case 'd': printf(fmt,Q*sqrt(StaVar(name))); break;
      case 'v': printf(fmt,Q*StaVar(name)); break;
      case 'n': printf(ifmt,StaN(name)); break;
      case 'N': printf(" %s",name); break;
      case 'q': printf(fmt,Q0); break;
      case 'Q': printf(fmt,Q); break;
      case 'f': printf(" %s",fn); break;
      case 't': putchar('\t'); break;
      case 'T': printf(fmt,StaLag(name,'t')); break;  
      case 'L': printf(ifmt,(int)StaLag(name,'l')); break;  
      case 'D': printf(fmt,StaLag(name,'d')); break;  
      default: putchar(*k); }

  putchar('\n');
}

int main(int narg,char **arg) /**************************************** main */
{
  initscroll(0);

  if (narg<2) {
    fprintf(stderr,"\
Prints the contents of .sta file(s).\n\
  staprt [OPTIONS] FILE.sta [FILE.sta ...]\n\
OPTIONS:\n\
  -a     more digits in table and line output (cf. -f,-i)\n\
  -c0    print standard statistics [default]\n\
  -c1    print time correlation functions (to FILE.NAME.tcf)\n\
  -c2    print time covariances (to FILE.NAME.cov)\n\
  -C     remove non-alphanum chars in NAME (default=/ changed to _ only)\n\
  -iFMT  number format for line output (KEY=n) [df.=-i\" %%7d\"], also -I\n\
  -k     (w/o  -n) only 1 line summary statistics for each recorded variable\n\
  -kKEY  (with -n) only 1 line of statistics, KEY = string of:\n\
         a=average\n\
         d=standard deviation of the data = sqrt(variance)\n\
         D=DT (time unit)\n\
         e=standard error of the average\n\
         f=file name\n\
         L=lag (integer)\n\
         n=number of data points\n\
         N=NAME (value of option -n or -N)\n\
         q=the value of Q (as set by options -q -Q -* -+ -* -/)\n\
         Q=the 1st number found in the filename\n\
         r=relative standard error of the average\n\
         s=sum of the data=a*n\n\
         t=tabulator\n\
         T=lag (in time units)\n\
         v=variance of the data\n\
         other=print as is\n\
  -fFMT  number format for line output (KEY=aed) [df.=-f\" %%.8g\"], also -F\n\
  -mMTD  0=conservative weighed std.error estimate, good for short series\n\
           blocked table based on tau=1+2*partial sum of c_k\n\
         1,2=based on autoregressive process of the 1st order=exp.decay of c_k\n\
           block-block (k=1) uses subaverage process formula\n\
           block-..-block (k>1) uses exponential extrapolation\n\
         1=final value based on the last blocked item w. c_1 (subaverage)\n\
         2=final value based on the last blocked item w. c_2 (exp.extrap.)\n\
  -nNAME print statistics for given name (this option may repeat)\n\
  -NNAME the same as -nNAME -kN:aef\n\
         NAME may contain wildcard * at end, but unique match is required\n\
  -p     the same as -a\n\
  -qQ    set scaling factor for line output (KEY=aedv, not .cov,.tcf) [df.=1]\n\
  -*# -/# --# -+# -^# modify the above factor Q\n\
  -Q#    override any -qQ -*# -/# --# -+# -^# by:\n\
         Q:=(1st number found in the filename)^#\n\
  -t#    the same as -c#\n\
Return value:\n\
  0 on success, 1 on error (incl. name not found)\n\
FILE:\n\
  if NAME is given, extension may be omitted; then .sta is assumed\n\
Examples:\n\
  staprt -N\"P [Pa]\" water256.sta water512.sta\n\
  staprt -n\"P*\" -kaenf -q1e-6 -F\" %%7.2f\" *.sta\n\
See also: autocorr showcp\n\
");
    exit(0); }

  loop (iarg,1,narg) if (arg[iarg][0]=='-') switch (arg[iarg][1]) {
    case 'p':
    case 'a': prec++; break;
    case 'f':
    case 'F': fmt=arg[iarg]+2; break;
    case 'i':
    case 'I': ifmt=arg[iarg]+2; break;
    case 'k': key=arg[iarg]+2; break;
    case 'm': StaErrMethod(atoi(arg[iarg]+2)); break;
    case 'N': key="N:aef";
    case 'n': isname++; break;
    case 'C': CHAR=1; break;
    case 'c':
    case 't': tcf=atoi(arg[iarg]+2); if (!tcf) tcf++; files++; break;
    case 'q': Q=atof(arg[iarg]+2); break;
    case 'Q': Qexp=atof(arg[iarg]+2); break;
    case '*': Q*=atof(arg[iarg]+2); break;
    case '/': Q/=atof(arg[iarg]+2); break;
    case '+': Q+=atof(arg[iarg]+2); break;
    case '-': Q-=atof(arg[iarg]+2); break;
    case '^': Q=pow(Q,atof(arg[iarg]+2)); break;
    default: ERROR(("%s is bad option",arg[iarg])); }

  if (isname) {
    loop (iarg,1,narg) if (arg[iarg][0]=='-' && tolower(arg[iarg][1])=='n')
      loop (iiarg,1,narg) if (arg[iiarg][0]!='-') {
	if (key)
	  one(arg[iiarg],arg[iarg]+2);
	else {
	  process(arg[iiarg]);
	  StaPrint(arg[iarg]+2,how[0]?how:""); }
	StaFree(); } }
  else
    loop (iiarg,1,narg) if (arg[iiarg][0]!='-') {
      process(arg[iiarg]);
      if (key) StaPrintAll(NULL);
      else StaPrintAll(how[0]?how:"");
      StaFree(); }

  return StaError;
}
