/* make autocorr
*/

#include "ground.h"
#include "varfile.h"
#include "statics.h"
#include "options.h"


#define x badoption
int optionlist[32] =
  /* ` a b c d e f g h i j k  l m  n o p q r s t u v w x y z { | } ~   */
    {x,x,2,0,x,1,0,x,0,x,x,x,29,0,15,0,0,0,x,0,0,0,x,x,x,x,x,x,x,x,x,x};
  /* @ A B C D E F G H I J K  L M  N O P Q R S T U V W X Y Z [ \ ] ^ _ */
#undef x

#define NLINE 32000

int main(int narg, char **arg) /*************************************** main */
{
  char line[NLINE];
  double x;
  int ln=0,i;
  char *fn="STDIN";
  char *tcffn;
  double M,m=0;
  int sum=0;
  double sumx=0,sumxx=0;
  double DT=1;
  char *STA=NULL;
  
  initscroll(0);
  
  if (narg<2) {
   fprintf(stderr,"\
Autocorrelation and error analysis of time-series. Call by:\n\
  autocorr {FILE|-} [OPTIONS]\n\
OPTIONS [default value]:\n\
  -b#   lag for calculations with data blocked by 2, 4, 8, ... items [2]\n\
  -c#   column of data [1]\n\
  -e#   missing field is: 1=error, 0=warning+skip, -1:silently ignored [1]\n\
  -f-1  follow trajectory: input data are periodic mod 2*PI\n\
  -f#   follow trajectory: input data are mod # (#>0 is integer)\n\
  -h#   set time step (often h*noint of cook) [1]\n\
  -l#   lag (# of autocorrelation coefficients) for statistics [29]\n\
  -l0   suppress autocorrelation calculations\n\
  -m#   method for stderr of corelated data\n\
  -n#   number of blocked calculations: the maximum block size is 2^# [12]\n\
  -o    one line of statistics instead of compact table\n\
  -p    higher precision (more digits) on output [off]\n\
  -sFILE.sta: create sta-file FILE.sta\n\
  -t1   write file of time correlation function (tcf(t=0)=1) [off]\n\
        (note that 1+2tau is w.r.t h=1, to be used in stderr calc.)\n\
  -t2   write file of time covariances (Cov(t=0)=Var) [off]\n\
  -t3   -t1+-t2\n\
  -u    uncorrelated data: average, stdev of aver., rel.err (most opts.ignored)\n\
  - instead of FILE means stdin; with -t also output to stdout\n\
See also:\n\
  tcov lr avfiles statfile sumfiles sumetc showcp runsum runint staprt spectrum\n"
);
    exit(0); }

  loop (i,1,narg)
    if (arg[i][0]=='-') 
      if (arg[i][1]) {
        getoption(arg[i])
        if (arg[i][1]=='h') DT=atof(arg[i]+2);
        if (arg[i][1]=='s') STA=arg[i]+2; }
      else
        in=stdin;
    else {
      in=fopen(fn=arg[i],"rt");
      if (!in) ERROR(("cannot open %s",arg[i])) }

  StaErrMethod(option('m'));
  
  StaSet(DT,option('l'),option('b'),option('n'));

  M=option('f'); if (M<0) M=2*PI;
  if (option('f')) fprintf(stderr,"WARNING: data modulo %g reconstructed\n",M);

  alloc(tcffn,strlen(fn)+8);
#ifndef DOS
  if (option('c'))
    sprintf(tcffn,"%s.%d",fn,option('c'));
  else
#endif
    strcpy(tcffn,fn);
  if (!option('c')) option('c')=1;

  while (fgets(line,NLINE,in)) {
    ln++;
    if (!strchr("!#",line[0])) {
      char *t=strtok(line," \t\n");

      loop (i,0,option('c')-1) t=strtok(NULL," \t\n");
      if (!t) {
        if (option('e')==-1) break;
        fprintf(stderr,"autocorr: line %d too short\n",ln);
        if (option('e')==1) exit(1);        
        break; }
      x=atof(t);
      if (option('f')) {
        while (x-m>=M/2) x-=M;
        while (x-m<-M/2) x+=M;
        m=x; }
      if (option('u')) { sum++; sumx+=x; sumxx+=x*x; }
      if (option('l')) StaAdd(tcffn,x); } }

  if (option('u')) {
    double av=sumx/sum;
    double std=sqrt((sumxx/sum-Sqr(av))/(sum-1));
    
    printf("%g %g %g %d av stderr rel.stderr n (uncorrelated)\n",av,std,std/av,sum); }

  if (option('t')&1) {
    if (!strcmp(fn,"STDIN")) StaPrint(tcffn,"-");
    else StaPrint(tcffn,"-@"); /* used to be "@" */ }

  if (option('t')&2) {
    if (!strcmp(fn,"STDIN")) StaPrint(tcffn,"-c-");
    else StaPrint(tcffn,"-c"); }

  if (!option('t') && !option('o') && option('l')) 
    StaPrint(tcffn,option('p')?"+":"");
  
  if (option('o') && option('l'))
    StaPrint(tcffn,NULL);

  if (STA) StaSave(STA);

  return 0;
}
