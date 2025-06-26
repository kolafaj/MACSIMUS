/* cc -O2 -o sumetc sumetc.c -lm
09/2021: FMT added
03/2011: first added
 */
#include "../gen/include.h"
#include "../gen/powi.c"

#define M 204

double sum[M],sumq[M],sumg[M],x[M],first[M];
/* note: first[M] = the first value, improves precision of Var of data with little relative fluctuation */
int N[M];
struct opt_s {
  struct opt_s *next;
  enum op_e { ISUM, XSUM, RSUM, NSUM, PROD, MIN, MAX, AV, GEOMAV, STDEV, ABSSTDEV, STDERR, VAR, ABSVAR, EQ,NE,LT,LE,GT,GE } op;
  double power; /* or relation */
  int ipower;
  char *fmt;
  double sum[M]; /* sum or other operation result */
  unsigned i[M];   /* for min/Max, the line # */
} *head,*last;
int maxn=0,printnline=0;
unsigned nline=0;

int main(int narg, char **arg) /*************************************** main */
{
  char line[M*16];
  char *FMT,*fn="-",*c;
  int n,i,iarg,sg=1;
  FILE *f;
  struct opt_s *one;
  char *SEP=" ";

  FMT=getenv("FMT");
  if (!FMT || !strchr(FMT,'%')) FMT="%.8g";

  if (narg<2) {
    fprintf(stderr,"\
Sum and other operations by columns of data. Call by:\n\
  sumetc [OPTIONS] [FILE|-]\n\
FILE = file of max. %d columns\n\
  Missing FILE or - = filter (input=stdin)\n\
  Input lines with non-numeric data are ignored beyond the non-numeric item\n\
OPTIONS:\n\
  -#       = sum of data^# (-2 = sum of squares)\n\
  -a       = arithmetic average, sum X_i/n\n\
  -d       = sample standard deviation of the data (aka Bessel corrected),\n\
             sqrt[sum(X_i-<X>)^2/(n-1)]\n\
  -D       = standard deviation of the data set, sqrt[sum(X_i-<X>)^2/n]\n\
  -e       = standard error of the average (assuming uncorrelated data)\n\
  -fFORMAT = C-style format for double, -f.9f=-f%%.9f [default=%s]\n\
  -g       = geometric average, (prod X_i)^(1/n)\n\
  -gt# -ge# -lt# -le# -ne# -eq#\n\
           = # of data > >= < <= != == given number #\n\
  -i       = print also the line number (1st instance) of min/max (-m,-M) in ()\n\
  -m       = minimum\n\
  -M       = maximum\n\
  -n       = number of data (the same as -0)\n\
  -p       = product\n\
  -s       = sum\n\
  -SSEP    = separator between outputs from 1 input column [1 space]\n\
  -x       = alternating sum (1st-2nd+3rd...)\n\
  -v       = sample variance (aka Bessel corrected), sum(X_i-<X>)^2/(n-1)\n\
  -V       = variance of the data set, sum(X_i-<X>)^2/n\n\
Environment:\n\
  FMT = format accepting double (will be overriden by -f)\n\
One line with results is printed. Options may repeat; then the results for\n\
the fist column are printed first, then for the second, etc.\n\
Examples\n\
  tab 1 10 | sumetc -f.2f --0.5                   # sum_i=1^10 1/sqrt(i)\n\
  tab 1 1000000 | tabproc \"rnd(0)\" | sumetc -a -e # check unformity of rnd\n\
See also:\n\
  tabproc mergetab stat difxmin norm sum filttab sumerr autocorr\n\
  sumfiles sumline avfiles avfiles2 avfiles3 difxmin\n\
",M,FMT);
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-') {
      if (arg[iarg][1]==0) fn="-";
      else if (arg[iarg][1]=='S')
        SEP=arg[iarg]+2;
      else if (arg[iarg][1]=='f') {
        FMT=arg[iarg]+2;
        c=strchr(FMT,'%');
        if (!c) {
          FMT=strdup(arg[iarg]+1);
          FMT[0]='%'; }
        else
          if (strchr(c+1,'%')) Error("sumetc: multiple %% in format"); }
      else if (arg[iarg][1]=='i')
        printnline++;
      else {
        alloconezero(one);
        if (head) last->next=one;
        else head=one;
        last=one;
        one->fmt=FMT;

        one->power=atof(arg[iarg]+3);

        switch (arg[iarg][1]) {
          case 's':
	    one->power=1;
	    one->op=ISUM;
	    break;
          case 'x':
            one->power=1;
            one->op=XSUM;
            break;
          case 'p':
	    one->op=PROD;
	    loop (i,0,M) one->sum[i]=1;
	    break;
          case 'e':
            one->op=arg[iarg][2]=='q'?EQ:STDERR; break;
          case 'd': one->op=STDEV; break;
          case 'D': one->op=ABSSTDEV; break;
          case 'v': one->op=VAR; break;
          case 'V': one->op=ABSVAR; break;
          case 'a': one->op=AV; break;
          case 'g':
            one->op=arg[iarg][2]=='e'?GE:arg[iarg][2]=='t'?GT:GEOMAV; break;
          case 'l':
            one->op=arg[iarg][2]=='e'?LE:LT; break;

          case 'M':
	    one->op=MAX;
	    loop (i,0,M) one->sum[i]=-9e99;
	    break;
          case 'm':
	    one->op=MIN;
	    loop (i,0,M) one->sum[i]=9e99;
	    break;
          case 'n':
            one->op=arg[iarg][2]=='e'?NE:NSUM; break;
          default: /* digit */
            if (strchr("0123456789.-+",arg[iarg][1])) {
              one->op=ISUM;
              one->power=atof(arg[iarg]+1); }
            else Error("sumetc: bad option"); }
        one->ipower=one->power;

        if (one->op==ISUM && one->ipower!=one->power) one->op=RSUM; } }
    else fn=arg[iarg];

  if (strcmp(fn,"-")) f=fopen(fn,"rt");
  else f=stdin;
  if (!f) Error("sumetc: sumetc: cannot open file");

  while (fgets(line,M*16,f))
    if (!strchr("!#",line[0])) {
      char *tok=strtok(line,"\t \r,");

      for (n=0;;n++) {
        if (n>=M) Error("sumetc: sumetc: too many columns");
        if (!tok || sscanf(tok,"%lf",&x[n])!=1) break;
        tok=strtok(NULL,"\t \r,"); }

      loop (i,0,n) {
        if (N[i]==0) first[i]=x[i];
        N[i]++;
        sum[i]+=x[i]-first[i];
        sumg[i]+=log(x[i]);
        sumq[i]+=Sqr(x[i]-first[i]); }

      looplist (one,head) switch (one->op) {
        case ISUM: loop (i,0,n) one->sum[i]+=powi(x[i],one->ipower); break;
        case XSUM: loop (i,0,n) one->sum[i]+=sg*powi(x[i],one->ipower); sg*=-1; break;
        case NSUM: loop (i,0,n) one->sum[i]+=1; break;
        case RSUM: loop (i,0,n) one->sum[i]+=pow(x[i],one->power); break;
        case PROD: loop (i,0,n) one->sum[i]*=x[i]; break;
        case MAX: loop (i,0,n) if (x[i]>one->sum[i]) one->sum[i]=x[i],one->i[i]=nline; break;
        case MIN: loop (i,0,n) if (x[i]<one->sum[i]) one->sum[i]=x[i],one->i[i]=nline; break;
        case EQ: loop (i,0,n) one->sum[i]+=(x[i]==one->power); break;
        case NE: loop (i,0,n) one->sum[i]+=(x[i]!=one->power); break;
        case GT: loop (i,0,n) one->sum[i]+=(x[i]> one->power); break;
        case GE: loop (i,0,n) one->sum[i]+=(x[i]>=one->power); break;
        case LT: loop (i,0,n) one->sum[i]+=(x[i]< one->power); break;
        case LE: loop (i,0,n) one->sum[i]+=(x[i]<=one->power); break;
        default: ; /* suppress compiler warning */
      }
      nline++;
      Max(maxn,n) }

  loop (i,0,maxn)
    looplist (one,head) {
      switch (one->op) {
        case AV: printf(one->fmt,sum[i]/N[i]+first[i]); break;
        case GEOMAV: printf(one->fmt,exp(sumg[i]/N[i])); break;
        case STDERR: printf(one->fmt,sqrt((sumq[i]/N[i]-Sqr(sum[i]/N[i]))/(N[i]-1))); break;
        case STDEV: printf(one->fmt,sqrt((sumq[i]-Sqr(sum[i])/N[i])/(N[i]-1))); break;
        case ABSSTDEV: printf(one->fmt,sqrt((sumq[i]-Sqr(sum[i])/N[i])/(N[i]))); break;
        case VAR: printf(one->fmt,(sumq[i]-Sqr(sum[i])/N[i])/(N[i]-1)); break;
        case ABSVAR: printf(one->fmt,(sumq[i]-Sqr(sum[i])/N[i])/(N[i])); break;
        case MAX:
        case MIN: printf(one->fmt,one->sum[i]);
          if (printnline) printf("(%u)",one->i[i]);
          break;
        default: printf(one->fmt,one->sum[i]); }

      if (one->next)
        putchar(' ');
      else {
        if (i==maxn-1) putchar('\n');
        else printf("%s",SEP); } }

  return 0;
}
