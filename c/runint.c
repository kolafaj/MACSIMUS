/* cc -O2 -o runint runint.c -lm
 */
#include "../gen/include.h"

#define M 1024

double sum[M],subsum[M],x[M],last[M];
char line[M*32];
int n;

/* see Ralston, p.141, and http://mathworld.wolfram.com/Newton-CotesFormulas.html */
double coeff[8][9]={
  {0.5,0.5}, // trapezoid
  {1./3,4./3,1./3}, // Simpson
  {3./8,9./8,9./8,3./8}, // Simson 3/8
  {14./45,64./45,24./45,64./45,14./45}, // Boole
  {95./288,375./288,250./288,250./288,375./288,95./288},
  {41./140,216./140,27./140,272./140,27./140,216./140,41./140},
  {5257./17280,25039./17280,9261./17280,20923./17280,20923./17280,9261./17280,25039./17280,5257./17280},
  {3956./14175,23552./14175,-3712./14175,41984./14175,-18160./14175,41984./14175,-3712./14175,23552./14175,3956./14175} };

int nextline(void)
/* reads x[], returns number of columns (=n) on success */
{
  char *tok;
  static char lastn=0;

  getsbufsize=M*32;
  
  do {
    if (!gets(line)) return 0;
  } while (strchr("!#",line[0]));

  tok=strtok(line,"\t \r,");

  for (n=0;;n++) {
    if (n>=M) Error("runint: too many columns");
    if (!tok || sscanf(tok,"%lf",&x[n])!=1) break;
    tok=strtok(NULL,"\t \r,"); }

  if (n==0) return 0;

  if (lastn && n!=lastn) n=lastn;
  lastn=n;

  return n;
}

int main(int narg, char **arg)
{
  char *FMT=getenv("FMT");
  double h=0,hlast=0,hnew,xlast=0,Q=1;
  int i,evenfunc=0,from=1;
  int method=1,phase,pass=0,err=0;
  int badstep=0,badsteps=0;

  if (!FMT || !strchr(FMT,'%')) FMT="%.8g";
  
  if (narg<2) {
    fprintf(stderr,"\
Newton-Cotes running integral from equidistant (or division-point matching)\n\
data. Call by:\n\
  runint OPTIONS < DATA\n\
DATA are lines of (if OPTION -h, [x] is missing)\n\
  [x] y1 [y2 ..]\n\
  [x] y1 stderr1 [y2 stderr2..] # if option -e\n\
  #- and !-lines are omitted on input\n\
OPTIONS:\n\
  -mMETHOD: number of subintervals, particularly:\n\
    -m1  2-point, 2nd order = trapezoid (default)\n\
    -m2  3-point, 4th order = Simpson\n\
    -m3  4-point, 4th order = Simpson 3/8, not recommended\n\
    -m4  5-point, 6th order = Boole\n\
    -m5  6-point, 6th order, not recommended\n\
    -m6  7-point, 8th order\n\
    -m7  8-point, 8th order, not recommended\n\
    -m8  9-point, 10th order, bad stability (negative coefficients)\n\
  -o the first interval is halved, METHOD must be even number\n\
     use if y(x-x0)=y(-x-x0) (x0-centered even function)\n\
  -e data in pairs with statistical errors\n\
  -q# multiply the results by real number #\n\
  -fFORMAT: C-style format (default=$FMT)\n\
  -hSTEP: there is no x column, x's are n*STEP (starts with n=0)\n\
Environment:\n\
  FMT=C-style format (now or default FMT=%s)\n\
Example:\n\
  runint -o -m2 -h.01 < spce256.Ptxz.cov | plot -\n\
See also:\n\
  runsum (both running sum and average)\n\
NB:\n\
  tab 0 3 | runsum         # -> 0 1 3 6\n\
  tab 0 3 | runint -h1 -m1 # -> 0.5 2 4.5\n\
",FMT);
    exit(0); }

  loop (i,1,narg) switch (tolower(arg[i][1])) {
    case 'h': from=0; h=atof(arg[i]+2); break;
    case 'o': evenfunc++; break;
    case 'e': err++; break;
    case 'm': method=atoi(arg[i]+2); break;
    case 'q': Q=atof(arg[i]+2); break;
    case 'f': FMT=arg[i]+2;
      if (!strchr(FMT,'%')) FMT=string("%%%s",FMT);
      break;
    default: Error("runint: bad option"); }

  if (method<1 || method>8) Error("runint: METHOD must be in interval [1,8]");

  //  if (evenfunc && err) Error("implementation limitation: cannot combine -e and -o");

  if (evenfunc && method&1) Error("even function (-o) requires even method (-m)");

  if (evenfunc) evenfunc=method/2;

  if (!nextline()) Error("no data");

  for (pass=0;;pass=1) {
    /* loop by (method+1)-chunks of data; 1st [0] assumed to have been read */

    loop (i,0,n) subsum[i]=0;
    if (from) {
      h=x[0];
      hlast=-3e33;
      badstep=0;
      xlast=x[0]; }

    loopto (phase,evenfunc,method) {
      if (phase!=evenfunc) {
        /* read new line and check step */
        if (!nextline()) goto TheEnd;
        if (from) {
          hnew=x[0]-xlast;
          xlast=x[0];
          if (hlast>-2e-33)
            if (hnew==0 || fabs((hlast-hnew)/(hlast+hnew))>1e-6) badstep++;
          hlast=hnew; }
      }
      //      printf("%g:%g\n",x[0],coeff[method-1][phase]);
      if (err) {
        loop (i,from,n)
          if ((i-from)&1) {
            subsum[i] += (pass&&phase==0?3:1)*Sqr(coeff[method-1][phase]*x[i]); }
          else
            subsum[i] += coeff[method-1][phase]*x[i];
      }
      else {
        if (evenfunc && phase==evenfunc)
          loop (i,from,n) subsum[i] += coeff[method-1][phase]*x[i]/2;
        else
          loop (i,from,n) subsum[i] += coeff[method-1][phase]*x[i]; } }

    if (from) {
      h=(x[0]-h)/method;
      if (evenfunc) h*=2;
      sum[0]=x[0]; }
    evenfunc=0;
    if (err) {
      if (from) { printf(FMT,sum[0]); putchar(' '); }
      for (i=from; i<n; i+=2) {
        sum[i] += h*subsum[i];
        printf(FMT,sum[i]);
        putchar(' ');
        sum[i+1] += h*h*subsum[i+1];
        printf(FMT,Q*sqrt(sum[i+1]));
        if (i<n-1) putchar(' '); }
      if (badstep) {
        badsteps+=badstep;
        printf(" nonequidistant steps detected!"); }
      putchar('\n'); }
    else {
      loop (i,from,n) sum[i] += h*subsum[i];
      loop (i,0,n) {
        if (i) printf(FMT,Q*sum[i]);
        else printf(FMT,sum[i]);
        if (i<n-1) putchar(' '); }
      if (badstep) {
        badsteps+=badstep;
        printf(" nonequidistant step(s)!"); }
      putchar('\n'); }
    badstep=0;
  } /* pass */


 TheEnd:
  if (badsteps) 
    fprintf(stderr,"WARNING %d nonequidistant step(s) detected (rel. threshold=1e-6)\n",
            badsteps);

  return 0;
}
