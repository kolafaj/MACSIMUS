/* make spectrum
  07/2021 -q fixed, -g added
  05/2021 extended: -I,-s,-C
 */
#include "ground.h"
#include "fft.h"
/* float if -DFLOAT, otherwise double
   CAVEAT: not conforming to -DPRECISION */
typedef fftcomplex complex;
typedef fftreal real;
#include "solve2.h"

int n,win,prec=8;

double window(double x) /******************************************** window */
{
  x=sin(x/2);
  return powi(x,win);
}

double eps=1e-7,gpass;
int npass;
double A,*ypass;

double f(double a,double s) /********************************************* f */
{
  int i;
  double e,p,ee=0,ye=0,eepps=0,yepps=0,eepss=0,yepss=0;

  loopto (i,0,npass) {
    p=i-a;
    e=exp(-Sqr(p*s));
    ee+=e*e;
    ye+=ypass[i]*e;
    eepps+=e*e*Sqr(p*p)*s;
    yepps+=ypass[i]*e*Sqr(p*p)*s;
    eepss+=e*e*p*Sqr(s);
    yepss+=ypass[i]*e*p*Sqr(s); }

  A=ye/ee;
  gpass=A*eepss-yepss;

  return A*eepps-yepps;
}

double g(double a,double s) /********************************************* g */
{
  return gpass;
}

double peak(double *y,int n) /***************************************** peak */
{
  double a=n/2,s=0.5;
  int err;

  npass=n;
  ypass=y;

  err=Solve2(f,g,
             &a,0,n,eps,
             &s,0,2,eps,-30,-2);

  if (err) return 0;
  else return a;
}

int isq=0,isverlet=0,iscm=1;
double h=0;

double conv(double f) /************************************************ conv */
{
  if (isq) f/=2;
  
  if (isverlet)
    f=sqrt(2-2*cos(2*PI*h/isverlet*f))/(2*PI*h/isverlet);

  if (iscm) f/=0.0299792458;

  return f;
}

int main(int narg,char **arg) /**************************************** main */
{
  complex *y,sum={0,0};
  real *sp;
  FILE *fpeak=NULL;
  char *fn=NULL,*Im=NULL,*F=NULL,*fmt;
  int iarg,ic,i,k=0,ik,col=0,Icol=-1,omit=0,plot=0,ndata=0,missing=0,isstep=0;
  int maxcol,Coutput=0,REfill=0,IMfill=0,copyin=0,direction=1,paddbyaverage=0;
  double shift=0,slowdown=3,minpeak=-1;

  if ( (fmt=getenv("FMT")) && fmt[0])
    fprintf(stderr,"WARNING: environment FMT specified but ignored, use option -d instead\n");
  
  if (narg<2) {
    fprintf(stderr,"\
MACSIMUS spectrum (discrete Fourier transform of real or complex data). Call by:\n\
  spectrum [OPTIONS] {FILE|-}\n\
FILES:\n\
  FILE = input file of white-separated data; empty, #-, !-lines are ignored\n\
  FILE.sp = output spectrum\n\
    * without -h: unnormalized by lines (see -C for shift):\n\
               n-1\n\
      A[k] = h SUM a[j]*exp[2*pi*i*(j+shift)*k/n], k=0..n-1\n\
               j=0\n\
    * with -h: 1st column = wavenumbers (w/o -T) or frequencies f (with -T),\n\
      conforming to formula:\n\
      A(f) = int a(t)*exp(2*pi*i*(t+tshift)*f) dt, tshift=shift/n*h\n\
  FILE.F.sp = as above with F given by option -f\n\
  FILE.pk = peaks analyzed (with -m) from sqrt(Re^2+Im^2)\n\
  FILE.F.pk = as above with F given by option -f\n\
  -: input=stdin, both outputs are concatenated to stdout\n\
OPTIONS:\n\
  -a   missing data on input are padded by the average [default=zeros]\n\
  -c[SG]RE[,[SG]IM] = column(s) of input data (default: Re=col 1, Im=zeros):\n\
       SG=+: filled to be an even function, 1st line -> t=h*shift\n\
       SG=-: filled to be an odd function, 1st line -> t=h*shift\n\
       Implies (-k,k) output. Both SGs must be defined, normally +- or -+.\n\
  -C#  complex output [default = sqrt(Re^2+Im^2)]\n\
       #=cyclic shift by # lines on input (line 0 becomes #), can be real\n\
         number (implemented after FFT by *exp[i*shift*k/n])\n\
       HACK for fast decaying functions for t-> +- infinity, with -k-\n\
         and half-shift: spectrum for k<0 is multiplied by -1\n\
       With even/odd filling (-cSG..), use:\n\
         -C or -C0: 1st line is for t=0 (velocity Verlet/trapezoid style)\n\
         -C0.5: 1st line is for t=h/2 (leap-frog/rectangle style)\n\
       NB: Re^2+Im^2 does not depend on the shift\n\
  -d#  decimal digits for output g-format [default=8] (NB: $FMT not used)\n\
  -e#  efficiency warning limit (estimated speed wrt n=2^m) [default=3]\n\
  -f   F in file names is column number (see -c)\n\
  -f@  @=F in file names (string)\n\
  -h#  timestep h [default=1]; usually interpreted in ps\n\
  -g#  grid (the same as h=1/grid)\n\
  -i   write a file (FILE.in or FILE.F.in) of preprocessed input data\n\
  -I   inverse Fourier transform (interpretation reverses)\n\
  -k#  # of output frequences [default=n/2]\n\
  -k-# as above, show in range (-k,k), where -1 == n-1; useful for even/odd\n\
       see also -C HACK\n\
  -K#  as -k in units of n (-n must precede -k) [default=0.5]\n\
  -m0  calculate all peaks [default=-1=do not calculate peaks]\n\
  -m#  calculate only peaks higher than #\n\
  -n#  # of valid data lines on input [no default]\n\
       more data ignored, missing data are padded according to option -a\n\
  -o#  omit # lines of data from start [default=0]\n\
  -p   plot the output spectrum (only if FILE given)\n\
  -q   quantity is quadratically dependent on amplitude (as kin.energy):\n\
       the frequency is halved, compatible with -v\n\
  -s   for step at t=0: the 1st value is halved, not good with -w,-C0.5,-cSG\n\
  -v#  apply Verlet frequency correction, h = #*Verlet timestep, -v=-v1 [none]\n\
  -w#  multiply input data by bell-like window [sin(PI*(i+1/2)/NDATA)]^#\n\
       accuracy of peaks increases, #=1 or 2 recommended, large # -> Gaussian\n\
       WARNING: with -w, only a few data may be padded\n\
  -T   with -h or -g: print frequencies in 1/h (this is THz for h in ps)\n\
       [default=wavenumbers in cm^-1 for -h in ps]\n\
METHODS:\n\
  * fast Fourier transform with arbitrary factoring\n\
  * may be slow if the number of data is not a product of small primes:\n\
    - a warning is printed (use -eBIG to suppress)\n\
    - change -n: either omit some data or padd by zeros\n\
  * the peak searcher uses a parabola of 4 points closest to a maximum\n\
EXAMPLES:\n\
  tab 1 1024 | tabproc \"abs(cos(A))\" | spectrum -n1024 - | plot -:0:1\n\
  showcp -a -b1 H2O; spectrum -n20000 -k5000 -h0.002 -c2 -p H2O.cpa\n\
");
    exit(0); }

  initscroll(0);

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-') {
      switch (arg[iarg][1]) {
        case 0: fn=NULL; break;
        case 'a': paddbyaverage++; break;
        case 'C':
          Coutput=1;
          shift=atof(arg[iarg]+2);
          break;
        case 'c':
          Im=arg[iarg]+2;
          if (*Im=='+') REfill=1,Im++;
          else if (*Im=='-') REfill=-1,Im++;
          col=atoi(Im)-1;
          Im=strchr(arg[iarg],',');
          if (Im) {
            Im++;
            if (*Im=='+') IMfill=1,Im++;
            else if (*Im=='-') IMfill=-1,Im++;
            Icol=atoi(Im)-1; }
          break;
        case 'd': prec=atoi(arg[iarg]+2); break;
        case 'f': F=arg[iarg+2]; break;
        case 'i': copyin=1; break;
        case 'I': direction=-1; break;
        case 'o': omit=atoi(arg[iarg]+2); break;
        case 'm': minpeak=atof(arg[iarg]+2); break;
        case 'n': n=atoi(arg[iarg]+2); break;
        case 'k': k=atoi(arg[iarg]+2); break;
        case 'K': k=(int)(n*atof(arg[iarg]+2)+0.5); break;
        case 'w': win=atoi(arg[iarg]+2); break;
        case 'h': h=atof(arg[iarg]+2); break;
        case 'g': h=1/atof(arg[iarg]+2); break;          
        case 'p': plot++; break;
        case 'q': isq++; break;
        case 's': isstep++; break;
        case 'v':
          isverlet=atoi(arg[iarg]+2);
          if (!isverlet) isverlet=1;
          break;
        case 'T': iscm=0; break;
        case 'e': slowdown=atof(arg[iarg]+2); break;
        default: ERROR(("%s is wrong argument",arg[iarg])); } }
    else
      fn=arg[iarg];

  //  if (Im && Icol<0) Icol=col+1;
  maxcol=max(col,Icol);

  if (REfill*IMfill==1) fprintf(stderr,"WARNING: odd+odd or even+even fill (option -c)\n");
  if (REfill && !IMfill || !REfill && IMfill) Error("spectrum: even/odd fill for both Re and Im must be specified");
  if (REfill*IMfill) {
    if (n&1) Error("spectrum: odd n not implemented with even/odd fill");
    if (shift!=0 && shift!=0.5) Error("even/odd fill requires shift=0 or 0.5");
    if (win) fprintf(stderr,"WARNING: window (-w) with REfill and IMfill dos not make sense\n"); }

  {
    int sum=0,blk=n,sqblk,d;
    double q;

    while (blk>1) {
      sqblk=(int)sqrt(blk+0.5);
      loopto (d,2,sqblk) if (blk % d == 0) break;

      if (blk % d != 0) d=blk; /* blk = prime */

      sum+=d;

      blk=blk/d; }

    q=sum/log(n)/3; /* should be /e ideally */

    if (q>slowdown) WARNING(("%d is a product of not so small primes,\n\
*** estimated slowdown = %.3g (-eBIG to suppress this message)",n,q))
  }

  if (fn) {
    in=fopen(fn,"rt");
    if (!in) ERROR(("%s not found",fn))
    if (F) {
      if (*F) out=fopen(string("%s.%s.sp",fn,F),"wt");
      else out=fopen(string("%s.%d.sp",fn,col+1),"wt"); }
    else out=fopen(string("%s.sp",fn),"wt");
    if (h && minpeak>=0) {
      if (F) {
        if (*F) fpeak=fopen(string("%s.%s.pk",fn,F),"wt");
        else out=fopen(string("%s.%d.pk",fn,col+1),"wt"); }
      else fpeak=fopen(string("%s.pk",fn),"wt"); } }
  else {
    in=stdin;
    out=fpeak=stdout; }

  if (n<=0) ERROR(("no or bad # of data (option -n)"))
  if (!k) k=n/2;
  if (k>n) ERROR(("k>n: change options -k (or -n)"))

  allocarrayzero(y,n);

  loop (i,-omit,n) {
    complex yy={0,0};
    char line[1024],*tok;

    if (!missing) {
      while ( (tok=fgets(line,1024,in)) && (line[0]==0 || strchr("!#",line[0])));
      if (!tok) {
        if (i<=0) Error("not enough data on input");
        fprintf(stderr,"WARNING: only %d data lines on input (%d expected)\n",i+omit,n+omit);
        missing++;
        sum.re/=ndata; sum.im/=ndata;
        if (!paddbyaverage) sum.re=sum.im=0;
        fprintf(stderr,"missing Re data will be padded by %.*g %.*g%s\n",
                prec,sum.re,prec,sum.im,
                REfill*IMfill?"\nand/or padded according to parity":"");
        if (win) fprintf(stderr,"WARNING: -w specified => the window may not be properly centered\n\
(only a few padded data are acceptable without a loss of precision)\n"); } }

    if (i>=0) {
      if (missing) {
        yy.re=sum.re; yy.im=sum.im; }
      else {
        tok=strtok(line," \t\n\r");

        loopto (ic,0,maxcol) {
          if (!tok) ERROR(("line %d not enough columns",i))
          if (ic==col) yy.re=atof(tok);
          if (ic==Icol) yy.im=atof(tok);
          tok=strtok(NULL," \t\n\r"); }

        sum.re+=yy.re; sum.im+=yy.im;
        ndata++; }

      y[i].re=yy.re;
      y[i].im=yy.im; } }

  /* correcting step at t=0 */
  if (isstep) y[0].re/=2,y[0].im/=2;

  if (win) loop (i,0,n) {
      double w=window(2*PI*(i+0.5)/n);
      y[i].re*=w;
      y[i].im*=w; }

  if (REfill*IMfill) {
    if (shift==0.5) {
      /* i=0 corresponds to t=h/2 */
      if (REfill) loop (i,0,n/2) y[n-1-i].re=REfill*y[i].re;
      if (IMfill) loop (i,0,n/2) y[n-1-i].im=IMfill*y[i].im; }
    else {
      /* i=0 corresponds to t=0 */
      if (REfill==-1 && y[0].re!=0) Error("with shift=0, odd function Re(t) implies y[0].re=0");
      if (IMfill==-1 && y[0].im!=0) Error("with shift=0, odd function Im(t) implies y[0].im=0");
      if (REfill) loop (i,1,n/2) y[n-i].re=REfill*y[i].re;
      if (IMfill) loop (i,1,n/2) y[n-i].im=IMfill*y[i].im; } }

  if (copyin) {
    FILE *aux;
    if (F) {
      if (*F) aux=fopen(string("%s.%s.in",fn,F),"wt");
      else aux=fopen(string("%s.%d.in",fn,col+1),"wt"); }
    else aux=fopen(string("%s.in",fn),"wt");
    loop (i,0,n)
      if (h) fprintf(aux,"%.*g %.*g %.*g\n",prec,(i+shift)*h,prec,y[i].re,prec,y[i].im);
      else fprintf(aux,"%.*g %.*g\n",prec,y[i].re,prec,y[i].im); }

  Fourier(y,direction*n);

  /* h=0: unnormalized (Fourier series)
     h defined: Fouries transform (integral) */
  if (h) loop (i,0,n) { y[i].re*=h; y[i].im*=h; }

  if (shift) {
    /* multiply y[j] by exp(i*shift*j/n) */
    complex q={direction*cos(2*PI*shift/n),sin(direction*2*PI*shift/n)},z=q,x;

    loop (i,1,n) {
      x.re=y[i].re*z.re-y[i].im*z.im;
      x.im=y[i].re*z.im+y[i].im*z.re;
      y[i]=x;

      if ((i&31)==31) {
        /* refresh for better precision */
        z.re=cos(direction*2*PI*shift/n*(i+1));
        z.im=sin(direction*2*PI*shift/n*(i+1)); }
      else {
        x.re=q.re*z.re-q.im*z.im;
        x.im=q.re*z.im+q.im*z.re;
        z=x; } } }

  if ((int)(shift*2)&1 && k<0)
    fprintf(stderr,"\
spectrum HACK! half-shift (\"leap-frog\") and negative -k\n\
  => spectrum for k<0 multiplied by -1\n\
function must have zero limits for t->+-infty\n\
");
  ik=0;
  if (k<0) k=abs(k),ik=-k+1;

  allocarrayzero(sp,k);


  loop (i,ik,k) {
    int ii=(i+n)%n;
    double spi,sg=1;

    if ((int)(shift*2)&1 && i<0) sg=-1;

    if (h) prt_("%.*g ",prec,conv(i/(h*n)));
    spi=sqrt(Sqr(y[ii].re)+Sqr(y[ii].im));
    if (Coutput) prt("%.*g %.*g",prec,y[ii].re*sg,prec,y[ii].im*sg);
    else prt("%.*g",prec,spi);
    if (i>=0) sp[i]=spi; }

  if (fn) {
    fclose(out);
    out=fpeak; }

  if (h && out)
    loop (i,1,k-5)
      if (sp[i]<sp[i+1] && sp[i+1]<sp[i+2] && sp[i+2]>sp[i+3] && sp[i+3]>sp[i+4]) {
        real p=peak(sp+i,4);
        if (p && A>minpeak) prt("\n%.*g 0\n%.*g %.*g",
                   prec,conv((i+p)/(h*n)),
                   prec,conv((i+p)/(h*n)),
                   prec,A); }

  if (out) fclose(out);

  if (plot) {
    if (h) i=system(string("plot %s.sp:1:2 %s.pk",fn,fn));
    else i=system(string("plot %s.sp:0:1",fn));
    if (i) {
      fprintf(stderr,"plot failed\n");
      return 1; } }

  return 0;
}
