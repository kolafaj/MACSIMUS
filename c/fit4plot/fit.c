/* \make fit
  This version is normally called from plot (see MACSIMUS/show/plot.c).
  for a stand-alone version, see myfit.c.

  Change 02/2022: small change of output

  Change 09/2021: higher derivatives, optionally 4th order formulas
  
  Change 12/2020: small change of output

  Change 09/2017: environment variable FIT =
     func X       : calculate the function value for X
     expr X       : calculate expression (y-column of the 3rd argument to plot)
     solve Y[,X0] : solve equation f(X)=Y, X0=initial iteration
     integ NINT[,XFROM,XTO] : calculate integral
       NINT = # if subintervals by 4th-order (2-point) Gauss quadrature
       XFROM = integrate from (default = 1st x)
       XTO = integrate to (default = last x)
     deriv X,DX : calculate the derivative at X
     (max 1 of each type, uncertainties included)

  Change 02/2016: environment variable INTEG (in 9/2017 replaced by FIT)

  getdata variables:
method: -1: steepest descent
         0: amoeba (default)
         1: conjugate gradient
         2: Newton-Raphson
         3: Monte Carlo
         9: exit getdata loop
D: step to calculate numerical derivatives (method=-1,1,2 only, default=1e-5)
eps: accuracy of the minimum (irrelevant for method=3, default=1e-8)
par: `nonsphericity' valley param. par>=1 (method=3 only, default=1)
     par<=0: repeats w. extrapolation (narrow valley, all methods)
maxit: (max) number of steps, negative = quiet
nerr: # of samples to calculate errors (default=0=do not calculate err.)
err: err=-1: uses stdev of data provided (assuming 1 if not provided)
     err=-2: as err=-1 if stdev of data available, err=0 otherwise (default)
     err=-3: stdev of dy_i is proportional to 3rd column,
             rescaled by stdev calculated from residual sum of squares
             (for data with known weights~1/dy_i^2 but unknown precision)
     err>=0: uses stdev calc. from resid. sum sq., multiplied by y^err:
     err=0: the same abs.err. (of stdev above) of all data
     err=1: the same rel.err. of all data
     err=0.5: error is proportional to y^0.5 (suited for histogram data)
*/

#include "ground.h"
#include "minimize.h"
#include "qmin.h"
#include "ms.h"

/* to be able to use pi */
#define pi PI

// created by plot, hot keys t T ^t
#include "fitinc.c"

REAL A[N],A0[N];

REAL Asum[N],Acsum[N][N];
volatile int sig;

#define MAXn 1048576
int n;
double x[MAXn],y[MAXn],*dy,dy0[MAXn];

REAL ff(REAL *A) /******************************************************* ff */
{
  REAL s=0,f;
  int i;

  loop (i,0,n) {
    f=func(A,x[i]);
    s+=Sqr((f-y[i])/dy[i]); }

  return s;
}

REAL xsqrt(REAL x) /************************************************** xsqrt */
{
  if (x<0) return 0;
  else return sqrt(x);
}

void dody(double err,double sigma) /*********************************** dody */
{
  int i;

  if (err==-3) loop (i,0,n) dy[i]=dy0[i]*sigma;
  else if (err<0) loop (i,0,n) dy[i]=dy0[i];
  else loop (i,0,n) dy[i]=pow(fabs(y[i]),err)*sigma;
}

int main(int narg,char **arg) /**************************************** main */
{
  char line[1024];
  FILE *f, *outf;
  REAL D=1e-5,eps=1e-8,par=1;
  int method=0,maxit=1000;
  int i,j,nerr=0,ierr,isstdev=0,order=1;
  const char *Order[5]={"","1st","2nd","3rd","4th"};
  double err=-2,sigma,sqsum;
  int nint=0;
  enum FIT_e { NONE,EXPR,INTEG,DERIV,SOLVE } fitkey=NONE;
  double xderiv=0,dxderiv=1e-5,xfunc=0,yy=0,x0=0;
  double xfrom=0,xto=0,integx=0,derivx=0,funcx=0,*y_0=NULL,*yfit=NULL;
  char *FIT=getenv("FIT");
  const char *feinfo;
  char *RNDSEED=getenv("RNDSEED");
  int rndseed=RNDSEED?atoi(RNDSEED):0;

  initscroll(0);

  if (getenv("INTEG"))
    ERROR(("INTEG no longer supported, use FIT=integ ... to integrate"))

  rndinit(rndseed,0);

  if (narg<2) {
    fprintf(stderr,"Call by (normally from plot):\n\
  fit INPUTFILE\n");
    exit(0); }

  f=fopen(arg[1],"rt");
  if (!f) ERROR(("open %s",arg[1]))

  while (fgets(line,1024,f)) {
    if (line[0]=='#') continue;
    if (n>=MAXn) ERROR(("tab ofl, increase MAXn and recompile"))
    i=sscanf(line,"%lf%lf%lf",x+n,y+n,dy0+n);
    if (i==1) ERROR(("%s bad format",line))
    if (i<=0) continue;
    if (i==2) {
      if (n) dy0[n]=dy0[n-1];
      else dy0[n]=1; }
    if (i==3) isstdev++;
    n++; }

  allocarrayzero(dy,n);
  if (FIT) {
    /* parsing FIT */
    char *comma,*c=FIT;

    switch (*c) {

      case 'e':
        while (!strchr(".0123456789+- ",*c)) { if (*c==0) goto endFIT; c++; }
        fitkey=EXPR;
        xfunc=atof(c);
        break;

      case 'i':
        while (!strchr(".0123456789+- ",*c)) { if (*c==0) goto endFIT; c++; }
        fitkey=INTEG;
        nint=atoi(c);
        xfrom=x[0]; /* default = first point */
        xto=x[n-1]; /* default = last point */
        if ( (comma=strchr(c,',')) ) {
          comma++;
          while (*comma==' ') comma++;
          if (*comma!=',') xfrom=atof(comma);
          if ( (comma=strchr(comma,',')) ) xto=atof(++comma); }
        break;

      case 'd':
        while (!strchr(".0123456789+- ",*c)) { if (*c==0) goto endFIT; c++; }
        fitkey=DERIV;
        xderiv=atof(c);
        if ( (comma=strchr(c,',')) ) dxderiv=atof(++comma);
        if ( (comma=strchr(comma,',')) ) order=atoi(++comma);
        else order=1;
        break;

      case 's':
        while (!strchr(".0123456789+- ",*c)) { if (*c==0) goto endFIT; c++; }
        fitkey=SOLVE;
        yy=atof(c);
        if ( (comma=strchr(c,',')) ) x0=atof(++comma);
        break;

      default:
        ERROR(("fit: %s is bad keyword in environment variable FIT\n\
*** use one of {expr,integ,deriv,solve}",c)) }
  endFIT:; }

  for (;;) {
    double dd;
    getdata
      getvec(A,,N)
      get(method)
      get(maxit) get(eps) get(par) get(D)
      get(nerr) get(err)
      checkdata
    enddata

    if (method<-1 || method>4) return 0;

    if (err==-2 && !isstdev) err=0;

    dody(err,1);
    if (nerr) {
      if (!yfit) allocarrayzero(yfit,n);
      if (!y_0) allocarrayzero(y_0,n); }

    Minimize(method,N,A,ff,maxit,eps,D,par);

    sqsum=Val(ff(A));
    sigma=sqrt(sqsum/(n-N));

    // prt("sigma=%g, adj. to deg.f.=%g",(double)(sqrt(sqsum/n)),(double)sigma);
    prt("! nu=%d ### s=%g ###",n-N,sigma);


    /*** file fit.fit ***/
    outf=fopen("fit.fit","wt");
    if (!outf) Error("write fit.fit");

    loop (i,0,N) fprintf(outf,"# %c=%.15g\n",i+'a',Val(A[i]));
    fprintf(outf,"#        x         y_data       y_fit       stdev_data   (y_fit-y_data)/stdev_data\n");

    loop (i,0,n) {
      fprintf(outf,"%12.6g %14.8g %14.8g %12.6g %10.4g\n",
              x[i],y[i],dd=Val(func(A,x[i])),dy[i],Val((func(A,x[i])-y[i])/dy[i]));
      if (nerr) {
        yfit[i]=dd;
        y_0[i]=y[i]; } }


    feinfo=funcexpr==func?"func":"expr";

    /*
      These calculations are repeated if error is to be calculated
      WARNING: codes are duplicated!
    */

    if (fitkey==EXPR) {
      REAL X=xfunc;

      funcx=Val(funcexpr(A,X));
      fprintf(outf,"# %s(%g) = %.9g\n",feinfo,xfunc,funcx);
      prt(           "%s(%g) = %.9g",  feinfo,xfunc,funcx); }

    if (fitkey==INTEG) {
      REAL sum=0;
      REAL d,X;
      int j;

      d=(xto-xfrom)/nint;
      loop (j,0,nint) {
        X=xfrom+d*(j+0.2113248654051871);
        sum+=funcexpr(A,X);
        X=xfrom+d*(j+0.7886751345948129);
        sum+=funcexpr(A,X); }

      sum*=d/2;
      integx=Val(sum);
      fprintf(outf,"# integral_%s[%g,%g]_%d = %.9g\n",feinfo,xfrom,xto,nint,integx);
      prt(           "integral_%s[%g,%g]_%d = %.9g",  feinfo,xfrom,xto,nint,integx); }

    if (fitkey==DERIV) {
      REAL XP=(REAL)xderiv+(REAL)dxderiv,XM=(REAL)xderiv-(REAL)dxderiv,XPP,XMM;
      REAL FP=funcexpr(A,XP),FM=funcexpr(A,XM);

      switch (order) {
        case 1:
          /* 2nd order 1st derivative */
          derivx=Val((FP-FM)/(REAL)(2.0*dxderiv));
          break;
        case -1:
          /* 4th order 1st derivative */
          XPP=(REAL)xderiv+(REAL)(2.0*dxderiv);
          XMM=(REAL)xderiv-(REAL)(2.0*dxderiv);
          derivx=Val((-funcexpr(A,XPP)+(FP-FM)*8.0+funcexpr(A,XMM))/(REAL)(12.0*dxderiv));
          break;
        case 2:
          /* 2nd order 2nd derivative */
          derivx=Val((FP-2.0*funcexpr(A,xderiv)+FM)/(REAL)(Sqr(dxderiv)));
          break;
        case -2:
          /* 4th order 2nd derivative */
          XPP=(REAL)xderiv+(REAL)(2.0*dxderiv);
          XMM=(REAL)xderiv-(REAL)(2.0*dxderiv);
          derivx=Val((-funcexpr(A,XPP)+(FP+FM)*16.0-funcexpr(A,xderiv)*30.0-funcexpr(A,XMM))/(REAL)(12.0*Sqr(dxderiv)));
          break;
        case 3:
          /* 2nd order 3rd derivative */
          XPP=(REAL)xderiv+(REAL)(2.0*dxderiv);
          XMM=(REAL)xderiv-(REAL)(2.0*dxderiv);
          derivx=Val((funcexpr(A,XPP)-funcexpr(A,XMM)-2.0*(FP-FM))/(REAL)(2.0*Cub(dxderiv)));
          break;
        case 4:
          /* 2nd order 4th derivative */
          XPP=(REAL)xderiv+(REAL)(2.0*dxderiv);
          XMM=(REAL)xderiv-(REAL)(2.0*dxderiv);
          derivx=Val((funcexpr(A,XPP)+funcexpr(A,XMM)-4.0*(FP+FM)+6.0*funcexpr(A,xderiv))/(REAL)Pow4(dxderiv));
          break;
        default:
          ERROR(("This combination of derivative and order is not supported, use:\n\
     1 = 1st derivative, central formula of 2nd order\n\
    -1 = 1st derivative, central formula of 4th order\n\
     2 = 2nd derivative, central formula of 2nd order\n\
    -2 = 2nd derivative, central formula of 4th order\n\
     3 = 3rd derivative, central formula of 2nd order\n\
     4 = 4th derivative, central formula of 2nd order"))
            }

      fprintf(outf,"# %s_derivative_%s(%g)_%g,%d  = %.9g\n",Order[abs(order)],feinfo,xderiv,dxderiv,order<0?4:2,derivx);
      prt(           "%s_derivative_%s(%g)_%g,%d  = %.9g",  Order[abs(order)],feinfo,xderiv,dxderiv,order<0?4:2,derivx); }

    if (fitkey==SOLVE) {
      REAL X=x0,Y=yy;
      int nit=100;

      MS_BEGIN(X,1e-10)
        MS_f=funcexpr(A,X)-Y;
        if (nit--<0) {
          fprintf(outf,"### solve failed\n");
          prt("*** solve failed");
          break; }
      MS_END(X,1)

      x0=Val(X);
      fprintf(outf,"# %g = %s(%.9g)\n",yy,feinfo,x0);
      prt(           "%g = %s(%.9g)\n",yy,feinfo,x0); }

    fclose(outf);

    /*** file fit.env ***/
    outf=fopen("fit.env","wt");
    loop (i,0,N) fprintf(outf,"export %c=%.15g\n",i+'a',Val(A[i]));
    fclose(outf);

    /*** standard output ***/
    loop (i,0,N) {
      prt("%c=%.15g",i+'a',Val(A[i]));
      A0[i]=A[i];
      Asum[i]=0;
      loop (j,0,N) Acsum[i][j]=0; }

    /* input for ev */
    loop (i,0,N) {
      //      char c=i+'a';
      //      if (strchr("cegh",c)) prt_("~%c;",c); // reserved constants in ev
      //      if (strchr("acegh",c)) prt_("~%c;",c); // reserved constants in evu
      prt_("%c=%.15g;",i+'a',Val(A[i])); }
    _n

    /*** error estimation by Gaussian sampling ***/
    if (nerr) {
      double *yav,*yavq;
      double SUMi=0,SUMiq=0;
      double SUMd=0,SUMdq=0;
      double SUMs=0,SUMsq=0;
      double SUMf=0,SUMfq=0;

      allocarrayzero(yav,n);
      allocarrayzero(yavq,n);

      dody(err,sigma);

      outf=fopen("fit.err","wt");

      loop (ierr,0,nerr) {
        loop (i,0,N) A[i]=A0[i];
        loop (i,0,n) y[i]=yfit[i]+rndgauss()*dy[i];

        Minimize(method,N,A,ff,-maxit,eps,D,par);
        loop (i,0,N) {
          Asum[i]+=A[i];
          loop (j,0,N) Acsum[i][j]+=A[i]*A[j]; }

        loop (i,0,n) {
          REAL fx=funcexpr(A,x[i]);
          yav[i]+=Val(fx);
          yavq[i]+=Sqr(Val(fx)); }

        if (fitkey==EXPR) {
          REAL X=xfunc;
          double f=Val(funcexpr(A,X));

          SUMf+=f;
          SUMfq+=Sqr(f); }

        if (fitkey==INTEG) {
          REAL sum=0;
          REAL d,X;
          int j;

          d=(xto-xfrom)/nint;
          loop (j,0,nint) {
            X=xfrom+d*(j+0.2113248654051871);
            sum+=funcexpr(A,X);
            X=xfrom+d*(j+0.7886751345948129);
            sum+=funcexpr(A,X); }

          fprintf(stderr,"%.16g\n",Val(sum));

          sum*=d/2;
          SUMi+=Val(sum);
          SUMiq+=Val(Sqr(sum)); }

        if (fitkey==DERIV) {
          REAL XP=xderiv+dxderiv,XM=xderiv-dxderiv,XPP,XMM,deriv;
          REAL FP=funcexpr(A,XP),FM=funcexpr(A,XM);

          switch (order) {
            case 1:
              /* 2nd order 1st derivative */
              deriv=(FP-FM)/(REAL)(2.0*dxderiv);
              break;
            case -1:
              /* 4th order 1st derivative */
              XPP=(REAL)xderiv+(REAL)(2.0*dxderiv);
              XMM=(REAL)xderiv-(REAL)(2.0*dxderiv);
              deriv=(-funcexpr(A,XPP)+(FP-FM)*8.0+funcexpr(A,XMM))/(REAL)(12.0*dxderiv);
              break;
            case 2:
              /* 2nd order 2nd derivative */
              deriv=(FP-2.0*funcexpr(A,xderiv)+FM)/(REAL)(Sqr(dxderiv));
              fprintf(stderr,"%.16g\n",Val(deriv));
              break;
            case -2:
              /* 4th order 2nd derivative */
              XPP=xderiv+dxderiv*2;
              XMM=xderiv-dxderiv*2;
              deriv=(-funcexpr(A,XPP)+(FP+FM)*16.0-funcexpr(A,xderiv)*30.0-funcexpr(A,XMM))/(REAL)(12.0*Sqr(dxderiv));
              break;
            case 3:
              /* 2nd order 3rd derivative */
              XPP=xderiv+dxderiv*2;
              XMM=xderiv-dxderiv*2;
              deriv=(funcexpr(A,XPP)-funcexpr(A,XMM)-2.0*(FP-FM))/(REAL)(2.0*Cub(dxderiv));
              break;
            case 4:
              /* 2nd order 4th derivative */
              XPP=xderiv+dxderiv*2;
              XMM=xderiv-dxderiv*2;
              deriv=(funcexpr(A,XPP)+funcexpr(A,XMM)-4.0*(FP+FM)+6.0*funcexpr(A,xderiv))/(REAL)Pow4(dxderiv);
              break;
            default: exit(1); }

          SUMd+=Val(deriv);
          SUMdq+=Val(Sqr(deriv)); }

        if (fitkey==SOLVE) {
          REAL X=x0,Y=yy;

          MS_BEGIN(X,1e-10)
            MS_f=funcexpr(A,X)-Y;
          MS_END(X,1)

          SUMs+=Val(X);
          SUMsq+=Sqr(Val(X)); }

      } /* ierr */

      loop (i,0,n) y[i]=y_0[i];

      prt("\n#     the_fit        stderr    <perturbed_fits> (nerr=%d)",nerr);
      fprintf(outf,"# i   parm[i]   stderr\n");

      loop (i,0,N) {
        A[i]=A0[i];
        if (nerr) {
          REAL Aerr=xsqrt((Acsum[i][i]-Sqr(Asum[i])/nerr)/(nerr-1));
          prt("%c=%14.11f %14.11f %14.11f",
                 i+'a', Val(A[i]), Val(Aerr), Val(Asum[i]/nerr));
          fprintf(outf,"# %d %14.11f %14.11f\n",
                  i,Val(A[i]),   Val(Aerr)); } }

      if (nerr) {
        prt_("\n# covariance matrix:");
        fprintf(outf,"# covariance matrix:");
        loop (i,0,N)
          loop (j,0,N) {
            REAL Aerr=(Acsum[i][j]-Asum[i]*Asum[j]/nerr)/(nerr-1);
            prt_("%c%17.11g", " \n"[j==0],Val(Aerr));
            if (j==0) fprintf(outf,"\n#%c", i+'a');
            fprintf(outf," %16.11g", Val(Aerr)); }
        _n _n
        fprintf(outf,"\n\n"); }

      fprintf(outf,"# nu=%d ### s=%g ###",n-N,(double)sigma);
      fprintf(outf,"# errors based on %d*sampled data\n",nerr);
      fprintf(outf,"\n#      x          y_fit           y_stdev       <y_sampled>\n");

      loop (i,0,n)
        fprintf(outf,"%10.6g  %13.9g %12.6g %13.9g\n",
                x[i],
                yfit[i],
                Val(sqrt((yavq[i]-Sqr(yav[i])/nerr)/(nerr-1))),
                yav[i]/nerr);

      if (fitkey==EXPR) {
        fprintf(outf,"# %s(%g) = %.9g %.6g %.9g\n",
                feinfo,
                xfunc,
                funcx,
                Val(sqrt((SUMfq-Sqr(SUMf)/nerr)/(nerr-1))),
                SUMf/nerr);
        prt("%s(%g) = %.9g %.6g %.9g",
            feinfo,
            xfunc,
            funcx,
            Val(sqrt((SUMfq-Sqr(SUMf)/nerr)/(nerr-1))),
            SUMf/nerr); }

      if (fitkey==INTEG) {
        fprintf(outf,"# integral_%s[%g,%g]_%d = %.9g %.6g %.9g\n",
                feinfo,
                xfrom,xto,nint,
                integx,
                /* Val because sqrt is cast to REAL: */
                Val(sqrt((SUMiq-Sqr(SUMi)/nerr)/(nerr-1))),
                SUMi/nerr);
        prt("integral_%s[%g,%g]_%d = %.9g %6g %.9g",
            feinfo,
            xfrom,xto,nint,
            integx,
            Val(sqrt((SUMiq-Sqr(SUMi)/nerr)/(nerr-1))),
            SUMi/nerr); }

      if (fitkey==DERIV) {
        fprintf(outf,"# %s_derivative_%s(%g)_%g,%d = %.9g %.6g %.9g\n",
                Order[abs(order)],
                feinfo,
                xderiv,dxderiv,order<0?4:2,
                derivx,
                Val(sqrt((SUMdq-Sqr(SUMd)/nerr)/(nerr-1))),
                SUMd/nerr);
        prt("%s_derivative_%s(%g)_%g,%d = %.9g %.6g %.9g",
            Order[abs(order)],
            feinfo,
            xderiv,dxderiv,order<0?4:2,
            derivx,
            Val(sqrt((SUMdq-Sqr(SUMd)/nerr)/(nerr-1))),
            SUMd/nerr); }

      if (fitkey==SOLVE) {
        fprintf(outf,"# %g = %s(%.9g %.6g %.9g)\n",
                yy, feinfo, x0,
                Val(sqrt((SUMdq-Sqr(SUMd)/nerr)/(nerr-1))),
                SUMd/nerr);
        prt("%g = %s(%.9g %.6g %.9g)",
            yy, feinfo, x0,
            Val(sqrt((SUMsq-Sqr(SUMs)/nerr)/(nerr-1))),
            SUMs/nerr); }

      fclose(outf);
    } // nerr

  }

  return 0;
}
