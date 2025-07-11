/*****************************************************************************

Fast (but memory/cache consuming) method to calculate the following functions:

  eru(x) = alpha*e(sqrt(x)*alpha) = erfc(sqrt(x)*alpha)/sqrt(x)
  erd(x) = -alpha^3*e'(sqrt(x)*alpha)/(sqrt(x)*alpha)

where

  e(y) = erfc(y)/y
  e'(y) = de(y)/dy

                         infinity
  erfc(y) = 2/sqrt(PI) * integral exp(-t**2) dt,
                            y

  alpha = Ewald separation parameter (el.alphar)

Method:
  hyperbolic spline (functions A+B/(x+C), no jumps in the derivative)

Interface:
  eru(x), erd(x) are macroized
  erud(x) = eru(x) and the value of erd(x) available as buyer
  use ermacrodcl do declare thinks these macros need (incl. byerd)
  #define FLOAT: DEPRECATED

Initialization:
  initerfc(int grid, double minr, double maxr, double cutoff, double alpha, int shift, double kappa, real *alphar)
where:
  grid is a number of grid points per 1 AA^2 unity; negative grid = less verbose
  [minr,maxr) defines the range of r, the range of the argument x of eru, erd,
    and erud is then [minr^2,maxr^2)
    r < minr leads to low-accuracy results
    r > maxr is out of table (not checked, the simulation will likely crash)
  cutoff and shift:
    shift&1: eru is shifted to avoid jump at r=cutoff
    shift&2: erd is shifted to avoid jump at r=cutoff
    shift&4: Hammonds-Heyes method
  cutoff<=rmax required; typically, maxr=box.cutoff+el.rplus
  alpha = separation parameter
  kappa = Ewald k-space parameter
  *alphar = alpha' if shift&4 (Hammonds-Heyes method ) returned (for r-space)

optimization of hyperbolic splines as well as TWODIM (2D Ewald support)
were removed, see cook V2.7l and older or sim/old+misc/erfc-eropt.c

*****************************************************************************/

/*
#define SUBGRID
  1 = use library function erfc() for calculations (recommended, C99 needed)
  >1 = use NUMERICAL INTEGRATION with SUBGRID division points in one spline 
       interval (minimum 2, recommended 4)
  Use SUBGRID>1 if erfc() is not in the math library or it is imprecise:
  This was the only option in V2.7l or older with Erfc.ndiv = SUBGRID,
  see sim/old+misc/erfc-eropt.c.
*/

#define SQRTPI 1.772453850905516

#if !(COULOMB==-1 || COULOMB==-2) 
#  error "COULOMB=-1 or -2 expected"
#endif /*# !(COULOMB==-1 || COULOMB==-2)  */

struct Erfc_s Erfc;

static void erfcs(int to,double *e,double *z) /********************** erfcs */
/*
  e=eru(h*to) and z=erd(h*to) are calculated - will be alpha-rescaled later
*/
{
#if SUBGRID==1
  /* library erfc() used */
  if (to) {
    double rr=to*Erfc.h,r=sqrt(rr),er=erfc(r);

    *e=er/r;
    *z=((2/SQRTPI)*exp(-rr)+er/r)/rr; }
  else
    *e=*z=9e9;
#else /*# SUBGRID==1 */
/*
  LEGACY: Numerical integration and asymptotic formula used to calculate erfc.
  (used to be also for TWODIM).
  Must be called with decreasing sequence of parameters  to , starting
  with a sufficiently high value because the integral is for this  to
  approximated by the asymptotic formula.
*/
#  define f(X) (sqrt(X)*exp(X))
  static long t= -1;
  static double sum;
  double
    w=Erfc.h/2, ha=Erfc.h*0.2113248654051871, hb=Erfc.h-ha,
    x=to*Erfc.h;

  if (to>t) {
    t=to; x=to*Erfc.h;
    if (x<1) ERROR(("erfcs: range too short"))
    /* the correct asymptotic formula is
       exp(-x)/sqrt(x+1-0.75/x+2/x^2-7.75/x^3+-... */
    sum=exp(-x)/sqrt(x+1-0.7/x);  }

  while (t>to) {
    x=Erfc.h*t;
    sum += w/f(x-ha)+w/f(x-hb);
    t--; }

  x=t*Erfc.h+1e-20;
  *e=sum/sqrt(x);
  *z=(*e+2/exp(x))/x/SQRTPI;
  *e /= SQRTPI;
#  undef f
#endif  /*#!SUBGRID==1 */
} /* erfcs */

#include "splines.c"

void initerfc(int grid, double minr, double maxr, double cutoff, double alpha,
              int shift, double kappa, erreal *alphar)           /* initerfc */
{
  int dump=alpha<0;
  int verbose=grid>0;
  long ifrom,ito;
  double x,z,ar;
  ermacrodcl

  grid=abs(grid);
  if (grid==0) ERROR(("spline grid=0"))
  ar=alpha=fabs(alpha);
                 
  if (shift&4) ar=alpha*erf(PI*kappa/alpha);
                 
  if (Erfc.grid==grid && Erfc.shift==shift
      && Erfc.minr==minr && Erfc.maxr==maxr
      && Erfc.alpha==ar) {
    prt(":::::: erfc already initialized ::::::");
    return; }

  if (shift&4) 
    prt("Hammonds-Heyes: alpha=%g -> alphar=%g (EXPERIMENTAL FEATURE)",alpha,*alphar);
  *alphar=ar;
  
  Erfc.shift=shift;
  Erfc.sgrid=(erreal)grid; /* scaled as for alpha=1 - will be rescaled below */
  Erfc.grid=grid;
  Erfc.minr=minr; Erfc.maxr=maxr;
  Erfc.alpha=*alphar;
  
  /* ifrom,ito are in units of 1/grid */
  ifrom=(int)(Sqr(minr*Erfc.alpha)*grid);
  ito=(int)(Sqr(maxr*Erfc.alpha)*grid)+1;

#ifdef erexact_eps
  Erfc.A=*alphar*(-2/SQRTPI);
  Erfc.alphaq=Sqr(Erfc.alpha);
  Erfc.B=Erfc.A*Erfc.alphaq*2;
#endif /*# erexact_eps */

  if (Erfc.tab!=NULL) free(Erfc.tab);
  allocarray(Erfc.tab,ito);

  prt("\n:::::: erfc :::::: tab 0..%ld..%ld = %ld B  grid=%d", ifrom,ito,ito*sizeof(ertab_t),grid);

  if (verbose) {
#if SUBGRID==1 
    prt("splines based on library erfc");
#else /*# SUBGRID==1  */
    prt("splines based on numerical integration (SUBGRID=%d)",SUBGRID);
#endif /*#!SUBGRID==1  */
    prt("alphar=%g  range=[%g,%g)  cutoff=%g  1st interval=[0,%g]",
        *alphar,minr,maxr,cutoff,sqrt(1/Erfc.sgrid)); }

  /* spline calculation */
  makesplines(ito,Erfc.alpha);

  z=Sqr(cutoff);
  
  if (shift&1)  {
    x=eru(z);
    loop (er_p,Erfc.tab,Erfc.tab+ito) er_p->Au -= x;
    prt("eru(r) shifted by %g to avoid jump in elst energy",x); }
  
  if (shift&2)  {
    x=erd(z);
    loop (er_p,Erfc.tab,Erfc.tab+ito) er_p->Ad -= x;
    prt("erd(r) shifted by %g to avoid jump in elst forces",x); }
  
#include "ertest.c"

} /* initerfc */

#ifdef erexact_eps
#  include "erexact.c"
#endif /*# erexact_eps */
