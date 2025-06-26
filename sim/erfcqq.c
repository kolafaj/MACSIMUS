/*****************************************************************************
  version of erfc.c with separate splines for all charge-charge pairs
  QQTAB must be #defined
  library erfc() used (even if SUBGRID>1, which is not recommended)

WARNING: 1-4 terms are incorrect
         hyperbolic splines wrong (unstable) for 1-4

*****************************************************************************/

#error "QQTAB with point charges: This module has not been tested recently. Please test, remove this message, and recompile."

#define SQRTPI 1.772453850905516

#if COULOMB>=0
#  error "COULOMB<0 expected"
#endif /*# COULOMB>=0 */

#ifndef QQTAB
#  error "QQTAB expected"
#endif

struct Erfc_s Erfc;

/* to be compatible with Gaussian charges module ... */
static struct ga_s {
  double qq;  /* product of charges */
  double qfactor; /* factor14; 0 means exclusion */
  double qfactor_1; /* factor14-1; -1 means exclusion */
  double alpha; /* Erfc.alpha = Ewald separation parameter */
} ga;

static void erfcs(int to,double *e,double *z) /*********************** erfcs */
/*
  e=qq*eru(rr) and z=qq*erd(rr) are calculated (no post-scaling as in erfc.c)
*/
{
  /* library erfc() used */
  if (to) {
    double rr=to*Erfc.h,r=sqrt(rr),e12=erfc(ga.alpha*r);

    *e=(ga.qq/r)*(e12+ga.qfactor_1);
    *z=(*e+ga.alpha*(2/SQRTPI)*ga.qq*exp(-Sqr(ga.alpha)*rr))/rr; }
  else
    *e=*z=9e9;
}

#include "splines.c"

double initerfc(ertab_p *tab, double qq, /************************* initerfc */
                double factor14,
                int grid, double minr, double maxr, double cutoff,
                double alpha, int shift, double kappa, erreal *alphar)
{
  int dump=alpha<0;
  int verbose=grid>0;
  /* Erfc.tab is allocated */
  long ifrom,ito;
  double x,z;
  ermacrodcl

  if (qq==0) {
    ito=1;
    rallocarrayzero(*tab,ito);
    //    (*tab)->Cd=(*tab)->Cu=1;
    return 1e-9; }

  grid=abs(grid);
  if (grid==0) ERROR(("spline grid=0"))

#if 0
  /* this is wrong for several q-q pairs - to check these pairs? */
  if (Erfc.grid==grid && Erfc.shift==shift
      && Erfc.minr==minr && Erfc.maxr==maxr
      && Erfc.alpha==*alphar) {
    prt(":::::: erfcqq already initialized ::::::");
    return; }
#endif

  alpha=fabs(alpha);
  if (shift&4) {
    *alphar=alpha*erf(PI*kappa/alpha);
    if (verbose) prt("WARNING: the Hammonds-Heyes method is experimental, be careful"); }
  else
    *alphar=alpha;

  Erfc.shift=shift;
  Erfc.sgrid=grid*Sqr(alpha);
  Erfc.grid=grid;
  Erfc.minr=minr; Erfc.maxr=maxr;
  Erfc.alpha=*alphar;

#ifdef erexact_eps
  Erfc.A=alpha*(-2/SQRTPI);
  Erfc.alphaq=Sqr(alpha);
  Erfc.B=Erfc.A*Erfc.alphaq*2;
#endif /*# erexact_eps */

  /* ifrom,ito are in units of 1/grid */
  ifrom=(int)(Sqr(minr*alpha)*grid);
  ito=(int)(Sqr(maxr*alpha)*grid)+1;

  /* cf. initss in simdef.c */
  rallocarray(*tab,ito);

  /* local macro using ERFC will be used below */
#undef ERFC
#  define ERFC Erfc
  Erfc.tab=*tab;

  prt("\n:::::: erfc(qq=%g p.u.) :::::: tab 0..%ld..%ld = %ld B  grid=%d", qq,ifrom,ito,ito*sizeof(ertab_t),grid);

  if (verbose) {
    prt("splines based on library erfc");
    prt("alphar=%g  range=[%g,%g)  cutoff=%g  1st interval=[0,%g]",
        *alphar,minr,maxr,cutoffsqrt(1/Erfc.sgrid)); }

  //  Erfc.from=Erfc.tab+ifrom; Erfc.to=Erfc.tab+ito;

  /* spline calculation */
  ga.qq=qq;
  ga.alpha=alpha;
  ga.qfactor=factor14;
  ga.qfactor_1=factor14-1;
  makesplines(ito,1);

  z=Sqr(cutoff);

  if (factor14!=1) shift&=0x7fffffc;
  /* note that Ewald corrections and 1-4 cannot be shifted */

  if (shift&1)  {
    x=eru(z);
    loop (er_p,Erfc.tab,Erfc.tab+ito) er_p->Au -= x;
    prt("eru(r) shifted by %g to avoid jump in elst energy",x); }

  if (shift&2)  {
    x=erd(z);
    loop (er_p,Erfc.tab,Erfc.tab+ito) er_p->Ad -= x;
    prt("erd(r) shifted by %g to avoid jump in elst forces",x); }

#include "ertest.c"

  return Erfc.sgrid;
} /* initerfc */

#ifdef erexact_eps
#  include "erexact.c"
#endif /*# erexact_eps */
