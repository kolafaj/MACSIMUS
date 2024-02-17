/*****************************************************************************
  support for Gaussian charges in the Ewald summation
  see also erfc.c and erfcqq.c
  requires library functions erfc()
*****************************************************************************/

/* compatibility WARNING: SUBGRID>1 is out of order */

#  define SUBGRID 1

#define SQRTPI 1.772453850905516

#if COULOMB>=-2
#  error "COULOMB<-2 expected"
#endif

#ifndef QQTAB
#  error "QQTAB expected"
#endif

#ifdef TWODIM
#  error TWODIM not supported
#endif

struct Erfc_s Erfc;

static struct ga_s {
  double qq;  /* product of charges */
  double qfactor; /* factor14; 1=full interaction, 0=exclusion */
  double qfactor_1; /* factor14-1; 0=full interaction, -1=exclusion */
  double alpha; /* Erfc.alpha = Ewald separation parameter */
  double g12; /* interaction of Gaussian charges GC1-GC2 */
  double sg1; /* screening Ewald charge and GC1 */
  double sg2; /* screening Ewald charge and GC2 */
  double maxg; /* max(g12,sg1,sg2) */
} ga;

static void erfcs(int to,double *e,double *z) /******************** erfcs */
{
  /* library erfc() used */
  double rr=to*Erfc.h,r=sqrt(rr);

  if (to<0) ERROR(("erfcs: negative rr (imaginary r) not supported"))

  if (r*ga.maxg>1) {
    double e12=erfc(ga.g12*r);
    double e1=erfc(ga.sg1*r);
    double e2=erfc(ga.sg2*r);
    *e=(ga.qq/r)*((e1+e2)/2-ga.qfactor*e12+ga.qfactor_1); }
  else if (r*ga.maxg>0.003) {
    double e12=erf(ga.g12*r);
    double e1=erf(ga.sg1*r);
    double e2=erf(ga.sg2*r);

    *e=(ga.qq/r)*(e12*ga.qfactor-(e1+e2)/2); }
  else {
    /* expansion, error ~ r^6; see gaussiancharges-limit.mw */
    *e=ga.qq/sqrt(PI)*(2*ga.g12*ga.qfactor-ga.sg1-ga.sg2
        +rr/3*(-2*Cub(ga.g12)*ga.qfactor+Cub(ga.sg1)+Cub(ga.sg2))
        +rr*rr/10*(2*Pow5(ga.g12)*ga.qfactor-Pow5(ga.sg1)-Pow5(ga.sg2)));
    *z=ga.qq/sqrt(PI)*((4./3)*((ga.g12*ga.g12)*ga.g12)*ga.qfactor
      -(2./3)*((ga.sg1*ga.sg1)*ga.sg1)
      -(2./3)*((ga.sg2*ga.sg2)*ga.sg2)
      +(-(4./5)*(((ga.g12*ga.g12)*(ga.g12*ga.g12))*ga.g12)*ga.qfactor
        +(2./5)*(((ga.sg1*ga.sg1)*(ga.sg1*ga.sg1))*ga.sg1)
        +(2./5)*(((ga.sg2*ga.sg2)*(ga.sg2*ga.sg2))*ga.sg2))*(r*r)
      +(-(1./7)*((((ga.sg1*ga.sg1)*ga.sg1)*((ga.sg1*ga.sg1)*ga.sg1))*ga.sg1)
        -(1./7)*((((ga.sg2*ga.sg2)*ga.sg2)*((ga.sg2*ga.sg2)*ga.sg2))*ga.sg2)
        +(2./7)*((((ga.g12*ga.g12)*ga.g12)*((ga.g12*ga.g12)*ga.g12))*ga.g12)*ga.qfactor)*((r*r)*(r*r)));
      return;
  }
  *z=*e/rr-ga.qq/rr*(1/SQRTPI)*(
       2*ga.g12*exp(-Sqr(ga.g12)*rr)*ga.qfactor
        -ga.sg1*exp(-Sqr(ga.sg1)*rr)
        -ga.sg2*exp(-Sqr(ga.sg2)*rr));
}

#include "splines.c"

double initerfc(ertab_p *tab, double qq,
                double sigma1, double sigma2,
                double factor14, int Grid, double maxr, double cutoff, double Alpha, int shift)
{
  /* Erfc.tab is allocated */
  double alpha=fabs(Alpha); /* refreshed meaning: -v4 verbose passed here */
  long ito;
  int grid=abs(Grid),verbose=Grid>0;
  double x,z;
  double minr=0; /* not argument */
  ermacrodcl

  if (qq==0) {
    ito=1;
    rallocarrayzero(*tab,ito);
    //    (*tab)->Cd=(*tab)->Cu=1;
    return 1e-9; }

  ito=(int)(Sqr(maxr*alpha)*grid)+1;
  /* ito is in units of 1/grid */

  /* to be solved: already initialized */
  Erfc.shift=shift;
  Erfc.sgrid=grid*Sqr(alpha);
  Erfc.grid=grid;
  Erfc.minr=0; Erfc.maxr=maxr;
  Erfc.alpha=alpha;

#ifdef erexact_eps
  Erfc.A=alpha*(-2/SQRTPI);
  Erfc.alphaq=Sqr(alpha);
  Erfc.B=Erfc.A*Erfc.alphaq*2;
#endif /*# erexact_eps */

  /* cf. initss in simdef.c */
  rallocarray(*tab,ito);

  // for macro compatibility
#  define ERFC Erfc
  Erfc.tab=*tab;

  prt("\n:::::: erfc(qq=%g p.u.) :::::: tab 0..%ld = %ld B  grid=%d", qq,ito,ito*sizeof(ertab_t),grid);

  if (verbose) {
    prt("splines based on library erfc");
    prt("alpha=%g  sigma=%g  range=[%g,%g)  cutoff=%g",
        alpha,sqrt(Sqr(sigma1)+Sqr(sigma2)), minr,maxr,cutoff); }

  //  Erfc.from=Erfc.tab+ifrom; Erfc.to=Erfc.tab+ito;

  /* spline calculation */
  ga.qq=qq;
  ga.g12=1/sqrt(2*(Sqr(sigma1)+Sqr(sigma2)));
  ga.alpha=alpha;
  ga.sg1=1/sqrt(1/Sqr(alpha)+2*Sqr(sigma1));
  ga.sg2=1/sqrt(1/Sqr(alpha)+2*Sqr(sigma2));
  ga.maxg=fmax(ga.sg2,ga.sg1);
  Max(ga.maxg,ga.g12)
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

  if (shift&4 || Alpha<0)
#include "ertest.c"

    if (Erfc.sgrid<3) WARNING(("Erfc.sgrid=%g<3: splines inaccurate at r=1 AA\n\
*** check Ewald parameters, particularly el.grid might be insufficient",Erfc.sgrid))
    
  return Erfc.sgrid;
} /* initerfc */

/* this irrelevant for QQTAB=0 because qq=0 */
double exacterud_sqrt(double rr,erreal *erd)
{
  *erd=0;
  return 0;
}

double exacterud_sqrt_1(double rr,erreal *erd)
{
  *erd=0;
  return 0;
}
