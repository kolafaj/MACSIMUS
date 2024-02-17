/*
  This is the Fennel, Gezelter variant of cutelst
  see:
    Fennell CJ, Gezelter JD
    Is the Ewald summation still necessary? Pairwise alternatives to the accepted standard for long-range electrostatics 
    JCP 124, 234104 (2006)
  Approximated by quadratic splines
*/

#if COULOMB!=3
#error COULOMB==3 expected
#endif

struct Erfc_s Erfc = { 1,4,0,NULL };

erreal byerd; /* by-product result of erud
  NOTE: if ermacro is not #defined then this byerd is used,
  if #ifdef ermacro then a local byerd may be faster, however,
  you must not declare local byerd if ermacro is not #defined ! */

static double r1;

static void erfcs(long to,double *e,double *z) /********************** erfcs */
/*
  e=eru(h*to) and z=erd(h*to) are computed.
  Must be called with decreasing sequence of parameters to, starting
  with a sufficiently high value because the integral is for this to
  approximated by the asymptotic formula
*/
{
  double
    h=1e0/(Erfc.ndiv*Erfc.sgrid),
    r=sqrt(h*to);

  if (r<r1) {
    if (r==0) r=sqrt(h)/2; /* patch for /0 */
    *e=erfc(Erfc.alpha*r)/r
      -erfc(Erfc.alpha*r1)/r1
      +(erfc(Erfc.alpha*r1)/Sqr(r1)+2*Erfc.alpha/sqrt(PI)
	*exp(-Sqr(Erfc.alpha*r1))/r1)*(r-r1);
    *z=((erfc(Erfc.alpha*r)/Sqr(r)+2*Erfc.alpha/sqrt(PI)*exp(-Sqr(Erfc.alpha*r))/r)
	-(erfc(Erfc.alpha*r1)/Sqr(r1)+2*Erfc.alpha/sqrt(PI)*exp(-Sqr(Erfc.alpha*r1))/r1))/r; }
  else
    *e=*z=0;
} /* erfcs */

#ifndef ermacro
#error not implemented

double eru(double x) /************************************************** eru */
{
  ertab_p er_p=Erfc.tab+(int)(x*Erfc.sgrid);
#ifdef erminmax
  if (x<Erfc.min) Erfc.min=x;
  if (x>Erfc.max) Erfc.max=x;
#endif
#ifdef errange
  /*.....if (er_p<Erfc.from || er_p>=Erfc.to) Error("eru:range");*/
  if (er_p>=Erfc.to) Error("eru:range");
#endif
  return er_p->Bu/(er_p->Cu+x) + er_p->Au;
} /* eru */

double erd(double x) /************************************************** erd */
{
  ertab_p er_p=Erfc.tab+(int)(x*Erfc.sgrid);
#ifdef erminmax
  if (x<Erfc.min) Erfc.min=x;
  if (x>Erfc.max) Erfc.max=x;
#endif
#ifdef errange
  /*.....if (er_p<Erfc.from || er_p>=Erfc.to) Error("erd:range");*/
  if (er_p>=Erfc.to) Error("erd:range");
#endif
  return er_p->Bd/(er_p->Cd+x) + er_p->Ad;
} /* erd */

double erud(double x) /************************************************ erud */
/* erud computes eru as above and gives erd as a by-product in variable byerd */
{
  ertab_p er_p=Erfc.tab+(int)(x*Erfc.sgrid);
#ifdef erminmax
  if (x<Erfc.min) Erfc.min=x;
  if (x>Erfc.max) Erfc.max=x;
#endif
#ifdef errange
  /*.....if (er_p<Erfc.from || er_p>=Erfc.to) Error("erud:range");*/
  if (er_p>=Erfc.to) Error("erud:range");
#endif
  byerd = er_p->Bd/(er_p->Cd+x) + er_p->Ad;
  return er_p->Bu/(er_p->Cu+x) + er_p->Au;
} /* eru */
#endif

#ifdef eropt
#error eropt not supported
#endif

static void eerfcs(double q, long ito) /***************************** eerfcs */
/* ito,ifrom in units of step */
{
  long x;
  double e,z,eh,zh,ee,zz,ed,zd,xh,xxh;
  double step=1e0/Erfc.sgrid, h=step/Erfc.ndiv;
  ertab_p p;

#ifdef eropt
  double sq;
  ertab_p er_p;
  long i;
#endif

  x=ito*Erfc.ndiv; /* in units of h */
  erfcs(x+1,&eh,&zh);
  erfcs(x  ,&e, &z);
  erfcs(x-1,&ee,&zz);
  /*.....ed=(eh-ee)/(2*h)*q;*/
  /*.....zd=(zh-zz)/(2*h)*q;*/
  ed=zd=0;
  xh=x*h;

  while (ito) { ito--;
    xxh=xh; ee=e; zz=z;
    erfcs(x=ito*Erfc.ndiv,&e,&z);
    xh=x*h;

    p=Erfc.tab+ito;

    p->Cu=(step*ed+e-ee)/Sqr(step);
    p->Bu=ed-2*xxh*p->Cu;
    p->Au=e-xh*p->Bu-xh*xh*p->Cu;
    ed= p->Bu+2*p->Cu*xh;
    
    p->Cd=(step*zd+z-zz)/Sqr(step);
    p->Bd=zd-2*xxh*p->Cd;
    p->Ad=z-xh*p->Bd-xh*xh*p->Cd;
    zd= p->Bd+2*p->Cd*xh;
  }

} /* eerfcs */

void initerfc(int Grid, double minr, double maxr, double cutoff, double alpha, double beta)
/* NEW: Grid<0 means suppressed info */
{
  alpha=fabs(alpha); /* negative = verbose */
{
  long ifrom,ito,i;
  int grid=abs(Grid),verbose=Grid>0;
#ifdef eropt
  double r0=0.0001,r,q,f,ff;
#endif
  double x,e,ee,z,zz;
  double h=1e0/(Erfc.ndiv*grid);
  double mine=0,maxe=0,minz=0,maxz=0,minre=0,maxre=0,minrz=0,maxrz=0;
  double Ad=beta, AuBd=beta, Bu=beta, Cud=1;
  ertab_p er_p;

  r1=cutoff;

  ifrom=(int)(Sqr(minr)*grid); ito=(int)(Sqr(maxr)*grid)+1;
  /* ifrom,ito are in units of 1/grid */

  if (Erfc.grid==grid && Erfc.minr==minr && Erfc.maxr>=maxr
      && Erfc.alpha==alpha && Erfc.beta==beta) {
    Erfc.to=Erfc.tab+ito; return; }

  Erfc.sgrid=(erreal)grid;
  Erfc.grid=grid; Erfc.minr=minr; Erfc.maxr=maxr;
  Erfc.alpha=alpha; Erfc.beta=beta;

  if (Erfc.tab!=NULL) free(Erfc.tab);
  alloc(Erfc.tab,i=ito*sizeof(ertab_t));

  prt("\n:::::: fgelst :::::: tab 0..%ld..%ld = %ld B  grid=%d", ifrom,ito,i,grid);

  if (verbose) prt("subgrid=%d  alpha=%f  beta=%f  range=[%f,%f)",
		   Erfc.ndiv,alpha,beta,minr,maxr);

  Erfc.from=Erfc.tab+ifrom; Erfc.to=Erfc.tab+ito;

  eerfcs(1.0,ito);

  /***** check accuracy *****/
  ifrom *= Erfc.ndiv; ito *= Erfc.ndiv;
  /* ifrom,ito are in units of h */
  for (i=ito-1; i>=ifrom; i--) {
    erfcs(i,&e,&z);
    x=i*h; ee=eru(x)-e; zz=erd(x)-z;
    if (i==ito-1) 
      if (verbose) prt("eru(%f)=%e  erd(%f)=%e",x,e*AuBd,x,z*Ad);
    if (ee<mine) mine=ee; if (ee>maxe) maxe=ee;
    if (zz<minz) minz=zz; if (zz>maxz) maxz=zz;
    if (e==0 || z==0) continue;
    ee /= e; zz /= z;
    if (ee<minre) minre=ee; if (ee>maxre) maxre=ee;
    if (zz<minrz) minrz=zz; if (zz>maxrz) maxrz=zz; }

  if (verbose) {
    prt("abs err:  eru in (%.3e,%.3e)  erd in (%.3e,%.3e)",
	mine*AuBd,maxe*AuBd, minz*Ad,maxz*Ad);
    prt("rel err:  eru in (%.3e,%.3e)  erd in (%.3e,%.3e)",
	minre,maxre, minrz,maxrz); }

#ifdef erminmax
  Erfc.min=1e30,Erfc.max=0;
#endif

/*
  if (fabs(minre) > maxre) maxre = fabs(minre);
  if (fabs(minrz) > maxre) maxre = fabs(minrz);
  if (     maxrz  > maxre) maxre = maxrz;
*/

/* rescaling */
  loop (er_p,Erfc.tab,Erfc.to) {
    er_p->Ad *= Ad;   er_p->Au *= AuBd;
    er_p->Bd *= AuBd; er_p->Bu *= Bu;
    er_p->Cd *= Cud;  er_p->Cu *= Cud; }

#include "ertest.c"
}
} /* initerfc */

#ifdef erexact_eps
#error erexact_eps (#undef EXACTERUD)
#endif
