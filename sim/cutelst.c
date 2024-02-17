/*
  Short-range electrostatic approximation, deprecated quadratic spline version
  (The k-space part should be turned off)
  Should conform elst.h, (ERFC not #defined)

  eru(r*r) = 1/r - shift,                             r<alpha*cutoff
           = (r-r1)^3 (A*FA(r)+B*FB(r)), alpha*cutoff<r<cutoff
           = 0,                                cutoff<r

  recommended alpha=0.7
  shift,A,B are so that the second derivative of eru is continuous
  FA and FB are appropriate base functions

  (old version cutelst0.c)
*/

#define FA(x) (1)
#define FAd(x) (0)  /* 1st derivative */
#define FAdd(x) (0) /* 2nd derivative */
#define FB(x) (x)
#define FBd(x) (1)  /* 1st derivative */
#define FBdd(x) (0) /* 2nd derivative */

#if COULOMB!=2
#  error COULOMB==2 expected
#endif

struct Erfc_s Erfc = { 0,0 };

erreal byerd; /* by-product result of erud
  NOTE: if ermacro is not #defined then this byerd is used,
  if #ifdef ermacro then a local byerd may be faster, however,
  you must not declare local byerd if ermacro is not #defined ! */

static double r0,r1,shift,A,B;
static int ir0,ir1;

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

  if (r<r0) {
    if (r==0) r=sqrt(h)/2; /* patch for /0 */
    *e=1/r-shift;
    *z=1/Cub(r); }
  else if (r<r1) {
    *e=Cub(r-r1)*(A*FA(r)+B*FB(r));
    *z=(3*Sqr(r-r1)*(A*FA(r)+B*FB(r))+Cub(r-r1)*(A*FAd(r)+B*FBd(r)))/-r; }
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

  alpha=fabs(alpha); /* negative = verbose */

  if (alpha==0) alpha=0.7;
  if (alpha<0 || alpha>=1) ERROR(("alpha=%g is impossible for cutelst",alpha))

  r1=cutoff;
  r0=cutoff*alpha;

  ifrom=(int)(Sqr(minr)*grid); ito=(int)(Sqr(maxr)*grid)+1;
  /* ifrom,ito are in units of 1/grid */

  ir0=(int)(Sqr(r0)*grid);
  ir1=(int)(Sqr(r1)*grid);
  r1=sqrt(ir1/(double)grid);
  r0=sqrt(ir0/(double)grid);
  prt("cutoff 0..%g..%g rounded to 0..%g..%g (int: %d %d)",cutoff*alpha,cutoff,r0,r1,ir0,ir1);

  if (Erfc.grid==grid && Erfc.minr==minr && Erfc.maxr>=maxr
      && Erfc.alpha==alpha && Erfc.beta==beta) {
    Erfc.to=Erfc.tab+ito; return; }

  Erfc.sgrid=(erreal)grid;
  Erfc.grid=grid; Erfc.minr=minr; Erfc.maxr=maxr;
  Erfc.alpha=alpha; Erfc.beta=beta;

  if (Erfc.tab!=NULL) free(Erfc.tab);
  alloc(Erfc.tab,i=ito*sizeof(ertab_t));

  prt("\n:::::: erfc :::::: tab 0..%ld..%ld = %ld B  grid=%d", ifrom,ito,i,grid);

  if (verbose) prt("subgrid=%d  alpha=%f  beta=%f  range=[%f,%f)",
		   Erfc.ndiv,alpha,beta,minr,maxr);
  
  Erfc.from=Erfc.tab+ifrom; Erfc.to=Erfc.tab+ito;

  {
    double d=r0-r1;
    double A1=3*d*d*FA(r0)+d*d*d*FAd(r0);
    double B1=3*d*d*FB(r0)+d*d*d*FBd(r0);
    double rhs1=-1/r0/r0;
    double A2=6*d*FA(r0)+6*d*d*FAd(r0)+d*d*d*FAdd(r0);
    double B2=6*d*FB(r0)+6*d*d*FBd(r0)+d*d*d*FBdd(r0);
    double rhs2=2/r0/r0/r0;
    double det=A1*B2-A2*B1;
    A=(rhs1*B2-rhs2*B1)/det;
    B=(A1*rhs2-A2*rhs1)/det;
    shift=1/r0-d*d*d*(A*FA(r0)+B*FB(r0));
    put3(r0,r1,shift) }

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

} /* initerfc */

#ifdef erexact_eps
#error erexact_eps (#undef EXACTERUD)
#endif
