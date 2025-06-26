/*
  Short-range electrostatic approximation provided in the same style as Ewald
  r-space functions via splines, cf. cutelstd.c
  The k-space part should be turned off.

  Should conform elst.h, (ERFC not #defined)

  eru(r*r) = 1/r - shift,                             r<alpha*cutoff
           = (r-r1)^3 (A*FA(r)+B*FB(r)), alpha*cutoff<r<cutoff
           = 0,                                cutoff<r

  recommended alpha=0.7
  shift,A,B are so that the second derivative of eru is continuous
  FA and FB are appropriate base functions

  see also erfc.c

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
#endif /*# COULOMB!=2 */

#define ERFC Erfc
struct Erfc_s Erfc;

static void erfcs(long to,double *e,double *z) /********************** erfcs */
/*
  e=eru(h*to) and z=erd(h*to) are computed.
  Must be called with decreasing sequence of parameters `to`, starting
  with a sufficiently high value because the integral is for this to
  approximated by the asymptotic formula
*/
{
  double r=sqrt(Erfc.h*to);

  if (r<Erfc.r0) {
    if (r==0) r=sqrt(Erfc.h)/2; /* patch for /0 */
    *e=1/r-Erfc.shift;
    *z=1/Cub(r); }
  else if (r<Erfc.r1) {
    *e=Cub(r-Erfc.r1)*(Erfc.A*FA(r)+Erfc.B*FB(r));
    *z=(3*Sqr(r-Erfc.r1)*(Erfc.A*FA(r)+Erfc.B*FB(r))+Cub(r-Erfc.r1)*(Erfc.A*FAd(r)+Erfc.B*FBd(r)))/(-r); }
  else
    *e=*z=0;
} /* erfcs */

#include "splines.c"

                                                                 /* initerfc */
void initerfc(int grid, double minr, double maxr, double cutoff, double alpha)
{
  int dump=alpha<0;
  int verbose=grid>0;
  long ifrom,ito;
  //  ertab_p er_p;

  grid=abs(grid);
  alpha=fabs(alpha);

  if (!grid) ERROR(("grid=0 in initerfc"))
  if (alpha==0) alpha=0.7;
  if (alpha<0 || alpha>=1) ERROR(("alpha=%g is impossible for cutelst",alpha))

  if (Erfc.grid==grid
      && Erfc.minr==minr && Erfc.maxr==maxr
      && Erfc.r1==cutoff
      && Erfc.alpha==alpha) {
    prt(":::::: cutelst (aka erfc) already initialized ::::::");
    return; }

  Erfc.r0=cutoff*alpha;
  Erfc.r1=cutoff;

  {
    double d=Erfc.r0-cutoff;
    double A1=3*d*d*FA(Erfc.r0)+d*d*d*FAd(Erfc.r0);
    double B1=3*d*d*FB(Erfc.r0)+d*d*d*FBd(Erfc.r0);
    double rhs1=-1/Sqr(Erfc.r0);
    double A2=6*d*FA(Erfc.r0)+6*d*d*FAd(Erfc.r0)+d*d*d*FAdd(Erfc.r0);
    double B2=6*d*FB(Erfc.r0)+6*d*d*FBd(Erfc.r0)+d*d*d*FBdd(Erfc.r0);
    double rhs2=2/Cub(Erfc.r0);
    double det=A1*B2-A2*B1;

    Erfc.A=(rhs1*B2-rhs2*B1)/det;
    Erfc.B=(A1*rhs2-A2*rhs1)/det;
    Erfc.A3=-3*Erfc.A;
    Erfc.B3=-3*Erfc.B;
    Erfc.shift=1/Erfc.r0-d*d*d*(Erfc.A*FA(Erfc.r0)+Erfc.B*FB(Erfc.r0)); }

  put3(minr,maxr,cutoff)
  put3(alpha,Erfc.r0,Erfc.r1)
  put3(Erfc.A,Erfc.B,Erfc.shift)

  Erfc.sgrid=(erreal)grid; // cf. makesplines(ito,SCALING);
  Erfc.grid=grid;
  Erfc.minr=minr; Erfc.maxr=maxr;
  Erfc.alpha=alpha;

  /* ifrom,ito are in units of 1/grid */
  ifrom=(int)(Sqr(minr)*grid);
  ito=(int)(Sqr(maxr)*grid)+1;

  if (Erfc.tab!=NULL) free(Erfc.tab);
  allocarray(Erfc.tab,ito);

  prt("\n:::::: cutelst by splines :::::: tab 0..%ld..%ld = %ld B  grid=%d", ifrom,ito,ito*sizeof(ertab_t),grid);

  if (verbose) {
    prt("alpha=%g  range=[%g,%g)  cutoff=%g  1st interval=[0,%g]",
        alpha,minr,maxr,cutoff,sqrt(1/Erfc.sgrid)); }

  makesplines(ito,1);

  if (dump) {
    int i;
    double r;

    prt("!r uref fref");
    for (i=ito; i>0; i--) {
      double rr=i*Erfc.h;
      double uu,dd;

      erfcs(i,&uu,&dd);
      prt("%.10g %.10g %.10g ELST1",sqrt(rr),uu,dd); }
    
    prt("\n!r uspline fspline");
    for (r=maxr-0.01; r>0.02; r-=0.01) {
      ermacrodcl
      double u=erud(r*r);

      prt("%.10g %.10g %.10g ELST2",r,u,byerd); } }

} /* initerfc */

#ifdef erexact_eps
#  error erexact_eps does not make sense for cutelst (#undef EXACTERUD)
#endif /*# erexact_eps */
