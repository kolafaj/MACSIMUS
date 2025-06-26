/*
  This is the Fennel, Gezelter variant of cutelst
  see:
    Fennell CJ, Gezelter JD
    Is the Ewald summation still necessary? Pairwise alternatives to the accepted standard for long-range electrostatics
    JCP 124, 234104 (2006)
  Approximated by quadratic splines
*/

#if COULOMB!=3
#  error COULOMB==3 expected
#endif

#define ERFC Erfc
struct Erfc_s Erfc;

erreal byerd; /* by-product result of erud
  NOTE: if ermacro is not #defined then this byerd is used,
  if #ifdef ermacro then a local byerd may be faster, however,
  you must not declare local byerd if ermacro is not #defined ! */

static void erfcs(long to,double *e,double *z) /********************** erfcs */
/*
  e=eru(h*to) and z=erd(h*to) are computed.
  Must be called with decreasing sequence of parameters `to`, starting
  with a sufficiently high value because the integral is for this to
  approximated by the asymptotic formula
*/
{
  double r=sqrt(Erfc.h*to);

  if (r<Erfc.r1) {
    if (r==0) r=sqrt(Erfc.h)/2; /* patch for /0 */
    *e=erfc(Erfc.alpha*r)/r
      -erfc(Erfc.alpha*Erfc.r1)/Erfc.r1
      +(erfc(Erfc.alpha*Erfc.r1)/Sqr(Erfc.r1)+2*Erfc.alpha/sqrt(PI)
	*exp(-Sqr(Erfc.alpha*Erfc.r1))/Erfc.r1)*(r-Erfc.r1);
    *z=((erfc(Erfc.alpha*r)/Sqr(r)+2*Erfc.alpha/sqrt(PI)*exp(-Sqr(Erfc.alpha*r))/r)
	-(erfc(Erfc.alpha*Erfc.r1)/Sqr(Erfc.r1)+2*Erfc.alpha/sqrt(PI)*exp(-Sqr(Erfc.alpha*Erfc.r1))/Erfc.r1))/r; }
  else
    *e=*z=0;
} /* erfcs */

#include "splines.c"

                                                                 /* initerfc */
void initerfc(int grid, double minr, double maxr, double cutoff, double alpha)
{
  int dump=alpha<0;
  int verbose=grid>0;
  long ifrom,ito,i;
  //  ertab_p er_p;

  grid=abs(grid);
  alpha=fabs(alpha);

  if (!grid) ERROR(("grid=0 in initerfc"))

  if (Erfc.grid==grid
      && Erfc.minr==minr && Erfc.maxr==maxr
      && Erfc.r1==cutoff
      && Erfc.alpha==alpha) {
    prt(":::::: Fennell-Gezelter/erfc already initialized ::::::");
    // Erfc.maxr>=maxr :  Erfc.to=Erfc.tab+ito;
    return; }

  Erfc.r1=cutoff;

  Erfc.sgrid=(erreal)grid; // cf. makesplines(ito,SCALING);
  Erfc.grid=grid;
  Erfc.minr=minr; Erfc.maxr=maxr;
  Erfc.alpha=alpha;

  /* ifrom,ito are in units of 1/grid */
  ifrom=(int)(Sqr(minr)*grid);
  ito=(int)(Sqr(maxr)*grid)+1;

  if (Erfc.tab!=NULL) free(Erfc.tab);
  alloc(Erfc.tab,i=ito*sizeof(ertab_t));

  prt("\n:::::: fgelst :::::: tab 0..%ld..%ld = %ld B  grid=%d", ifrom,ito,i,grid);
  if (verbose) prt("alpha=%g range=[%f,%f)",alpha,Erfc.minr,Erfc.maxr);

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

#if 0
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

  prt("abs err:  eru in (%.3e,%.3e)  erd in (%.3e,%.3e)",
      mine*AuBd,maxe*AuBd, minz*Ad,maxz*Ad);
  prt("rel err:  eru in (%.3e,%.3e)  erd in (%.3e,%.3e)",
      minre,maxre, minrz,maxrz);
#endif

} /* initerfc */
