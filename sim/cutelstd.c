/*
  Short-range electrostatic approximation, DIRECT version (no splines)
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

#if COULOMB!=0
#  error COULOMB==0 expected
#endif

struct Erfc_s Erfc = { 0,0 };

erreal byerd; /* by-product result of erud
  NOTE: if ermacro is not #defined then this byerd is used,
  if #ifdef ermacro then a local byerd may be faster, however,
  you must not declare local byerd if ermacro is not #defined ! */

void initerfc(int Grid, double minr, double maxr, double cutoff, double alpha, int shifted)
/* NEW: Grid<0 means suppressed info */
{
  int verbose=alpha<0;

  alpha=fabs(alpha);

  if (alpha==0) alpha=0.7;
  if (alpha>=1) ERROR(("alpha=%g is impossible for cutelst",alpha))

  Erfc.r0=cutoff*alpha;

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

  put2(Erfc.r0,cutoff)
  put3(Erfc.A,Erfc.B,Erfc.shift)

  if (verbose) {
#define xput(X) prt("%s = %.15g",#X,X)
    xput(cutoff);
    xput(Erfc.shift);
    xput(Erfc.A3);
    xput(Erfc.B3);
    xput(Erfc.A);
    xput(Erfc.B);
    xput(Erfc.r0);
#undef xput
  }

#if 0
    {
      ermacrodcl
      double r;
      for (r=1; r<cutoff; r+=0.001) {
	double byerd,rr=r*r;
	double u=erud(rr);
	prt("%.15g %.15g %.15g",er_x,u,byerd*er_x); }
      exit(0); 
    }
#endif

} /* initerfc */

#ifdef erexact_eps
#error erexact_eps (#undef EXACTERUD)
#endif
