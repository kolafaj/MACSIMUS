/*
initialization of site-site potentials
site-site potential interfaced through sitesite.c and sitesite.h

this module depends on #defines WORM, CUT, TWODIM, POLAR, PARALLEL

NOTE: initcombrule() and combrule() are in sitesite.c

Do not change the header of setss !!!

NEW: rewritten for internal use of sigvdW,EvdW
  sigvdW=potential minimum distance (van der Waals diameter)
  epsvdW=potential minimum, should be NEGATIVE
For LJ, it holds sigvdW=LJsigma*2^(1/6), epsvdW=-4*LJeps
*/

#include "ground.h"
#include "simglob.h"

#include "interpot.h"
#ifdef TWODIM
#  include "inter2d.h"
#else /*# TWODIM */
#  include "intermac.h"
#endif /*#!TWODIM */
#include "simdef.h"
#include "units.h"

#include "setss.h"

double poteps=-1e-5;
/*
  Internal test of the potential and forces, incl. smoothed cutoff,
  the value of poteps is some relative error.
  Occasionally (functions close to zero), an error may be falsely indicated.
  Use option -v4 to get more info!
  former switch: #define SS_DEBUG 1e-5
*/

#ifndef SS_MEASURE_rep
#  define SS_MEASURE_rep /**/
#endif /*# SS_MEASURE_rep */

#ifndef SS_NOMEASURE_rep
#  define SS_NOMEASURE_rep /**/
#endif /*# SS_NOMEASURE_rep */

/* for calculating cutoff corrections */
#define NGauss4 64

#ifdef WIDOM
static struct widomrdf_s {
  struct widomrdf_s *next;
  rdf_t **addr;
  rdf_t *rdf;
} *widomrdfhead;

void widomrdf(int mode) /****************************************** widomrdf */
/*
  mode=0: rdf = NULL (no calculation of rdf in site-site functions)
  mode=1: reset back (calculate rdf in site-site functions)
  mode=-1: purge the list (ralloc+release cannot be used - why???)
*/
{
  static struct widomrdf_s *wr,*wrnext;

  for (wr=widomrdfhead; wr; wr=wrnext) {
    wrnext=wr->next;
    switch (mode) {
      case -1: free(wr); break;
      case 0: *(wr->addr)=NULL; break;
      case 1: *(wr->addr)=wr->rdf; break;
      default: ERROR(("widomrdf: mode")) } }

  if (mode==-1) widomrdfhead=NULL;
}
#endif /*# WIDOM */

#ifdef METAL
static double rho1=1,rho2=1,rho12=0; /* for potential test */
#endif /*# METAL */

/*
   SS_U and/or may be #defined in XXX/sitesite.h as simple
   (unoptimized) variants of the site-site potential energy/force as a
   function of distance, using the tabulated parameters (pairparm_t).
   If not #defined, they are translated from (optimized) versions here.
*/

#ifdef POLAR
#  if POLAR&1
#    ifndef SS_Urep
#      define SS_Urep(R) ss_urep(ss,R)
static double ss_urep(sitesite_t *ss,double r)
{
  double x,y,z,U=0,f=0,rr=r*r;
  double Urep,frep;
  SS_MEASURE_rep
  return Urep;
}
#    endif /*# SS_Urep */
#  endif /*# POLAR&1 */
#endif /*# POLAR */

#ifndef SS_U
#  define SS_U(R) ss_u(ss,R)
static double ss_u(sitesite_t *ss,double r)
{
  double x,y,z,U=0,f=0,rr=r*r;
#  ifdef POLAR
#    if POLAR&1
  double Urep,frep;
  SS_MEASURE_rep
#    endif /*# POLAR&1 */
#  endif /*# POLAR */
  SS_MEASURE
  return U;
}
#endif /*# SS_U */

#ifndef SS_F
#  define SS_F(R) ss_f(ss,R)
static double ss_f(sitesite_t *ss,double r)
{
  double x,y,z,U=0,f=0,rr=r*r;
#  ifdef POLAR
#    if POLAR&1
  double Urep,frep;
  SS_MEASURE_rep
#    endif /*# POLAR&1 */
#  endif /*# POLAR */
  SS_MEASURE
  return f*r;
}
#endif /*# SS_F */

void setss(sitesite_t *ss, int i,int j, double C2,int onefour) /****** setss */
/***
  Table *ss of potential constants for a pair (i,j) (where i and j are
  site types) is initialized, incl. the normalized cutoff correction.

  Normally, C2=LJcutoff as given in the input data.
  If C2<0 then C2=abs(C2)*sigvdW is used (i.e., -C2 is in units of
  potential minimum sigvdW)

  The pair potential U(r) is approximated by:
    u(r) = U(r)              for r<C1
    u(r) = A*(r^2-C2^2)^2    for C1<r<C2
    u(r) = 0                 for C2<r
  The pair force F(r)/r is approximated by:
    f(r)/r = F(r)/r          for r<C1
    f(r)/r = -4*A*(r^2-C2^2) for C1<r<C2
    f(r)/r = 0               for C2<r
  (Note that the vector of force = f(r)/r*(vector r), so f(r)/r is needed)
  A,C1 are calculated from the max. cutoff C2 so that the forces are continuous.

  The lowest possible value for LJ is C2/sigvdW=1.59 (C2/LJsigma=1.7819),
  if shorter C2 is requested, the calculations are inaccurate and a WARNING
  is printed

  The cutoff corrections are calculated using the assumption that
  g(r)=1 for r>C1.  The normalized cutoff correction corr for a pair
  is returned, then the cutoff correction of total energy is sum over
  pairs of corr/V and the cutoff correction of pressure is sum over
  pairs of corr/V^2.  (It follows from integration by parts that the
  factor corr is the same for both pressure and energy.)

  The code is mostly potential-independent (the dependency by macros
  U,F only).

  onefour & 2 suppresses verbose debugging
***/
{
  int n,isfixed=0;
  double C1, oldC1,r,h,corr0,corr1,x;
  siteparm_t *lj0,*lj1;
  pairparm_t pp;
  nbfix_t *fix;
  int nerr=0,ninconsis=0;
  int warnonly=onefour&2;
  int verbose=!warnonly && option('v')&4;

  box.V=PROD(box.L); /* should be elsewhere... */

  onefour=onefour&1;

  lj0=&sitedef[i].LJ[onefour];
  lj1=&sitedef[j].LJ[onefour];
  /* try NBFIX first */
  looplist (fix,nbfix0)
    if ( (fix->indx[0]==i && fix->indx[1]==j)
     ||  (fix->indx[0]==j && fix->indx[1]==i) ) {
      pp=fix->onefour[onefour];
      lj0=lj1=NULL; /* combrule() just prints info */
      isfixed++;
      goto fixed; }

 fixed:;
  /* verbose level changed 4->2; NEW: lj0==NULL just prints info */
  combrule(&pp,lj0,lj1,
           option('v')&2 ? string("%s%s %s-%s",lj0?"comb":"fix",onefour?"14":"",sitedef[i].name,sitedef[j].name) : NULL);

  /* site-site specific support incl. equivalent wall-atom potential */
  initssaux(&ss->a,&pp);

  /* if (!rdf) Error("rdf not initialized"); */
  if (rdf) ss->rdf = i<j ? rdf[i][j] : rdf[j][i];
  else ss->rdf=NULL;
  if (onefour && !onefourinrdf) ss->rdf=NULL;

#ifdef WIDOM
  if (ss->rdf) {
    struct widomrdf_s *wr;

    allocone(wr);
    wr->addr=&ss->rdf;
    wr->rdf=ss->rdf;
    wr->next=widomrdfhead;
    widomrdfhead=wr; }
#endif /*# WIDOM */

  if (option('v')&4)
    prt("setss: %i %i  epsvdW=%.3f (%.5f kcal/mol) sigvdW=%.5f cutoff=%.4f %s",
        i,j,pp.eps,pp.eps*(Eunit/kcal),pp.sig,C2,isfixed?"nbfix":"comb.r.");

#ifndef NIBC
  if (pp.eps==0 || C2==0)
    C1=C2=ss->A=0;
  else {
    /* it is assumed that sigvdW has meaning of (a sort of) atom diameter */
    if (C2<0) C2= -pp.sig*C2;
    if (box.cutoff>0 && C2>box.cutoff)
      WARNING(("%d-%d: LJcutoff=%g > cutoff=%g",i,j,C2,box.cutoff))

#  ifdef SHARPCUTOFF
    C1=C2; ss->A=0;
    // extremely dirty: very special project only
    globalLJshift=Pow6(pp.sig/C2)/2;
    globalLJshift=globalLJshift-Sqr(globalLJshift);
#  else /*# SHARPCUTOFF */
    /* smooth cutoff
       calculates C1 from C2 - this part is potential-independent */

    C1=0.775*C2;
    n=0;
    do {
      double den=SS_F(C1);

      if (fabs(den)<1e-50) {
        prt("! %d-%d: sigvdW=%g C1=%g force %g, smooth cutoff off",
            i,j,pp.sig,C1,den);
        C1=-C1; break; }
      oldC1=C1;
      C1=1+4*SS_U(C1)/C1/den;
      if (n++>2000 || C1<0) {
        C1=0.775*C2;
        ninconsis++;
#    ifdef GOLD
        if (i!=nsites-1 || j!=nsites-1)
#    endif       /*# GOLD */
          WARNING(("%d-%d: sigvdW=%g C2=LJcutoff=%g\n\
*** Cannot determine consistent smooth cutoff, C1=0.775*C2 will be used:\n\
*** - The energy conservation precision will decrease.\n\
*** - (This component to) the homogeneous cutoff correction in pressure will\n\
***   be by less than %.4g %% off the calculated one.\n\
*** - The internal energy/force test below will not pass, but will be reported\n\
***   as a warning only.\n\
*** DO NOT IGNORE THIS MESSAGE unless this interaction is never calculated.\n\
*** NOTE: For LJ, the minimum LJcutoff=%.6f",
                   i,j,pp.sig,C2,
                   (-SS_F(C1)/(C1*SS_U(C1)*(C2*C2-C1*C1)))-1,
                   pp.sig*1.58741609+5e-7))
/*.....      epsvdW=0; goto again;*/
        break; }
      C1=C2/sqrt(C1);
    } while (fabs(1-C1/oldC1) > No.eps*n);

    if (C1>0) if (C1>0.999*C2 || C1<0.5*C2)
#    ifdef GOLD
      if (i!=nsites-1 || j!=nsites-1)
#    endif /*# GOLD */
        ERROR(("%d-%d: sigvdW=%f C1=%f C2=%f bad cutoff",i,j,pp.sig,C1,C2))

/* potential-independent constants */
      ss->A=SS_U(C1)/Sqr(C2*C2-C1*C1);
#  endif /*#!SHARPCUTOFF */
  }

  if (C1<0) { C1=-C1; ss->A=ss->A4=0; }
  else ss->A4= -SS_F(C1)/(C2*C2-C1*C1)/C1;

  ss->C2q=C2*C2; ss->C1q=C1*C1;
#endif /*# NIBC */

#ifdef POTCONST
  prt("ss->A = %.15g",ss->A);
  prt("ss->A4 = %.15g",ss->A4);
  prt("C1 = %.15g",C1);
  prt("C2 = %.15g",C2);
#endif /*# POTCONST */

  if (poteps) {
    /* test of macros defining site-site energies and forces */
    if (pp.sig!=0 && pp.eps!=0) {
      double R,rr,y,z,f,q,D,Fmac;
      int k,nplus;
      static double d[3]={-5e-6,5e-6,0}; /* delta r for numerical derivatives */
      double Umac[3],Umeas[3],fnmeas,dUmac,dUmeas,minU;
      double UX=0,fX=0;
#if POLAR&1
      double frep,Urep,Umeasrep[3],Umacrep[3];
#endif /*# POLAR&1 */

      if (verbose)
        prt("\nsetss: C1=%.14g C2=%.14g epsvdW=%.14g sigvdW=%.14g\n\
smoothing u=A*(r^2-C2^2)^2: A=%.14g",C1,C2,pp.eps,pp.sig,ss->A);
      q=pow(pp.sig/C2,0.03125);
      if (q>=0.99) q=0.99; /* patch... */

#if POLAR&1
      if (!onefour) {
        static int pass=0;

        if (verbose) {
          rr=ss->C1q;
          SS_MEASURE_rep
          SS_MEASURE
          prt("%d-%d repulsive cutoff at C1=%f: Urep=%g frep=%g",
               i, j,                     C1,    Urep,   frep); }
        rr=ss->C2q;
        SS_MEASURE_rep
        SS_MEASURE
        if (verbose || pass<4*Sqr(nsites)) {
          pass++;
          prt("%d-%d repulsive cutoff at C2=%f: Urep=%g frep=%g",
               i, j,                     C2,    Urep,   frep); } }
#endif /*# POLAR&1 */

      if (verbose)
        header("   r           SS_U     U(MEASURE)     SS_F     f(MEASURE)   f(NOMEAS)    -dSS_U/dr    -dU/dr  "
#if POLAR&1
"       frep_AB       SS_Urep   Urep(MEASURE)"
#endif /*# POLAR&1 */
             );

      //    for (R=pp.sig*pow(q,-8.5); R<C2*Sqr(q); R*=q) // orig with q>1
      nplus=0; minU=0;
      for (R=C2*pow(q,1.5); R>0.1 && nplus<10; R*=q) {
// version for the ucut plot:
//      for (R=C2*pow(q,1.5); R>0.1 && nplus<10; R*=q) {

        loop (k,0,3) {
          double U=0;

          r=R*(d[k]+1);
          rr=r*r;

          Umac[k]=SS_U(r);
#if POLAR&1
          Umacrep[k]=SS_Urep(r);
          Umeasrep[k]=0;
#endif /*# POLAR&1 */

#ifdef NIBC
          SS_NOMEASURE_rep
          SS_NOMEASURE(=) fnmeas=f;
          SS_MEASURE_rep
          SS_MEASURE
#else /*# NIBC */
          if (rr>=ss->C2q)
            U=f=fnmeas=0;
          else if (rr < ss->C1q) {
            SS_NOMEASURE_rep
            SS_NOMEASURE(=) fnmeas=f;
            SS_MEASURE_rep
            SS_MEASURE
#  if POLAR&1
            Umeasrep[k]=Urep;
#  endif /*# POLAR&1 */
          } else {
            x=rr-ss->C2q; U=x*x*ss->A; fnmeas=f=ss->A4*x; }
#endif /*#!NIBC */
          Umeas[k]=U; } /* k */

        D=R*(d[1]-d[0]);
        dUmac=(Umac[0]-Umac[1])/D;
        dUmeas=(Umeas[0]-Umeas[1])/D;
        Fmac=SS_F(r);
        f*=R;
        fnmeas*=R;

        Min(minU,Umac[2])
        Min(minU,Umeas[2])

        if (minU<0) if (Umac[2]>-minU && Umeas[2]>-minU) nplus++;
        if (Umac[2]<0 && Umeas[2]<0) nplus=0;

        if (verbose) prt_("\
%8.5f  %11.5f %11.5f  %11.5f %11.5f %11.5f  %11.5f %11.5f"
#if POLAR&1
            "  %11.5f %11.5f"
#endif /*# POLAR&1 */
,r,    Umac[2],Umeas[2],Fmac,  f,     fnmeas, dUmac, dUmeas
#if POLAR&1
                ,Umeasrep[2],Umacrep[2]
#endif /*# POLAR&1 */
            );

        UX=fabs(UX/3)+fabs(Umac[2])+fabs(Umeas[2]);
        fX=fabs(fX/3)+fabs(Fmac)+fabs(f)+fabs(fnmeas)+fabs(dUmac)+fabs(dUmeas);

        if ( (R<C1*sqrt(q)
              && (fabs(Umac[2]-Umeas[2])>Sqr(poteps)*UX
               || fabs(Fmac-f)>Sqr(poteps)*fX))
             || fabs(fnmeas-f)>Sqr(poteps)*fX
             || fabs(dUmac-Fmac)>fabs(poteps)*fX
             || !isfinite(Umeas[2]) || !isfinite(f) || !isfinite(Fmac) || !isfinite(dUmac)
             || fabs(dUmeas-f)>fabs(poteps)*fX ) {
          nerr++;
          if (verbose) prt(" %dPOT PROBLEM",onefour); }
        else
          if (verbose) prt(" %dPOT OK",onefour);
      } /* r */

      if (verbose) header("");

      if (nerr) {
        if (warnonly || nerr<4 || ninconsis) WARNING(("setss: %d potential/force errors detected\n\
*** try option -v4, see manual on variable poteps",nerr))
        else ERROR(("setss: %d errors detected\n\
*** try option -v4, see manual on variable poteps",nerr)) }
    else
      if (verbose) prt("setss: site-site tests passed"); }
  }

  /* The normalized cuttoff correction is corr0+corr1, where
                        C2        2   2 2  2
       corr0 = -4 PI integral A (r -C2 )  r  dr
                        C1
  and
                    infinity       2
       corr1 = 4 PI integral U(r) r  dr
                       C1
  corr0 is given analytically while corr1 is calculated by 4-order Gauss
  formula with NGauss4*2 points after substition to r=1/x.

  The code is potential-independent */

#ifdef NIBC
  /* in future, real cutoff correction should be calculated !!! */
  /* no longer the same correction for P and U !!! */
  /* see program outcube.c : */
  corr0=corr1=0;
#elif defined(SHARPCUTOFF)
  /* by #definition, no cutoff corrections */
  corr0=corr1=0;
#else /*? NIBC */ /*#!NIBC!defined(SHARPCUTOFF) */
  if (pp.eps==0)
    corr0=0;
  else {

    corr0=4*PI*ss->A*
      (C1*(Pow6(C1)/7-0.4*Pow4(C1)*Sqr(C2)+Sqr(C1)*Pow4(C2)/3)-8.0/105*powi(C2,7));
    /* original form
       corr0=4*PI*ss->A*( (powi(C1,7)-powi(C2,7))/7
       -0.4*Sqr(C2)*(powi(C1,5)-powi(C2,5))
       +Pow4(C2)*(powi(C1,3)-powi(C2,3))/3);
    */

    /* integral of r^2 U(r) dr numerically */
#  define Gauss4 0.2113248654051871

    h=1.0/NGauss4/C1;
    corr1=0;
    loop (n,0,NGauss4) {
      r=1/((n+Gauss4)*h); corr1+=Pow4(r)*SS_U(r);
      r=1/((n+(1-Gauss4))*h); corr1+=Pow4(r)*SS_U(r); }
    corr1*=2*PI*h;
#  undef Gauss4

    corr0+=corr1; }
#endif /*#!NIBC!defined(SHARPCUTOFF) */

#ifdef SLAB
#  define NINT 200
  if (ss->Skk) free(ss->Skk); /* NB: sstab must have been initialized to zeros */
  if (slab.K && !onefour) {
    int k;
    int NINTRG=(int)(NINT*slab.range+0.5);

    if (ss->A) {

      allocarrayzero(ss->Skk,slab.K);

      loop (k,0,slab.K) {
        int idz;

        loop (idz,0,NINTRG) {
          double dz=(idz+0.5)/NINT*box.L[2];
          double tou=1/fmax(Sqr(C1),Sqr(dz));
          int iu,nu=1+(int)(tou*NINT*Sqr(C1));
          double h=tou/nu; /* global h masked */
          double a=2*PI*k*dz/box.L[2];
          double cah=cos(a)*h;
          double sah=sin(a)*h*a;

          loop (iu,0,nu) {
            double u=(iu+0.5)*h;
            double rr=1/u;
            double y,z,f,U=0; /* for SS_MEASURE */
#  if POLAR&1
            ERROR(("not implemented (?)"))
#  endif /*# POLAR&1 */

#  ifdef NIBC
            ERROR(("not implemented"))
#  else /*# NIBC */
            if (rr<Sqr(C1)) ERROR(("internal"))
            SS_MEASURE /* f not used, but must be declared */
            if (rr<Sqr(C2)) {
              x=rr-Sqr(C2); U-=x*x*ss->A; }
#  endif /*#!NIBC */
            ss->Skk[k].E+=cah*U/Sqr(u);
            ss->Skk[k].D+=sah*U/Sqr(u); } }

        /* divided by V (rho in slabcutcor is multiplied by V instead)
           and multiplied by (1+(k>0))/2 */
        ss->Skk[k].E*=(1+(k>0))*PI/(box.L[0]*box.L[1]*NINT);
        ss->Skk[k].D*=(1+(k>0))*PI/(box.L[0]*box.L[1]*NINT);
        if (option('v')&4) prt("ss->Skk[%d] = %g %g",k,
                               ss->Skk[k].E*box.V,
                               ss->Skk[k].D*box.V); } /* slab.K */
    } /* if (ss->A) */
    else 
      ss->Skk=NULL;
  }
#endif /*# SLAB */

  ss->corr=corr0;

} /* setss */
