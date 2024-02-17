/* support for optimized water models */

#include "ground.h"
#include "simglob.h"
#include "interpot.h"
#include "intermac.h"
#include "water.h"

#ifdef COULOMB
#  include "elst.h"
#else /*# COULOMB */
#  define ermacrodcl /* empty */
#  undef ermacro
#endif /*#!COULOMB */

#define DEBUGU(X)

struct waters_s waters;

#define ENVIR En.vir
#define ENEL En.el
#define ENPVIR En.Pvir

#if PARALLEL==3
/* this module is NOT parallelized but still must be compatible */
#  define STARTPAR pll_global_t glob; memset(&glob,0,sizeof(glob));
#  define ENDPAR En.vir+=glob.Envir; En.el+=glob.Enel;
#else /*# PARALLEL==3 */
#  define STARTPAR
#  define ENDPAR
#endif /*#!PARALLEL==3 */

double TIP3P(vector f1[],vector f2[],molecule_t *m1,molecule_t *m2,vector *rp)
/* order of sites = HHO */
{
#ifdef POLAR
  ERROR(("sorry, optimized polar TIP3P not supported"))
  return 0; /* to avoid "control reaches end of non-void function" */
#else /*# POLAR */
  vector *r1=rof(m1,rp), *r2=rof(m2,rp);
  double U=0;
  STARTPAR

  if (m1!=m2) {
    /*** intermolecular interactions ***/
    if (measure)
#  if PARALLEL==3
#    define GLOB ,&glob
#  else /*# PARALLEL==3 */
#    define GLOB
#  endif /*#!PARALLEL==3 */
#  include "cookmeas.c"
#  include "tip3pmm.c"
#  include "cookundf.c"
#  undef GLOB
    else
#  define GLOB
#  include "cooknmea.c"
#  include "tip3pmm.c"
#  include "cookundf.c"
#  undef GLOB
  } else {
    /*** intramolecular interactions ***/
    /* all site-site interactions */

#  define I0 0
#  define I1 ns1

    if (measure)
#  if PARALLEL==3
#    define GLOB ,&glob
#  else /*# PARALLEL==3 */
#    define GLOB
#  endif /*#!PARALLEL==3 */
#  include "cookmeas.c"
#  include "tip3pm.c"
#  include "cookundf.c"
#  undef GLOB
    else
#  define GLOB
#  include "cooknmea.c"
#  include "tip3pm.c"
#  include "cookundf.c"
#  undef GLOB
/*
   removed in 2.7b: flexible TIP3P-like models no longer supported as optimized
   #include "bonded.c"
*/
  }

  ENDPAR

  return U;
#endif /*#!POLAR */
} /* TIP3P */

double TIP4P(vector f1[],vector f2[],molecule_t *m1,molecule_t *m2,vector *rp)
/* order of sites = HOMH */
{
#ifdef POLAR
  ERROR(("sorry, optimized polar TIP4P not supported"))
  return 0; /* to avoid "control reaches end of non-void function" */
#else /*# POLAR */
  vector *r1=rof(m1,rp), *r2=rof(m2,rp);
  double U=0;
  STARTPAR

  if (m1!=m2) {
  /*** intermolecular interactions ***/
    if (measure) {
#  if PARALLEL==3
#    define GLOB ,&glob
#  else /*# PARALLEL==3 */
#    define GLOB
#  endif /*#!PARALLEL==3 */
#  include "cookmeas.c"
#  include "tip4pmm.c"
#  include "cookundf.c"
#  undef GLOB
      /* measuring OH rdf */
      if (waters.ss.OH.rdf) {
	/* for sure, waters.ss.OH.rdf checked in rdfM, too */
	rdfM(&waters.ss.OH,r1[1],r2[0]);
	rdfM(&waters.ss.OH,r1[1],r2[3]);
	rdfM(&waters.ss.OH,r1[3],r2[1]);
	rdfM(&waters.ss.OH,r1[3],r2[1]); } }
    else
#  define GLOB
#  include "cooknmea.c"
#  include "tip4pmm.c"
#  include "cookundf.c"
#  undef GLOB
  } else {
    /*** intramolecular interactions ***/
    /* all site-site interactions */

#  define I0 0
#  define I1 ns1

    if (measure)
#  if PARALLEL==3
#    define GLOB ,&glob
#  else /*# PARALLEL==3 */
#    define GLOB
#  endif /*#!PARALLEL==3 */
#  include "cookmeas.c"
#  include "tip4pm.c"
#  include "cookundf.c"
#  undef GLOB
      else
#  define GLOB
#  include "cooknmea.c"
#  include "tip4pm.c"
#  include "cookundf.c"
#  undef GLOB

/*
  removed in 2.7b: flexible TIP4P-like models no longer supported as optimized
  #include "bonded.c"
*/
  }

  ENDPAR

  return U;
#endif /*#!POLAR */
} /* TIP4P */

/* ST2 */
#define rL 2.0160
#define rU 3.1287

#define CUB Cub(rU-rL)

/* this is switch-1 function because we SUBTRACT from Ewald */
static double Switch_1(double r) /********************************* Switch_1 */
{
  if (r<rL) return -1;
  else if (r>rU) return 0;
  else return Sqr(r-rL)*( (3*rU-rL)/CUB - (2/CUB)*r) - 1;
}

static double dSwitch_1(double r) /******************************* dSwitch_1 */
{
  if (r<rL) return 0;
  else if (r>rU) return 0;
  else return (6/CUB)*(r-rL)*(r-rU);
}

#undef CUB

#include "waterni.c"

double ST2(vector f1[],vector f2[],molecule_t *m1,molecule_t *m2,vector *rp)
/* order of sites=OLHLH */
/* only rigid model */
{
  double U=0;

#if PARALLEL==3 || defined(POLAR)
  ERROR(("sorry, ST2 not supported"))
#else /*# PARALLEL==3 || defined(POLAR) */

  int i,j;
  double fdr,qq=waters.qq.HH;
#  ifdef GOLD
  double xxyy;
#  endif /*# GOLD */
  vector *R1=rof(m1,rp), *R2=rof(m2,rp);
  double *r1=R1[0],*r2=R2[0];
  sitesite_t *ss=&waters.ss.OO;

  if (m1==m2)
    if (measure) {
      static double U0;
#  ifdef COULOMB
      static real alpha0=-1;

      if (fabs(el.alpha-alpha0)>1e-12) {
	double tau=1.910633236249019,r;
	ermacrodcl

        r=2*sin(tau/2); U0=eru(r*r)-1/r;
	r*=0.8; U0+=eru(r*r)-1/r;
	r=sqrt(sqr(1-0.8*cos(tau))+sqr(0.8*sin(tau))); U0-=4*(eru(r*r)-1/r);
	U0*=qq;
	alpha0=el.alpha; }
#  endif /*# COULOMB */
      U=U0; }
    else;
  else /* m1!=m2 */
    if (measure) {
      /*** intermolecular interactions ***/
      double rr,f,x,y,z,sw,dsw,srr;
      rdf_t *ssrdf;
      ermacrodcl
      vector dr;

      beginCellImage(0,box.cq)
	f=0;
	srr=sqrt(rr);
	sw=Switch_1(srr);
	dsw=dSwitch_1(srr)/srr;
	if ( (ssrdf=ss->rdf) ) {
	  int i=(int)(srr*ssrdf->grid);
	  if (i<ssrdf->nhist) (ssrdf->hist[i])++; }
#  ifdef NIBC
	SS_MEASURE
        ENVIR+=rr*f;
#  else /*# NIBC */
	if (rr < ss->C2q) {
	  if (rr < ss->C1q)
	    SS_MEASURE
          else {
	    x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
	  ENVIR+=rr*f; }
#  endif /*#!NIBC */

	/* O..H rdf must be measured separately */
	if ( (ssrdf=waters.ss.OH.rdf) ) {
	  rdfM(&waters.ss.OH,R1[0],R2[2]);
	  rdfM(&waters.ss.OH,R1[0],R2[4]);
	  rdfM(&waters.ss.OH,R1[2],R2[0]);
	  rdfM(&waters.ss.OH,R1[4],R2[0]); }

	loop (i,1,5) { qq=-qq; loop (j,1,5) {
	  vector ddr;
	  double ff,RR,sRR,xx;

	  qq=-qq;
	  VVVV(ddr,=R1[i],-R2[j],+dni)
          RR=SQR(ddr);
	  sRR=sqrt(RR);
	  if ( (ssrdf=waters.ss.HH.rdf) ) if ((i&1)+(j&1)==0) {
	    int i=(int)(sRR*ssrdf->grid);
	    if (i<ssrdf->nhist) (ssrdf->hist[i])++; }
	  ff=qq/sRR;
	  ENEL+=ff*sw;
	  f+=dsw*ff;
	  ff*=sw/RR; /* subtracted switch part */
#  ifdef COULOMB
	  ENEL+=erud(RR)*qq; ff+=byerd*qq;
#  else /*# COULOMB */
	  ENEL+=xx=qq/sRR; ff+=xx/RR;
#  endif /*#!COULOMB */
	  MADDFORCES(f1[i],f2[j],ff,ddr) } }

	MADDFORCES(f1[0],f2[0],f,dr)
      endCellImage
    }
    else { /* !measure */
      /*** intermolecular interactions ***/
      double rr,f,x,y, /*z,*/ sw,dsw,srr;
      ermacrodcl
      vector dr;

      beginCellImage(0,box.cq)
	srr=sqrt(rr);
	sw=Switch_1(srr);
	dsw=dSwitch_1(srr)/srr;
#  ifdef NIBC
	SS_NOMEASURE(=)
#  else /*# NIBC */
        f=0;
	if (rr < ss->C2q) {
	  if (rr < ss->C1q)
	    SS_NOMEASURE(+=)
        else {
          x=rr-ss->C2q; f+=ss->A4*x; } }
#  endif /*#!NIBC */

	loop (i,1,5) { qq=-qq; loop (j,1,5) {
	  vector ddr;
	  double ff,RR,sRR;

	  qq=-qq;
	  VVVV(ddr,=R1[i],-R2[j],+dni)
          RR=SQR(ddr);
	  sRR=sqrt(RR);
	  ff=qq/sRR;
	  f+=dsw*ff;
	  ff*=sw/RR; /* subtracted switch part */
#  ifdef COULOMB
	  ff+=erd(RR)*qq;
#  else /*# COULOMB */
	  ff+=qq/(RR*sRR);
#  endif /*#!COULOMB */
	  ADDFORCES(f1[i],f2[j],ff,ddr) } }

	ADDFORCES(f1[0],f2[0],f,dr)
      endCellImage
    }

#endif /*#!PARALLEL==3 || defined(POLAR) */
  return U;
}

#ifdef WATERPLUS
#  if WATERPLUS==2
/* optimized codes for MeOH, acetone, CO2
   WARNING: undocumented feature */
#    include "xxxplus.c"
#  endif /*# WATERPLUS==2 */
#endif /*# WATERPLUS */
