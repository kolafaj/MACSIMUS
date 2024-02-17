/*
  polar version of internp.c

  WARNING: #if POLAR&1 then SS_NOMEASURE_rep, SS_MEASURE_rep required
           (#defined in sitesite.h)
           #if not POLAR&1 then they (if not #defined) are considered empty

  REMAINS TO CHECK: 
    zero charge and nonzero Drude charge in some versions (OK with GAUSSIANCHARGES)

#define DEBUG // brute force - elst energy WRONG!
*/

#ifndef POLAR
#  error "POLAR not #defined"
#endif /*# POLAR */

#if POLAR&32
#  error "interpol: POLAR&32 not expected"
#endif /*# POLAR&32 */

#if !defined(SS_MEASURE_rep) || !defined(SS_NOMEASURE_rep)
#  if POLAR&1
#    error "missing SS_MEASURE_rep or SS_NOMEASURE_rep for POLAR&1"
#  else /*# POLAR&1 */
#    define SS_MEASURE_rep /*empty*/
#    define SS_NOMEASURE_rep /*empty*/
#  endif /*#!POLAR&1 */
#endif /*# !defined(SS_MEASURE_rep) || !defined(SS_NOMEASURE_rep) */

#ifdef STARS
#  error "STARS not supported for POLAR"
#endif /*# STARS */

double polarLJQQM(sitesite_t *ss,                          /***** polarLJQQM */
                  vector r1,vector r2,vector f1,vector f2,
                  siteinfo_t *si1, siteinfo_t *si2 GLOB)
/* basic nonbonded version with measurements */
{
  double U=0,rr,f,x,y,z,qq;
  rdf_t *ssrdf;
  ermacrodcl
  vector dr,dr0;
  real
    fdr,
    *r1pol=(real*)((char*)r1+polar_off),
    *r2pol=(real*)((char*)r2+polar_off),
    *f1pol=(real*)((char*)f1+polar_off),
    *f2pol=(real*)((char*)f2+polar_off);

#ifdef DEBUG
  En.el=0;
#endif /*# DEBUG */

  beginCellImage(0,box.cq)

#if POLAR&1
    double frep=0,Urep=0;

    int anion1=si1->qtype&1 && si2->qtype&2;
    int anion2=si2->qtype&1 && si1->qtype&2;
#endif /*# POLAR&1 */

    if ( (ssrdf=ss->rdf) ) {
      int i=(int)(sqrt(rr)*ssrdf->grid);
      if (i<ssrdf->nhist) (ssrdf->hist[i])++; }

    f=0;
    if (rr < ss->C2q) {
      SS_MEASURE_rep
      if (rr < ss->C1q)
        SS_MEASURE
      else {
        x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
      ENVIR+=rr*f; }

#if POLAR&4
    /* should be checked/reconsidered */
    if (si1->qtype+si2->qtype==0) {
      /* uncharged nonpolarizable */
      MADDFORCES(f1,f2,f,dr) /* LJ */
      return U; }
#endif /*# POLAR&4 */

#ifdef QQTAB
#  include "drudedrudeqqtabm.c"
#else /*# QQTAB */
#  include "drudedrudem.c"
#endif /*#!QQTAB */

  endCellImage

#ifdef DEBUG
  prt("%d %d r=%16.12f U=%.8f El=%.8f ljqqm%d%d ~",
      (vector*)r1-cfg[0]->rp,(vector*)r2-cfg[0]->rp,
      sqrt(rr),
      //      (U+En.el)*0.0019872036,
      U,En.el,
      si1->qtype&1 && si2->qtype&2,
      si2->qtype&1 && si1->qtype&2);
#endif /*# DEBUG */

  return U;
} /* polarLJQQM */


void polarLJQQ(sitesite_t *ss,                              /***** polarLJQQ */
               vector r1,vector r2,vector f1,vector f2,
               siteinfo_t *si1, siteinfo_t *si2)
/* basic nonbonded version without measurements */
{
  double rr,f,x,y,qq,z;
  ermacrodcl  /* byerd not needed */
  vector dr,dr0;
  real
    fdr,
    *r1pol=(real*)((char*)r1+polar_off),
    *r2pol=(real*)((char*)r2+polar_off),
    *f1pol=(real*)((char*)f1+polar_off),
    *f2pol=(real*)((char*)f2+polar_off);

  beginCellImage(,box.cq)

#if POLAR&1
    double frep=0;

    int anion1=si1->qtype&1 && si2->qtype&2;
    int anion2=si2->qtype&1 && si1->qtype&2;
#endif /*# POLAR&1 */

    f=0;
    if (rr < ss->C2q) {
      SS_NOMEASURE_rep
      if (rr < ss->C1q)
        SS_NOMEASURE(=)
      else {
        x=rr-ss->C2q; f=ss->A4*x; } }

#if POLAR&4
    /* should be checked/reconsidered */
    if (si1->qtype+si2->qtype==0) {
      /* uncharged nonpolarizable */
      ADDFORCES(f1,f2,f,dr) /* LJ */
      return; }
#endif /*# POLAR&4 */

#ifdef QQTAB
#  include "drudedrudeqqtab.c"
#else /*# QQTAB */
#  include "drudedrude.c"
#endif /*#!QQTAB */

  endCellImage

  return;
} /* polarLJQQ */


double polarXQQM(sitesite_t *ss,                         /******* polarXQQM */
                 vector r1,vector r2,vector f1,vector f2,
                 siteinfo_t *si1, siteinfo_t *si2 GLOB)
/*
  1-2 and 1-3 exclusions, both measure+nomeasure
  #if POLAR&64: intramolecular induced dipole-induced dipole interactions
*/
{
#if !(POLAR&64) && (defined(FREEBC) || defined(NIBC) || COULOMB>=0)
  /* free b.c. and cutoff electrostatics: no 1-2 and 1-3 correction */
  return 0;
#else /*# !(POLAR&64) && (defined(FREEBC) || defined(NIBC) || COULOMB>=0) */
  /* now:
     - either POLAR&64 (intramolecular induced dipole-induced dipole)
     - or Ewald corrections 
     - or both
  */
  double U=0;
  double rr,f,x,y,z,qq,qqpp;
#  if (defined(COULOMB) && COULOMB>=-1) || defined(QQTAB)
  ermacrodcl
#  endif
  vector dr,dr0;
  real
    fdr,
    *r1pol=(real*)((char*)r1+polar_off),
    *r2pol=(real*)((char*)r2+polar_off),
    *f1pol=(real*)((char*)f1+polar_off),
    *f2pol=(real*)((char*)f2+polar_off);

  beginDistance(0,box.cq)

#  ifdef QQTAB
#    include "drudedrudeqqtabx.c"
#  else /*# QQTAB */
#    include "drudedrudex.c"
#  endif /*#!QQTAB */

  endDistance

  return U;
#endif /*#!!(POLAR&64) && (defined(FREEBC) || defined(NIBC) || COULOMB>=0) */
} /* polarXQQM */


double polarLJQQ14M(sitesite_t *ss,                      /***** polarLJQQ14M */
                    vector r1,vector r2,vector f1,vector f2,
                    siteinfo_t *si1, siteinfo_t *si2,
                    int ifLJ, double qfactor_1 GLOB)
/*
  1-4 interactions, scaled qfactor times (LJ-like from nbfixes)
  #if POLAR&64: (ADIM)
    intramolecular induced dipole-induced dipole interactions in full
  #if POLAR&1:
    shell core terms 1-4 are not allowed (silently ignored)
*/
{
  double U=0;
  double rr,f,x,y,z,qq,qqpp;
  rdf_t *ssrdf;
  ermacrodcl
  vector dr,dr0;
  real
    fdr,
    *r1pol=(real*)((char*)r1+polar_off),
    *r2pol=(real*)((char*)r2+polar_off),
    *f1pol=(real*)((char*)f1+polar_off),
    *f2pol=(real*)((char*)f2+polar_off);

  double qfactor=qfactor_1+1; /* not ALWAYS needed */

#ifdef DEBUG
  En.el=0;
#endif /*# DEBUG */

  beginDistance(0,box.cq) 

#if POLAR&1
    /* shell core terms 1-4 are not allowed (silently ignored) */
    double frep=0,Urep=0;
#endif /*# POLAR&1 */

    if (measure)
      if ( (ssrdf=ss->rdf) ) {
        int i=(int)(sqrt(rr)*ssrdf->grid);
        if (i<ssrdf->nhist) (ssrdf->hist[i])++; }

    f=0;
    if (ifLJ && rr < ss->C2q) {
      SS_MEASURE_rep
      if (rr < ss->C1q)
        SS_MEASURE
      else {
        x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
      ENVIR+=rr*f; }

#if POLAR&4
    if (si1->qtype+si2->qtype==0) {
      /* uncharged nonpolarizable */
      MADDFORCES(f1,f2,f,dr) /* LJ */
      return U; }
#endif /*# POLAR&4 */

#ifdef QQTAB
    ERROR(("14 forces not supported with QQTAB"))
#else /*# QQTAB */
#  include "drudedrude14.c"
#endif /*#!QQTAB */

  endDistance

#ifdef DEBUG
  prt("%d %d r=%8.5f  U=%.4f 14~",(vector*)r1-cfg[0]->rp,(vector*)r2-cfg[0]->rp,sqrt(rr),(U+En.el)*0.0019872065);
#endif /*# DEBUG */
  return U;
} /* polarLJQQ14M */
