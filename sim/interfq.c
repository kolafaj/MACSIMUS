/*
  fluctuating charge version of internp.c/interpol.c

  see also interpol.c
*/

//#define DEBUG // brute force - elst energy WRONG!

#ifndef POLAR
#  error "POLAR not #defined"
#endif /*# POLAR */

#if !(POLAR&32)
#  error "interfq: POLAR&32 expected"
#endif /*# !(POLAR&32) */

#if POLAR&64
#  error "interfq: POLAR&64 not supported for FQ (probably OK for Drude - to be checked)"
#endif /*# POLAR&64 */

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
    /* POLAR&4 just changes the order of the test */
    if (si1->qtype+si2->qtype==0) {
      /* uncharged nonpolarizable */
      MADDFORCES(f1,f2,f,dr) /* LJ */
      return U; }
#endif /*# POLAR&4 */

    /* to be optimized by a ss table */
    if (si1->qtype&DRUDE && si2->qtype&DRUDE) {
#include "drudedrudem.c"
    } 
    else if (si1->qtype&FQ && si2->qtype&FQ) {
      qq=r1pol[0]*r2pol[0];
#ifdef COULOMB
      y=erud(rr);
      f1pol[0]+=y*r2pol[0];
      f2pol[0]+=y*r1pol[0];
      ENEL+=y*qq; f+=byerd*qq;
#else /*# COULOMB */
      y=sqrt(rr);
      f1pol[0]+=r2pol[0]/y;
      f2pol[0]+=r1pol[0]/y;
      ENEL+=y=qq/y; f+=y/rr;
#endif /*#!COULOMB */
      MADDFORCES(f1,f2,f,dr) }
#if !(POLAR&4)
      else if (si1->qtype+si2->qtype==0) {
      /* uncharged nonpolarizable */
      MADDFORCES(f1,f2,f,dr) /* LJ */ }
#endif /*# !(POLAR&4) */
    else {
      // this includes, e.g., 36 vs. 0 , which is zero interaction
      //      ERROR(("qtypes %d vs. %d not implemented",si1->qtype, si2->qtype))
    }
  endCellImage

#ifdef DEBUG
  prt("%d %d r=%8.5f  U=%.4f ljqqm%d%d ~",(vector*)r1-a[0]->rp,(vector*)r2-a[0]->rp,sqrt(rr),(U+En.el)*0.0019872065,
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
    if (si1->qtype+si2->qtype==0) {
      /* uncharged nonpolarizable */
      ADDFORCES(f1,f2,f,dr) /* LJ */
      return; }
#endif /*# POLAR&4 */

    /* could be optimized to the ss table */
    if (si1->qtype&DRUDE && si2->qtype&DRUDE) {
#include "drudedrude.c"
    } 
    else if (si1->qtype&FQ && si2->qtype&FQ) {
      qq=r1pol[0]*r2pol[0];
#ifdef COULOMB
      y=erud(rr);
      f1pol[0]+=y*r2pol[0];
      f2pol[0]+=y*r1pol[0];
      f+=byerd*qq;
#else /*# COULOMB */
      y=sqrt(rr);
      f1pol[0]+=r2pol[0]/y;
      f2pol[0]+=r1pol[0]/y;
      f+=qq/(y*rr);
#endif /*#!COULOMB */
      MADDFORCES(f1,f2,f,dr) /* charge-charge + LJ */ }
#if !(POLAR&4)
    else if (si1->qtype+si2->qtype==0) {
      /* uncharged nonpolarizable */
      ADDFORCES(f1,f2,f,dr) /* LJ */ }
#endif /*? POLAR&4 */ /*# !(POLAR&4) */
    else {
      // this includes, e.g., 36 vs. 0 , which is zero interaction
      //      ERROR(("qtypes %d vs. %d not implemented",si1->qtype, si2->qtype))
    }
  endCellImage

  return;
} /* polarLJQQ */


double polarXQQM(sitesite_t *ss,                         /******* polarXQQM */
                 vector r1,vector r2,vector f1,vector f2,
                 siteinfo_t *si1, siteinfo_t *si2 GLOB)
/*
  1-2 and 1-3 exclusions, both measure+nomeasure
  #if POLAR&64: ?
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
#  if defined(COULOMB) && COULOMB>=-1
  ermacrodcl
#  endif /*# defined(COULOMB) && COULOMB>=-1 */
  vector dr,dr0;
  real
    fdr,
    *r1pol=(real*)((char*)r1+polar_off),
    *r2pol=(real*)((char*)r2+polar_off),
    *f1pol=(real*)((char*)f1+polar_off),
    *f2pol=(real*)((char*)f2+polar_off);

  beginCellImage(0,box.cq)

#  if POLAR&4
    /* POLAR&4 just changes the order of the test */
    if (si1->qtype+si2->qtype==0) {
      /* uncharged nonpolarizable  - no term */
      return U; }
#  endif /*# POLAR&4 */

    /* to be optimized by a ss table */
    if (si1->qtype&DRUDE && si2->qtype&DRUDE) {
#  include "drudedrudex.c"
    } 
    else if (si1->qtype&FQ && si2->qtype&FQ) {
#  if defined(FREEBC) || defined(NIBC) || COULOMB>=0
#    error internal
      /* void */
      return 0;
#  else /*# defined(FREEBC) || defined(NIBC) || COULOMB>=0 */
      qq=r1pol[0]*r2pol[0];
#    if COULOMB==-2
      x=exacterud_sqrt(rr,&y);
      f=y*qq;
#    else /*# COULOMB==-2 */
      y=1/sqrt(rr);
      x=erud(rr)-y;
      f=(byerd-y/rr)*qq;
#    endif /*#!COULOMB==-2 */
      ENEL+=x*qq; 
      f1pol[0]+=r2pol[0]*x;
      f2pol[0]+=r1pol[0]*x;
      MADDFORCES(f1,f2,f,dr) 
#  endif /*#!defined(FREEBC) || defined(NIBC) || COULOMB>=0 */
    }
#  if !(POLAR&4)
    else if (si1->qtype+si2->qtype==0) {
      /* uncharged nonpolarizable */
      return 0;
     }
#  endif /*# !(POLAR&4) */
    else {
      // this includes, e.g., 36 vs. 0 , which is zero interaction
      //      ERROR(("qtypes %d vs. %d not implemented",si1->qtype, si2->qtype))
    }

  endCellImage

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
  ERROR(("1-4 for FQ not implemented"))
  return 0;
} /* polarLJQQ14M */
