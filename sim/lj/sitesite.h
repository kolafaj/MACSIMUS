/*
  Lennard-Jones support for BLEND and COOK

    uLJ(r) = -eps*[(sig/r)^12 - 2*(sig/r)^6]

  NB: eps<0, sig=2^(1/6)*LJsigma = van der Waals diameter

  Combining rules:

  comb_rule&1:
    0 = Kirkwood combining rule (with polarizabilities) for energy terms
    1 = square-root combining rule for energy terms:
        eps_ij = -sqrt(epsvdW_i*epsvdW_j)

  comb_rule&2:
    0 = additive (arithmetic mean) diameters
        sig_ij = RvdW_i + RvdW_j
    2 = geometric mean for diameters
        sig_ij = sqrt (4*RvdW_i*RvdW_j)

  (in COOK before V2.8a, epsvdW included kcal->K conversion)
*/

#ifndef SITE_SITE
#  define SITE_SITE

/* refers to atom (site) table name in the par- and ble- files */
#  define SS_TABLE "Lennard-Jones"
#  include "simopt.h"

/* enable optimized water models (unless turned off later by LINKCELL/POLAR) */
#  define WATER

/*** LJ-> tables and combining rules ***/
#  define SS_PARMS 0    /* no extra parms for LJ */
/* these parameters will be set by initcombrule() and then used by combrule() */
#  define DECLARE_INITCOMBRULE \
  real Esqrt; /* =sqrt(2*Emin), if alpha==0 => sqrt(Emin1*Emin2) comb. rule */ \
  real gamma; /* =alpha^2/EvdW/RvdW^6/256 */

/* constants for site-site interactions (former SS_DECLARE):
   will be accessible via  sitesite_t ss->a.*
   will be calculated by initssaux() after combrule() */
typedef struct pairaux_s {
  /* there parameters refer to form
    uLJ(r) = 4*LJeps*[(LJsig/r)^12 - (LJsig/r)^6]
  */
  real E4; /* 4*LJeps   NB: LJeps>0, LJeps=-eps=-Emin */
  real E24;/* 24*LJeps                          */
  real E48;/* 48*LJeps                          */
  real Sq; /* LJsig^2   NB: LJsig = sig/2^(1/6) */
} pairaux_t;

/* common declarations */
#  include "ssdecl.h"

/* basic macro for the LJ potential and forces, measurements included
   to be used by both cook and blend
   x,y,z are declared in advance */
#  define SS_MEASURE {     \
  x=ss->a.Sq/rr; y=x*x*x;  \
  U+=ss->a.E4*((z=y*y)-y); \
  f=ss->a.E24*(z+z-y)*x; }

/* LJ forces only (no measurements), for cook only */
#  define SS_NOMEASURE(OP) {    \
  x=ss->a.Sq/rr; y=x*x*x;       \
  f OP ss->a.E48*y*(y-0.5)*x; }

/* 
   Optional macros calculating the potential and force directly from
   pairparm_t pp.  These are used by `cook' to set up the smooth cutoff
   and calculate cutoff corrections.  All versions are tested for
   consistency.
   If missing, potential and force will be calculated from SS_MEASURE
   (using the auxiliary parameters).
   Real x,y,z are declared in advance and parameters taken from pp.
*/
#  define SS_U(r) (x=Pow6(pp.sig/r),-pp.eps*(x*x-2*x))
#  define SS_F(r) (x=Pow6(pp.sig/r),-12*pp.eps*(x*x-x)/r)

#endif /*# SITE_SITE */
