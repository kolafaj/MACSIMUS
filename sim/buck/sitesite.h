/* Buckingham potential with optional +D/r^12 or -D/r^8 terms
   exp-6 aka Buckingham form

   support for Gaussian charges  (#define GAUSSIANCHARGES)

    u_rep(r) = A_ij exp(-B_ij r)

    u(r) = u_rep(r) - C_ij/r^6
    u(r) = u_rep(r) - C_ij/r^6 - D_ij/r^8      (#define POW 8)
    u(r) = u_rep(r) - C_ij/r^6 + D_ij/r^12     (#define POW 12)

  The site parameters (X_i, X={A,B,C,D}) appearing in the "buck" table
  are the same as the site-site parameters X_ii.

  Compatible with:
    not POLAR
    POLAR&1==0 (standard linear polarizability)
    POLAR&1==1 (shell-core model, cf. Busing version)
    POLAR + GAUSSIANCHARGES

  Parameters:
    alpha   -- unused
    EvdW    -- A_ii  [kcal/mol], converted to K in simdef.c
    RvdW    -- B_ii  [1/AA]
    parm[0] -- C_ii  [kcal/mol.AA^6], converted to K AA^6 in simdef.c
    parm[1] -- D_ii  [kcal/mol.AA^12], converted to K AA^12 in simdef.c (REP12)
    parm[SS_PARMS-1] - sigma_i, Gaussian charge standard width (GAUSSIANCHARGES)
                       setss and charge profile assume that sigma_i is the last parameter
  Combining rules:

  comb_rule&3:
    0: B_ij=[2/(1/B_ii^6+1/B_jj^6)]^(1/6) (Waldman-Hagler)
    1: B_ij=(B_ii B_jj)^1/2               (sqrt = geometric)
    2: B_ij=(B_ii+B_jj)/2                 (additive repulsive parameters B)
    3: B_ij=2/(1/B_ii+1/B_jj)             (additive radii 1/B; Busing)
   0,1,3: B_ij=0 if B_ii=0 or B_jj=0

  comb_rule&12:
    0: A_ij=(A_ii A_jj)^1/2 B_ij^6/(B_ii B_jj)^3  (Waldman-Hagler)
    4: A_ij=(A_ii A_jj)^1/2                       (sqrt = geometric)
    8: A_ij=[ A_ii*(A_jj*B_jj/A_ii/B_ii)^qi + A_jj*(A_ii*B_ii/A_jj/B_jj)^qj ]/2 (Kong)
       qi=B_ii/(B_ii+B_jj)
       qj=B_jj/(B_ii+B_jj)
   12: A_ij=(A_ii B_ii)^qi (A_jj B_jj)^qj / B_ij  (Busing)
         qi=1/B_ii/(1/B_ii+1/B_jj)
         qj=1/B_jj/(1/B_ii+1/B_jj)
   A_ij=0 if A_ii=0 or A_jj=0


   NOTES:
     for additive radii 1/B, Busing = Kong
     other rule for B is possible, then Busing != Kong (none meaningful)
     V2.7q and older: comb_rule&12==8: A_ij=(A_ii+A_jj)/2 (arithmetic, not useful) 

  r^-6, r^-8 and r^-12 terms:
    C_ij=(C_i C_j)^1/2
    D_ij=(D_i D_j)^1/2
    NOTE: some other modules may use C_ij=C_i C_j and/or D_ij=D_i D_j

  Gaussian charges:
    sigmaij=sqrt(sigmai^2+sigmaj^2) [in the gcelst module, as info here]
    WARNING: cannot be changed in table nbfixes (any value given will be ignored)

  See also:
    sim/busing (The Busing combining rule directly)
*/

#ifndef SITE_SITE
#  define SITE_SITE

#  ifndef POW
#    error "POW not #defined, one of {0,8,12} expected"
#  endif /*# POW */

#  if POW==12
#    define SS_TABLE "Buckingham+12"
#  elif POW==8
#    define SS_TABLE "Buckingham-8"
#  elif POW==0
#    define SS_TABLE "Buckingham"
#  else /*#!POW==12!POW==8!POW==0 */
#    error "POW is not one of {0,8,12}"
#  endif /*#!POW==12!POW==8!POW==0 */

#  include "simopt.h"
#  include "units.h"
#  include "int4.h"

#  ifdef GAUSSIANCHARGES
#    define PLUS 1
#  else /*# GAUSSIANCHARGES */
#    define PLUS 0
#  endif /*#!GAUSSIANCHARGES */

/*** exp-6 tables and combining rules ***/
#  if POW
#    define SS_PARMS      (2+PLUS)
#    define DECLARE_INITCOMBRULE \
  double sqA,sqB,Bi6,Bi,sqABi3,AB,sqC, sqD;

typedef struct pairaux_s {
  real Aij,Bij,Bijq,C,C6,D,DPOW;
#    ifdef SLAB
    /* The parameters refer to the equivalent LJ potential (the same minimum).
       The LJ potential is in the form:
         uLJ(r) = 4*LJeps*[(LJsig/r)^12 - (LJsig/r)^6]
       Sq=LJsig^2, E4=4*LJeps */
  real LJeps,LJsig,Sq,E4;
#    endif /*# SLAB */
} pairaux_t;

#  else /*# POW */

#    define SS_PARMS      (1+PLUS)
#    define DECLARE_INITCOMBRULE \
  double sqA,sqB,Bi6,Bi,sqABi3,AB,sqC;

typedef struct pairaux_s {
  real Aij,Bij,C,C6;
#    ifdef SLAB
    /* The parameters refer to the equivalent LJ potential (the same minimum).
       The LJ potential is in the form:
         uLJ(r) = 4*LJeps*[(LJsig/r)^12 - (LJsig/r)^6]
       Sq=LJsig^2, E4=4*LJeps>0 */
  real LJeps,LJsig,Sq,E4;
#    endif /*# SLAB */
} pairaux_t;
#  endif /*#!POW */

/* auxiliary constants for site-site interactions */
/* will be accessible via sitesite_t->a */

#  include "ssdecl.h"

#  if defined(POLAR) && (POLAR&1)

/* Urep'', used by blend, r=sqrt(rr) must be passed */
#    define SS_ffrep { \
  if (ss->a.Aij==0) ffrep=0; \
  else { \
    y=ss->a.Aij*exp(-ss->a.Bij*r); \
    ffrep=ss->a.Bijq*y; } }

/* site-site repulsive energy and forces */
#    define SS_MEASURE_rep { \
  if (ss->a.Aij==0) Urep=frep=0; \
  else { \
    x=sqrt(rr); /* double eval? */ \
    y=ss->a.Aij*exp(-ss->a.Bij*x); \
    Urep=y;  \
    frep=ss->a.Bij*y/x; }

/* site-site repulsive forces only */
#    define SS_NOMEASURE_rep  {  \
  x=sqrt(rr); /* double eval? */ \
  y=ss->a.Aij*exp(-ss->a.Bij*x); \
  frep=ss->a.Bij*y/x; }

#    if POW==12

/* site-site energy and forces (called after SS_[NO]MEASURE_rep) */
#      define SS_MEASURE { double zz; \
  z=Cub(rr); \
  zz=ss->a.D/Sqr(z); \
  z=ss->a.C/z; \
  U+=Urep-z+zz; \
  f=frep+(12*zz-6*z)/rr; }

/* site-site forces only (called after SS_NOMEASURE_rep) */
#      define SS_NOMEASURE(OP) { \
  z=Cub(rr); \
  f OP frep+(ss->a.DPOW/z-ss->a.C6)/(rr*z); }

#    elif POW==8

/* site-site energy and forces (called after SS_MEASURE_rep) */
#      define SS_MEASURE { double zz; \
  z=Cub(rr); \
  zz=ss->a.D/(z*rr); \
  z=ss->a.C/z; \
  U+=Urep-z-zz; \
  f=frep-(8*zz+6*z)/rr; }

/* site-site forces only (called after SS_NOMEASURE_rep) */
#      define SS_NOMEASURE(OP) {              \
  z=Cub(rr);                                  \
  f OP frep-(ss->a.DPOW/z+ss->a.C6)/(rr*z); }

#    else /*#!POW==12!POW==8 */

#      define SS_MEASURE { \
  z=ss->a.C/Cub(rr); \
  U+=Urep-z; \
  f=frep-6*z/rr; }
#      define SS_NOMEASURE(OP) { \
  f OP frep-ss->a.C6/Pow4(rr); }
#    endif /*?!REP12 */ /*#!POW==12!POW==8 */

#  else /*# defined(POLAR) && (POLAR&1) */

/* standard (POLAR not shell-core or nonPOLAR) */
#    if POW==12

/* site-site energy and forces */
#      define SS_MEASURE { double zz; \
  x=sqrt(rr); /* double eval? */  \
  y=ss->a.Aij*exp(-ss->a.Bij*x);  \
  z=Cub(rr); \
  zz=ss->a.D/Sqr(z); \
  z=ss->a.C/z; \
  U+=y-z+zz; \
  f=ss->a.Bij*y/x+(12*zz-6*z)/rr; }

/* site-site forces only */
#      define SS_NOMEASURE(OP) { \
  x=sqrt(rr); /* double eval? */ \
  y=ss->a.Aij*exp(-ss->a.Bij*x); \
  frep=ss->a.Bij*y/x; \
  z=Cub(rr); \
  f OP frep+(ss->a.D12/z-ss->a.C6)/(rr*z); }

#    elif POW==8

/* site-site energy and forces */
#      define SS_MEASURE { double zz; \
  x=sqrt(rr); /* double eval? */  \
  y=ss->a.Aij*exp(-ss->a.Bij*x);  \
  z=Cub(rr); \
  zz=ss->a.D/(z*rr); \
  z=ss->a.C/z; \
  U+=y-z-zz; \
  f=ss->a.Bij*y/x-(8*zz+6*z)/rr; }

/* site-site forces only */
#      define SS_NOMEASURE(OP) { \
  x=sqrt(rr); /* double eval? */ \
  y=ss->a.Aij*exp(-ss->a.Bij*x); \
  frep=ss->a.Bij*y/x; \
  z=Cub(rr); \
  f OP frep+(ss->a.DPOW/rr-ss->a.C6)/(rr*z); }

#    else /*#!POW==12!POW==8 */

/* site-site energy and forces */
#      define SS_MEASURE { \
  x=sqrt(rr); /* double eval? */ \
  y=ss->a.Aij*exp(-ss->a.Bij*x); \
  z=ss->a.C/Cub(rr); \
  U+=y-z; \
  f=ss->a.Bij*y/x-6*z/rr; }

/* site-site forces only */
#      define SS_NOMEASURE(OP) { \
  x=sqrt(rr); /* double eval? */ \
  y=ss->a.Aij*exp(-ss->a.Bij*x); \
  f OP ss->a.Bij*y/x-ss->a.C6/Pow4(rr); }

#    endif /*#!POW==12!POW==8 */
#  endif  /*#!defined(POLAR) && (POLAR&1) */
#endif /*# SITE_SITE */
