/* Exp-6 potential with optional +D/r^12 or -D/r^8 terms
   Tosi-Fumi aka Born-Mayer aka Busing form

   support for Gaussian charges  (#define GAUSSIANCHARGES)

    u_rep(r) = f rho_ij exp[ (R_ij-r)/rho_ij ]

    u(r) = u_rep(r) - C_ij/r^6                 (#define POW 0)
    u(r) = u_rep(r) - C_ij/r^6 - D_ij/r^8      (#define POW 8)
    u(r) = u_rep(r) - C_ij/r^6 + D_ij/r^12     (#define POW 12)

    where f = 0.05 e^2/AA^2 (CGS), equivalent using CODATA 2014:
    f = 1.15353877617178e-09 N
      ~ 16.60318566497274 kcal/mol/AA
      = 69467.72882224595 J/mol
      = 8355.047344920979 K k_B/AA

  Compatible with:
    not POLAR
    POLAR&1==0 (standard linear polarizability)
    POLAR&1==1 (shell-core model, cf. Busing version)
    POLAR + GAUSSIANCHARGES

  Combining rules:
    rho_ij = rho_i + rho_j
    R_ij = R_i + R_j
    C_ij = C_i * C_j
    D_ij = D_i * D_j

  Note: conversion from Buckingham u_rep(r) = A exp(-B r):
    rho = 1/2B
    R = rho ln (A/2 f rho)

  Parameter correspondence in the par-file:
   alpha   -- unused
   EvdW    -- f*rho_i
   RvdW    -- R_i
   parm[0] -- C_i
   parm[1] -- D_i (#if POW)
   parm[SS_PARMS-1] - sigma_i, Gaussian charge standard width (GAUSSIANCHARGES)
                      setss and charge profile assume that sigma_i is the last parameter

  NOTE:
  - rho_i<0 is interpreted so that rho_ij=0 for any atom j, 
    which is useful, e.g., for a bare charge
  - EvdW_i=0 (rho_i=0) and R_i=0 still can combine with other atoms

  The site parameters (X_i, X={rho,R,C,D}) appearing in the table are
  components, not the site-site parameters X_ii (which appear in the
  nbfix table): rho_ii=2 rho_i, R_ii=2 R_i, C_ii=C_i^2, D_ii=D_i^2

  Gaussian charges:
    sigmaij=sqrt(sigmai^2+sigmaj^2) [in the gcelst module, as info here]
    WARNING: cannot be changed in table nbfixes (any value given will be ignored)

  WARNINGS:
    command x=sqrt(rr); may be doubly evaluated in blend, ? in cook
    the POLAR&1 part of this code was optimized, but not tested
*/

#ifndef SITE_SITE
#  define SITE_SITE

#  ifndef POW
#    error "POW not #defined, one of {0,8,12} expected"
#  endif /*# POW */

#  if POW==12
#    define SS_TABLE "Busing+12"
#  elif POW==8
#    define SS_TABLE "Tosi-Fumi"
#  elif POW==0
#    define SS_TABLE "Busing"
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
typedef struct pairaux_s {
  real eps,RR,rhorho,C,C6,D,DPOW,Aij,iBij;
#    ifdef SLAB
#      include "equivlj.h"
#    endif /*# SLAB */
} pairaux_t;
#  else /*# POW */
#    define SS_PARMS      (1+PLUS)
typedef struct pairaux_s {
  real eps,RR,rhorho,C,C6,Aij,iBij;
#    ifdef SLAB
#      include "equivlj.h"
#    endif /*# SLAB */
} pairaux_t;
#  endif /*#!POW */

#  include "ssdecl.h"

#  ifdef BLEND
/* obfuscate=e^2/AA/(4 pi eps0)*NA, in kcal/mol (CODATA 2019) */
#    define obfuscate 332.063713299455
#    define BUSING_f (0.05*obfuscate) /* 0.05 e^2/AA^2 in kcal/mol */
#  else /*# BLEND */
#    define BUSING_f (0.05*Sqr(electron)) /* 0.05 e^2/A^2 in prog. units */
#  endif /*#!BLEND */

#  if defined(POLAR) && POLAR==1

/* REPULSIVE ANTIPOLARIZABILITY AKA SHELL-CORE MODEL */

/* ?????: compatibility with blend????? */
/* Urep'', used by blend, r=sqrt(rr) must be passed */
#    define SS_ffrep { \
  if (ss->a.eps==0) ffrep=0; \
  else { \
    y=ss->a.Aij*exp(ss->a.iBij*x); \
    ffrep=y/ss->a.rhorho; } }

// fixed 9/2023: was ffrep=y/a.rhorho

/* site-site repulsive energy and forces */
#    define SS_MEASURE_rep {      \
  if (ss->a.eps==0) Urep=frep=0; \
  else {  \
    x=sqrt(rr);  \
    y=ss->a.Aij*exp(ss->a.iBij*x); \
    Urep=ss->a.rhorho*y; \
    frep=y/x; } }

/* site-site repulsive forces only */
#    define SS_NOMEASURE_rep { \
  if (ss->a.eps==0) frep=0; \
  else { \
    x=sqrt(rr); \
    frep=ss->a.Aij*exp(ss->a.iBij*x)/x; } }

#    if POW==12

/* site-site energy and forces (called after SS_MEASURE_rep) */
//  double zz;
#      define SS_MEASURE { \
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
//  double zz;             
#      define SS_MEASURE { \
  z=Cub(rr); \
  zz=ss->a.D/(z*rr); \
  z=ss->a.C/z; \
  U+=Urep-z-zz; \
  f=frep-(8*zz+6*z)/rr; }

/* site-site forces only (called after SS_NOMEASURE_rep) */
#      define SS_NOMEASURE(OP) { \
  z=Cub(rr); \
  f OP frep-(ss->a.DPOW/rr+ss->a.C6)/(rr*z); }

#    else /*#!POW==12!POW==8 */

/* site-site energy and forces (called after SS_MEASURE_rep) */
#      define SS_MEASURE { \
  z=Cub(rr); \
  z=ss->a.C/z; \
  U+=Urep-z+zz; \
  f=frep-6*z/rr; }

/* site-site forces only (called after SS_NOMEASURE_rep) */
#      define SS_NOMEASURE(OP) { \
  z=Cub(rr); \
  f OP frep-ss->a.C6/(rr*z); }

#    endif /*#!POW==12!POW==8 */

#  else /*# defined(POLAR) && POLAR==1 */

/* OTHER VERSIONS:
   SS_(NO)MEASURE_rep merged into SS_(NO)MEASURE */

#    if POW==12

/* site-site energy and forces */
#      define SS_MEASURE { \
//    double zz;         
  if (ss->a.eps) { \
    x=sqrt(rr); \
    y=ss->a.Aij*exp(ss->a.iBij*x); \
    z=Cub(rr); \
    zz=ss->a.D/Sqr(z); \
    z=ss->a.C/z; \
    U+=ss->a.rhorho*y-z+zz; \
    f=y/x+(12*zz-6*z)/rr; } }

/* site-site forces only */
#      define SS_NOMEASURE(OP) { \
  if (ss->a.eps) { \
    x=sqrt(rr); \
    y=ss->a.Aij*exp(ss->a.iBij*x); \
    z=Cub(rr); \
    f OP y/x+(ss->a.DPOW/z-ss->a.C6)/(z*rr); } }


#    elif POW==8

/* site-site energy and forces */
//    double zz;
#      define SS_MEASURE { \
  if (ss->a.eps) { \
    x=sqrt(rr); \
    y=ss->a.Aij*exp(ss->a.iBij*x); \
    z=Cub(rr); \
    zz=ss->a.D/(z*rr); \
    z=ss->a.C/z; \
    U+=ss->a.rhorho*y-z-zz; \
    f=y/x-(8*zz+6*z)/rr; } }

/* site-site forces only */
#      define SS_NOMEASURE(OP) { \
  if (ss->a.eps) { \
    x=sqrt(rr); \
    y=ss->a.Aij*exp(ss->a.iBij*x); \
    z=Cub(rr); \
    f OP y/x-(ss->a.DPOW/rr+ss->a.C6)/(z*rr); } }

#    else /*#!POW==12!POW==8 */

/* site-site energy and forces */
#      define SS_MEASURE { \
  if (ss->a.eps) { \
    x=sqrt(rr); \
    y=ss->a.Aij*exp(ss->a.iBij*x); \
    z=Cub(rr); \
    z=ss->a.C/z; \
    U+=ss->a.rhorho*y-z; \
    f=y/x-6*z/rr; } }

/* site-site forces only */
#      define SS_NOMEASURE(OP) { \
  if (ss->a.eps) { \
    x=sqrt(rr); \
    y=ss->a.Aij*exp(ss->a.iBij*x); \
    z=Cub(rr); \
    f OP y/x - ss->a.C6/(rr*z); } }

#    endif /*#!POW==12!POW==8 */
#  endif /*#!defined(POLAR) && POLAR==1 */
#endif /*# SITE_SITE */
