/*

  Exp-6 (Buckingham) potential for COOK and BLEND
  with Lennard-Jones-like combining rules

                 /  6                       a           6 \
    u(r) = eps * | --- exp[a(1-r/sigma)] - --- (sigma/r)  |
                 \ a-6                     a-6            /

  (a is parm[0], eps = |E_min|)

  Combining rules:

  combrule&1:
    0 = n.a.
    1 = square-root combining rule for energy terms
        eps_ij = -sqrt(eps_i*eps_j)

  combrule&2:
    0 = additive (arithmetic mean) diameters
        sig_ij = RvdW_i + RvdW_j
    2 = geometric mean for diameters
        sig_ij = sqrt (4*RvdW_i*RvdW_j)

  COOK: epsvdW includes kcal->K conversion

  #ifdef REPULSIVE then the r^-6 term is turned off (brute force)

  this module modifies interpot.c and partly cook/simdef.c 
*/

#ifndef SITE_SITE
#  define SITE_SITE

#  ifdef REPULSIVE
#    define SS_TABLE "repulsive-exp-6"
#    define TURNVDW 0* /* just multiply by 0 */
#  else /*# REPULSIVE */
#    define SS_TABLE "exp-6"
#    define TURNVDW /*empty*/
#  endif /*#!REPULSIVE */

#  include "simopt.h"

#  define SS_PARMS      1      /* parameter a (alpha) */
#  define DECLARE_INITCOMBRULE /* */

/* constants for site-site interactions (former SS_DECLARE):
   will be accessible via  sitesite_t ss->a.*
   will be calculated by initssaux() after combrule() */
typedef struct pairaux_s {
  real tau,eu,ef,ru,rf;
} pairaux_t;

/* common declarations */
#  include "ssdecl.h"

/* (old) warning: fix of negative u(r) at very small r removed, fingers crossed */

/* LJ/exp-6 in site-site functions with measurements included */
#  define SS_MEASURE {                          \
    x=sqrt(rr); /* double eval? */              \
    y=exp(ss->a.tau*x);                         \
    z=Cub(rr);                                  \
    U+=ss->a.eu*y-TURNVDW ss->a.ru/z;            \
    f=ss->a.ef*y/x-TURNVDW ss->a.rf/z/rr; }

/* no measurements */
#  define SS_NOMEASURE(OP) {                    \
    double z;                                   \
    x=sqrt(rr);                                 \
    y=exp(ss->a.tau*x);                         \
    z=Cub(rr);                                  \
    f OP ss->a.ef*y/x-TURNVDW ss->a.rf/z/rr; }

#endif /*# SITE_SITE */
