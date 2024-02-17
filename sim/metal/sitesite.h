#ifndef SITE_SITE
#  define SITE_SITE

/* METAL: RGL potential, tight binding potential (2nd order) 
   NEW version transparently implemented into MACSIMUS
   see also goldTB and goldErcan
*/

#  define SHARPCUTOFF /* no smooth cutoff, no cutoff corrections */

/* refers to atom (site) table name in the par- and ble- files */
#  define SS_TABLE "Metal"
#  include "simopt.h"

/*** parm tables and combining rules ***/
#  define SS_PARMS      3
#  define DECLARE_INITCOMBRULE \
  double dummy;
typedef struct pairaux_s {
  real E4,Sq; /* required by initwall() #ifdef SLAB, but wall likely cannot be used anyway */
  real r0,q,p,xi,A,f1,f2;
} pairaux_t;

/* auxiliary constants for site-site interactions */
/* will be accessible via sitesite_t->a */

/* common declarations */
#  include "ssdecl.h"

extern double globalLJshift;

#  if 0
/* old non-optimized version */
#    define SS_MEASURE { x=sqrt(rr); \
  y=exp(-ss->a.p*(x/ss->a.r0-1)); \
  U+=(ss->a.A*2)*y; \
  f=   ( (2*ss->a.A*ss->a.p/ss->a.r0)*y \
       -(ss->a.q/ss->a.r0)*(1/sqrt(rho1)+1/sqrt(rho2))*rho12 )/x; }
#    define SS_NOMEASURE(OP) { x=sqrt(rr); \
  f OP ( (2*ss->a.A*ss->a.p/ss->a.r0)*exp(-ss->a.p*(x/ss->a.r0-1)) \
       -(ss->a.q/ss->a.r0)*(1/sqrt(rho1)+1/sqrt(rho2))*rho12 )/x;}
#  else /*# 0 */
/* new optimized */
#    define SS_MEASURE { x=sqrt(rr); \
  y=exp(-ss->a.p*(x/ss->a.r0-1)); \
  U+=(ss->a.A*2)*y; \
  f=   ( ss->a.f1*y \
       -ss->a.f2*(1/sqrt(rho1)+1/sqrt(rho2))*rho12 )/x; }
#    define SS_NOMEASURE(OP) { x=sqrt(rr); \
  f OP ( ss->a.f1*exp(-ss->a.p*(x/ss->a.r0-1)) \
       -ss->a.f2*(1/sqrt(rho1)+1/sqrt(rho2))*rho12 )/x;}
#  endif /*#!0 */
#endif /*# SITE_SITE */
