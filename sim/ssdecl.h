/* 
   Common structures and types to be #included in XXX/sitesite.h 
   and used by both blend and cook.

   site-site strategy
   ^^^^^^^^^^^^^^^^^^

   Single atom parameters in the par-files are:
     alpha = auxiliary parameter for a combining rule:
            LJ, exp6: polarizability (for the Kirkwood-Slater comb. rule)
            other: unused
     EvdW = any energy-like parameter:
            LJ: energy minimum
            Busing: f*rho_i
            Buckingham: A_i=A_ii
     RvdW = any size-like parameter:
            LJ: van der Waals radius = sigma/2^2^(1/6)
            Busing: R_i
            Buckingham: B_i=B_ii (actually it has dimension A^-1)
   additional parameters (their number is given by SS_PARMS):
     parm[0] = additional parameter:
            Busing: C_i=sqrt(C_ii)       
            Buckingham: C_ii (for -C_ii/r^6 term)
     parm[1] = additional parameter (cf. REP12):
            Busing: D_i=sqrt(D_ii)       
            Buckingham: D_ii (for D_ii/r^12 term)

   Function initcombrule() calculates (optional) auxiliary single atom
   parameters, see macro DECLARE_INITCOMBRULE in structure siteparm_t.
   These parameters are used to optimize the combining rule calculations.

   Function combrule() calculates pair atom-atom (site-site) parameters
   (from single atom parameters and the auxiliary single atom parameters).
   These atom-atom parameters may be overridden by data in the "nbfixes"
   table.

   The pair atom-atom parameters are:
     eps = energy-like parameter: 
           LJ,exp6: energy minimum
           Busing: f*(rho_i+rho_j)
           Buckingham: A_ij
     sig = size-like parameter:
           LJ: van der Waals diameter (sig=sigma*2^(1/6))
           Busing: R_i+R_j
           Buckingham: B_ij (actually 1/size)
   additional parameters:
     parm[0] = additional parameter:
           Busing, Buckingham: C_ij (for -C_ij/r^6)
     parm[1] = additional parameter:
           Busing, Buckingham: D_ij (for D_ij/r^12)

   NOTE: the single atom parameters are not generally the pair
   atom-atom ones. E.g., for LJ: sig_ii=2*RvdW_i whereas eps_ii=EvdW_i.

   Function initssaux() calculates auxiliary constants needed for pair
   atom-atom interaction (see pairaux_t) and used by SS_MEASURE,
   SS_NOMEASURE, etc.

   The pair atom-atom interactions come in two flavors, with
   measurements of energy (macro SS_MEASURE) and without it (macro
   SS_NOMEASURE).  The forces are provided as f=-dU/dr/r.

   The shell-core model needs a separate repulsive part, which then
   does not enter SS_MEASURE and SS_NOMEASURE.  In addition, blend
   (using point dipoles) needs the derivative of force, SS_ffrep.

   Cook contains a module checking the macros for consistency.  In
   addition, if separate unoptimized "fool-proof" versions of energy
   and forces are give, they are checked, too (and used in cutoff
   error calculations).
*/

#ifdef FLOAT
/* FLOAT implies IFLOAT */
#  define IFLOAT
#endif /*# FLOAT */

#ifndef IREAL
#  define IREAL
#  ifdef IFLOAT
typedef float ireal;
#    define irealfmt "%f"
#  else /*# IFLOAT */
typedef double ireal;
#    define irealfmt "%lf"
#  endif /*#!IFLOAT */
#endif /*# IREAL */

#include "int4.h"

/* 
   SS_PARMS is #defined in XXX/sitesite.h and means the number of parameters 
   in addition to the obligatory alpha,EvdW,RvdW (see siteparm_t)
*/
#if SS_PARMS==0
/* e.g., Lennard-Jones potential - no additional parameter */
#  define DECLARE_PARMS /* empty */
#  define PASS_PARMS    /* empty */
/* exponents of the energy unit; 0 for lengths (in Angstroem) etc.; 
   [0] is for sites, [1] for nbfixes */
extern struct Epow_s {
  double EvdW,RvdW;
} Epow[2];
#else /*# SS_PARMS==0 */
/* e.g., exp-6 and its variants */
/* new: sigma of Gaussian charges is the last parm, used in setqq() */
#  define DECLARE_PARMS double parm[SS_PARMS];
#  define PASS_PARMS    ,parm 
/* exponents of the energy unit; 0 for lengths (in Angstroem) etc.; 
   [0] is for sites, [1] for nbfixes */
extern struct Epow_s {
  double EvdW,RvdW;
  double parm[SS_PARMS];
} Epow[2];
#endif /*#!SS_PARMS==0 */

/* 
   parameters of atoms (sites), as read from par-files or ble-files from 
   table SS_TABLE (#defined in the respective XXX/sitesite.h)
*/
typedef struct siteparm_s {
  /* obligatory parameters for any force field (even if alpha is not used): */
  real alpha;   /* polarizability in A^3, for combining rules only */
  real EvdW;    /* energy parameter (e.g., the minimum energy; LJ: EvdW<0) */
  real RvdW;    /* atom-size parameter (e.g., the van der Waals radius) */
  DECLARE_PARMS /* more parameters (e.g., empty for LJ) */
#ifdef DECLARE_INITCOMBRULE
  DECLARE_INITCOMBRULE  /* optional auxiliary parameters set in initcombrule()
                           and used by combrule() */
#endif /*# DECLARE_INITCOMBRULE */
#ifdef BLEND
  /* to get architecture-portable files */
  int4 comb_rule; /* the combining rule, see combrule() */
  int4 paddto8B;
#else /*# BLEND */
  int comb_rule;
#endif /*#!BLEND */
} siteparm_t;

/* parameters for the pair potential, as calculated by combrule() or
   taken from table nbfixes (old name NBFIX) */
typedef struct pairparm_s {
  real eps;  /* min energy (generally any energy-like parameter; A) */
  real sig;  /* van der Waals diameter (atom-size-like parameter; B) */
  DECLARE_PARMS /* more parameters (empty for LJ) */
} pairparm_t;

/* the nbfixes list */
typedef struct nbfix_s {
  struct nbfix_s *next;
  int indx[2]; /* site (atom) numbers in the table of atoms */
  pairparm_t onefour[2]; /* [0]: basic non-bonded, [1]: for 1-4 interactions */
} nbfix_t;

/* calculate auxiliary atom parameters in DECLARE_INITCOMBRULE */
void initcombrule(siteparm_t *lj,int comb_rule);

/* the combining rule: calculate pp from lj0,lj1 */
void combrule(pairparm_t *pp,siteparm_t *lj0,siteparm_t *lj1,char *info);

/* 
   calculate auxiliary parameters for site-site macros SS_MEASURE and
   SS_NOMEASURE (former macro SS_DEF) 
   pairaux_t is declated in XXX/sitesite.h and should appear inside sitesite_t
*/
void initssaux(pairaux_t *a,pairparm_t *pp);

/* print verbose info on the combining rules */
void combruleinfo(int comb_rule);
