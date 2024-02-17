#include "ground.h"
#include "vector.h" /* simglob.h not needed */
#include "sitesite.h"
#ifdef COOK
#  include "units.h"
#  define PROMPT
struct Epow_s Epow[2]={
       {1,0,{0.5
#  if POW
              ,0.5
#  endif /*# POW */
    }},{1,0,{1
#  if POW
              ,1
#  endif /*# POW */
         }}}; /* EvdW is energy, RvdW is other (length) */
#else /*# COOK */
#  define PROMPT "! "
#endif /*#!COOK */

void combruleinfo(int comb_rule) /**************************** combruleinfo */
{
#if POW==12
  prt(PROMPT "Exp-6-12 potential:\n"
PROMPT "  u_rep(r) = f rho_ij exp[ (R_ij-r)/rho_ij ]\n"
PROMPT "  u(r) = u_rep(r) - C_ij/r^6 + D_ij/r^12");
#elif POW==8
  prt(PROMPT "Exp-6-8 (Tosi-Fumi) potential:\n"
PROMPT "  u_rep(r) = f rho_ij exp[ (R_ij-r)/rho_ij ]\n"
PROMPT "  u(r) = u_rep(r) - C_ij/r^6 - D_ij/r^8");
#else /*? POW==8 */ /*#!POW==12!POW==8 */
  prt(PROMPT "Exp-6 potential:\n"
PROMPT "  u_rep(r) = f rho_ij exp[ (R_ij-r)/rho_ij ]\n"
PROMPT "  u(r) = u_rep(r) - C_ij/r^6");
#endif /*?!POW==8 */ /*#!POW==12!POW==8 */
  prt(PROMPT "where f = 0.05 e^2/AA^2 (CGS) = 16.603186 kcal/mol/AA");
  if (comb_rule) ERROR(("Busing combining rule (comb_rule=0) expected"))
  prt(PROMPT "The Busing (aka Kong) combining rule for u_rep:\n"
PROMPT "  rho_ij = rho_i + rho_j, R_ij = R_i + R_j\n"
PROMPT "Geometric mean combing rule for r^n:\n"
PROMPT "  C_ij=C_i*C_j");
#if POW
  prt(PROMPT "  D_ij=D_i*D_j");
#endif /*# POW */
}

void initcombrule(siteparm_t *sp,int comb_rule) /************** initcombrule */
/* unit conversion: note that parm[0] is C */
{
  /* sp->alpha unused */
  // since V2.8a the conversion is in cook/simdef.c
  // #ifdef COOK
  //   sp->EvdW*=(kcal/Eunit);
  //   /* sp->R unchanged */
  //   sp->parm[0]*=sqrt(kcal/Eunit); // C
  // #  if POW==12
  //   sp->parm[1]*=sqrt(kcal/Eunit); // D
  // #  endif
  // #endif
  sp->comb_rule=comb_rule;
}

void combrule(pairparm_t *pp,siteparm_t *lj0,siteparm_t *lj1,char *info)
{
  if (lj0 && lj1) {
    /* C_ij=C_i C_j */
    pp->parm[0]=lj0->parm[0]*lj1->parm[0];
#if POW
    pp->parm[1]=lj0->parm[1]*lj1->parm[1];
#endif /*# POW */
    pp->sig = lj0->RvdW+lj1->RvdW;
    if (lj0->EvdW<0 || lj1->EvdW<0) pp->eps=0;
    else pp->eps = lj0->EvdW+lj1->EvdW; }

  if (info) {
#if POW
    prt_("%s: eps=%.9g sig=%.9g C=%.9g D=%.9g",info,pp->eps,pp->sig,pp->parm[0],pp->parm[1]);
#else /*# POW */
    prt_("%s: eps=%.9g sig=%.9g C=%.9g",info,pp->eps,pp->sig,pp->parm[0]);
#endif /*#!POW */
    if (pp->eps) prt(" => A=%.9g [K], B=%.9g [AA-1]", pp->eps*exp(BUSING_f*pp->sig/pp->eps), BUSING_f/pp->eps);
    else prt(" (turned off)"); }
  /* NOTE: #ifdef GAUSSIANCHARGES, sigmaij=sqrt(sigmai^2+sigmaj^2) is calculated in gcelst */
}

/* calculate auxiliary constants for site-site interactions */
void initssaux(pairaux_t *a,pairparm_t *pp)
{
  a->eps=pp->eps;
  a->RR=pp->sig;
  a->rhorho=a->eps/BUSING_f;
  //  a->f_rhorho=BUSING_f/a->rhorho;
  a->C=pp->parm[0];
  a->C6=pp->parm[0]*6;
  /* optimization: */
  if (a->rhorho) {
    a->Aij=BUSING_f*exp(a->RR/a->rhorho);
    a->iBij=-BUSING_f/a->eps; }
  else
    a->Aij=a->iBij=0;
#if POW
  a->D=pp->parm[1];
  a->DPOW=pp->parm[1]*POW;
#endif /*# POW */
#ifdef SLAB
#  include "equivlj.c"
#endif /*# SLAB */
}
