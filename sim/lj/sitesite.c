/* 
  non-bonded (atom-atom) potential support + combining rules
  see sitesite.h for more info
*/

#include "ground.h"
#include "vector.h" /* simglob.h not needed */
#include "sitesite.h"
#ifdef COOK
#  include "units.h"
#  define PROMPT
/* NEW in cook V2.8a: exponents of the energy unit */
struct Epow_s Epow[2]={{1,0},{1,0}}; /* EvdW is energy, RvdW is other (length) */
#else
#  define PROMPT "! "
#endif

/*** LJ combining rules ***/

void combruleinfo(int comb_rule) /***************************** combruleinfo */
{
  prt(PROMPT "Lennard-Jones potential:\n"
PROMPT "  uLJ(r) = -eps*[(sig/r)^12 - 2*(sig/r)^6]  (where eps<0)");
  if (comb_rule==1) prt(PROMPT "The Lorentz-Berthelot combing rule:");
  if (comb_rule==0) prt(PROMPT "The CHARMM 21 combing rule:");
  prt(PROMPT "comb_rule&2: atom-atom diameters (sig) are %s",
      comb_rule&2?"given by the geometric mean":"additive");
  prt(PROMPT "comb_rule&1: atom-atom min. energies (eps) are given by the\n"
PROMPT "             %s",  comb_rule&1?"geometric mean":"Kirkwood-Slater formula using atomic polarizabilities\n"
PROMPT "             (or geometric mean if either polarizability is zero)");
}

void initcombrule(siteparm_t *LJ,int comb_rule) /************** initcombrule */
/***
  pre-calculates some constants (like gamma) for the LJ combining rule
  comb_rule&1:
    0 = Kirkwood combining rule (with polarisabilities) for energy terms
    1 = square-root combining rule for energy terms
  comb_rule&6:
    0 = additive (arithmetic mean) diameters
    2 = geometric mean for diameters
  comb_rule for a pair is obtained by `or' operator, i.e. geometric means have
  higher priority
  (if LJ->alpha<=0 then also comb_rule applies for energy terms)
  COOK: includes kcal->K conversion
***/
{
  LJ->comb_rule=comb_rule | (LJ->alpha<=0);

  /* NB: BLEND uses kcal/mol, COOK uses K as energy units;
     since V2.8a, energy conversion is in cook/simdef.c */
  /* #IF COOK: LJ->EvdW*=(kcal/Eunit); */

  if (LJ->EvdW>0) {
    WARNING(("EvdW=%g is positive.\n\
*** (EvdW is energy minimum and should not be positive)\n\
*** continuing with EvdW=%g",LJ->EvdW,-LJ->EvdW))
    LJ->EvdW=-LJ->EvdW; }
  LJ->Esqrt=sqrt(-2*LJ->EvdW);

  /* Kirkwood combining rule */
  if (LJ->EvdW!=0) LJ->gamma=Sqr(LJ->alpha)/LJ->EvdW/Pow6(LJ->RvdW)/-256;
  else LJ->gamma=3e37; /* undefined -- solved later */
}

void combrule(pairparm_t *pp,siteparm_t *lj0,siteparm_t *lj1,char *info)
/*
  the combining rule: calculate van der Waals sig,eps from lj0,lj1
  note that the LJ potential is written as
    u(r)=-eps*[(sig/r)^12-2*(sig/r)^6]
  and that eps<0
  (a general version may have parm[] appended to the list of arguments)
*/
{
  if (lj0 && lj1) {
    double r,A,B;
    int comb_rule=lj0->comb_rule | lj1->comb_rule;

    if (comb_rule&2) r=Cub(lj0->RvdW*lj1->RvdW*4);
    else r=Pow6(lj0->RvdW + lj1->RvdW);

    if (comb_rule&1)
      B = lj0->Esqrt*lj1->Esqrt*r;
    else
      B = lj0->alpha*lj1->alpha / (lj0->gamma + lj1->gamma);
    A=B*r;

    if (B==0) pp->sig=0; else pp->sig=pow(A/B,1./6);
    if (A==0) pp->eps=0; else pp->eps=B*B/(-2*A); }

  if (info)
    prt("%s: eps=%.9g sig=%.9g",info,pp->eps,pp->sig);
}

/* calculate auxiliary constants for site-site interactions (former SS_DEF) */
void initssaux(pairaux_t *a,pairparm_t *pp)
{
  a->Sq=Sqr(pp->sig)*0.7937005259840998;
  a->E4=-4*pp->eps;
  if (a->Sq==0) {
    a->E24=a->E48=0; }
  else {
    a->E24=-24*pp->eps/a->Sq;
    a->E48=-48*pp->eps/a->Sq;  }

#ifdef POTCONST
    prt("a->Sq = %.15g",a->Sq);
    prt("a->E4 = %.15g",a->E4);
    prt("a->E24 = %.15g",a->E24);
    prt("a->E48 = %.15g",a->E48); 
#endif
}
