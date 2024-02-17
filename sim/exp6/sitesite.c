#include "ground.h"
#include "vector.h" /* simglob.h not needed */
#include "sitesite.h"
#ifdef COOK
#  include "units.h"
#  define PROMPT
struct Epow_s Epow[2]={ {1,0,{0}}, {1,0,{0}}};
#else
#  define PROMPT "! "
#endif

/* exp-6 potential, LJ-like combining rules, cf. sim/buck and sim/busing */

void combruleinfo(int comb_rule) /***************************** combruleinfo */
{
  prt(PROMPT "Exp-6 potential, Lennard-Jones-like form:\n"
PROMPT "               /  6                       a           6 \\\n"
PROMPT "  u(r) = eps * | --- exp[a(1-r/sigma)] - --- (sigma/r)  |\n"
PROMPT "               \\ a-6                     a-6            /\n"
PROMPT "  (a is parm[0], eps = |E_min|)");
  if (comb_rule==1) prt(PROMPT "The Lorentz-Berthelot combing rule:");
  prt(PROMPT "comb_rule&2: atom-atom diameters (sigma) are %s",
      comb_rule&2?"given by the geometric mean":"additive");
  prt(PROMPT "comb_rule&1: atom-atom min. energies (epsilon) are given by the\n\
             %s",  comb_rule&1?"geometric mean":"NOT IMPLEMENTED");
  prt(PROMPT "the stiffness parameters a are given by the geometric mean");
}

void initcombrule(siteparm_t *LJ,int comb_rule) /************** initcombrule */
/***
  pre-calculates some constants (like gamma) for the combining rule
  comb_rule&1:
    0 = n.a.
    1 = square-root combining rule for energy terms
  comb_rule&6:
    0 = additive (arithmetic mean) diameters
    2 = geometric mean for diameters
  comb_rule for a pair is obtained by `or' operator, i.e., the geometric 
  means have higher priority
  (if LJ->alpha<=0 then also comb_rule applies for energy terms)
  COOK: includes kcal->K conversion
***/
{

  if ((comb_rule&1) == 0) {
    ERROR(("comb_rule&1 (Kirkwood combining rule) not implemented"))
    comb_rule|=1; }
  LJ->comb_rule=comb_rule | (LJ->alpha<=0);
#ifdef COOK
  LJ->EvdW*=(kcal/Eunit);
#endif
  if (LJ->EvdW>0) {
    WARNING(("EvdW=%g is positive.\n\
*** (EvdW is energy minimum and should not be positive)\n\
*** continuing with EvdW=%g",LJ->EvdW,-LJ->EvdW))
    LJ->EvdW=-LJ->EvdW; }
}

void combrule(pairparm_t *pp,siteparm_t *lj0,siteparm_t *lj1,char *info)
{
  if (lj0 && lj1) {
    int comb_rule=lj0->comb_rule | lj1->comb_rule;

    if (comb_rule&2) pp->sig=sqrt(lj0->RvdW*lj1->RvdW*4);
    else pp->sig=lj0->RvdW + lj1->RvdW;

    if (comb_rule&1)
      pp->eps = sqrt(lj0->EvdW*lj1->EvdW);
    else {
      ERROR(("comb_rule&1 (Kirkwood combining rule) not implemented"))
      pp->eps=0; }

    pp->parm[0]=sqrt(lj0->parm[0]*lj1->parm[0]); }

  if (info)
    prt("%s: eps=%.9g sig=%.9g alpha=%.9g",info,pp->eps,pp->sig,pp->parm[0]);
}

/* calculate auxiliary constants for site-site interactions */
void initssaux(pairaux_t *a,pairparm_t *pp)
{
  a->tau=-(pp->parm[0]/pp->sig);
  a->eu=6*pp->eps/(pp->parm[0]-6)*exp(pp->parm[0]);
  a->ef=a->eu*(pp->parm[0]/pp->sig);
  a->ru=TURNVDW pp->eps*pp->parm[0]*Pow6(pp->sig)/(pp->parm[0]-6);
  a->rf=TURNVDW 6*a->ru;
}
