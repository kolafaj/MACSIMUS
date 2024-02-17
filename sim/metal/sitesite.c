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
struct Epow_s Epow[2]={{1,0,{1}},{1,0,{1}}}; /* EvdW is energy, RvdW is other (length) */
#else
#  define PROMPT "! "
#endif

/*** combining rules ***/

void combruleinfo(int comb_rule) /***************************** combruleinfo */
{
  if (comb_rule!=3)
    ERROR(("comb_rule=%d is not supported (now only 3 is supported)",comb_rule))
}

void initcombrule(siteparm_t *LJ,int comb_rule) /************** initcombrule */
{
  LJ->comb_rule=comb_rule;
  //#ifdef COOK
  //  LJ->EvdW*=(kcal/Eunit);
  //  LJ->parm[0]*=(kcal/Eunit); /* xi */
  //#endif
}

double globalLJshift;

void combrule(pairparm_t *pp,siteparm_t *lj0,siteparm_t *lj1,char *info)
/*
  the combining rule
  calculate sigvdW,epsvdW from lj0,lj1
  note that the LJ potential is u(r)=-epsvdW*[(sigvdW/r)^12-2*(sigvdW/r)^6]
  and that epsvdW is negative
  (a general version may have parm[] appended to the list of arguments)
*/
{
  if (lj0 && lj1) {
    int k;

    pp->sig=sqrt(lj0->RvdW*lj1->RvdW);
    pp->eps=sqrt(lj0->EvdW*lj1->EvdW);
    loop (k,0,SS_PARMS) pp->parm[k]=sqrt(lj0->parm[k]*lj1->parm[k]); }

  if (info) prt("%s: eps=%.8g sig=%.8g parm = %.8g %.8g %.8g",
                info,pp->eps,pp->sig,
                pp->parm[0],pp->parm[1],pp->parm[2]);
}

/* calculate auxiliary constants for site-site interactions (former SS_DEF) */
void initssaux(pairaux_t *a,pairparm_t *pp)
{
  a->r0=pp->sig;
  a->A=pp->eps;
  a->xi=pp->parm[0];
  a->p=pp->parm[1];
  a->q=pp->parm[2];
  a->f1=2*a->A*a->p/a->r0;
  a->f2=a->q/a->r0;
}
