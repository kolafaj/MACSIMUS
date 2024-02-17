/*
  non-bonded (atom-atom) potential support + combining rules
  see sitesite.h for more info
  see also sim/busing/
*/

#include "ground.h"
#include "vector.h" /* simglob.h not needed */
#include "sitesite.h"

#ifdef COOK

#  include "units.h"
#  define PROMPT
/* NEW in cook V2.8a: exponents of the energy unit */
struct Epow_s Epow[2]={
       {1,0,{1
#  if POW
              ,1
#  endif /*# POW */
    }},{1,0,{1
#  if POW
              ,1
#  endif /*# POW */
         }}}; /* EvdW is energy, RvdW is other (length) */

#else /*# COOK */
#  define PROMPT "! "
#endif /*#!COOK */

/*** combining rules ***/

void combruleinfo(int comb_rule) /***************************** combruleinfo */
{
#if POW==12
  prt(PROMPT "Exp-6+12 potential:\n"
PROMPT "  u_rep(r) = A_ij exp(-B_ij r)\n"
PROMPT "  u(r) = u_rep(r) - C_ij/r^6 + D_ij/r^12");
#elif POW==8
  prt(PROMPT "Exp-6-8 potential:\n"
PROMPT "  u_rep(r) = A_ij exp(-B_ij r)\n"
PROMPT "  u(r) = u_rep(r) - C_ij/r^6 - D_ij/r^8");
#else /*#!POW==12!POW==8 */
  prt(PROMPT "Exp-6 potential:\n"
PROMPT "  u_rep(r) = A_ij exp(-B_ij r)\n"
PROMPT "  u(r) = u_rep(r) - C_ij/r^6");
#endif /*#!POW==12!POW==8 */
  if (comb_rule==0) prt(PROMPT "The Waldman-Hagler combing rule:");
  if (comb_rule==3+8 || comb_rule==3+12) prt(PROMPT "The Busing alias Kong combing rule:");
  prt_(PROMPT "sqrt_rule&3=%d: B_ij=",comb_rule&3);
  switch (comb_rule&3) {
    case  0: prt("[2/(1/B_ii^6+1/B_jj^6)]^(1/6)"); break;
    case  1: prt("(B_ii B_jj)^1/2 (geometric)"); break;
    case  2: prt("(B_ii+B_jj)/2 (additive repulsion parameters B)"); break;
    case  3: prt("2/(1/B_ii+1/B_jj) (additive radii 1/B)"); break;
  }

  prt_(PROMPT "sqrt_rule&12=%d: A_ij=",comb_rule&12);
  switch (comb_rule&12) {
    case  0: prt("(A_ii A_jj)^1/2 B_ij^6/(B_ii B_jj)^3"); break;
    case  4: prt("(A_ii A_jj)^1/2 (geometric)"); break;
      /*    case  8: prt("(A_ii A_jj)^1/2 (arithmetic)"); break; */
    case  8: prt("[ A_ii*(A_jj*B_jj/A_ii/B_ii)^qi\n\
               + A_jj*(A_ii*B_ii/A_jj/B_jj)^qj ]/2\n\
               qi=B_ii/(B_ii+B_jj), qj=B_jj/(B_ii+B_jj) (Kong)"); break;
    case 12: prt("A_ij=(A_ii B_ii)^qi (A_jj B_jj)^qj / B_ij (Busing)\n\
             qi=1/B_ii/(1/B_ii+1/B_jj), qj=1/B_jj/(1/B_ii+1/B_jj)"); break;
  }
  prt(PROMPT "C_ij/r^6: geometric combining rule: C_ij=(C_i C_j)^1/2");
#if POW
  prt(PROMPT "D_ij/r^12: geometric combining rule: D_ij=(D_i D_j)^1/2");
#endif /*# POW */
}

void initcombrule(siteparm_t *sp,int comb_rule) /************** initcombrule */
/*
   for historical reasons, site parameters have the following names:
     sp->EvdW = A_ii, assumed in kcal/mol, converted to K (in simdef.c)
     sp->RvdW = B_ii, assumed in 1/AA
     sp->parm[0] = C_ii, assumed in kcal/mol AA^6, converted to K AA^6
     sp->parm[1] = D_ii, assumed in kcal/mol AA^12, converted to K AA^12
  */
{
  if ( ((comb_rule&12)==0) != ((comb_rule&3)==0) )
    ERROR(("Illegal combining rule comb_rule=%d\n\
*** Waldman-Hagler must be set either for both A and B or for none",comb_rule))

// FIXED in V2.8a:
// #ifdef COOK
//   sp->EvdW*=(kcal/Eunit);    /* A_ii */
//   sp->parm[0]*=(kcal/Eunit); /* C_ii */
//   /* WARNING: no conversion when reading nbfixes */
// #  if POW
//   sp->parm[1]*=(kcal/Eunit); /* D_ii */
//   /* WARNING: no conversion when reading nbfixes */
// #  endif
// #endif

  sp->sqA=sqrt(sp->EvdW);
  sp->sqB=sqrt(sp->RvdW);
  sp->Bi6=0.5/Pow6(sp->RvdW);
  sp->Bi=0.5/sp->RvdW;
  sp->sqABi3=sp->sqA/Cub(sp->RvdW);
  sp->AB=sp->EvdW*sp->RvdW;
  sp->sqC=sqrt(sp->parm[0]);
#if POW
  sp->sqD=sqrt(sp->parm[1]);
#endif /*# POW */
  sp->comb_rule=comb_rule;
}

void combrule(pairparm_t *pp,siteparm_t *lj0,siteparm_t *lj1,char *info)
/*
  for historical reasons, site-site parameters have the following names:
  pp->eps = A_ij
  pp->sig = B_ij
  pp->parm[0] = C_ij
  pp->parm[1] = D_ij
*/
{
  if (lj0 && lj1) {
    double aux=0,qi;
    if (lj0->comb_rule!=lj1->comb_rule)
      ERROR(("atom combining rules differ (%d vs. %d)",lj0->comb_rule,lj1->comb_rule))

    if ((lj0->comb_rule&3)!=2 && lj0->RvdW*lj1->RvdW==0)
      pp->sig=0;
    else switch (lj0->comb_rule&3) {
      case 0:
        aux=lj0->Bi6+lj1->Bi6;
        pp->sig=pow(aux,-1./6);
      break;
      case 1:
        pp->sig=lj0->sqB*lj1->sqB;
        break;
      case 2:
        pp->sig=(lj0->RvdW+lj1->RvdW)/2; /* NB: lj0->RvdW = B_i (not "diameter) */
        break;
      case 3:
        pp->sig=1/(lj0->Bi+lj1->Bi); /* Bi=0.5/B_i */
        break; }

    if (lj0->EvdW*lj1->EvdW==0)
      pp->eps=0;
    else switch (lj0->comb_rule&12) {
      case 0:
        pp->eps=lj0->sqABi3*lj1->sqABi3/aux;
        break;
      case 4:
        pp->eps=lj0->sqA*lj1->sqA;
        break;
        /*    case 8: pp->eps=(lj0->EvdW+lj1->EvdW)/2; break; VERSION <= 2.7q */
      case 8: /* Kong */
        qi=lj0->RvdW/(lj0->RvdW+lj1->RvdW);
        aux=lj1->AB/lj0->AB;
        pp->eps=( lj0->EvdW*pow(aux,qi) + lj1->EvdW*pow(aux,qi-1) )/2;
        break;
      case 12: /* Busing */
        qi=lj0->Bi/(lj0->Bi+lj1->Bi);
        pp->eps=pow(lj0->AB,qi)*pow(lj1->AB,1-qi)/pp->sig;
        break; }

#if POW
    pp->parm[1]=lj0->sqD*lj1->sqD;
#endif /*# POW */
    pp->parm[0]=lj0->sqC*lj1->sqC; }

  if (info)
#if POW
    prt("%s: eps=A=%.12g sig=B=%.12g C=%.12g D=%.12g",info,pp->eps,pp->sig,pp->parm[0],pp->parm[1]);
#else /*# POW */
    prt("%s: eps=A=%.12g sig=B=%.12g C=%.12g",info,pp->eps,pp->sig,pp->parm[0]);
#endif /*#!POW */
  /* NOTE: #ifdef GAUSSIANCHARGES, sigmaij=sqrt(sigmai^2+sigmaj^2) is calculated in gcelst */
}

/* calculate auxiliary constants for site-site interactions */
void initssaux(pairaux_t *a,pairparm_t *pp)
{
  a->Aij=pp->eps;
  a->Bij=pp->sig;
  a->C=pp->parm[0];
  a->C6=pp->parm[0]*6;
#if POW
  a->Bijq=Sqr(a->Bij);
  a->D=pp->parm[1];
  a->DPOW=pp->parm[1]*POW;
#endif /*# POW */
#ifdef SLAB
#  include "equivlj.c"
#endif /*# SLAB */
}
