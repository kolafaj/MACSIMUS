/* function forces()

   common code for LINKCELL (formerly in lcforce.c) and all-pair

   LINKCELL: add some tests for internal consistency:
   #define TEST

   CALLING SCHEME (as in constrd.c, rhs.c = el.eps-controlled):

   depend_r(A,0);
 #ifdef POLAR
   scf.nit=0;
   do {
 #endif
     zeroEn();
     forces(B,A);
 #ifdef POLAR
   } while (!selffield(B,A,scf.eps,scf.omega,1));
 #endif
   depend_f(A,B);

   FULL ITERATED:
 #ifdef POLAR
     scforces(locB,locA); // see scf.c
 #else
     zeroEn();
     forces(locB,locA);
 #endif
*/

#include "ground.h"
#include "sds.h"
#include "simglob.h"
#include "norm.h"
#include "ewald.h"
#include "simmeas.h"

#include "forces.h"
#include "xforces.h"

#include "siminit.h"  /* because of rhs .. move??? */
#include "intermac.h" /* because of fixsites */
#include "interpot.h" /* needed? */
#include "simdef.h"   /* needed? */
#include "simpot.h"   /* needed? */
#include "units.h"    /* needed? */
#include "simils.h"   /* dumpall(B,A) */

#include <time.h>
#include "cputime.h"

#ifdef METAL
#  include "eldens.c"
#endif /*# METAL */

#if PARALLEL==2
struct par2_data_s {
  ToIntPtr A,B,pair;
  clock_t (*times)[2];
} par2_data;

static void pairforces(ToIntPtr B, ToIntPtr A);

void *par2_func(void *arg) /************************************** par2_func */
{
  int ith=(pthread_t*)arg-No.thread;
  ToIntPtr B=par2_data.B;

  /* partimes.kspace and partimes.rspace always allocated at once */
  if (partimes.kspace) par2_data.times[ith][0]=clock();
  if (!ith) {
    Ewald(1,par2_data.B->rp,par2_data.A->rp);
    Ewald(2,par2_data.B->rp,par2_data.A->rp); }
  else
    pairforces(par2_data.pair,par2_data.A);
  if (partimes.kspace) par2_data.times[ith][1]=clock();

  return arg;
}
#endif /*# PARALLEL==2 */

#ifdef GAUSSIAN
#  include "gaussian.c"
#else /*# GAUSSIAN */

#  ifdef LINKCELL
#    include "lc.c"
#  else /*# LINKCELL */
#    include "pair.c"
#  endif /*#!LINKCELL */

/* standard MACSIMUS forces */
void forces(ToIntPtr B, ToIntPtr A) /******************************** forces */
/*
   B=forces, A=configuration
   if POLAR, should be finished by selffield to pass forces to the central atom
*/
{
  /* this should be somehow unified... */
  VV(box.Lh,=0.5*box.L)
  box.cq=Sqr(box.cutoff);

#  if PARALLEL==1
#    ifndef LINKCELL
#      error "PARALLEL==1 version requires LINKCELL which is not #defined"
#    endif /*# LINKCELL */

#    ifdef FREEBC
#      error "FREEBC cannot be combined with PARALLEL==1"
#    endif /*# FREEBC */

  /*** clear forces ***/
  sdszero(B);

  /*** Ewald k-forces (parallelized) ***/
  Ewald(1,B->rp,A->rp);
  En.el=Ewald(2,B->rp,A->rp);
  if (option('v')&4) prt("after Ewald k-space: En.el=%.12g (P1)",En.el);

                                                      CPUtime("Ewald k-space");
  /* pair forces (x-slab-parallelized) + intramolecular forces */
  pairforces(B,A);
                                                         CPUtime("pairforces");

#  elif PARALLEL==2
/* parallel version for 2 processors:
   Ewald and pair forces run in parallel */

#    ifdef FREEBC
#      error "FREEBC cannot be combined with PARALLEL==2"
#    endif /*# FREEBC */

#    ifndef COULOMB
#      error "PARALLEL==2 does not make sense without COULOMB"
#    endif /*# COULOMB */

  int i,n;
  ToIntPtr pair;
  double EnEwald;
  clock_t times[2][2];
  real *rB,*rpair;

  sdsalloczero(pair,B->size);
  sdszero(B);

  /* does not ... thread for Ewald must be 0 */
  if (!pll_ewald) allocone(pll_ewald);
  // if (!pll_ewald) allocarrayzero(pll_ewald,2);

  par2_data.A=A;
  par2_data.B=B;
  par2_data.pair=pair;
  par2_data.times=times;
  PARALLELIZE(par2_func,2)

  CPUtime("pair||k-space");

  if (measure) {
    EnEwald=pll_ewald->Ereturned;

#    if PRESSURETENSOR&PT_VIR
    VO(En.Pvir,+=pll_ewald->E)
    loop (i,0,PT_DIM) En.Pvir[i]-=0.5/PI*pll_ewald->Pvir[i];
#    endif /*# PRESSURETENSOR&PT_VIR */
  }

  /* summing up all forces */
  rB=(real*)B->rp;
  rpair=(real*)pair->rp;

  loop (i,0,No.nreal) rB[i]+=rpair[i];
  if (measure) En.el+=EnEwald;
  free(pair);

  if (partimes.kspace) partimes.kspace[0]+=(times[0][1]-times[0][0])/(double)CLOCKS_PER_SEC;
  if (partimes.rspace) partimes.rspace[1]+=(times[1][1]-times[1][0])/(double)CLOCKS_PER_SEC;

  if (option('v')&4) prt("after Ewald k-space: En.el=%.12g (P2, from stored value)",EnEwald);

#  else /*#!PARALLEL==1!PARALLEL==2 */
  /* serial version (and PARALLEL==3 - no longer active) */

  sdszero(B);

#    ifdef COULOMB
  /*** Ewald k-forces ***/
  Ewald(1,B->rp,A->rp);
  En.el=Ewald(2,B->rp,A->rp);
                                                      CPUtime("Ewald k-space");

  if (option('v')&4) prt("after Ewald k-space: En.el=%.12g (P0)",En.el);
#    else /*# COULOMB */
  En.el=0;
#    endif /*#!COULOMB */

#    ifdef METAL
  eldensities(A);
#    endif /*# METAL */
  pairforces(B,A);
  CPUtime("pair+in forces");

  /* end of serial version */
#  endif /*#!PARALLEL==1!PARALLEL==2 */

  /* these forces not parallelized...*/
#  ifdef SLAB
  wallforces(B,A);
  if (slab.K) slabcutcor(B,A);
#  endif /*# SLAB */

  fixforces(B,A);
  centerforces(B,A);
  elstforces(B,A);
  /* NOTES:
   * En.pot += En.el is AFTER forces() and optionally selffield() is called
                     (because selffield() modifies both En.el and En.self)
   * En.vir += En.el is in measureconstraints()
   */
  userforces(B,A);
  advancerdf();

  if (option('v')&128) dumpall(B,A);

} /* forces */
#endif /*#!GAUSSIAN */

#ifdef WIDOM
#  include "widom.c"
#endif /*# WIDOM */
