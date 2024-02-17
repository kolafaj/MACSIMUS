/* C2c lc.C ; cppnest lc.c

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! DO NOT EDIT lc.c               !!!
  !!! * edit lc.C and run C2c lc.C   !!!
  !!! * optionally, run cppnest lc.c !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   pair forces by the linked-cell list method
   #included by cook/forces.c #ifdef LINKCELL
*/

#include "elst.h"
#include "linklist.h"

#ifdef GOLD
#  error GOLD for LINKCELL not implemented
#endif /*# GOLD */

#ifdef POLAR
#  if !defined(SS_MEASURE_rep) || !defined(SS_NOMEASURE_rep)
#    if POLAR&1
#      error "missing SS_MEASURE_rep or SS_NOMEASURE_rep for POLAR&1"
#    else /*# POLAR&1 */
#      define SS_MEASURE_rep /*empty*/
#      define SS_NOMEASURE_rep /*empty*/
#    endif /*#!POLAR&1 */
#  endif /*# !defined(SS_MEASURE_rep) || !defined(SS_NOMEASURE_rep) */
#endif /*# POLAR */

#if PARALLEL==1

#  ifdef POLAR
#    define DOM(X) polarLJqqm(X,&lsloc,ssi,qi,qipol,glob);
#    define DO(X)  polarLJqq (X,&lsloc,ssi,qi,qipol,glob);
#  else /*# POLAR */
#    define DOM(X) LJqqm(X,&lsloc,ssi,qi,glob); /* now U is in glob */
#    define DO(X)  LJqq (X,&lsloc,ssi,qi,glob);
#  endif /*#!POLAR */
#  define ENEL glob->Enel
#  define ENVIR glob->Envir
#  define ENPVIR glob->Pvir
#  define MAX14 glob->max14

#else /*# PARALLEL==1 */

#  ifdef POLAR
#    define DOM(X) polarLJqqm(X,&lsloc,ssi,qi,qipol,&U);
#    define DO(X)  polarLJqq (X,&lsloc,ssi,qi,qipol);
#  else /*# POLAR */
#    define DOM(X) LJqqm(X,&lsloc,ssi,qi,&U);
#    define DO(X)  LJqq (X,&lsloc,ssi,qi);
#  endif /*#!POLAR */
#  define ENEL En.el
#  define ENVIR En.vir
#  define ENPVIR En.Pvir
#  define MAX14 max14
#endif /*#!PARALLEL==1 */

static real
  excrrlimit, /* excrlimit^2, sure distance limit for 1-2,1-3,1-4 exceptions */
  max14;      /* max 1-2,1-3,1-4 distance found in the configuration */

#ifdef POLAR
#  if POLAR&32
#    define MEASURE
#    include "ljqqfq.c"
#    undef MEASURE
#    include "ljqqfq.c"
#  elif defined(QQTAB)
#    define MEASURE
#    include "ljqqtabpol.c"
#    undef MEASURE
#    include "ljqqtabpol.c"
#  else /*#!POLAR&32!defined(QQTAB) */
#    define MEASURE
#    include "ljqqpol.c"
#    undef MEASURE
#    include "ljqqpol.c"
#  endif /*?!POLAR&32 */ /*#!POLAR&32!defined(QQTAB) */
#endif /*# POLAR */

#if !defined(POLAR) || defined(WIDOM)
/* POLAR+WIDOM need both modules, but this is problematic anyway */
#  ifdef QQTAB
#    define MEASURE
#    include "ljqqtabnp.c"
#    undef MEASURE
#    include "ljqqtabnp.c"
#  else /*# QQTAB */
#    define MEASURE
#    include "ljqqnp.c"
#    undef MEASURE
#    include "ljqqnp.c"
#  endif /*#!QQTAB */
#endif /*# !defined(POLAR) || defined(WIDOM) */

static void update14(void) /*************************************** update14 */
/*
   update the 1-2,1-3,1-4 limit (excrlimit)
*/
{
  real v14=0,damp14;
  static double oldt=-99999;
  static unsigned pass=0,recalculate=4;
  static real oldmax14; /* previous max14 */
  static real eps14=1e-3; /* default addeed to excrlimit += box.over14*eps14 */
  double safety=0;

  if (box.max14==0 && box.over14==0) {
    box.max14=box.cutoff;
    WARNING(("box.max14 set to the cutoff %g (because box.over14=0)",box.cutoff)) }

  v14=t-oldt;
  oldt=t;
  if (fabs(v14)<h/20) return;

  pass++;

  max14=sqrt(max14); /* max 1-4 (or 1-3,1-2) distance found */
  Max(box.Max14,max14) /* global maximum */

  if (pass==1) {
    if (fabs(epsc)<1) eps14=1e-3+epsc*8;
    else eps14=1.2e-3;

    /* excrlimit calculated in simdef.c checked here */
    //    put2(excrlimit,box.cutoff)
    if (box.threebonds>box.cutoff)
      WARNING(("The pessimistic estimate of the max 1-4 distance (3*max_bond) exceeds the\n\
*** cutoff. The 1-4 interactions might be incorrect.\n\
*** Note: 1st step max 1-4 distance=%g",max14))
    /* ... the above just check, 2nd step still with box.cutoff */
    oldmax14=max14;
    return; }

  if (box.max14) {
    box.excrlimit=box.max14;
    return; }

  /* now pass>1 and oldmax14 is known */
  if (max14>box.cutoff*0.99)
    WARNING(("The maximum 1-2,1-3,1-4 distance found (%g) >~ the cutoff (%g).\n\
*** I can't handle models with cutoff>1-2,1-3,1-4 distance!",max14,box.cutoff))

  safety=(box.excrpred-max14)/(box.excrpred-oldmax14);
  if (box.excrlimit) {
    if (box.excrlimit>box.excrpred) {
      if (box.over14<0)
        prt("recalculate: max14=%g box.excrpred=%g pass=%d",max14,box.excrpred,pass);
      if (max14>box.excrpred)
        WARNING(("\
The restart of the max 1-2, 1-3, 1-4 distance predictor at pass=%d\n\
*** found a leak (= missed interaction): max14=%g > expected max14=%g.\n\
*** Unless this happens during unequilibrated start/change, stop the simulation\n\
*** and increase box.max14 or box.over14.\n\
*** I will recalculate the max distance.",pass,max14,box.excrpred)) }
    else
      if (safety<0.25) {
        WARNING(("\
The max 1-2, 1-3, 1-4 distance is close to the safety limit, the\n\
*** 1-2, 1-3, 1-4 forces might have been omitted!\n\
*** Ignore this message if the system is unequilibrated, e.g, swelling,\n\
*** or if vibrating bonds have been changed to fixed, etc.)\n\
*** I will recalculate the max distance.")) } }

  damp14=1-0.005/Sqr(box.over14); /* damping of old max14v */

  /* v14 = speed with which max14 changes
     note that diminishing max14 is still considered as
     an estimate of v14 because of time reversibility */

  v14=fabs(max14-oldmax14)/h;
  //  if (pass<4) v14*=1+1./pass; changed after leak found for EtOH
  if (pass<9) v14*=1.05+1./pass;
  box.max14v=fmax(damp14*box.max14v,v14);

  box.excrlimit=box.excrpred=max14+(box.max14v*h+eps14)*fabs(box.over14);
  if (pass>=recalculate) {
    box.excrlimit=box.cutoff;
    recalculate+=4+recalculate/4; }

  if (box.over14<0)
    prt("max14=%.8f v14=%g box.max14v=%g excrlimit=%.4f (pred=%.4f) safety=%.4f pass=%d",
        max14,v14,box.max14v,box.excrlimit,box.excrpred,safety,pass);

  oldmax14=max14;
} /* update14 */

linklist_t ****lastlist; /* last list */

#if PARALLEL==1
struct par1_lc_s {
  int by;     /* stride of the 1st slab loop, normally by=No.th */
#  ifndef BARRIER
  int dix;
#  endif /*# BARRIER */
  /* dix=0: inter-slab interactions, dix=1: neighboring slabs, etc. */
  int maxdix; /* maximum dix (distance of the 2 slabs processed) */
  linklist_t ****list;
  pll_global_t *glob0;
  int forceserial; /* NOTE: No.measureserial is leftover from PAR==3 and is no longer supported */
  pthread_barrier_t barrier; /* FOR THE barrier VERSION */
} par1_lc;

void *pll_linkcell(void *arg) /************************ pll_linkcell:threads */
/*
   Scans x-slabs |ixs|:
     for (ixs=ixs0; ixs<No.cell[0]; ixs+=by)
   where ixs0=ith=thread number, by=No.th=number of threads,
   and calculates the interactions with x-slab |ixs+dix|.
   The forces to slab |ixs| are summed in linklist_t->f1,
   the forces to slab |ixs+dix| are summed in linklist_t->f2
   (thus, no concurrent memory clash can occur).
*/
#else /*# PARALLEL==1 */
static double linkcell(linklist_t ****list) /************** linkcell: serial */
#endif /*#!PARALLEL==1 */
{
  linklist_t ***lxs, **ly,**lys,*ls;
  linklist_t lsloc;
  sitesite_t *ssi;
  real DX,dy,DY,dyq,dz;
  double qi;
#ifdef POLAR
  double qipol;
#endif /*# POLAR */
  vector cellr;
  int ix,ixm,ixs,ix1;
  int iy,iym,iys,iy0,iy1;
  int iz,izm,izs,iz0,iz1;
  vector cell,rL; /* cell size and rL=1/cell */
#if PARALLEL!=1
  double U=0;
#else /*# PARALLEL!=1 */
  clock_t clock0;
  int by=par1_lc.by;
  linklist_t ****list=par1_lc.list;
  int dix,ith=(pthread_t*)arg-No.thread;
  pll_global_t *glob;

#  ifdef BARRIER
  dix=ith/No.th;
  ith%=No.th;
#  else /*# BARRIER */
  dix=par1_lc.dix;
#  endif /*#!BARRIER */
  glob=par1_lc.glob0+ith;
#endif /*#!PARALLEL!=1 */

  VVV(cell,=box.L,/No.cell)
  VVV(rL,=No.cell,/box.L)
#if PARALLEL==1
  if (partimes.rspace) clock0=clock();
  for (ixs=ith; ixs<No.cell[0]; ixs+=by) {
#else /*# PARALLEL==1 */
  if (No.bonded) {
    max14=0;
    if (box.excrlimit==0) excrrlimit=Sqr(box.cutoff);
    else excrrlimit=Sqr(box.excrlimit); }
  loop (ixs,0,No.cell[0]) {
#endif /*#!PARALLEL==1 */

    /* "core" of the linked-cell list method, serial+parallel */

    lxs=list[ixs];
    loop (iys,0,No.cell[1]) {
      lys=lxs[iys];
      loop (izs,0,No.cell[2]) {
        for (ls=lys[izs]; ls; ls=ls->next) {
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
          /* central site *ls selected:
             draw sphere of radius cutoff around it and find matching cells */
          ssi=sstab[ls->si->st];
          qi=ls->si->charge;
#ifdef POLAR
          qipol=ls->si->chargepol;
#endif /*# POLAR */
          /* lsloc is the local copy of the central cell, *ls
             - needed if PARALLEL because we cannot change ls->r
               (could be used by other threads)
             - not needed in the serial version (cf. lcforce-nolsloc.c);
               for code simplicity however, we use it, too */
          lsloc=*ls;
          VVV(cellr,=rL,*lsloc.r) /* in cell units */
          ix1=cellr[0]+box.cutoff*rL[0];

          loopto (ix,ixm=ixs,ix1) {
            if (ix==ixs) DX=0;
            else DX=ix-cellr[0];
            if (ixm==No.cell[0]) {
              lsloc.r[0]-=box.L[0]; ixm-=No.cell[0]; }
            dyq=box.cq-Sqr(DX*cell[0]);

            if (dyq<0) {
              if (dyq/box.cq<-0.1)
                ERROR(("dyq=%g<-0.1*cq (atom with negative coordinate?)\n\
*** check box.rmin (maybe %g is not enough)",dyq,box.r0))
              if (dyq/box.cq<-1e-13)
                WARNING(("dyq=%g<0 (atom with negative coordinate?)\n\
*** check box.rmin (maybe %g is not enough)",dyq,box.r0))
              dyq=0; }

            dy=sqrt(dyq)*rL[1];

            if (ixm==ixs)
              iy0 = iym = iys;
            else {
              iym = iy0 = (int)(cellr[1]+No.cell[1]-dy)-No.cell[1];
              if (iym<0) {
                lsloc.r[1]+=box.L[1]; iym+=No.cell[1]; } }

            if (iym<0 || iym>No.cell[1])
              ERROR(("linked-cell core: cutoff too large or internal error"))

#if PARALLEL==1
            if (ix==ixs+dix) { /* (((((((((((((((( */
#endif /*# PARALLEL==1 */

            iy1 = (int)(cellr[1]+dy);

            loopto (iy,iy0,iy1) {
              if (iym==No.cell[1]) {
                lsloc.r[1]-=box.L[1]; iym-=No.cell[1]; }
              ly=list[ixm][iym];
              if (iy>iys) DY=iy-cellr[1];
              else if (iy<iys) DY=iy+1-cellr[1];
              else DY=0;
              dz=dyq-Sqr(DY*cell[1]);

              if (dz<0) {
                if (dz/box.cq<-0.1)
                  ERROR(("dz=%g<-0.1*cq (atom with negative coordinate?)\n\
*** check box.rmin (maybe %g is not enough)",dz,box.r0))
                if (dz/box.cq<-1e-13)
                  WARNING(("dz=%g<0 (atom with negative coordinate?)\n\
*** check box.rmin (maybe %g is not enough)",dz,box.r0))
                dz=0; }

              dz=sqrt(dz)*rL[2];
              iz1 = (int)(cellr[2]+dz);

#ifdef SLIT
              /* wall or slit pore: not periodic in z */
              if (iz1>=No.cell[2]) iz1=No.cell[2]-1;

              if (ixm==ixs && iym==iys) {
                iz0 = izs;
                if (measure) loopto (izm,iz0,iz1) {
                  if (izm==izs) /* the same cell */
                    DOM(lsloc.next)
                  else
                    DOM(ly[izm]) }
                else loopto (izm,iz0,iz1) {
                  if (izm==izs) /* the same cell */
                    DO(lsloc.next)
                  else
                    DO(ly[izm]) } }
              else {
                iz0 = (int)(cellr[2]+No.cell[2]-dz)-No.cell[2];
                if (iz0<0) iz0=0;

                if (measure) loopto (izm,iz0,iz1) {
                  DOM(ly[izm]) }
                else loopto (izm,iz0,iz1) {
                  DO(ly[izm]) } }
#else /*# SLIT */
              /* periodic b.c. in z */
              if (ixm==ixs && iym==iys) {
                izm = iz0 = izs;
                if (measure) loopto (iz,iz0,iz1) {
                  if (izm==No.cell[2]) {
                    lsloc.r[2]-=box.L[2]; izm-=No.cell[2]; }
                  if (izm==izs) /* the same cell */
                    DOM(lsloc.next)
                  else
                    DOM(ly[izm])
                  izm++; }
                else loopto (iz,iz0,iz1) {
                  if (izm==No.cell[2]) {
                    lsloc.r[2]-=box.L[2]; izm-=No.cell[2]; }
                  if (izm==izs) /* the same cell */
                    DO(lsloc.next)
                  else
                    DO(ly[izm])
                  izm++; } }
              else {
                izm = iz0 = (int)(cellr[2]+No.cell[2]-dz)-No.cell[2];
                if (izm<0) {
                  lsloc.r[2]+=box.L[2]; izm+=No.cell[2]; }

                if (measure) loopto (iz,iz0,iz1) {
                  if (izm==No.cell[2]) {
                    lsloc.r[2]-=box.L[2]; izm-=No.cell[2]; }
                  DOM(ly[izm])
                  izm++; }
                else loopto (iz,iz0,iz1) {
                  if (izm==No.cell[2]) {
                    lsloc.r[2]-=box.L[2]; izm-=No.cell[2]; }
                  DO(ly[izm])
                  izm++; } }
#endif /*#!SLIT */
              lsloc.r[2]=ls->r[2]; /* restore original */

              iym++; } /*iy*/

#if PARALLEL==1
            } /* ))))))))) */
#endif /*# PARALLEL==1 */

            lsloc.r[1]=ls->r[1];

            ixm++; } /*ix*/
          lsloc.r[0]=ls->r[0]; /* to omit this? (lsloc.r no longer needed...) */

          /* forces in lsloc must be copied back (if not indirected) */

#if PARALLEL==1
          /* lsloc.f2 not changed */
#  if !(LINKCELL&1)
	VV(ls->f1,=lsloc.f1)
#  endif /*# !(LINKCELL&1) */
          //#  IF !(LINKCELL&2)   VV(ls->f2,=lsloc.f2)
#  ifdef POLAR
#    if !(LINKCELL&1)
	VV(ls->f1pol,=lsloc.f1pol)
#    endif /*# !(LINKCELL&1) */
            //#    IF !(LINKCELL&2)   VV(ls->f2pol,=lsloc.f2pol)
#  endif /*# POLAR */

#else /*# PARALLEL==1 */

#  if !(LINKCELL&1)
	VV(ls->f,=lsloc.f)
#  endif /*# !(LINKCELL&1) */
#  ifdef POLAR
#    if !(LINKCELL&1)
	VV(ls->fpol,=lsloc.fpol)
#    endif /*# !(LINKCELL&1) */
#  endif /*# POLAR */
#endif /*#!PARALLEL==1 */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
        } /*ls*/
      } /*izs*/
    } /*iys*/

#if PARALLEL==1
  } /* ixs: selected c-slabs */
  if (partimes.rspace) glob->t+=clock()-clock0;
  //fprintf(stderr," r-space: thread %d clock=%d\n",ith,glob->t/10000);
  return arg;
#else /*# PARALLEL==1 */
//  prt("%d: U=%12.3f el=%12.3f JK",ixs,U,En.el);
  } /* ixs: all slabs */

  if (No.bonded) update14();

  return U;
#endif /*#!PARALLEL==1 */
}

#if PARALLEL==1
#  include "lcpar.c"
#endif /*# PARALLEL==1 */

static void pairforces(ToIntPtr B, ToIntPtr A) /***************** pairforces */
{
  int n,sp;
  vector *f,*rp=A->rp;
  molecule_t *mn;
  double Enel=En.el; /* V2.7c: to suppress rounding errors (En.el
                        is huge for POLAR) */

  En.el=0;

  /*** intramolecular forces ***/
  loop (n,0,No.N) if ((sp=(mn=molec+n)->sp)>=0) {
    f=rof(mn,B->rp);
    /* this is not correct but the total sum should be OK */
    En.pot+=(*pot1[sp])(f,mn,rp); }

#if PARALLEL!=2
                                                              CPUtime("intra");
#endif /*# PARALLEL!=2 */

#if PARALLEL==1
  /*** pair forces: parallel version of the linked-cell list called  ***/


  {
    /*
      to be able to run LJqq* in parallel:
      forces from atom1 summed to normal forces (or local copies), f1
      forces from atom2 summed to forces2
    */
    ToIntPtr forces2;
    alloczero(forces2,cfg[0]->size);

    lastlist=linklist(B->rp,forces2->rp,A->rp);
#  if LOG
	En.pot += En.intra+linkcell(lastlist);
#  else /*# LOG */
	En.pot += linkcell(lastlist);
#  endif /*#!LOG */

#  if LINKCELL==0
  {
    int ix,iy,iz;
    linklist_t *l;

    loop (ix,0,No.cell[0])
      loop (iy,0,No.cell[1])
        loop (iz,0,No.cell[2])
          looplist (l,lastlist[ix][iy][iz]) {
            VVV(l->fptr,+=l->f1,+l->f2)
#    ifdef POLAR
            VVV(((real*)((char*)(l->fptr)+polar_off)),+=l->f1pol,+l->f2pol)
              //            VVV(l->fptrpol,+=l->f1pol,+l->f2pol)
#    endif /*# POLAR */
          }
  }
#  elif LINKCELL==1
  {
    int ix,iy,iz;
    linklist_t *l;

    loop (ix,0,No.cell[0])
      loop (iy,0,No.cell[1])
        loop (iz,0,No.cell[2])
          looplist (l,lastlist[ix][iy][iz]) {
            VV(l->f1,+=l->f2) /* f2 is secondary local copy */
#    ifdef POLAR
            VV(((real*)((char*)(l->f1)+polar_off)),+=l->f2pol)
#    endif /*# POLAR */
          }
  }
#  endif /*#!LINKCELL==0!LINKCELL==1 */

    /* add forces2 to the forces */
#  if 0
    /* this is not good for POLAR */
    loop (n,0,No.N) if ((sp=(mn=molec+n)->sp)>=0) {
      vector *f1=rof(mn,B->rp);
      vector *f2=rof(mn,forces2->rp);
      ns=mn->ns;
      loop (i,0,ns) VV(f1[i],+=f2[i]) }
#  else /*# 0 */
    /* POLAR-proof */
    loop (n,0,No.nreal) ((double*)B->rp)[n] += ((real*)forces2->rp)[n];
#  endif /*#!0 */
    free(forces2);
  }

#else /*# PARALLEL==1 */
  /*** pair forces: serial version of the linked-cell list called  ***/
  lastlist=linklist(B->rp,A->rp);
#  if LOG
	En.pot += En.intra+linkcell(lastlist);
#  else /*# LOG */
	En.pot += linkcell(lastlist);
#  endif /*#!LOG */

#  if (LINKCELL&1)==0
  {
    int ix,iy,iz;
    linklist_t *l;

    loop (ix,0,No.cell[0])
      loop (iy,0,No.cell[1])
        loop (iz,0,No.cell[2])
          looplist (l,lastlist[ix][iy][iz]) {
            VV(l->fptr,+=l->f)
#    ifdef POLAR
            VV(((real*)((char*)(l->fptr)+polar_off)),+=l->fpol)
              //            VV(l->fptrpol,+=l->fpol)
#    endif /*# POLAR */
          }
  }
#  endif /*#!(LINKCELL&1)==0 */
#endif /*#!PARALLEL==1 */

  En.el+=Enel;
}
