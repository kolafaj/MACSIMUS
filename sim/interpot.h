#include "sitesite.h" /* normally in a subdirectory, e.g., lj/sitesite.h */

/*
  control module for definitions of site-site potentials
  (#included from sim/DIR/sitesite.c and sim/DIR/sitesite.h)

  this module depends on the following #defines:
     WORM, CUT?, TWODIM?, POLAR, PARALLEL, STARS

  NOTE: initcombrule() and combrule() are in sitesite.c
*/

#ifdef QQTAB 
#  include "elst.h"
#  define DOUBLEQQ /* empty - parameter qq not used */
#  ifdef POLAR
typedef struct qqtabpol_s {
  ertab_p tab;      /* table offset */
  erreal sgrid;     /* spline scaling */
  double qq;        /* passed to calculate tab */
  //  ertab_p from,to;  /* valid range of the table */
} qqtabpol_t;  
#  endif /*# POLAR */
#else /*# QQTAB  */
#  define DOUBLEQQ ,double qq
#endif /*#!QQTAB  */


typedef struct sitesite_s {
  real C1q,C2q;  /* squares of cutoffs */
  real A,A4;     /* for smooth cutoff */
  pairaux_t a;   /* former SS_DECLARE */
  real corr;     /* homogeneous cutoff correction term (sometimes not necessary) */
#ifdef SLAB
  struct Skk_s { real E,D; } *Skk;    
                 /* [slab.K] for slab cutoff correction, E and deriv: k-component */
#endif /*# SLAB */
#ifdef QQTAB
#  ifdef POLAR
  qqtabpol_t qqtab[2][2],qqtabx[2][2]; /* to be extended for qqtabx14... */
#  else /*# POLAR */
  ertab_p tab;      /* table offset */
  //  ertab_p from,to;  /* valid range of the table */
  erreal sgrid;     /* spline scaling */
  double qq;        /* passed to calculate tab */
#  endif /*#!POLAR */
#endif /*# QQTAB */
#if POLAR&32
  // TO BE OPTIMIZED LIKE HERE:
  //  enum fq_e { UNCHARGED, DRUDE_DRUDE, DRUDE_FQ, FQ_DRUDE, FQ_FQ } fq;
#endif /*# POLAR&32 */
  rdf_t *rdf;    /* pointer to the site-site rdf histogram */
#if PARALLEL==1
  rdf_t **rdfs;  /* [No.thread] per/thread instances */
#endif /*# PARALLEL==1 */
} sitesite_t;

/* NOTES:
   sstab[i][j] and sstab[j][i] are identical pointers 
   nsites is (inconsistently?) defined in cook/simdef.h
*/
extern sitesite_t **sstab,**sstab14; /* [nsites][nsites] */
extern real factor14,factor14_1;

#if PARALLEL==3
/* NOW NOT ACTIVE */
#  define GLOB ,pll_global_t *glob
#else /*# PARALLEL==3 */
#  define GLOB /*empty*/
#endif /*#!PARALLEL==3 */

#if PARALLEL==1 || PARALLEL==3
/* different instances of max14,En.el,En.vir for par threads */
typedef struct pll_global_s {
  double U;        /* to be returned from pll_linkcell */
  double Enel,Envir;
  double max14;
  int t;           /* exec time (clock ticks) */
  int ith;         /* thread (!SHM only); padd to 8-byte boundary */
#  if (PRESSURETENSOR&PT_OFF) || CACHELINE==128
  double Pvir[11]; /* incl. padding to 128 */
#  else /*# (PRESSURETENSOR&PT_OFF) || CACHELINE==128 */
  double Pvir[3];  /* may not be used, but padding needed anyway */
#  endif /*#!(PRESSURETENSOR&PT_OFF) || CACHELINE==128 */
} pll_global_t;
#endif /*# PARALLEL==1 || PARALLEL==3 */

#ifndef POLAR

/* common fo all versions (LINKCELL or not) */
double LJQQM   (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ GLOB);
void   QQM     (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ GLOB);
double LJM     (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 GLOB);
void   rdfM    (sitesite_t *ss,vector r1,vector r2);

#  ifdef PERSUM
double LJQQM_MM(sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ GLOB);
void   QQM_MM  (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ GLOB);
double LJM_MM  (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 GLOB);
void   rdfM_MM (sitesite_t *ss,vector r1,vector r2);
#  endif /*# PERSUM */
#  ifndef LINKCELL

void   LJQQ    (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ);
void   LJ      (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2);

#    ifdef PERSUM
void   LJQQ_MM (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ);
void   LJ_MM   (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2);
#    endif /*# PERSUM */

#    ifdef QQTAB
void   QQ      (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2);
#    else /*# QQTAB */
void   QQ      (               vector r1,vector r2,vector f1,vector f2,double qq);
#    endif /*#!QQTAB */

double LJQQ14M (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ GLOB);
void   LJQQ14  (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ);
/* XQQM, XQQ do not use qq-based splines even if QQTAB */
double XQQM    (vector r1,vector r2,vector f1,vector f2,double qq GLOB);
void   XQQ     (vector r1,vector r2,vector f1,vector f2,double qq);

#  endif /*# LINKCELL */

#  ifdef DEBUG
double LJQQtest(sitesite_t *ss,double qq,int noqq);
#  endif /*# DEBUG */

#else /*# POLAR */
double polarLJQQM(sitesite_t *ss,                          /***** polarLJQQM */
                  vector r1,vector r2,vector f1,vector f2,
                  siteinfo_t *si1, siteinfo_t *si2 GLOB);
void polarLJQQ(sitesite_t *ss,                              /***** polarLJQQ */
               vector r1,vector r2,vector f1,vector f2,
               siteinfo_t *si1, siteinfo_t *si2);
double polarXQQM(sitesite_t *ss,                            /***** polarXQQM */
                  vector r1,vector r2,vector f1,vector f2,
                  siteinfo_t *si1, siteinfo_t *si2 GLOB);
double polarLJQQ14M(sitesite_t *ss,                      /***** polarLJQQ14M */
                  vector r1,vector r2,vector f1,vector f2,
                  siteinfo_t *si1, siteinfo_t *si2,
                  int ifLJ,double qfactor GLOB);
#endif /*#!POLAR */

#undef GLOB

#ifdef GOLD
void goldQQ(vector f,double z, double qq);
#endif /*# GOLD */

#ifdef METAL
extern struct eldens_s {
  int ns;
  vector *rref;
  real *rho;
  real **rhorho; 
} eldens;
#endif /*# METAL */
