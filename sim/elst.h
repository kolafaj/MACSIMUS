#ifndef ELST_H
#  define ELST_H

/*
Interface for r-space COULOMB energy and force calculations in
the periodic boundary conditions.
Macros eru,erd,erud are defined here.

Compile-time switch COULOMB
-3 Gaussian charges Ewald summation, r-space terms are splined
   for all charge-charge terms separately.
-2 The same as -1 except that splines are not used for charges close together.
-1 Point charges Ewald summation, one set of hyperbolic splines for r-space
   energy and one for r-space forces. Splines are used for all charge pairs.
   Bad if charges close together are present.
 0 MACSIMUS cut-off electrostatics (directly - no splines).
 1 (reserved)
 2: MACSIMUS cut-off electrostatics by quadratic splines (not recommended).
 3: Fennell and Gezelter [JCP 124, 234104 (2006)] cut-and-shift
    potential by quadratic splines.
*/

#  include "simopt.h"

#  ifndef COULOMB
#    error "module elst.h requires #define COULOMB"
#  endif /*# COULOMB */
#  ifdef FREEBC
#    error "module elst.h incompatible with FREEBC"
#  endif /*# FREEBC */
#  ifdef NIBC
#    error "NIBC probably wrong here"
#  endif /*# NIBC */

#  ifdef QQTAB
#    ifdef POLAR
/* ERFC is #defined in drudedrudeqqtab*.c as needed */
#    else /*# POLAR */
#      define ERFC (*ss)
#    endif /*#!POLAR */
#  else /*# QQTAB */
#    define ERFC Erfc
#  endif /*#!QQTAB */

#  ifdef FLOAT
typedef float erreal;
#  else /*# FLOAT */
typedef double erreal;
#  endif /*#!FLOAT */

#  if COULOMB<=-2
#    ifndef erexact_eps
#      ifdef POLAR
#        define erexact_eps 1e-15
#      else /*# POLAR */
#        define erexact_eps 1e-10
#      endif /*#!POLAR */
#    endif /*# erexact_eps */
#  endif /*# COULOMB<=-2 */

#  if COULOMB!=0
typedef struct ertab_s {
#    if SPLINE==3
  erreal Au,Bu,Cu,Du, Ad,Bd,Cd,Dd;
#    else /*# SPLINE==3 */
  erreal Au,Bu,Cu, Ad,Bd,Cd;
#    endif /*#!SPLINE==3 */
} ertab_t,*ertab_p;
#  endif  /*# COULOMB!=0 */

#  ifndef SUBGRID
#    define SUBGRID 1
#  endif /*# SUBGRID */

extern struct Erfc_s {
  /* data for calculating Ewald r-space sums. Not used in k-space (sim/ewald.c) */
  int grid;          /* d(r^2)=1/grid */
  erreal h;          /* 1/grid/SUBGRID */
#  ifdef erexact_eps
  erreal alphaq,A,B; /* derived from alphar */
#  endif /*# erexact_eps */
  erreal minr,maxr,alpha; /* alpha=alphar for the Hammonds-Heyes method */
#  if COULOMB>=0
  erreal shift,A,B,A3,B3,r0,r1; /* MACSIMUS short-range */
#    if COULOMB>=2
  ertab_p tab;       /* spline table; #ifdef QQTAB only passed */
  erreal sgrid;      /* spline scaling */
#    endif /*# COULOMB>=2 */
#  else /*? COULOMB>=2 */ /*# COULOMB>=0 */
  ertab_p tab;       /* spline table; #ifdef QQTAB only passed */
  erreal sgrid;      /* spline scaling */
  int shift;         /* splines: 1 = shifted to be continuous (tiny error) */
#  endif /*?!COULOMB==0 || COULOMB==2 */ /*?!COULOMB>=2 */ /*#!COULOMB>=0 */
} Erfc;

#  if COULOMB==0
/* cut-and-smooth 1/r directly (no splines) */

#    define ermacrodcl erreal byerd,er_x;
#    define eru(x) (er_x=sqrt(x), \
  er_x<ERFC.r0 ? \
  1/er_x-ERFC.shift : \
  Cub(er_x-box.cutoff)*(ERFC.A+ERFC.B*er_x))
#    define erd(x) (er_x=sqrt(x), \
  er_x<ERFC.r0 ? \
  er_x/Sqr(x) : \
  Sqr(er_x-box.cutoff)*(ERFC.A3+ERFC.B3*er_x-(er_x-box.cutoff)*ERFC.B)/er_x)
#    define erud(x) (er_x=sqrt(x), \
  er_x<ERFC.r0 ? \
  (byerd=er_x/Sqr(x),1/er_x-ERFC.shift) : \
  (byerd=Sqr(er_x-box.cutoff)*(ERFC.A3+ERFC.B3*er_x-(er_x-box.cutoff)*ERFC.B)/er_x,\
    Cub(er_x-box.cutoff)*(ERFC.A+ERFC.B*er_x)))

/* initerfc() defined in cutelstd.c */
void initerfc(int grid, double minr, double maxr, double cutoff, double alpha);

#  else /*# COULOMB==0 */

/* splines for Ewald r-space terms and Fennell-Gezelter */
#    define ermacrodcl ertab_p er_p; erreal byerd;

#    if SPLINE==-2
/* hyperbolic splines, good for erfc */
#      define eru(x) (er_p=ERFC.tab+(int)((x)*ERFC.sgrid), \
                      er_p->Bu/(er_p->Cu+(x))+er_p->Au)
#      define erd(x) (er_p=ERFC.tab+(int)((x)*ERFC.sgrid), \
                      er_p->Bd/(er_p->Cd+(x))+er_p->Ad)
#      define erud(x) (er_p=ERFC.tab+(int)((x)*ERFC.sgrid), \
                       byerd = er_p->Bd/(er_p->Cd+(x)) + er_p->Ad, \
                       er_p->Bu/(er_p->Cu+(x)) + er_p->Au)

#    elif SPLINE==2
/* quadratic splines, with cutelst.c (deprecated) or fgelst.c */
#      define eru(x) (er_p=ERFC.tab+(int)((x)*ERFC.sgrid), \
                      er_p->Au+(x)*(er_p->Bu+er_p->Cu*(x)))
#      define erd(x) (er_p=ERFC.tab+(int)((x)*ERFC.sgrid), \
                      er_p->Ad+(x)*(er_p->Bd+er_p->Cd*(x)))
#      define erud(x) (er_p=ERFC.tab+(int)((x)*ERFC.sgrid), \
                       byerd=er_p->Ad+(x)*(er_p->Bd+er_p->Cd*(x)), \
                       er_p->Au+(x)*(er_p->Bu+er_p->Cu*(x)))

#    elif SPLINE==3
/* cubic splines */
#      define eru(x) (er_p=ERFC.tab+(int)((x)*ERFC.sgrid), \
                      er_p->Au+(x)*(er_p->Bu+(x)*(er_p->Cu+(x)*er_p->Du)))
#      define erd(x) (er_p=ERFC.tab+(int)((x)*ERFC.sgrid), \
                      er_p->Ad+(x)*(er_p->Bd+(x)*(er_p->Cd+(x)*er_p->Dd)))
#      define erud(x) (er_p=ERFC.tab+(int)((x)*ERFC.sgrid), \
                       byerd=er_p->Ad+(x)*(er_p->Bd+(x)*(er_p->Cd+(x)*er_p->Dd)), \
                       er_p->Au+(x)*(er_p->Bu+(x)*(er_p->Cu+(x)*er_p->Du)))
#    endif /*#!SPLINE==-2!SPLINE==2!SPLINE==3 */

#    ifdef QQTAB
/*
   separate tables for all QQ pairs: Gaussian charges only (to be extended?)
   1-4 forces not implemented
*/
#      define erud14(x) 0
#      define erd14(x) 0
#      if COULOMB<-2
#        ifndef GAUSSIANCHARGES
#          error GAUSSIANCHARGES expected
#        endif /*# GAUSSIANCHARGES */
/* Gaussian charges: initerfc() defined in gcelst.c */
double initerfc(ertab_p *tab, double qq,
                double sigma1, double sigma2,
                double factor14,
                int grid, double maxr, double cutoff,
                double alpha, int shift, double kappa, erreal *alphar);
#      else /*# COULOMB<-2 */
/* initerfc() defined in erfcqq.c -- NOT TESTED */
double initerfc(ertab_p *tab, double qq,
                double factor14,
                int grid, double minr, double maxr, double cutoff,
                double alpha, int shift, double kappa, erreal *alphar);
#      endif /*#!COULOMB<-2 */
#    else /*# QQTAB */
#      if COULOMB>=0
/* initerfc() for no-Ewald (and no QQTAB) */
void initerfc(int Grid, double minr, double maxr, double cutoff, double alpha);
#      else /*# COULOMB>=0 */
/* initerfc() for Ewald (no QQTAB) */
void initerfc(int Grid, double minr, double maxr, double cutoff, double alpha, int shift,double kappa,erreal *alphar);
#      endif /*#!COULOMB>=0 */
#    endif /*#!QQTAB */

#  endif /*#!COULOMB==0 */

#  ifdef erexact_eps
double exacterud_sqrt(double xx,erreal *erd);
double exacterud_sqrt_1(double xx,erreal *erd);
#  endif /*# erexact_eps */

#endif  /*? ELST_H */ /*? COULOMB>=0 */ /*# ELST_H */
