#ifndef _H_INTRAPOT_
#  define _H_INTRAPOT_

#  include "simopt.h"

/* this module contains intrapotential support + LJ combining rules */

#  ifndef TINY
#    define UREY_BRADLEY
/* i.e., no CHARMM22 in TINY version */
#  endif /*# TINY */

#  ifdef FLOAT
/* FLOAT implies IFLOAT */
#    define IFLOAT
#  endif /*# FLOAT */

#  ifndef IREAL
#    define IREAL
#    ifdef IFLOAT
typedef float ireal;
#      define irealfmt "%f"
#    else /*# IFLOAT */
typedef double ireal;
#      define irealfmt "%lf"
#    endif /*#!IFLOAT */
#  endif /*# IREAL */

#  include "int4.h"

typedef struct {
  ireal length; /* equilibrium bond length r0 */
  ireal K;      /* force constant:  U=K*(r-r0)^2 */
  real  Ki2;    /* -2*K */
#  ifdef MORSE
#    ifdef BLEND
#      error expanded MORSE not implemented
#    endif /*# BLEND */
  real a,Ka,Ka2,fa1,fa2; /* a, K*a, K*7/12*a^2, 3*a, 7*a^2/3 */
#  endif /*# MORSE */
} bondparm_t;

typedef struct {
  ireal angle;  /* equilibrium angle;
                    angle<0 means angle=PI and special algorithm */
  ireal K;      /* force constant:  U=K*(phi-angle)^2 */
  real K2;      /* 2*K */
#  ifdef UREY_BRADLEY
  ireal Kub;    /* Urey-Bradley force constant: U += K*(r13-length)^2 */
  ireal Kubi2;  /* -2*Kub */
  ireal length; /* Urey-Bradley equilibrium length */
#  endif /*# UREY_BRADLEY */
} angleparm_t;

#  ifdef DIHHIST
typedef struct dihhist_s {
  struct dihhist_s *next;
  unsigned4 *hist;   /* [grid] - total histogram */
  char type[4][6];   /* site types - must be 1st after *next ! */
  int4 sp;           /* species - valid only for DIHHIST==-1 */
  int4 grid;
  int4 gauche,trans; /* one-measurement values to be put to StaAdd */
  ireal gauche_trans;/* actually, gauche may be cis... */
  real q;            /* factor for angle->index conversion */
} dihhist_t;

extern FILE *dihcp;
#  endif /*# DIHHIST */

typedef struct {
  int n;         /* n=0: U=K[0]*(phi-K[1])^2
                    n={1,2,3,4,6}: U=|K[0]|+K[0]*cos(n*phi) [1 term]
                    n<0: U=SUM{0<=i<=-n} K[i] cos^i(phi)    [more terms] */
#  ifdef DIHHIST
  dihhist_t *dihhist;
#  endif /*# DIHHIST */
  ireal K[1];   /* force constant(s) (variable length array) */
} torsionparm_t;

extern double phi; /* if measure==2 then
                      angles and bond length are returned here */

double bondpot(
  vector r0,vector r1,
  vector f0,vector f1,
  bondparm_t *parm);

double anglepot(
  vector r0,vector r1,vector r2,
  vector f0,vector f1,vector f2,
  angleparm_t *parm);

typedef double torsionpot_t(
  vector r0,vector r1,vector r2,vector r3,
  vector f0,vector f1,vector f2,vector f3,
  torsionparm_t *parm);
extern torsionpot_t dihedralpot,improperpot;

/***
the following typedefs have been moved here because they are common for
blend and cook
***/

typedef struct any_s {
  struct any_s *next;
  int indx[2];
} any_t;

typedef struct bond_s {
  struct bond_s *next;
  int indx[2];
  bondparm_t parm;
} bond_t;

typedef struct angle_s {
  struct angle_s *next;
  int indx[3];
  angleparm_t parm;
} angle_t;

typedef struct torsion_s { /* both for dihedrals and impropers */
  struct torsion_s *next;
  int indx[4];
  char X;   /* # of X matches */
  torsionparm_t parm; /* ! may have variable length ! */
} torsion_t;

#endif /*# _H_INTRAPOT_ */
