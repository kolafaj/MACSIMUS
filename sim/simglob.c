/*
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                        some global variables                          %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include <stdio.h>
#include "simglob.h"

double optionscaling;/* set to option('_') */

pot2_t ***pot;       /* array [nspec][nspec] of pointers to functions computing
                        INTERmolecular interactions */

#ifdef LINKCELL
pot1_t  **pot1;      /* array [nspec] of pointers to functions computing
                     INTRAmolecular interactions */
#endif /*# LINKCELL */


specinfo_t **spec;   /* array [nspec] of pointers to species definitions */
M_t **Ms;            /* array [nspec or N] of pointers to matrices M */

molecule_t *molec;   /* array [N] of molecule descriptors */
ToIntPtr cfg[MAXGEARORDER+3]; /* [0]=positions, [1]=h*velocities, [2]=forces .. */
ToIntPtr forDUMP;    /* for debugging - see internp.C */

/* warning: even #ifdef FLOAT, the constraint calculations are in double */
double **gs;         /* array [No.pred]:
                        tables of Lagrange multipliers g to be predicted
                        NOTE: pointers are NOT constant because array swaping
                        is implemented by pointer interchange */
double *binom;       /* coefficients of this predictor */

rdf_t ***rdf;        /* table of pointers to rdf's
                        NOTE: also in each sitesite_t */

int init=0;          /* INITIAL KEY :
__init____.CFG___.CP____.STA___t____logs
    0     read  append  read  read  read
    1     read  append  new    0    read
    2     read  rewrite new    0     0
    3     init  rewrite new    0     0  (sequential version only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int init_append; /* = init but :=0 after one step (to control appending) */

No_t No;             /* # of sites, constraints, degrees of freedom etc. */

double t;            /* MD time */
real   h=0.001;      /* time step */
int virial;          /* virial of forces for final pressure:
                        0 = automatic setup
                        1 = uses elst_virial = -elst)energy
                        2 = uses virtual volume change (if dV and rescale set)
                        3 = uses the trace of the pressure tensor */

/* equalization factors and species */
struct equalize_s equalize={0,0,-0x7fffffff};

double eps=1e-7;     /* accuracy for computing Lagrange multipliers
                        if >=1 then # of iteration steps
                        (was 1e-5 prior V 2.0d, 1e-6 prior 3.6k) */
double omega=0;      /* ? */

double epsc=1e-7;    /* accuracy for correcting constraints
                        if epsc>=1, then # of iteration steps
                        (the default is changed for -c4);
                        accuracy for Shake (was 1e-5 prior V 2.0d, 1e-6 prior 3.6k) */
double omegac=-1.2;  /* relaxation parameter for 'correct' or
                        relaxation parameter for Shake */

#ifdef POLAR
             /* eps,epsq, omega, epsx, omegax, maxit,domega margin */
struct scf_s scf={1,0.8, -9,     3e33,   0.95, 30,   0,     -0.03};
#endif /*# POLAR */

double T=300;        /* Nose|friction temperature (if tau.T) */
double P=1e5;        /* pressure (friction, if tau.P) */

En_t En;             /* energies, virial, pressure tensor, constraint errors.. */
int rescale=RESCALE_XYZ|RESCALE_CM;

struct lag_s lag={32,15
#ifdef RGYR
                  ,3,0x7fffffff
#endif /*# RGYR */
}; /* err,nblocks,[dim,nv,v]... */

double E=0;          /* total energy in Eunit to be kept constant */
double T_tr_in=1;    /* ratio T_tr/T_in for thermostat=13 */
int thermostat=0;    /* see simglob.h */
tau_t tau;           /* thermostat, volume change,... */

int measure;         /* 0: only forces required
                        1: also energy etc.
                        2: special for calculating # of pairs (DEBUG only) */

#ifdef LINKCELL
struct box_s box = { /* default over14 = */ 8 };
#else /*# LINKCELL */
struct box_s box;
#endif /*#!LINKCELL */

#ifdef COULOMB
#  ifndef  SPLINE
#    define SPLINE -2
#  endif /*#  SPLINE */

#  if SPLINE==3
#    define GRID 32+32*(COULOMB==-1)
#  elif SPLINE==-2
#    define GRID 256+256*(COULOMB==-1)
#  else  /*#!SPLINE==3!SPLINE==-2 */
#    define GRID 512
#  endif /*#!SPLINE==3!SPLINE==-2 */

/*** Ewald related global variables - see  comments in ewald.c and erfc.c ***/
        /* double...................................................  int............................... */
        /* alpha,kappa, epsr,epsk, minqq,epsinf,diff,Perr,rplus, sat, rshift, grid,test,diag,centroid,bg */
#  if COULOMB==2 || COULOMB==0
/* cut-off electrostatics, MACSIMUS style */
struct el_s el={0.7,0,   0,   0,   1,    3e33,  0.05,1e7, 2.5,   0,   0,      GRID,  0, 0,   1,       0};
#  elif COULOMB==3
/* Fennel-Gezelter */
struct el_s el={0.2,0,   0,   0,   1,    3e33,  0.05,1e7, 2.5,   0,   3,      GRID,  0, 0,   1,       0};
#  elif COULOMB<=-1 && COULOMB>=-3
struct el_s el={0.2,0.2, 0.05,0.5, 1,    3e33,  0.05,1e6, 2.5,   0,   3,      GRID,-10, 0,   1,       0};
#  else  /*#!COULOMB==2 || COULOMB==0!COULOMB==3!COULOMB<=-1 && COULOMB>=-3 */
#    error "unsupported COULOMB (one of -3,-2,-1,0,2,3 expected)"
#  endif /*#!COULOMB==2 || COULOMB==0!COULOMB==3!COULOMB<=-1 && COULOMB>=-3 */
#  undef GRID
#else /*# COULOMB */
/* FREEBC, NIBC(not tested recently) */
struct el_s el={0.2,0.2, 0.05,0.5, 1,    3e33,  0.05,1e6, 2.5,   0,   3,      32,  -10, 0,   1,       1};
#endif /*#!COULOMB */

int cache=64;        /* block size for pair sums:
                        cache=1: simplest pair forces in 2 nested loops
                           (good for small systems)
                        cache>0: blocked (4 nested loops total)
                           (large systems of small molecules)
                        cache<=0: recursive triangulation
                           (large systems of polyatomic molecules) */

/* Q_t *Q; is in ewald.h and ewald.c */

#ifdef POLAR
int polar_off;       /* offset of auxiliary charges in the units of bytes */
#endif /*# POLAR */

#ifdef SHEAR
double shear=0;      /* max gradient of acceleration for cos(z) method */
#endif /*# SHEAR */

fixsites_t *fixsites0;
ToIntPtr fixa;

int drift=DRIFT_AUTO; /* controls how numerical drift of CM, momentum, ang.mom. is fixed */
int conserved=-1;     /* number of conserved degrees of freedom; -1=auto */

/* int norm; OLD: 0: default according to b.c. and fixsites0
             1: always remove drift
            -1: never remove drift of momentum
            0,2,-2: the same for rotations and angular momentum */

struct center_s center;

#ifdef SLAB
/*                    grid max  out, outc, torem, geom, sp         sym,       mode prt  T  Tz0  Tz1   range K */
struct slab_s slab = { 0,  0,   0,   0,    -1,    0,    0x7fffffff,0x7fffffff,0,   31, -1, 0.25,0.75, 8,    0 };

#  if SLAB & 2
/*                        init,n,i */
struct cleave_s cleave = { 0,  0,-0x7fffffff };
#  endif /*# SLAB & 2 */
#endif /*# SLAB */

struct diff_s diff;

#ifdef XSECTION
/*              grid  RvdW Rscale             mode freq maxs  sizelimit     maxsize  */
struct xs_s xs={10.0, 1.4, 0.8908987181403393,0,   0,   -1,   256*1024*1024,0};
#endif /*# XSECTION */

struct Eext_s Eext;

#ifdef ANCHOR
struct anchor_s anchor;
#endif /*# ANCHOR */

int onefourinrdf=0;   /* control whether 1-4 included in RDF's */

int measuredrift=0; /* see simglob.h for info */

#  ifdef SHARPCUTOFF
double globalLJshift;
#  endif
