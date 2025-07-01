#ifndef SIMGLOB_INCLUDED
#  define SIMGLOB_INCLUDED

#  ifndef SHAKE
#    define SHAKE 1
#  endif /*# SHAKE */

/*
#  define SERIAL
for debugging only - serial loop instead of starting threads in parallel */

extern double optionscaling;
#  include "simopt.h"

#  ifdef GAUSSIANCHARGES
#    define COULOMB -3
#    ifndef QQTAB
#      error "#define GAUSSIANCHARGES requires QQTAB  or 1"
#    endif /*# QQTAB */
#  endif /*# GAUSSIANCHARGES */

#  if defined(FREEBC)+defined(NIBC)+defined(COULOMB) != 1
#    error Exactly one of {FREEBC,NIBC,COULOMB} must be #defined. Check simopt.h!
#  endif /*# defined(FREEBC)+defined(NIBC)+defined(COULOMB) != 1 */

#  ifdef GOLD
#    error GOLD version not tested recently and must be reconsidered
#    define SLIT
#  endif /*# GOLD */

/* do not forget to #include "simopt.h" into simulation and project-specific
modules that do not #include "simglob.h" ! */

#  include "vector.h"

#  ifdef FLOAT
#    define realfmt "%f"
#  else /*# FLOAT */
#    define realfmt "%lf"
#  endif /*#!FLOAT */

#  ifndef PARALLEL
#    define PARALLEL 0
#  endif /*# PARALLEL */

#  if PARALLEL
#    include <pthread.h>
/* NB: include "sys4par.h" removed */
#  endif /*# PARALLEL */

#  if defined(WIDOM) && defined(FREEBC)
#    error WIDOM requires periodic boundary conditions but FREEBC selected
#  endif /*# defined(WIDOM) && defined(FREEBC) */

/* PRESSURETENSOR: sum of flags: */
#  define PT_VIR 1 /* diagonal virial (=configurational + constraint) part */
#  define PT_KIN 2 /* diagonal kinetic part */
#  define PT_ANY 3 /* diagonal full (both virial and kinetic) */
#  define PT_OFF 4 /* also off-diagonal terms (with PT_VIR and PT_KIN)
                 7 = PT_OFF|PT_ANY = full pressure tensor */
#  define PT_MOL 8 /* also molecular-based (center-of-mass), modifies PT_KIN */
#  define PT_MOM 16 /* SPECIAL: 2nd and 3rd moments of center-of-mass energy */
#  define PT_OFF_TR 32 /* DEBUG: with PT_OFF: transposed off-diag dependant term */
#  define PT_OFF_AV 64 /* DEBUG: with PT_OFF: average of transpose and standard
                          (for debugging dependants, see norm.c) */
#  ifndef PRESSURETENSOR
#    define PRESSURETENSOR PT_ANY
#  endif /*# PRESSURETENSOR */
#  if PRESSURETENSOR == PT_OFF
#    error PRESSURETENSOR&PT_OFF (off-diagonal) and none of PRESSURETENSOR&PT_VIR (virial) and PRESSURETENSOR&PT_KIN (kinetic)
#  endif /*# PRESSURETENSOR == PT_OFF */
#  if PRESSURETENSOR == PT_MOL
#    error PRESSURETENSOR&PT_MOL (off-diagonal) and not PRESSURETENSOR&PT_KIN (kinetic)
#  endif /*# PRESSURETENSOR == PT_MOL */
#  if PRESSURETENSOR&PT_OFF
#    define PT_DIM (2*DIM)
#  else /*# PRESSURETENSOR&PT_OFF */
#    define PT_DIM (DIM)
#  endif /*#!PRESSURETENSOR&PT_OFF */
#  if (PRESSURETENSOR & PT_MOM) && !(PRESSURETENSOR&PT_MOL)
#    error PRESSURETENSOR: higher moments of CM velocity requested without CM velocity (PT_KIN)
#  endif /*# (PRESSURETENSOR & PT_MOM) && !(PRESSURETENSOR&PT_MOL) */

#  include "intrapot.h"

/*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     basic control structures                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Species are numbered by integers from zero.
Information on each species is stored in a structure of type `specinfo_t' which
contains number of sites, masses and charges of sites, number of constraints
and corresponding sites and lengths.  These structures are of variable length
and there is an array of pointers 'spec' to these structures.
Example:  spec[1]->si[2].charge  is charge of site 2 of species 1

Any of N molecules (numbered from 0 to N-1) is described by a structure of
type `molecule_t' which contains the species number, number of sites and
constraints, information where to find the physical configuration.  There is
an array of these molecule descriptors called `cfg'.

The physical configuration is stored in a special contiguous array of vectors
which is a part of a self-defined structure (`SDS', see "alloc.c") of type
`ToInt'.  There are several such self-defined structures containing higher
derivatives etc. as required by the Gear predictor-corrector method for solving
the differential equations.

The input information on interaction sites is stored in array of type
`sitedef_t' - see "cook/simdef.h"

9/96: dependants:
   these are auxiliary sites whose positions are defined as linear
   combinations of other sites.
   Let
     depend_t *d;
     struct depitem_s *dep=d->dep;

   then
     r[d] =       SUM       r[dep[j].i] * dep[j].w
                j<d->n

   and analogously for the forces
9/2010: dependants extended to out-of-plane and code changed
*/

/* WARNING: this differs from blend: NBFIX solved there */
enum exc_e { END=-1,EXCL=0,ONEFOUR=1 };

typedef struct exception_s {
  int indx;   /* atom index */
  enum exc_e type;
} exception_t;

typedef int pair_t[2];

typedef struct siteinfo_s {
  real mass;      /* mass of site [extended! - more efficient, more memory] */
  real imass;     /* 1/mass of site */
  real charge;    /* charge of site */
  real bondq;     /* squared lengths of constraints (bonds) */
  pair_t pair;    /* site1 and site2 of constraint */
#  if COULOMB<-2 && defined(SLAB)
  real esig;      /* GAUSSIANCHARGES: sqrt(2)*sigma (POLAR: the same sigmas) */
#  endif           /*# COULOMB<-2 && defined(SLAB) */
#  ifdef POLAR
  real chargepol; /* charge of auxiliary site to mimic polarisation */
  real qqpol;     /* chargepol^2 */
  real alpha_qq;  /* alpha/chargepol^2 where alpha=polarisability */
#    if POLAR&2
  real Esat;      /* saturation energy */
#    endif /*# POLAR&2 */
#    if POLAR&1
  real qkappa;    /* chargepol*kappa */
#    endif /*# POLAR&1 */
#    define QTYPE_SHELL1  1 /* (POLAR&1) 1st of the shell-core pair (anion,kappa!=0), must be 1 */
#    define QTYPE_SHELL2  2 /* (POLAR&1) 2nd of the pair (cation, repulsive); param. rep, must be 2 */
#    define QTYPE_CENTER  4 /* nonzero center charge - COS models have this 0 */
#    define QTYPE_DRUDE   8 /* nonzero pol (Drude) charge */
#    define QTYPE_CHARGE 16 /* nonzero total charge (=atom charge), sum of center+pol */
#    define QTYPE_FQ     32 /* fluctuating charge */
  int qtype;  /* 0 = uncharged nonpolarizable atom, otherwise sum of keys */
#    if POLAR&8
  real alphazz_qq;/* alphazz/chargepol^2 (z direction) */
  int tozz;       /* -1: no axial polarizability,
                     >=0: atom defining the direction */
#    endif /*# POLAR&8 */
#  endif /*# POLAR */
  int st;         /* site type - not in old versions! */
  int isLJ;       /* 1 if nonzero LJ or similar term (with any other atom, incl. nbfixes, not 1-4) */
  exception_t *exc;  /* pointer to the list of exceptions */
#  ifndef FREEBC
  real sfweight;  /* atom weight for structure factor, see also sitedef_t */
#  endif /*# FREEBC */
} siteinfo_t;

struct depitem_s {
  int i;   /* independent site in species */
  real w;  /* weight */
};

typedef struct depend_s {
  struct depend_s *next;   /* next in the list */
  int indx;                /* the dependent site */
  int n;                   /* number of reference sites (parents) */
  enum deptype_e { DEP_OLDM, DEP_M, DEP_R, DEP_L , DEP_N} type;
  /* DEP_OLDM = old syntax of DEP_M; No.depend[0] = total # of all dependants
     DEP_M = "Middle" (linear combination) dependants
     DEP_R = "Rowlinson", site perpendicular to a (generally flexible) triangle
     DEP_L = "Lone", general rigid triangle, dependant out of plane
     NB: order significant, operators < > used */
  struct depitem_s dep[3]; /* [n]
                              indices and weights of the parents
                              assumed: SUM dep[].w = 1 */
  /* DEP_L, for n=3 only (wz for DEP_R):
     local coord. system = (X,Y,Z), Z=X x Y (x=cross product)
     f_x=X.f, etc. */
  vector x;                /* X from 3 sites: X=sum_i x[i] r[i] */
  vector y;                /* Y from 3 sites: Y=sum_i y[i] r[i] */
  real wz;                 /* z-component of the site; DEP_R: LO distance */
  vector tx;               /* tx[a]*f_x added to f[a]_z */
  vector ty;               /* ty[a]*f_y added to f[a]_z */
} depend_t;

#  if defined(POLAR) && POLAR&32
/* fluctuating charge for water */
typedef struct fq4_s {
  struct fq4_s *next;
  int indx[4]; /* sites in order: H H' L L' */
  double AHH,ALL,AHL; /* interaction constants */
  double AHL2; /* 2 AHL */
  double H1H1,H1H2,L1L1,L1L2,HL,LH; /* Phi->q constants */
} fq4_t;
#  endif /*# defined(POLAR) && POLAR&32 */

typedef struct {
  int ns;           /* # of sites */
  int nc;           /* # of constraints */
  int N;            /* # of molecules */
  int config;       /* set if abs init config (from blend) */
  int group;        /* group number */
  int intra;        /* if the molecule has at least one bond,angle,torsion */
  int pot;          /* 0=default, 3=TIP3P, 4=TIP4P, 5=ST2 (optimized code) */
  bond_t *bond;     /* bonds [external constraints] table head */
  angle_t *angle;   /* angle table head */
  torsion_t *dihedral,*improper,*aromatics; /* heads of torsion tables */
  depend_t *dependants; /* table of lin. dependent sites */
#  if defined(POLAR) && POLAR&32
  fq4_t *fq4;       /* fluctuating charge */
#  endif /*# defined(POLAR) && POLAR&32 */
  double zero_energy; /* const to be added to the pot energy */
  double charge;    /* molecular charge */
  double mass;      /* molecular mass */
  char *name;       /* species name */
  siteinfo_t si[1]; /* [max{ns,nc}] - the last item */
} specinfo_t;

typedef struct {
  int sp;    /* species */
  int ns;    /* # of sites */
  int nc;    /* # of constraints */
  int ir;    /* relative address of 1st site in rp - in bytes (use rof) */
  int ig;    /* relative address of 1st g[amma] in gs */
#  ifdef ANCHOR
  int anchor;  /* OR of the following bits: */
                         /* MEASURE and PRINT (no constraint): */
#    define ANCHOR_r   0x01  /* position of center-of-mass */
#    define ANCHOR_v   0x02  /* velocity of center-of-mass */
#    define ANCHOR_x   0x04  /* acceleration of center-of-mass */
#    define ANCHOR_f   0x08  /* force to center-of-mass */
#    define ANCHOR_m   0x10  /* momentum of force (torque) to center-of-mass */

                         /* CONSTRAIN ("anchor" = keep constant): */
#    define ANCHOR_s  0x100  /* site (at r0) */
#    define ANCHOR_c  0x200  /* center-of-mass */
#    define ANCHOR_i  0x400  /* principial axis of the inertia tensor
                                (parallel to vector axis) */
#    define ANCHOR_a  0x800  /* axis given by center-of-mass and site */
#    define ANCHOR_p 0x1000  /* axis given by a pair of sites */
#    define ANCHOR_t 0x2000  /* axis perpendicular to a plane of 3 sites */
#    define ANCHOR_g 0x4000  /* center-of-mass of group(s) */
#    define ANCHOR_pos (ANCHOR_s|ANCHOR_c)
                         /* constrain any position BUT group CM */
#    define ANCHOR_axis (ANCHOR_i|ANCHOR_a|ANCHOR_p|ANCHOR_t)
                         /* constrain any axis */
  struct anchorgroup_s {
    struct anchorgroup_s *next;
    vector r0; /* the anchor point */
    vector cf; /* constrain force measured here */
    int ns;
    int site[1]; /* [ns] */ } *group;
  vector r0;   /* the anchor point for ANCHOR_{s,c} */
  vector axis; /* the axis direction (normalized; ANCHOR_{i,a,p,t} only) */
  int iaxis[3];/* site numbers for ANCHOR_{i,a,p,t}
                  iaxis[0] = site to keep for ANCHOR_s */
  int xyz;     /* ANCHOR_c applies to coordinates (OR): 1=x,2=y,4=z */
#  endif /*# ANCHOR */
} molecule_t;

/***
    Type of function computing forces between m1 and m2.
    It should return energy.
    Forces are summed to f1 and f2.
***/
typedef double pot2_t(
  vector f1[/*m1->ns*/], vector f2[/*m2->ns*/],
  molecule_t *m1, molecule_t *m2,
  vector *rp);

/***
    Type of function computing intramolecular forces of molecule m and eventually
    correcting for wrong site-site forces.  It should return energy.
    Forces are summed to f.
    Used if pairs are calculated by the link-cell method.
***/
typedef double pot1_t(
  vector f[/*m1->ns*/],
  molecule_t *m,
  vector *rp);

#  ifndef LINKCELL
#    define pot_t pot2_t
#  endif /*# LINKCELL */

#  define RPOFFSET 7 /* @(rp) - @(logs), or 1(logs) + 3(lambda) +3(aux) */

typedef struct /*SDS*/ {
  int  size;     /* whole struct in bytes */
  int  dep;      /* set if dependants calculated, unset if cfg changes */
  real logs;     /* log of the Nose variable s */
  vector lambda; /* log(L): active for derivatives,
                    cfg[0].lambda is always derived from box.L */
  vector offdiag;/* reserved for future non-rectangular box */
  vector rp[1];  /* contiguous array of No.s vectors (POLAR: 2*No.s) */
} *ToIntPtr,ToInt;
/* WARNING: better DO NOT change this struct, dirty tricks used inside cook
   - logs is the 1st real variable to be integrated
   - reals must follow until the end of the structure
   - MUST begin with 2 ints
   lambda,offdiag are new since V2.7a */

/***
    Macro  rof(M,RP)  returns the address of array of vectors containing
    positions of sites of molecule *M in array RP of all sites and momenta
    (typically RP=a[i]->rp).  In other words,  rof(M,RP)[i]  is the vector
    of 3 Cartesian coordinates of the i-th site in molecule *M.
    POLAR: see also polarrof
***/
#  define rof(M,RP) ((vector*)((char*)(RP)+(M)->ir))

/* sparse matrix M for one molecule: see simcg.c */
typedef struct M_s {
  int b;         /* 2nd subscript, b<=a */
  real imass;    /* 1/mass of the common site (for a!=b)
                    1/mass_i+1/mass_j (for a==b) */
  real Mab;      /* values of M[a,b] = SUM_i grad[i]c[a].grad[i]c[b] / m[i] */
} M_t;

typedef struct No_s {
  double mass;    /* total mass */
  double free_mass; /* mass in all movable molecules */
  double minmass; /* minimum mass (after equalization) */
  double charge;  /* total charge */
  double eps;     /* machine precision */
  double t0;      /* time at sweep start (for .cp) */
  double dt;      /* h*noint */
  double NPT;     /* 1+3/No.f for MTK thermostat+barostat; 1 for HooverNPT */
  double invf;    /* 1/No.f for MTK thermostat+barostat; 0 for HooverNPT */
  double Pkinq;   /* kinetic pressure correction factor, for NVT, site-based */
  double PkinqCM; /* kinetic pressure correction factor, for NVT, CoM-based */
  double M_T ;    /* mass of thermostat (Nose+MTK) */
  double M_Th;    /* 1/2 of the mass of thermostat (Nose+MTK) */
  double M_P;     /* the mass of barostat (Nose+MTK) */
  double M_Ph;    /* 1/2 of the mass of barostat (Nose+MTK) */
  double P;       /* P in p.u., for barostat */
  double bulkmodulus; /* bulkmodulus in p.u., for barostat */
#  ifdef PERSUM
  double molspan; /* max. molecule site-site distance in any rotation */
  int nimg[3];    /* number of images considered */
#  endif /*# PERSUM */
#  ifdef LINKCELL
  double percell; /* av. # of sites/cell */
  int cell[DIM];  /* # of cells per box */
  int occup;      /* report cell occupancies up to occup (occup=1 -> guess) */
#  endif /*# LINKCELL */
  int rotatefrom; /* do not rotate molecules < rotatefrom in initialization */
  double A;       /* sum of all polarizabilities (polarizability volumes), in AA^3 */
                  /* POLAR or ECC */
#  ifdef POLAR
  int pol;        /* # of dipole sites (used for Car-Parrinello) */
  int free_pol;   /* as above, in movable molecules
                     BUG: fixed molecules are not polarizable */
#  endif /*# POLAR */
#  if PARALLEL
  int th;         /* number of threads */
  pthread_t *thread; /* [No.th]: pthread-style threads */
  pthread_attr_t thread_attr; /* thread attributes */
  int measureserial; /* force measurements in serial (void now) */
#  endif /*# PARALLEL */
  int N;          /* total # of molecules */
  int s;          /* total number of sites (sum over all molecules) */
  int massy_s;    /* total number of massy sites (i.e., excl. dependants) */
  int c;          /* total # of constraints */
  int bonded;     /* if there are bonded terms: for 1-4 logic of linked-cell list */
  int free_N;     /* # of movable molecules */
  int free_s;     /* # of sites in movable molecules */
  int free_c;     /* # of constraints in movable molecules */
  int free_depend;/* total # of dependants in movable molecules */
  int maxc;       /* max # of constraints per molecule */
  int maxs;       /* max # of sites per molecule */
  int f;          /* effective # of degrees of freedom */
  int fx;         /* NPT: f+extra degrees of freedom */
  int conserved;  /* # of conserved degrees of freedom excl. Etot (i.e., momentum + angular momentum components) */
  int ncoord;     /* # of rescalable coordinates */
  int f0;         /* # of degrees of freedom excl. conserved quantities */
  int f_tr;       /* effective # of intermolecular degrees of freedom */
  int f0_tr;      /* # of intermolecular degrees of freedom excl. conserved quantities */
  int f_in;       /* effective # of intramolecular degrees of freedom */
  int nreal;      /* # of reals in vectors (POLAR: 2*) */
  int eq;         /* # of equations of motion */
#  ifdef LOG
  int first;      /* # of molecules to be included in En.first */
#  endif /*# LOG */
  int pred;       /* option('p')%10 : predictor of constraints or velocity */
  int depend[DEP_N]; /* total # of dependants of given type */
  int ion;        /* number of ions */
} No_t;
extern No_t No;   /* # of sites, constraints, degrees of freedom etc. */

typedef struct tau_s {
  double T;      /* Nose correlation time or friction thermostat time */
  double E;      /* correlation time to keep energy constant */
  double rho;    /* typical time for density changes to reach given density */
  double P;      /* isobaric friction-like correlation time */
  double CM;     /* center of mass thermostat correlation time */
#  ifdef POLAR
  double dip;    /* for Lagrange (Car-Parrinello) mechanical dipoles */
#  endif /*# POLAR */
#  ifdef SLAB
#    if SLAB & 1
  double L;      /* for NPzzE ensemble */
#    endif /*# SLAB & 1 */
#  endif /*# SLAB */
  double sat;    /* for automatic saturation calculation */
  double sig;    /* correlation time to change size parameter in nbfixes table */
  int i,j;       /* sites for i-j sigvdW parameter */
} tau_t;

/* to control different sources of pressure */
typedef struct P_s {
  double n;        /* without homogeneous cutoff corrections */
  double c;        /* with homog. cutoff corrections En.corr/V^2 (if active) */
} P_t;

/*** En_t: summary of cycle statistics ***/
typedef struct En_s {
  double pot;      /* potential energy, as sum of bond,LJ,el etc. */
  double kin;      /* kinetic energy (first calculated as 2*kin.en., later
                      halved to get the true kin.en.) */
  double tot;      /* total Hamiltonian (pot+kin+zero_energy+Nose degrees+corr) */
  double el;       /* electrostatic energy; for Ewald both k- and r-space
                      POLAR: excl. En.self */
  double usr;      /* see userforces() */
  double corr;     /* LJ cutoff correction: Ecorr=corr/V, Pcorr=corr/V^2 */
  double virc;     /* virial of constraint forces */
  double bonded;   /* sum of bond+angle+torsion terms */
#  ifdef LOG
  double LJmn;     /* mol-mol LJ energy, passed from charmm2 */
  double LJ;       /* Lennard-Jones energy (sum over all pairs) */
  double intra;    /* site-site intramolecular terms */
  double pot0,el0,LJ0,bonded0;  /* as above - No.first molecules */
  double potX,elX,LJX;          /* as above - all but No.first molecules */
#  endif /*# LOG */
  double fix;      /* extra terms to fix some atom positions */
  double kin_tr;   /* 2*translational kin. energy (from center-of-mass) */
  double kin_in;   /* 2*intramolecular+rotatioal energy */
#  ifdef POLAR
  double self;     /* self energy */
  double chi_hf;   /* optical frequency susceptibility by elst field sampling */
#    if POLAR&32
  double fqself;   /* self energy - fluctuating charge (NB: not in virial) */
#    endif /*# POLAR&32 */
  double kinpol,Tpol; /* mech.dip.: kin.energy and temperature */
  double polmaxerr,polstderr,polmaxdr;
                   /* SCF errors during predictor/corrector or iterations */
  double pollaststderr,pollastmaxerr;
                   /* SCF error of the last iteration (useful in the re-read mode) */
  double sumpolmaxerr,sumpolstderr;
                   /* sums over a cycle */
  double Polmaxerr,Polstderr;
                   /* SCF errors testSCF */
  double Fstderr,Fmaxerr;
                   /* errors of forces on Drude charges (vs. full SCF) */
  double polnit;   /* = scf.nit, in double */
#  endif /*# POLAR */
  double vir;      /* virial of force, using formula elstvirial = -elstenergy */
                   /* pressures are in 2 versions (see P_t)
                      all pressures are in [p.u] */
  double U;        /* internal energy = Epot + Ekin + Ecorr (no Nose etc.) [J/mol] */
  double Unc;      /* internal energy = Epot + Ekin (no cutoff corr, Nose etc.) */
  double H;        /* enthalpy, with cutoff correction (if active) */
  //  REMOVED double Hnc;      /* enthalpy, without cutoff correction */
  double T;        /* kinetic temperature */
  double kinT;     /* kin. energy of the thermostat degree of freedom (sometimes) */
  double kinP;     /* kin. energy of the barostat degree of freedom */
#  if PRESSURETENSOR&PT_VIR
  double Pvir[PT_DIM];  /* xx yy zz [yz zx xy] (virial part) [Pa] */
#  endif /*# PRESSURETENSOR&PT_VIR */
#  if PRESSURETENSOR&PT_KIN
  double Pkin[PT_DIM];  /* xx yy zz [yz zx xy] (kinetic part - site-based) [Pa] */
#  endif /*# PRESSURETENSOR&PT_KIN */
#  if PRESSURETENSOR&PT_MOL
  double PKin[PT_DIM];  /* xx yy zz [yz zx xy] (kinetic part - CM-based) [Pa] */
#  endif /*# PRESSURETENSOR&PT_MOL */
#  if PRESSURETENSOR&PT_MOM
  double PKin2[DIM];  /* SPECIAL: PKin^2 (to calc. 2nd moment = kurtosis of v) [Pa] */
  double PKin3[DIM];  /* SPECIAL: PKin^3 (to calc. 3rd moment) [Pa] */
#  endif /*# PRESSURETENSOR&PT_MOM */
#  if PRESSURETENSOR
  double Ptens[PT_DIM];  /* xx yy zz [yz zx xy] total [Pa] */
#  endif /*# PRESSURETENSOR */
  double Pref;     /* reference pressure [p.u.] monitored in col 5 of .cp
                      pressure used in a barostat (if active)
                      for thermostat="MTK" calculated every step */
  P_t Pevir;       /* pressure [p.u.] using elst virial (En.vir); virial=1
                      NB: inaccurate for CUTELST, wrong for Gaussian charges */
  P_t PdV;         /* pressure from virtual V change; virial=2 */

#  if (PRESSURETENSOR&PT_ANY) == PT_ANY
  P_t Ptr;         /* from pressure tensor = SUM(En.Ptens)/3; virial=3 */
  /*  double trPt;    SUM(En.Ptens)/3+En.corr/V/V
                      =P with correction, to put to .cp
                      NB: for cutoff electrostatics, the standard pressure uses Eelst=virial  */
  double Hz;       /* uncorrected enthalpy with uncorrected Pzz */
#  endif /*# (PRESSURETENSOR&PT_ANY) == PT_ANY */
  double T_in;     /* kin. T from intramolecular motions and rotations*/
  double T_tr;     /* kin. T from intermolecular motions (translations) */
  double r1;       /* error of constraints of predictor */
  double r2;       /* error of constraints after integration step */
  double r3;       /* error of constraints (of true configuration) */
  double v1,v2,v3; /* errors of velocity constraints */
#  ifdef SHEAR
  double Cv;       /* SUM v_i[xy] cos(2PIz/L) / SUM cos^2(2PIz/L) */
  double kinshear; /* kin E of shear flow */
#  endif /*# SHEAR */
  /* drift to be removed (calculated from the selected coordinates only) */
  double CMdrift;   /* CM (center of mass from En.center) drift */
  double vdrift;    /* velocity drift */
  double AMdrift;   /* angular momentum drift */
  double omegadrift;/* angular velocity drift */
  /* measurement (independent from above, code duplicated) */
  vector CM;        /* CM (center of mass from En.center), cf. lag.CM */
  vector LM;        /* linear momentum, cf. lag.LM */
  double TLM;       /* corresponding translational temperature */
  vector AM;        /* angular momentum, cf. lag.AM */
  double TAM;       /* corresponding rotational temperature */
  vector J;         /* current density in SI [A/m2] */
  vector Jg[3];     /* as above, by groups */
  vector M;         /* dipole moment of the box (determined outside Ewald) */
  vector Mg[3];     /* as above, by groups */
  double logs,dlogs; /* a[0]->logs, a[1]->logs/h */
  vector lambda,dlambda; /* a[0]->lambda, a[1]->lambda/h */
  double ext;       /* extended degrees of freedom (Nose etc.) kin+pot energy */

                    /* problems here: unify: SLAB/here */
#  ifdef RGYR
  vector RGm;       /* inertia tensor, diagonalized */
#  endif /*# RGYR */
#  ifdef SLAB
                    /* see also slab_s below */
  vector sym;       /* (x,y,z of slab center from sp<=slab.sym) - (x,y,z of slab center) */
  double Pt;        /* tangential pressure tensor (Ptxx+Ptyy-2*Ptzz)/3 [Pa] */
  double Ptv;       /* virial part of tangential pressure tensor (Pvxx+Pvyy-2*Pvzz)/3 [Pa] */
  double Pwall[2];  /* pressure on the walls (see wall.*) */
#    if SLAB & 2
  double WPhi;      /* PROBABLY NOT NEEDED... */
  double Wcleave;   /* dG (rev. work) per d sigma and unit area */
  double zcleave;   /* actual positions, in the units of box.L[2] */
  double fcleave[2];/* cleaving force */
#    endif /*# SLAB & 2 */
#  endif /*# SLAB */
#  if (COULOMB<=-1) && (PRESSURETENSOR&PT_VIR) && defined(SLAB)
  double YB[3];     /* Yeh-Berkowitz pressure correction */
#  endif /*# (COULOMB<=-1) && (PRESSURETENSOR&PT_VIR) && defined(SLAB) */
#  ifdef XSECTION
  double xsec[3];   /* for xyz, mode=3 */
#  endif /*# XSECTION */
#  ifdef ECC
  double ECC_U2;     /* ECC correction to Coulomb interaction */
  double ECC_U3;     /* charge or dipole self-energy */
  double ECC_virnel; /* virial of force without electrostatics */
  double ECC_Pcorr;  /* P2corr+P3corr [p.u.] */
  double ECC_P3ref;  /* U3/3V [Pa] */
  double ECC_Pvir;   /* Pvir - see the paper */
  double ECC_Pscaled;/* conventional pressure w/o any ECC term */
#  endif /*# ECC */
} En_t;

#  include "rdf.h"

/****** miscellaneous ******/
/* for some strange compatibility reasons,
   we do not #include "options.h" here */
#  ifndef option
extern int optionlist[32];
#    define option(X) optionlist[X & 31]
#  endif /*# option */

#  if defined(LOG)
#    define FROM No.first
#  else /*# defined(LOG) */
#    define FROM 0
#  endif /*#!defined(LOG) */

/*
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                     declarations of global variables                  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

extern
pot2_t ***pot;       /* array [nspec][nspec] of pointers to functions computing
                        INTERmolecular interactions */

#  ifdef LINKCELL
extern
pot1_t  **pot1;      /* array [nspec] of pointers to functions computing
                        INTRAmolecular interactions */
#  endif /*# LINKCELL */

extern
specinfo_t **spec;   /* array [nspec] of pointers to species definitions */

extern
M_t **Ms;            /* array [nspec or N] of pointers to matrices M */

/* max length of the Gear predictor, change initialization in simgear.c */
#  define MAXGEARORDER 9 /* affects: simgear.c, gear2.c, gear2pol.c */

extern
molecule_t *molec;   /* array [N] of molecule descriptors */
extern
ToIntPtr cfg[MAXGEARORDER+3]; /* [0]=positions, [1]=h*velocities, [2]=forces .. */
extern
ToIntPtr fixa;       /* fixed sites (by harmonic springs) */
extern
ToIntPtr forDUMP;    /* for debugging - see internp.C */

/* warning: even #ifdef FLOAT, the constraint calculations are in double */
extern
double **gs;         /* [option('p')] tables of Lagrange multipliers g to
                        predict their values
                     NOTE: pointers are NOT constant because swaping of cor-
                     responding arrays is implemented by pointer interchange */
extern
double *binom;       /* [option('p')] coefficients of this predictor */

extern
rdf_t ***rdf;        /* table of pointers to rdf's
                        see rdf.h for rdf_t
                        note: also in each sitesite_t sstab[][] (see setss) */

extern int init;        /* INITIAL KEY :
__init____.cfg___.cp____.sta .rdf .dih__logs
    0     read  append      read       read
    1     read  append      new        read
    2     read  rewrite     new         0
   >=3    init  rewrite     new         0  (sequential version only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
extern int init_append; /* = init but :=0 after one step (to control appending) */

extern double t;     /* MD time */
extern real h;       /* time step */
extern int virial;   /* virial of forces for final pressure:
                        0 = automatic setup
                        1 = uses elst_virial = -elst_energy
                        2 = uses virtual volume change (if dV and rescale set)
                        3 = uses the trace of the pressure tensor */

extern struct equalize_s {
  double mol; /* masses equalization factor for all molecules separately */
  double cfg; /* global equalization factor for given species */
  int sp;     /* specied to equalize; for sp<0, range 0..|sp| (incl.) */
} equalize;

extern double eps;   /* accuracy for computing Lagrange multipliers
                        if >=1 then # of iteration steps */
extern double omega; /* ? */

extern double epsc;  /* accuracy for correcting constraints
                        if epsc>=1, then # of iteration steps
                        (default is changed for -c4) */
extern double omegac;/* relaxation parameter for 'correct' */

extern int measuredrift; /* used by norm.c:removedrifts(),
                            set in simmeas.c:initCP() etc.:
                            1 = measure CoM (center of mass)
                            2 = measure linear momentum
                            4 = measure angular momentum */

#  ifdef POLAR
extern struct scf_s {
  /* public (in data): */
  double eps;   /* eps>0: accuracy for self-field (induced dipole moments)
                        in program units (p.u.= 0.00245 eA = 0.01175 D)
                   eps<0: |eps| iterations */
  double epsq;  /* eps enhancement factor for simulation (re)start
                   (close to 1:more, epsq<=0:off) */
  double omega; /* relaxation paramater for self-field iterations,
                   see also option -^ (default = max. stability or 0.95) */
  double epsx;  /* as above, for P by virtual volume change and similar
                   testing purposes (default 9e9 => epsp*(1-epspq)/10 */
  double omegax;/* as above, for P by virtual volume change and similar
                   testing purposes; default = 0.95  */
  int maxit;    /* max # of self-field iterations until WARNING
                   (for positive eps,epsx) */
  double domega;/* increment for automatic setup of scf.omega */
  double margin;/* scf.omega = scf.omega(DIVERGENCE) + scf.margin */
  double E;     /* for determination of epspol, in V/m (see also Eext) */
  int Estride;  /* epspol sampled once a Estride cycles only (0=1=every) */
  int test;     /* SCF (field, forces) accuracy vs. fully iterated solution */
  /* private: (see also En.*/
  double maxdr; /* Drude: max amplitude */
  double eps0;  /* for start (see epsq) */
  double maxerr;/* max dipole error passed from selffield; orig. maxerrold */
  double rate;  /* convergence rate (if defined) */
  int nit;      /* running number of iterations */
  int cycmaxit; /* for scfautoset(): max nit in a cycle */
} scf;
#  endif /*# POLAR */

extern double T;     /* thermostat temperature (if tau.T) [K] */

extern En_t En;      /* energies, virial, pressure tensor, constraint errors.. */

#  define RESCALE_X         1 /* rescale x */
#  define RESCALE_Y         2 /* rescale y */
#  define RESCALE_Z         4 /* rescale z */
#  define RESCALE_XYZ       7 /* rescale x,y,z */
#  define RESCALE_CM        8 /* center-of-mass based rescaling */
#  define RESCALE_PT       16 /* use PRESSURETENSOR (components) */
#  define RESCALE_XisY     32 /* modifier of RESCALE_PT (xscale=yscale) */
#  define RESCALE_XisYisZ  64 /* uniform (xscale=yscale=zscale) */
#  define RESCALE_H       128 /* rescale also logs,lambda */
#  define RESCALE_L       256 /* rescale also box */
                              /* RESCALE_WALL removed */
#  define RESCALE_SLAB   1024 /* SLAB version only (virtual area change) */
#  define RESCALE_CLEAVE 2048 /* z-box barostat w. cleaving boundaries */
extern int rescale;   /* what is allowed to rescale (with tau.P,tau.rho,dV):
                         sum of the flags above */

extern struct lag_s {
  int err;           /* lag for error analysis of important variables */
  int n;             /* number of blocked calculations: max.block size is 2^n */
#  ifdef RGYR
  int dim;           /* how many coordinates recorded (3 for vx,vy,vz) */
  int nv;            /* how many molecules */
  int v;             /* lag for velocity-velocity time autocorr. functions */
#  endif /*# RGYR */
  int ierr;          /* lag.err*noint */
  int in;            /* int(lag.n+log2(noint)+0.5) */
  int J;             /* lag for current time autocorr.f. (conductivity) */
  int M;             /* lag for dipole moment autocorr.f. */
  int CM;            /* lag for center of mass */
  int LM;            /* lag for linear momentum */
  int AM;            /* lag for angular momentum */
#  if PRESSURETENSOR&PT_OFF
  int Pt;            /* lag for Ptxy,... time autocorr.f. (viscosity) */
#  endif /*# PRESSURETENSOR&PT_OFF */
#  ifdef SPCTCF
  int tcf;           /* for special project SPC time correlation functions: individual */
  int TCF;           /* for special project SPC time correlation functions: sums */
#  endif
} lag;

extern double E;     /* total energy in [J/mol] (=Eunit) to be kept constant */
extern tau_t tau;    /* relaxation times for thermostat, volume change,... */
extern double T_tr_in; /* ratio T_tr/T_in for thermostat=13 */
extern int thermostat;
/* NOTE: do not change - relations < and > used, strchr used */
#  define T_BERENDSEN   1 /* Berendsen (simple friction) thermostat */
#  define T_NOSE        2 /* Nose-Hoover; with Verlet it uses predicted velocities */
#  define T_ANDERSEN    3 /* randomly selected atoms updated */
#  define T_MAXWELL     4 /* periodically all atoms updated at once */
#  define T_ANDERSEN_CM 5 /* randomly selected molecules (center-of-mass) updated */
#  define T_MAXWELL_CM  6 /* periodically all molecules updated at once */
#  define T_LANGEVIN    7 /* atom-based Langevin thermostat */
#  define T_LANGEVIN_CM 8 /* CM-based Langevin thermostat */
#  define T_BUSSI       9 /* canonical sampling through velocity rescaling */
#  define T_TR         11 /* friction only translations */
#  define T_IN         12 /* friction only intramolecular velocities */
#  define T_FRICTIONS  13 /* combined 11,12: see T_tr_in */
#  define T_DUMMY      19 /* (for inequalities only) */
#  define T_NPT        22 /* MTK NPT; WARNING all >= T_NPT are NPT*/
#  define T_HOOVERNPT  23 /* Hoover NPT */

extern int measure;  /* 0: only forces required
                        1: also energy etc.
                        2: special for calculating # of pairs (DEBUG only) */

extern int corr; /* made global in V3.6v */
/* 1=cutoff corrections included in final results, not in P for the barostat
   2=cutoff corrections included (also in barostat)
   4=kinetic pressure finite size correction (atom-based)
   8=kinetic pressure finite size correction (molecule-based)
  16=4+8 used if not NPT/MTK
  32=# of pairs is N^2/2, not N(N-1)/2
  64=subtract 1 from degrees of freedom for NVE,Berendsen
     (reversed in V3.1m, default in older versions)
 128=suppress warning */

extern int onefourinrdf; /* 1 if 1-4 pairs are included in RDFs */

/*** Ewald related global variables - see  comments in ewald.c and erfc.c ***/
/* they are used also for not-Ewald sums, so they are here... */
/* WARNING: change initialization as well! */
extern struct el_s {
  real alpha;     /* Ewald: r-space separation parameter in units of 1/AA;
                     cutoff electrostatics: sew point parameter */
  real kappa;     /* k-space cutoff in form kappa=K[i]/L[i] */
  real epsr,epsk; /* expected r- and k-space errors in forces, in K/AA */
  real minqq;     /* minimum expected charge-charge separation
                     (to test Drude and to report error of the erfc approximation) */
  real epsinf;    /* dielectric constant at infinity
                     (must be very big if free charges are present */
  real diff;      /* max rel.change of L before a warning about "unchanged K" */
  real Perr;      /* pressure tensor threshold (basic; WARNING for 10*Perr) */
  real rplus;     /* additional range to erfc splines (former ERFCPLUS)
                     needed for optimized water models,
                     does not hurt if unnecessarily large */
  real sat;       /* target saturation (see tau.sat) */
  int rshift;     /* &1: eru(r)=erfc(r)/r spline shifted to avoid jump at cutoff
                     &2: as above for the derivative, erd(r)
                     &4: Hammonds-Heyes method (normally with 1 and 2) */
  int grid;       /* # of grid points per unity of squared distance */
  int test;       /* test module switch (former etest) */
  int diag;       /* whether to include diagonal correction */
  int centroid;   /* 1=dipole moments (of ions) are w.r.t charge centroid */
  int bg;         /* 0: charged system is error
                     1: charged system is OK, Ewald: background, cut.elst=OK
                     2: charged system is OK, periodic b.c.: redistribution */
                  /* ---------- end of initialization in simglob.c ---------- */
  real alphar;    /* normally =alpha; may differ for the Hammonds-Heyes method */
  int sf;         /* structure factor switch:
                     0 = standard k-space Ewald sum
                     1 = sphericalized (cube only)
                     3 = 3D */
  int corr;       /* extended Yeh-Berkowitz correction, former el.slab:
                     full slab correction: corr=4
                     pressure and energy correction only: corr=4+8
                     [In-Chul Yeh and Max L. Berkowitz, JCP 111, 3155 (1999)] */
  vector E;       /* elst. field [V/m], was Eelst[] in main.c before V3.4g */
  vector phase;   /* phase: E(t)=E*cos(2*PI*(f*t-phase)), x,y,z separately */
  real f;         /* frequency [Hz]; f=0 is static field */
  vector B;       /* static magnetic field [T] */
  struct m_t {
    double m;       /* magnetic charge value */
    int sp;         /* species containing the magnetic dipole moment */
    int plus;       /* site number of positive monopole */
    int minus;      /* site number of negative monopole */
  } m;
#  ifdef ECC
  int ecc;        /* EXPERIMENTAL: Electronic Continuum Correction (ECC)
                     (aka Fast Dielectric Approximation)
                     0=no ECC, 1=ions, 2=dipoles -1,-2 = use ble charges*/
  double epsf;    /* for reference V */
#    if defined(POLAR) || defined(FREEBC)
#      error ECC with POLAR or FREEBC
#    endif /*# defined(POLAR) || defined(FREEBC) */
#  endif /*# ECC */
  int sfsize;     /* to check overflow */
  struct sfr_s {  /* for sphericalized (radial) structure factor */
    int nk;          /* # of independent vectors */
    real q;          /* squared structure factor */
  } *sfr;
  struct sf3d_s { /* for 3D structure factor */
    int k[DIM];      /* k-vector */
    real q;          /* squared structure factor */
  } *sf3d;
  double sumM;    /* sum of molecule dipole moments = max.cell moment
                     for saturation autoset, one cfg value (cf. Eext.sumM)
                     POLAR: includes induced moments */
  double xinf;    /* 1/(el.epsinf*2+1) */
  double epsq;    /* for neutrality test if automatic setup does not work */
} el;

extern struct box_s {
#  ifdef LINKCELL
  real over14;    /* For automatic excrlimit setup (if box.max14=0):
                     Factor by which the estimate of the change of the max
                     1-2,1-3,1-4 distance is multiplied.
                     Must be >1; the higher, the safer
                     negative: debug mode (prints distances and speeds)
                     1.5 is pessimistic for proteins,
                     but just enough for water */
  real max14;     /* max distance 14 = excrlimit set by user (not automatic) */
  real excrlimit; /* current radius of 1-2,1-3,1-4 interactions */
  real excrpred;  /* as above, predicted */
  real threebonds;/* for init check only */
  real max14v;    /* max (safe) speed of max14 change from step to step */
  real Max14;     /* max 1-4 distance found during the sweep */
  real rmin;      /* minimum x,y,z of (any site of) a molecule - user value */
  real r0;        /* minimum x,y,z of (any site of) a molecule - used */
  real r0limit;   /* to warn only once */
#  endif /*# LINKCELL */
  vector L;       /* size of the periodic box;
                     it's important that L is real, not double, for automatic
                     float<->double conversion in loadcfg() */
  vector center;    /* FREEBC: 0 (or CoM?), p.b.c.: Lh, SLAB: may be set to CoM */
                    /* droplet (x,y,z), cylinder (y,z), slab(z) position by autocenter */
#  if 0
  /* obsolete: */
  double Lx[4];   /* if tau.rho=-3: L[0] as a poly of t */
  double Ly[4];   /* if tau.rho=-3: L[1] as a poly of t */
  double Lz[4];   /* if tau.rho=-3: L[2] as a poly of t */
  double rho[4];  /* if tau.rho=-1 rho as a poly of t */
#  endif /*# 0 */
  vector Lh;      /* Lh=L/2 */
  real V;         /* Lx*Ly*Lz */
  real cutoff;    /* r-space cutoff */
  real cq;        /* =cutoff^2 */
} box;

extern int cache;    /* molecules in block for r-space sum */

/* Q_t *Q; is in ewald.h and ewald.c */

#  ifdef POLAR
extern
int polar_off; /* offset of auxiliary charges in the units of bytes */
#    define polarrof(M,RP) ((vector*)((char*)(RP)+(M)->ir+polar_off))
#  endif /*# POLAR */

#  ifdef SHEAR
extern double shear; /* max gradient of acceleration for cos(z) method */
#  endif /*# SHEAR */

typedef struct fixsites_s {
  struct fixsites_s *next;
  int from,to; /* C-style range of sites: [from,to) */
  int isr; /* whether r applies */
  vector r; /* position */
} fixsites_t;
extern fixsites_t *fixsites0;

extern int drift;
#  define DRIFT_X         1 /* center-of-mass */
#  define DRIFT_Y         2
#  define DRIFT_Z         4
#  define DRIFT_VX        8 /* momentum (velocity) */
#  define DRIFT_VY       16
#  define DRIFT_VZ       32
#  define DRIFT_AX       64 /* angular momentum drift */
#  define DRIFT_AY      128
#  define DRIFT_AZ      256
#  define DRIFT_WHEN   1536 /* mask for when remove drift */
#  define DRIFT_CYCLE     0 /* init/load and every cycle [default] */
#  define DRIFT_START   512 /* init/load only */
#  define DRIFT_STEP   1024 /* init/load and every step */
#  define DRIFT_SAVE   1536 /* init/load and before final save */
// #  define DRIFT_DEPEND 2048 /* dependants every step (former norm&4) */
                               /* since V3.6k set to always */
#  define DRIFT_AUTO   4096 /* auto set */
/* WARNING: do not change numbers - coded sometimes directly by numbers */

extern int conserved; /* number of conserved degrees of freedom */

extern struct center_s {
  int on;              /* flag: 2=to box center, 1=z-pos */
  vector K,K2;         /* force constant to center-of-box */
  vector r0;           /* range (force starts at center+-r0) */
  vector cmK;          /* global force of the center-of-mass to (box) center */
  double cmmass;       /* mass of cmn molecules (from FROM) */
  int cmn;             /* # of molecules affected; 0=all but FROM */
} center;

#  ifdef SLAB
extern struct slab_s {
  double grid; /* histogram bins/AA (if slab.max>0) or bins/range */
  double max;  /* range is [0,max] if given, or according to box sizes */
  double out;  /* remove molecule if > out from slab/trickle/drop center */
  double outc; /* "cluster size" to set outx (cheaply solve evaporated dimers) */
  int torem;   /* molecule to be removed (cf. slab.out or removemol.n) */
  int geom;    /* geometry: 2=slab (default); new in V3.5h: 1=cylinder, 0=droplet */
  int sp;      /* species to autocenter the slab:
                  negative: all species 0..|sp| (incl.)
                  non-negative: species sp
                  0x7fffffff: no autocenter */
  int sym;     /* species for slab asymmetry, key as slab.sp */
  int mode;    /* 1: calculate surface tension incl. cutoff corr.
                  2: with dV, forces const-volume scaling (test area method)
                  4: POLAR only: Drude dipoles redistributed to a grid */
  int prt;     /* output files, sum of flags (default=31=all)
                  1 = .cm.A-3.z
                  2 = .cm.kgm-3.z
                  4 = .site.A-3.z
                  8 = .site.kgm-3.z
                 16 = .q.eA-3.z */
  /* slab thermostat */
  double T;       /* temperature in [Tz0,Tz1]: negative = off */
  double Tz0,Tz1; /* range for T */
  /* slab cutoff corrections */
  double range;   /* slab cutoff correction: integration range in Lz, default = 8 */
  int K;          /* max k-vector (not incl.) for slab cutoff correction */
  /* -------- end of initialized variables --------- */
  double outx; /* effective out; outx=out-outc at start */
  double qT;      /* private: sqrt(slab.T/T); 0=off */
  /* stacking correction */
  struct slabext_s { /* see extenddpr() */
    int zero;     /* gas slab extended */
    int center;   /* liquid slab extended */
    int span;     /* length of repeating pattern */
  } ext;
  /* forces to organize the slab: */
#    define NCENTER 10 /* number of force terms */
  int n[NCENTER];      /* how many molecules is affected (input) */
  int nn[NCENTER];     /* molecules affected [range) */
  int ns[NCENTER];     /* affected sites: 0: all, >0: i<ns only, <0: i>=|ns| */
  double z[NCENTER];   /* center relatively to Lh[2] */
  double z0[NCENTER];  /* radius of zero-pot. range */
  double Kz[NCENTER];  /* z-force const, for molecules up to n[] */
  double Kz2[NCENTER]; /* z1=0: 2*Kz, z1!=0: step */
  double z1[NCENTER];  /* for slab/outside bias function (input) */
  double dz[NCENTER];  /* z1-z0, for z1!=0 */
#    if SLAB & 1
  double Lx[4]; /* Lx for NPzzE zone melting, see mainresc.c */
  double Ly[4]; /* Ly for NPzzE zone melting, see mainresc.c */
#    endif /*? SLAB & 8 */ /*# SLAB & 1 */
#    if SLAB & 4
  struct {
    double sig,epsrho;   /* integrated LJ, from z=0; [A], [K/AA^3] DEPRECATED */
    double A,B;
  } wall;
#    endif /*# SLAB & 4 */
} slab;

#    if SLAB & 2
extern struct cleave_s {
  /* cleaving potential (is a bell-like function): */
  int init;       /* if set then copy z0,z1 to z[] */
  int n;          /* number of cleaving forces: 0,1,2 */
  int i;          /* site susceptible to the cleaving potential;
                     i<0: all sites 0..|i| (incl.) are susceptible */
  /* -------- end of initialized variables --------- */
  double sigma;   /* half-width of the cleaving potential, in AA */
  double K;       /* 1/2 of force, in K/AA */
  double Ksigmah; /* cleave*sigma/2 */
  double z[2];    /* z-positions (used + stored in .cfg), in units of Lz */
  double z0,z1;   /* z-positions (initial only: see cleave.init) */
} cleave;
#    endif /*# SLAB & 2 */
#  endif /*# SLAB */

/* move to simmeasx.h ? */
extern struct diff_s {
  int mode; /* sum of flags: 1=MSD + diffusivities, 2=MSCD+conductivities */
} diff;

#  ifdef XSECTION
/* move to simmeasx.h ? */
extern struct xs_s {
  /* public */
  double grid;   /* grid of the mesh, in points per AA */
  double RvdW;   /* van der Waals radius of the test particle, additive */
  double Rscale; /* multiply van der Waals radii; default=2^(-1/6) */
  int mode;      /* 0..4; 16 for x-section, 8 for whole cfg: see simmeasx.c */
  int freq;      /* how often to calculate, in # of cycles (1=every cycle) */
  int maxs;      /* max. # of species !(mode&8) or sites (mode&8)*/
  int sizelimit; /* max size of the grid array */
  int maxsize;   /* that one actually encountered */
#    ifdef CLUSTERS
  int mincluster;/* ignore all smaller for X-section calculations */
  /* private */
  int *color;    /* bad style: copy of color from CLUSTERS, if also active */
  ToIntPtr A;    /* configuration with clusters not breaking periodic b.c. */
#    endif /*# CLUSTERS */
} xs;
#  endif /*# XSECTION */

extern struct Eext_s {
  int isE;        /* 1 if there is general nonzero E, -3=x, -2=y, -1=z */
  vector E;       /* in program units */
  vector phase;   /* in units of 2PI (phase=0->cos, phase=0.25->sin) */
  vector arg;     /* arg=2*PI*(f*t-phase); E(t)=E*cos(arg) */
  double f;       /* frequency, in program units */
  double Eepspol; /* for epspol, in prog. units */
  double sumM;    /* sum of molecule dipole moments passed (incl. induced)
                     = max. dipole moments of the cell, averaged value
                     cf. el.sumM (for 1 cfg) */
  int isB;        /* 1 if nonzero el.B */
  vector B;       /* B in the program units */
  vector Bh;      /* h*B in the program units */
  double m;       /* magnetic dipole in p.u. */
  int ism;        /* 1 if magnetic dipole and field present (cf. el.m.m) */
} Eext;

#  ifdef ANCHOR
extern struct anchor_s {
  FILE *f;        /* output file of anchor forces .anc */
  double rr;      /* max displacement, squared */
  double cos,sin; /* max cos,sin of rot angle */
  vector r0;      /* positions to anchor the CM of all remaining molecules */
  int xyz;        /* cordinates to anchor the CM of all remaining molecules (OR): 1=x,2=y,4=z */
  int col;        /* number of recorded variables */
  int i;          /* counter, up to col */
  struct anchorrec_s {
    double sum;     /* [col] variable accumulated here over noint measurements */
    double var;     /* sum/noint */
    char info[8];   /* for StaAdd */
  } *rec;
} anchor;
#  endif /*# ANCHOR */

#  if PARALLEL
  /* usage of the pthread-based parallel macro:

     PARALLELIZE(FUNCTION)

     will call option('\\')-times function FUNCTION, which is passed
     the thread number in the form of void* casted int*
     (thread 0 is the caller)

     Example:

     void *FUNCTION(void *arg)
     {
       int ithread=*(int*)arg;
       fprintf(stderr,"thread %d started\n",ithread);
       ...
       return arg;
     }
  */

#    ifdef SERIAL
/* SERIAL VERSION FOR DEBUGGING - "parallelized" code run serially */
#      define PARALLELIZE(FUNC,NTHREADS) { \
  int ith; \
  loop (ith,0,NTHREADS) FUNC((void*)(No.thread+ith)); }

#    else /*# SERIAL */
/* true parallel code run in threads */
#      define PARALLELIZE(FUNC,NTHREADS) { \
  int ith; void *rc; \
  loop (ith,1,NTHREADS) \
    if (pthread_create(No.thread+ith,&No.thread_attr,FUNC,(void*)(No.thread+ith))) \
      ERROR(("cannot create thread %d",ith)) \
  ith=0; FUNC((void*)No.thread); \
  loop (ith,1,NTHREADS) \
    if (pthread_join(No.thread[ith],(void**)&rc)) \
      ERROR(("thread %d failed (rc=%p)",ith,rc)) }
#    endif /*#!SERIAL */

#  endif /*# PARALLEL */

#endif /*# SIMGLOB_INCLUDED */
