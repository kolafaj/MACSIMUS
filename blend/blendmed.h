/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLENDMED.H

This module reads the `Molekyle Editor' definition of a molecule (species),
defines internal representation of species, and performs some checks
of molecule definiton consistency.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* "vector.h" must #include first */

char *tempname(void);

#define UNDEFATOM 9e5     /* marks undefined (WANTED) atoms */

enum exc_e {
  EXCL=0,      /* excluded pairs (normally 1-2 and 1-3) */
  ONEFOUR=1,   /* 1-4 interaction (if distance14!=4 then 1-distance14) */
  EXCDIM=2 };  /* to define array sizes */

enum keep_e {
  FREE=0,   /* atom can move while minimizing */
  INJAIL=1, /* position of atom is kept fixed while minimizing */
  WANTED=2, /* atoms with unknown positions that are to be calculated */
  FILLED=8, /* for recursive filling */
  POLICE=4  /* was special for `probe' - not used with the optimized version */
};

typedef struct exception_s {
  struct exception_s *next;
  int indx;   /* atom index */
  enum exc_e type;
} exception_t;

typedef struct probe_s {
  double d;     /* dx,dy,dz */
  double shell; /* shell around molecule */
  int i;        /* index of site */
  int ns;       /* number of sites of the probe */
  int show;     /* how many molecules/atoms to show */
} probe_t;

#ifdef POLAR
typedef struct axial_s {
  /* kappa must be 1st in both axial_t and isotropicparm_t (see blendss.c) */
  ireal kappa;       /* repulsive core parameter */
  axialparm_t parm;  /* the parameters */
  int tozz;          /* AXIAL: atom definig the z-direction */
} axial_t;
#endif

typedef struct site_s {
  char *id;          /* symbolic atom name */
  vector r;          /* coordinates of the site */
  double charge;     /* (partial) charge in e (? ireal: see water_s)*/
#ifdef POLAR
#  if POLAR==2 // or GAUSSIANCHARGES ?
#  error NOT FINISHED
  double COS;        /* charge on spring - for GAUSSIANCHARGES */
#  endif
  float mu;          /* for hotkey 'd' only */
  float sat;         /* for hotkey 's' only */
#define POL_ISO 1    /* isotropic polarizability */
#define POL_AXI 2    /* axial polarizability */
#define POL_SAT 4    /* saturated polarizability */
#define POL_REP 8    /* shell-core: repelling atom (usu not pol or cation) */
#define POL_SHL 16   /* shell-core: deformable dipole (usu anion) */
#define POL_NZISO 128 /* nonzero ans POL_ISO (temporary) */
  int polar;         /* type of polarizability, as above */
  void *pol;         /* POL_ISO: isotropicparm_t*, POL_AXI: axial_t* */
#endif
  enum keep_e keep;  /* status as regards keeping positions and filling H */
  int type;          /* site type refers to the `atom' table */
  signed char chir;  /* chirality (with respect to numbering of atoms) */
  signed char pch;   /* partial charge replacement status */
  signed char backbone; /* 1,2,3 for the backbone N CA C; 4 for carbonyl O */
  signed char nest;  /* tree nest (from CA) */
  signed char count; /* 1) multiplicity, like CG1 CG2 
                        2) 1=atom marked (click from GUI), 0=not marked */
  signed char rea;   /* rea status, see react() */
  int rsd;           /* residue number -- used by writePDB() etc. */
  int key;           /* residue key, see writePDB1() */
  int chain;         /* chain number */
  int clust;         /* cluster (submolecule) number */
  int nnbr;          /* # of neighbors connected by a bond */
  int nbr[MAXVAL];   /* list of neighbors connected by a bond */
  exception_t *exc;  /* pointer to the list of exceptions */
} site_t;

/* loop over all all neighbors *NBR around site SITE in the molecule */
#define loopnbr(NBR,SITE) \
  for (NBR=site[SITE].nbr; NBR<site[SITE].nbr+site[SITE].nnbr; NBR++)

extern int nspec;    /* number of species (molecules) as specified */
extern int newnspec; /* number of species after possible splitting */
extern int chargewarned; /* 0: will warn on fractional charge */

#define PAIRINFOLEN 72

typedef struct maxpair_s {
  double U;               /* max energy */
  char info[PAIRINFOLEN]; /* info string */
} maxpair_t;

extern double Uanglemax;
extern int ianglemax;
extern int PDBstyle;

typedef struct Xopt_s { /* see also initialization in blendmed.c !!!!! */
  ireal dr,Jeps,amplitude,T,core;
  ireal bdihedral[3]; /* dihedrals (phi,psi,omega), to build the chain */
  int ropt; /* set if more than 1 param in -rX:X: (e.g., toframe,byframe) */
  int toframe,byframe;
  int A,C,D,E,F,G,N,I,P,S,V,W; /* -Options */
  char *fn;
} Xopt_t;
extern Xopt_t Xopt;

typedef struct species_s {
  char *fn;               /* file name without extension (must be writable) */
  char *ext;              /* points to `.' in fn */
  char *info;             /* info header */
  int ns;                 /* number of sites (size of site[]) */
  int nclust;             /* number of clusters (submolecules) */
  int nchains;            /* number of (peptide) chains (backbones) */
  int wrmol;              /* whether to write *.mol */
  int edit;               /* effective (= if not 2D) option -e (edit) [int] */
  int *newindex;          /* renum. table if indices have been renumbered */
  double Emin;            /* the result of minimization */
  double zero_energy;     /* energy if set (if chem. reactions have occured) */
                          /*** values of options active for the species : ***/
  int opt_h;              /* option -h (bending angle mass limit) [int] */
  int opt_j;              /* option -j (jet - additional constraints) [0:1] */
  int opt_k;              /* option -k (keep in place key) [-3:3] */
  int opt_m;              /* option -m (minimize steps) [int] */
  int opt_p;              /* option -p (playback output) [int] */
  int opt_u;              /* option -u (print pair pot. energy limit) [int] */
  int opt_w;              /* option -w (write cfg) [-1:#] */
  int opt_y;              /* option -y (center-of-mass key) [0:3] */
  int reordw;             /* reorder water flag */
  Xopt_t Xopt;
  probe_t probe;          /* for probing */
  int frame;              /* frame for option -r4; toframe, byframe: Xopt */
  float L[3];             /* box from .plb; new in V2.1a */
  int N;                  /* option -n (# of molecules) [int] */
  int C1,C2;              /* cutoff values in A (options -t -[ -] ) */
  struct species_s *next; /* next in the list */
  site_t *site;           /* array[ns] of sites */
  maxpair_t *maxpair;     /* array[maxpairs] */
  int chargewarned;       /* not to repeat warning "charge fractional" */
} species_t;

extern species_t *spec0;  /* first in the list */

int chirality(int j,int k,int l); /* permutation sign */
species_t *readMOL(char *fn);     /* reads *.mol file in the MEDIT format */
species_t *readCHE(char *fn);     /* reads *.che file in the chemical format */
void randomcfg(species_t *spec);
int read3D(species_t *spec,int key);     /* read cfg */
int read2D(species_t *spec,const char *mode);  /* read 2D screen cfg */
void chiralities(species_t *spec,int assign); 
                          /* calc./check chiralities from config. */
extern char colors[17];
int atomcolor(int type);
void write3D(species_t *spec,int playback); /* write cfg in binary */
void perturb(species_t *spec,vector *r,int wave,double ran);
                          /* adds a perturbation to the molecule */
int checkneighbors(species_t *spec,int ns); /* returns 1 if bad */
double totalcharge(species_t *spec,int clust); /* charge of the molecule */
double roundedcharge(species_t *spec,int clust); /* to fix charge rounding */
void writeMOL(species_t *spec); /* write .mol */
void writeMARKED(species_t *spec,int key); /* write mark/keep info */
void readMARKED(species_t *spec,int key); /* read mark/keep info */

#ifndef TINY
void alphahelixold(species_t *spec); /* -r5 */
void alphahelix(species_t *spec);    /* -r6 */
void makechain(species_t *spec);    /* -r7 */
void writePDB1(species_t *spec); /* write .pdb, generate */
void writePDB2(species_t *spec); /* write .pdb, use id */
void averagecharges(species_t *spec, int *center); /* average charges on marked atoms */
void writeATM(species_t *spec); /* write .atm */
void writeCFG(species_t *spec); /* write .cfg, V2.7a and newer */
void writeCFGV26(species_t *spec); /* write .cfg, old format */
#endif

double masscenter(species_t *spec,int centermask);

void partialcharges(species_t *spec);
void cluster(void);

int connectsites(site_t *site, int n1,int n2);
void newemptyspec(void);
void appendspec(species_t *to,species_t *app);
