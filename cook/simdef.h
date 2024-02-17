/*                                 cookdef.h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "sitesite.h"

#define LINELEN 512 /* for line[LINELEN], not quite widely used... (increased) */

extern int nsites,nspec,ngroups,nbfixes,n14;

/* site type group (support for distinguishing site types for rdf) */
typedef struct stitem_s {
  int sp; /* species */
  int i;  /* site */
} stitem_t;

typedef struct stgroup_s {
  int n; /* # of atoms in the group */
  stitem_t stitem[1]; /* [n] */
} stgroup_t;

typedef struct sitedef_s {
  char name[8];    /* symbolic site name (CHARMM atom) */
  stgroup_t *g;    /* NULL if not a group */
  real M;          /* molar mass of the site in g/mol before equalization */
  real charge;     /* charge of the first site of this type */
  int ncharges;    /* 0: no site of this type found
                      1: all sites of this type have the same charge
                      2: sites of this type have different charges */
  siteparm_t LJ[2];/* LJ[0]: potential parameters, LJ[1]: for 1-4 interactions 
                      GAUSSIANCHARGES: last parameter = Gauss charge sigma */
#ifdef POLAR
  real alphapol;   /* polarisability in A^3 - to be used in self-field */
  real chargepol;  /* auxiliary charge in e - to be used in self-field */
#  if POLAR&1
  real kappa;      /* repulsive antipolarization (in 1/e) */
  int rep;         /* 2=2nd in the pair shell-core pair, 0=none */
#  endif /*# POLAR&1 */
#  if POLAR&2
  real Esat;       /* saturation energy for hyperpolarizability */
#  endif /*# POLAR&2 */
#endif /*# POLAR */
#ifndef FREEBC
  real sfweight;   /* atom weight for structure factor, see also siteinfo_t */
#endif /*# FREEBC */
} sitedef_t;

extern nbfix_t *nbfix0;
extern double *sigvdWptr; /* =fixij(tau.i,tau.j)->onefour[0].sig */
nbfix_t *fixij(int i,int j);

/* the energy scaling factor in the ble-file, in p.u. = k_B * K 
   added to the ble-file in V2.8a 
   NB: names Eunit, energyunit are already blocked - see units.h 
*/
extern double eunit;

extern sitedef_t *sitedef;

vector *initr(int sp);

int atomn(char *name,int sp,int i);
void readblend(void);

void initss(double LJcutoff);

#ifndef NIBC
#  ifdef QQTAB
#include "elst.h"
double setqq(ertab_p *tab,double qq,
#  if COULOMB<-2
             double sigma1,double sigma2,
#endif
             double factor14);
#endif
void initrspace(void);
#endif /*# NIBC */

#ifdef SLAB
extern struct wall_s {
  double rho;  /* density of LJ atoms [kg/m3] */
  double z[2]; /* wall position [0]=bottom wall, [1]=top wall */
  double g;    /* gravity */
  struct {
    /* wall-atom Lennar-Jones parameters (to replace the automatic setup
       based on the last site in the ble-file): */
    double sig;  /* LJ sigma in AA */
    double neps; /* product (number density)*(LJ epsilon) in K/AA3 */
  } *LJ; /* [nsites] */
  double scalez,minz; /* (PRIVATE) */
  double numden;      /* number density (PRIVATE) */
  double Punit_A;     /* constant for Pwall calculations */
  int n;       /* 1: wall 0, 2: wall 1, 3: both; n<0: atractive, n>0: repulsive */
  int is;      /* 1 if any wall or g (PRIVATE) */
} wall;

extern struct wsstab_s {
  double A,B;   /* atom-wall energy parameters */
  double AA,BB; /* atom-wall force parameters */
} *wsstab; /* [nsites] */

void initwall(void);
#endif /*# SLAB */
