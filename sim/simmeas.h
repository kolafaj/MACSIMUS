/* measuring CP, rdf and dihedral angle distribution (#ifdef DIHHIST) */

extern int NCP;      /* # of items in a record - to grow */
extern float *CPbuf; /* as array of float[NCP] */

extern struct ssd_s {
  struct ssd_s *next;
  int indx[2]; /* two site indices for a distance (u.dist will be recorded)
		  indx[0]==-1: energy variable pointed to by u.q
                               will be recorded with unit conversion (if SI)
                  indx[0]==-2: variable pointed to by u.q will be recorded */
  int stat;    /* 1 if also statistics (+ in front of line in SIMNAME.cpi) */
  union distq_u {
    real dist;   /* the distance */
    real *q;     /* pointer to the quantity */
  } u;
  char name[24];
} *ssd0;

void initCP(int no,int nbit,char *col4,char *col5,int corr);
void saveCP(int icyc,time_t stoptime);

void initrdf(double grid,double cutoff,char *fn);
void advancerdf(void);
void rdfgraph(void);
void loadrdf(double rdfgrid);
void saverdf(double rdfgrid);
#if PARALLEL==1
void parallelizerdf(int replicate);
#endif /*# PARALLEL==1 */

#ifdef DIHHIST
extern dihhist_t *dihhead; /* dihhist_t is defined in intrapot.h */
void includedihhist(torsion_t *t, char type[4][6], int sp);
void initdihhist(int dihgrid);
double gauchetrans(int key);
void printdihhist(void);
void loaddih(int dihgrid);
void savedih(int dihgrid);
#endif /*# DIHHIST */

void prtdielconst(double V,int ewald); /* dielectric constant analysis */
void prtviscosity(void);


#ifdef RGYR
/* rad of gyration, end-to-end distance + dipolemoments */
double measureplus(double *endtoend,int end[2],int I);
#else /*# RGYR */
/* dipolemoments only */
void measureplus(void);
#endif /*#!RGYR */

void ssdistance(void);

#ifdef HARMONICS
#  include "harmonic.h"
void initharmonics(double gr,double cutoff,double T,double rho);
void measureharmonics(void);
void loadharmonics(double gr,double cutoff);
void saveharmonics(void);
#endif /*# HARMONICS */

#ifdef WIDOM
double XWidom(int sp,int spvirt);
#  ifdef SLAB
double Widom(int sp,int ncyc,double z0,double z1,double dz,int sym);
#  else /*# SLAB */
double Widom(int sp,int ncyc,double corr);
#  endif /*#!SLAB */
#endif /*# WIDOM */

#ifdef SLAB
void initdpr(double grid,double zmax,char *fn);
void loaddpr(double dprgrid);
void savedpr(double dprgrid);
void measuredpr(int slabmode); // slab.mode&8 needed for POLAR not Gaussian
void printdpr(int slabmode,int slabprt,double dV,char *PdVname);
double slabcorr(int isP);
void extenddpr(int slabmode,int slabprt,double dV,char *PdVname,
               int nzero,int ncenter,int nspan);
int toremove(double rr);
#endif /*# SLAB */

/* fool-proof dipole moment (for debugging, POLAR+DRUDE) */
void dipmom(ToIntPtr A,char *info);

#ifdef POLAR
/* because polarrof(mn,cfg[1]->rp) is NOT velocity */
extern vector *lastrpols;
#endif
