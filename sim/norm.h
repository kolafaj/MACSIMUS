#define COR_BONDS 0 /* correct bonds (Lagrange), also SHAKE */
#define COR_VELOC 1 /* correct velocity constraints (Lagrange) */
#define LAGR_MULT 2 /* Lagrange multipliers */
#define COR_BOTH  3 /* Lagrange: -c2, both (?) */
#define COR_LAST  4 /* for SHAKE autoset of omegac */
#define NCONSTRIT 5 /* size */

typedef struct {
  int nit[NCONSTRIT]; /* numbers of iterations, indices see above */
  double omegac;      /* active omegac (SHAKE) */
  double omega;       /* active omegac/2 (SHAKE) */
  double omegaclast;  /* last (refers to nit[COR_LAST]) */
  double d;           /* step for omegac optimization */
  double weight;      /* weight for omegac optimization */
} constrit_t;
extern constrit_t *constrit; /* [nspec] */  
void initconstrit(double omegac);
void doconstrit(double omegac);

extern double vconstrainterror;
double constrainterror(ToIntPtr A,ToIntPtr V);
int Scorrect(molecule_t *m, vector *r, vector *rpvel, double eps);
int shake_r(ToIntPtr A);
void Lcorrect(ToIntPtr B, ToIntPtr V);
int normalize(int m);
void unsplitpbc(int xyz);
void homogeneity(int corr);
void sortmolecules(int sort);
void bounddr(double drmax);
void setdrift(void);

char *prtxyz(int xyz);

void CoM(vector CM,ToIntPtr A);
void removedrifts(int pr); /* see also measuredrift() ? */

double rescalecfg(ToIntPtr A,int mode, double q1,double *q);
int rhorescale(int rescale,double tau,vector finalL,double halftcyc);
double scaling(double tau,int noint,double factor,double maxscale);

void addEnnit(double *En_nit,int n,double add);
void distancecheck(void);
void zeroEn(void);

void depend_r(ToIntPtr A,int always);
void depend_f(ToIntPtr A,ToIntPtr B);

int iscube(void);

int centergroups(ToIntPtr A,int valence);
void addshift(int key,int nshift,vector shift);

#if defined(SLAB) && SLAB & 2
void cleavescaling(double *Pscale,double tau,int noint,double factor,double maxscale);
#endif
