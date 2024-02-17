/*** Lagrange equations of motion ***/
void constraintdynamics(ToIntPtr B, ToIntPtr A, ToIntPtr V);

/* Verlet integration + SHAKE; prob is for Maxwell/Andersen thermostat */
void Shake(double eps, double prob);

#ifdef POLAR
int selffield(ToIntPtr B,ToIntPtr A,double eps,double omega,int run);
void mechpolar(ToIntPtr B,ToIntPtr A);
int scforces(ToIntPtr B,ToIntPtr A);
void measureepspol(void);
void testSCF(void);
int scfautoset(int icyc,int noint);
#else
double calculateepsf(int rescale);
#endif /*# POLAR */

#ifdef SHEAR
void Shear(ToIntPtr B, ToIntPtr A, ToIntPtr VH,double H);
#endif /*# SHEAR */

extern struct constrd_s {
  double dV;     /* for pressure measurements via <dU/dV>, in units of V */
  int mode;      /* see rescalecfg() in norm.c */
  //  int ncoord;    /* # of coordinates scaled */ now No.ncoord
  char *PdVname; /* "PdVmol [Pa]" or "PdVatom [Pa]" */
} constrd;

void measureP(int pass);

/* now always included */
double Jacobi(int n,double **A,double **R,double eps);
double Jacobins(int n,double **A,double **R,double eps,int key);

extern struct nm_s {
  double dr;   /* step for numerical derivative (negative OK) */
  double eps;  /* precision of Jacobi */
  double ampl; /* amplitude for visualization, = blend -M */
  double mem;  /* max memory, in GiB (only terms ~O(N^2) counted) */
  int modes;   /* number of output plbs (modes), = blend -P */
  int frames;  /* frames for output plbs to mimic the motion, = blend -F */
  int zero;    /* limit in cm-1 to include */
  int key;     /* for model with constraints:
                  key<0: interrupt on a negative diagonal element
                  key>0: interrupt if key consecutive iterations do not converge
               */
  int method;  /* 0: for models without constraints, Jacobi method
                  1: for models with constraints, generalized Jacobi
                  2: n.a.
                  3: for models with constraints, calls octave to calculate eigenvalues
               */
} nm;

void normalmodes(void);
void normalmodesc(void);
