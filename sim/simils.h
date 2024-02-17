extern char *lastFn;
char *Fn(char *suffix);
const char *getext(const char *fn);
char *stripext(const char *fn);
void backup(char *ext);
void waitfordiskspace(int blocks);

extern int MC;
extern int mirror;
void rndorientation(vector *c,int ns);

extern struct simils_s {
  int pc[DIM];     /* replicate factors, initialized to 1,1,1 */
  int changeL[3];  /* = load.L */
  char *sysname;   /* system file name: for .ble and .def */
  char *simname;   /* simulation name, basename for most file names */
  char *plbname;   /* name for the playback mode (option -m1) */
  char *cfgname;   /* name for the playback mode (option -m0) */
  ToIntPtr cfg[9]; /* if changecfg then allocated, otherwise = cfg[9] */
  int N;           /* loaded: no of molecules */
  int Neq;         /* loaded: no of equations (incl. logs, lambda, r) */
  int Ns;          /* loaded: no of sites */
  int size;        /* size of a[] loaded */
  int changecfg;   /* 1 if replicate/omit/add molecules */
  int to;          /* max index of a[*] (not incl.) */
  int frommol;     /* 1st molecule to insert if incomplete configuration is loaded */
  int specsize;    /* size of allocated spec[] */
  struct spec_s {
    int N;           /* number of molecules */
    int ns;          /* number of sites */ 
  } *spec;         /* [specsize] */
} simils;

void initcfg(double pins,double Emax,int nplb,int slab_geom);
void initcryst(int npins,double Emax);

extern struct load_s {
  int n[DIM]; /* replicate factors, initialized to 1,1,1 */
  int L[DIM]; /* 0=L loaded, 1=replace by data if larger, 2=if smaller, 3=always */
  int N;      /* as above */
  int tr;     /* transformation on load */
  int zero;   /* coordinates and velocities :=0 on load */
} load;

int initreplicate(void);
void replicatecfg(void);
void loadcfg(int n,double *RvdW);
void savecfg(int n,int4 savekey,double *RvdW);
void replaceL(vector L);

void openplayback(int nm,int head);
void writeplayback(void);
void closeplayback(void);
void readplayback(int frame,int cfg);

void writeasc(void);
void readasc(int order,int cfg);

void Maxwell(int from,int to,double prob);
void MaxwellCM(double prob);

void makefixa(void);
void loadfixa(void);
void savefixa(void);
void dumpall(ToIntPtr B,ToIntPtr A);


#if defined(POLAR) && POLAR&32
void initfqcharges(void);
#endif

void zerocfg(int zero);
double Lfromfile(int k,double t);
void remove1mol(int n);
