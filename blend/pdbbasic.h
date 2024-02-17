/* basic structures and #defines for pdb */
#define VERSION "1.5b"

// DOS removed in 1.5a

#define SLASH "/"
#define ATOMBUFLEN 128 /* pdb ATOMs allocate buffer */
#define STRLEN 6     /* for atom id */
#define LINELEN 1024   /* max length of one line in PDB, mol, and par files */
#define FNLEN 1024 /* max length of fully qualified file name */

#define UNDEFATOM 1e6     /* marks undefined (WANTED) atoms */
/* #define UNDEFATOM 9e5  used by blend: a bit lower */

typedef float vector[3];
extern vector undefvector;

/* see sim/vector.h */
#define VO(A,O) { A[0] O; A[1] O; A[2] O; }
#define VV(A,B) { A[0] B[0]; A[1] B[1]; A[2] B[2]; }
#define VVV(A,B,C) { A[0] B[0] C[0]; A[1] B[1] C[1]; A[2] B[2] C[2]; }
#define SQR(A) (A[0]*A[0]+A[1]*A[1]+A[2]*A[2])

enum moltype_e {
  UNKNOWN,   /*0*/
  AMINOACID, /*1 aminoacid residue */
  WATER,     /*2*/
  MOLECULE,  /*3 any other molecule */
  PATCH,     /*4 deprecated */
  NTER,      /*5 N-terminus patch, specified by -h (or the default nter) */
  CTER       /*6 C-terminus patch, specified by -e (or the default cter) */
};
/* WARNING: the end of this list must be PATCH,NTER,CTER (?- probably not) */

/* list of sites as read from *.rsd files */

struct id_s {
  struct id_s *next;
  char id[STRLEN];
};  
/* WARNING: don't use typedef ... id_t - reserved on Digital (sys/types.h) */

typedef struct site_s {
  struct site_s *next;/* next in the list */
  int n;              /* site number */
  float *r;           /* indirection to the coordinates of the site (vector) */
  struct id_s *ids;   /* symbolic atom name - list of several names (aliases) */
  char patch[STRLEN]; /* atom name to be replaced if *.rsd is patch */
  char type[5];       /* force field atom type */
  signed char chir;   /* chirality (with respect to numbering of atoms) */
  signed char unique; /* set if only 1 such element in residue */
  float charge;       /* (partial) charge in e */
  int nnbr;           /* # of neighbors connected by a bond */
  int *nbr;           /* [option('l')] list of neighbors connected by a bond */
  int line;           /* pdb-file line # */
} site_t;

/* information from ATOM and HETATM command from the PDB file */

typedef struct atom_s {
  struct atom_s *next; /* next in the list */
  int    atomno;       /* PDB atom number in sequence */
  char   id[5];        /* atom type in PDB convention */
  char   altloc;       /* alternative location */
  vector r;            /* cartesian coordinates */
  float  occup;        /* occupancy */
  int    line;         /* pdb-file line # */
} atom_t;

/* one residue from the PDB file */

typedef struct residue_s {
  struct residue_s *next;/* next in the list */
  atom_t *atom;          /* -> list of pdb ATOMs (read from *.pdb) */
  site_t *site;          /* -> list of mol sites (read from *.rsd) */
  int    resno;          /* PDB residue number in sequence */
  char   resnm[4];       /* residue type in 3-letter code */
  char   chain;          /* chain identity letter */
  char   resins;         /* residue insertions (only partly implemented!) */
  signed char ter;       /* terminate chain flag (irrelevant for !AMINOACID) */
                         /* 0: no term, -1: nter, 1: cter */
  enum moltype_e moltype;/* AMINOACID, WATER, ... */ 
} residue_t;

extern residue_t *reshead;

/* SSBOND command from the PDB file */

typedef struct ssbond_s {
  struct ssbond_s *next;
  char chain[2];      /* chain letters (array[2] of char, not string) */
  int resno[2];       /* residue numbers */
  site_t *site[2];    /* sites (=SG atoms) to connect */
} ssbond_t;

extern ssbond_t *sshead;

void includeSS(ssbond_t *ss); /* includes ss into list if not already there */

/* CONECT command from the PDB file */

typedef struct connect_s {
  struct connect_s *next;
  int atomno[2];       /* atom numbers to connect (PDB numbering) */
  site_t *site[2];     /* sites to connect */
} connect_t;

typedef struct reconnect_s { /* to reconnect patch: site list */
  struct reconnect_s *next;
  struct site_s *site[2]; 
} reconnect_t;

extern connect_t *connecthead;
extern reconnect_t *reconnecthead;

/* info lines from the pdb file */
extern char pdbheader[82],pdbcompnd[82],optioninfo[82];

extern char *molname; /* mol-file name without extensions */
extern char *pdbname; /* PDB file name without extensions */

FILE *openfile(char *ext,char *mode,int chk); /* open molname.ext with mode */

char *removedigits(char *id); /* ret. id without leading and trailing digits */

/* omitted sites (are in *.pdb but not supported by *.rsd), see also -i */

typedef struct omit_s {
  struct omit_s *next;
  int atomno;          /* pdb atom number */
} omit_t;

extern omit_t *omithead;

void prtres(residue_t *res); /* prt residue info */

site_t *findsite(int n); /* finds site with given number n */

void checkneighbors(void); /* checks consistency of table of neighbors */

/* DEBUGGING TOOLS - may be turned off */
void prtsitelist(char *key,site_t *s0);

extern char *blendpath;
extern FILE *file;
extern char line[LINELEN];
extern char fn[FNLEN];
extern const char *rsdext;
extern char *rsddir;
extern char *parameter_set;

char *getrsddir(void); /* get rsddir form the parameter file */

void freeid(site_t *s); /* free the s->ids list */
void appendid(site_t *s, char *id); /* append id to site->ids */
