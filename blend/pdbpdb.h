/* one ATOM or HETATM line just read from the PDB file */

struct A_s {     /* field length */
  int    atomno;     /*   5  atom number in sequence */
                     /*   1  space */
  char   id[5];      /*   4  atom type in BDB convention */
  char   altloc;     /*   1  alternative locations */
  char   resnm[4];   /*   3  residue type in 3-letter code */
                     /*   1  space */
  char   chain;      /*   1  chain identity letter */
  int    resno;      /*   4  residue number in sequence */
  char   resins;     /*   1  residue insertions */
                     /*   3  spaces */
  vector r;          /* 3*8  cartesian coordinates */
  float  occup;      /*   6  occupancy */
  float  bvl;        /*   6  fourth ( B- ) value */
                     /*   1  space */
  int    footnote;   /*   3  footnote number */
                     /*   2  spaces */
  char   ydent[5];   /*   4  BDB identity number (alphanumeric) */
  int    recno;      /*   4  original record number in BDB list */
  int    hetatm;     /*   -  1 if read by HETATM statement */
  int    line;       /*   -  line # of pdb file */
  };

extern struct A_s A;

extern FILE *pdb;

int includeA(residue_t *r); /* incl. A in atom list of residue r */
residue_t *appendres(residue_t *res); /* new resid. appended, 1st atom incl. */
residue_t *newpatchres(residue_t *res,char *ter); /* new residue, to hold a patch */
void readPDB(void); /* reads the PDB file into program structures */
void pastePDB(char *pdbname,int iframe); /* pastes coord to a pdb-file */
