extern int ndependants;

extern struct dependant_s {
  struct dependant_s *next;
  int nnbr;         /* # of parents */
  int nbr[MAXVAL];  /* parents */
  int ndep;         /* # of dependants */
  int dep[1/*ndep*/];/* list of dependants */
  } *dep0;

extern int removeb,removea,removei;
extern bond_t *b00;

extern int anymass;
int readdependants(species_t *spec);
void removebond(int i,int j);
void removeangle(int i,int j,int k);
void removetorsion(torsion_t **t0,int i,int j,int k,int l);
int countdependants(void);
void prtdependants(void);

