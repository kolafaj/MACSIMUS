/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% residue statistics %%%%%%%%%%%%%%%%%%%%%%%%%%
*/

typedef struct omitrsd_s {
  struct omitrsd_s *next;
  int n; /* # of residues omitted */
  int l; /* # of lines omitted */
  char resnm[4];
  } omitrsd_t;

extern omitrsd_t *omitrsd0;

void count(char *resnm); /* counts residue names */
void prtcount(void);     /* prints statistics */
void prtomitted(void);   /* omitted (-oRSD) statistics */
void help(void);         /* prints help and exits */
