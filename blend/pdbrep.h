void replacepattern(void); /* replaces patterns according *.rep */

typedef struct select_s {
  struct select_s *next;
  char resnm[4];
  int resno;
  char chain;
  char resfn[1]; /* var len */
  } select_t;
  
extern select_t *select0;

void readSEL(void);
