#define N 1000 /* max # of columns; see also tabproc.c */
int drawerrbar(char *fn,
	       char *colx,char *coly,
	       char *coldy1,char *coldy2,
	       int maxcol,int color,int style,int errstyle);
int drawfile(char *fn,
             char *colx,char *coly,
             int maxcol,int color,int style);
int minmaxfile(char *fn,
               char *colx,char *coly,
               int maxcol,
               double *X0,double *X1,double *Y0,double *Y1);

extern double *count; /* needed in tabinc.c */

extern struct plotopt_s { /* passed from plot */
  int slowdraw; /* delay in us */
  char *key;    /* plot only lines containing string key */
} plotopt;
  
extern struct _Idlist_s *id;
