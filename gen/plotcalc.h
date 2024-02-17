#define N 27 /* max # of columns; see also tabproc.c */
int drawfile(char *fn,char *colx,char *coly,int maxcol,int color,int style);
int minmaxfile(char *fn,int cont,char *colx,char *coly,int maxcol,
               double *X0,double *X1,double *Y0,double *Y1);
extern double *count;
extern struct _Idlist_s *id;
