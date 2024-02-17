/* to be #included AFTER simopt.h */
/* serial version: time-measurements functions */

double mytime(void); /* returns real time in s, resolution: see CHEAPTIME */
char *myctime(double t); /* =ctime(&(int)t) */

void printtime(void);

double inittime(int maxitems); /* reset tables; returns real time in s */
void CPUtime(char *n);

#if PARALLEL
extern struct partimes_s {
  int nrspace,nkspace;
  double *rspace; /* [nrspace] */
  double *kspace; /* [nkspace] */
} partimes;
void printpartimes(void);
/* note: for PARALLEL==1 partimes->t[0]=k-space, partimes->t[1]=r-space */
#endif
