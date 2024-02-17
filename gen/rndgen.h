/*
  Common header file for random number generators rndgen.c, rndgen0.c, rndgen1.c, ...

  Recommended generator = rndgen.c (it #includes rndgeni.c)

  NB: rnd() and rndcos() may give, albeit with tiny probablity ~1e-16, the
  interval end values. See rndgen.c for more info! Particular version may be
  guaranteed to give rnd() in [0,1) or (0,1), etc., though.
  
  Similarly, rndball and rnddisk may give points by a tiny value outside the
  sphere or disk.
*/

#include "int4.h"

unsigned4 rndinit(int4 tablen,int4 seed);
                             /* tablen<0: query whether already initialized
                                tablen=0: initialize with seed (default tablen)
                                tablen>0: initialize with seed and shuffle
                                table length (if applicable) */
unsigned4 Irnd(void);        /* function returning random integer upto max */
double rnd(void);            /* random number in [0,1] */
unsigned4 irnd(unsigned4 n); /* random number in {0,...,n-1} */
double rndcos(void);         /* random number in [-1,1] */

/* see rndetc.c: */
double rndball(double *r);    /* random vector in 1-sphere */
double rnddisk(double *r);   /* random vector in 1-disk */
void rndsphere(double *v);    /* random vector on 1-sphere */
double rndgauss(void);       /* normalized Gausian distribution */
