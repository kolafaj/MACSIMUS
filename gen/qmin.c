#include <math.h>
#include "qmin.h"

#ifdef __cplusplus
extern "C" { int prt_(const char *format, ...); }
extern "C" { int prt(const char *format, ...); }
#else
int prt_(const char *format, ...);
int prt(const char *format, ...);
#define Val(X) (double)(X)
#endif

int qmin_prt=1;

REAL qmin_end(REAL *xxx, struct QMIN_s QMIN[3], int phase)
{
  REAL fx,mx,Mx,det,D[3];
  int i,swapped,repl;
  struct QMIN_s Q;

/*.....for (i=0; i<3; i++) prt_("%.16g",Val(1-QMIN[i].x)),RPut(QMIN[i].f);*/

  det=QMIN[0].f*(D[0]=QMIN[2].x-QMIN[1].x)
     +QMIN[1].f*(D[1]=QMIN[0].x-QMIN[2].x)
     +QMIN[2].f*(D[2]=QMIN[1].x-QMIN[0].x);
/*.....RPut(det);*/
  if (det>0)

/*** original ill-conditioned formula: ***
    *xxx=(QMIN[0].f*(QMIN[2].x-QMIN[1].x)*(QMIN[2].x+QMIN[1].x)
         +QMIN[1].f*(QMIN[0].x-QMIN[2].x)*(QMIN[0].x+QMIN[2].x)
         +QMIN[2].f*(QMIN[1].x-QMIN[0].x)*(QMIN[1].x+QMIN[0].x))/(det*2);
*** rewritten as difference of x1: ***/

    *xxx=( QMIN[0].f*(D[0]*D[0]) 
         + QMIN[1].f*D[1]*(D[0]-D[2]) 
         - QMIN[2].f*(D[2]*D[2]) ) / (det*2) + QMIN[1].x;

  mx=9e9;
  Mx=0;

  if (det<=0 || *xxx<2*QMIN[0].x-QMIN[2].x || *xxx>2*QMIN[2].x-QMIN[0].x)
    /* either not convex or extrapolated too much:  just go in the direction 
       of gradient by a step (if phase) or double-step (if !phase) */
    if (QMIN[0].f<QMIN[2].f)
      if (phase) {
	*xxx=2*QMIN[0].x-QMIN[1].x;
	repl=2; }
      else {
	*xxx=2*QMIN[0].x-QMIN[2].x;
	repl=1; }
    else
      if (phase) {
	*xxx=2*QMIN[2].x-QMIN[1].x;
	repl=0; }
      else {
	*xxx=2*QMIN[2].x-QMIN[0].x;
	repl=1; }

  else 
    /* new approx to the minimum: looking for the farthest old to replace */
    for (i=0; i<3; i++) {
      fx=fabs(*xxx-QMIN[i].x);
      if (fx>Mx) { Mx=fx; repl=i; }
      if (fx<mx) mx=fx; }
  
  /* replacing by new approx */
  QMIN[repl].x=*xxx; QMIN[repl].valid=0;

  /* to order by x again */
  do {
    swapped=0;
    for (i=0; i<2; i++) if (QMIN[i].x>QMIN[i+1].x) {
      swapped=1;
      Q=QMIN[i]; QMIN[i]=QMIN[i+1]; QMIN[i+1]=Q; }
  } while (swapped);

#ifdef STDERR
  if (qmin_prt) fprintf(stderr,"qmin: det=%9.2e",Val(det));
#else
  if (qmin_prt) prt_("qmin: det=%9.2e",Val(det));
#endif

#if 0 /* debug */
#ifdef STDERR
  {
    fprintf(stderr," x=");
    for (i=0; i<3; i++) fprintf(stderr,"%15.13f ",Val(QMIN[i].x));
    fprintf(stderr," min=%.4g",Val(Mx));
  }
#else
  {
    prt_(" x=");
    for (i=0; i<3; i++) prt_("%15.13f ",Val(QMIN[i].x));
    prt_(" min=%.4g",Val(Mx));
  }
#endif
#endif

#ifdef STDERR
  if (qmin_prt) fprintf(stderr," dx=%g\n",Val(mx));
#else
  if (qmin_prt) prt(" dx=%g",Val(mx));
#endif

  return mx;
}

REAL qmin_seteps(REAL eps)
{
  REAL x=1.999999,xx=1.9,obf=sqrt((REAL)7);
  REAL *p=&x;
  int i;

  if (eps>0) return eps;
  
  /* determining machine precision */
  while (x>xx) {
    eps=x-xx;
    xx=x;
    obf+=1e-6;
    i=Int(obf)-2; /* i=0 : the optimizer is terribly intelligent... */
    x=(p[i]*2+2)/3; }
 
  return sqrt(eps)*16;
}
