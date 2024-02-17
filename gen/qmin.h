#include "prec.h"

/*** Minimization of a function of one variable ***

???, parm=1 does not work?

Quadratic algorithm if close to minimum, otherwise `ad hoc' search 

USAGE:
^^^^^^
  QMIN_BEGIN(x,eps,minpar) // x must be REAL variable
                           // eps>0: precision=eps
                           //        inital triplet=(x-D*eps,x,x+D*eps)
                           // eps=0: precision=sqrt(machine precision)*16
                           //        inital triplet=(x-D*precision,x,x+D*precision)
                           // eps<0: precision=sqrt(machine precision)*16
                           //        inital triplet=(x-D*eps,x,x+D*eps)
                           // minpar=1: faster search, D=11
                           // minpar=3: safer but slower, D=3
                           // minpar>3: even safer and slower
    if (QMIN_it>MAXITER*(1+minpar)) Error("too many iterations");
                                  // optional convergence test
    QMIN_f=FUNCTION_TO_MINIMIZE(x);
  QMIN_END(x)                     // x must be the same as in QMIN_BEGIN

NOTES:
^^^^^^
  Can be nested.
  Recommended accuracy = several*(machine_precision)^0.5 ~ 1e-7 for REAL=double

WARNING:
^^^^^^^^
  No check is made for improper use like not matching x, impossible
  values of arguments etc.

EXAMPLE: 
^^^^^^^^
  find minimum of f(x)=x^-12-x^-6
  REAL x=10,xx;

  QMIN_BEGIN(x,1e-7,5)
    xx=1./(x*x*x); xx*=xx;
    QMIN_f=xx*xx-xx;
  QMIN_END(x)

***/

struct QMIN_s {
  REAL x,f;
  int valid; };

extern int qmin_prt;

REAL qmin_end(REAL *xxx, struct QMIN_s QMIN[3], int phase);
REAL qmin_seteps(REAL eps);

#define QMIN_BEGIN(QMIN_x,QMIN_deps,QMIN_minpar) { \
  int QMIN_i,QMIN_it=0,QMIN_par=QMIN_minpar; \
  struct QMIN_s QMIN[3]; \
  REAL QMIN_f,QMIN_eps=qmin_seteps(QMIN_deps),QMIN_d=QMIN_deps; \
 \
  if (QMIN_par<1) QMIN_par=3; \
  if (QMIN_deps>=0) QMIN_d=QMIN_eps*(2+9/Sqr(QMIN_par)); \
  for (QMIN_i=0; QMIN_i<3; QMIN_i++) { \
    QMIN[QMIN_i].valid=0; \
    QMIN[QMIN_i].x=(QMIN_x)+(QMIN_i-1)*QMIN_d; } \
  do { \
    for (QMIN_i=0; QMIN_i<3; QMIN_i++) \
      if (!QMIN[QMIN_i].valid) { \
        (QMIN_x)=QMIN[QMIN_i].x;

#define QMIN_END(QMIN_x) \
        QMIN[QMIN_i].valid++; \
        QMIN[QMIN_i].f=QMIN_f; } \
  } while (qmin_end(&(QMIN_x),QMIN,QMIN_it++%QMIN_par)>QMIN_eps); }

