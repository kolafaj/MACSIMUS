#ifndef REAL
#define REAL double
#endif

/* Method of Secants 

USAGE:

  MS_BEGIN(x,accuracy)            // x must be real-type variable
                                  // initial x = x-accuracy,x+accuracy
                                  // accuracy must be positive
    if (MS_it>MAXITER) Error("too many iterations"); 
                                  // optional convergence test
    MS_f=EQUATION(x);
    MS_END(x,extra_iter)          // x must be the same as in MS_BEGIN
                                  // extra_iter must be 0,1,...

NOTES:

  Can be nested to solve set of equations.
  Recommended accuracy = (machine_precision)^0.618 = 1e-10 for REAL=double
  Recommended extra_iter = 0 if only resulting x is of interest
  Recommended extra_iter = 1 if any value depending on x and calculated
    inside MS is used.  This happens e.g. if MS are nested: all but
    innermost MS should have extra_iter = 1 (if all roots are of interest)
  Recommended extra_iter > 1 if EQUATION(x) is noisy.  In this case,
    the value of eps should be comparable with the noise

WARNING:

  No check is made for improper use like not matching x, impossible
  values of arguments etc.

EXAMPLE 1: solve x**2==2

  REAL sqrt2=2;

  MS_BEGIN(sqrt2,1e-10)
    MS_f=sqrt2*sqrt2-2;
    MS_END(sqrt2,1)

EXAMPLE 2: solve (x*y==1,x**2==y**2)

  REAL x=2,y=2,eps=1e-10;

  MS_BEGIN(x,eps)

    MS_BEGIN(y,eps)
      MS_f=x*y-1;
      MS_END(y,0)
    MS_f=x*x-y*y;

    MS_END(x,1)

HISTORY:
  1/2005: MS_lastd added to solve very improbable case of x^3-x=0 with x0=0.5

*/

#define MS_BEGIN(MS_x,MS_deps) {                 \
  REAL MS_eps=(MS_deps), MS_d=MS_eps*2, MS_f, MS_g=0, MS_der, MS_lastd=1000*MS_eps;   \
  int MS_it=0, MS_n=-2;                          \
  (MS_x)+=MS_eps;                                \
  do {

#define MS_END(MS_x,MS_OK)                       \
  if (MS_g==MS_f) {                              \
    MS_d*=-1.618; (MS_x)+=MS_d; }                \
  else {                                         \
    if (MS_it<2) MS_der=MS_d/(MS_g-MS_f);        \
    else MS_der=(MS_d/(MS_g-MS_f)*(Sqr(MS_d)+Sqr(MS_lastd))+MS_der*Sqr(MS_eps)) \
               /(Sqr(MS_d)+Sqr(MS_lastd)+Sqr(MS_eps));         \
    MS_lastd=MS_d;                               \
    (MS_x)+=MS_d=MS_f*MS_der; }                  \
  MS_g=MS_f;                                     \
  if (fabs(MS_d)>MS_eps) MS_n=-1; else MS_n++;   \
  MS_it++; } while (MS_n<(MS_OK)); }

