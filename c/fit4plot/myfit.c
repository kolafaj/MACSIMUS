/* \make myfit.double
  see also fit.c
  
  getdata:
method: -1: steepest descent
         0: amoeba (default)
         1: conjugate gradient
         2: Newton-Raphson
         3: Monte Carlo 
         9: exit getdata loop
D: step to calculate numerical derivatives (method=-1,1,2 only, default=1e-5)
eps: accuracy of the minimum (irrelevant for method=3, default=1e-8)
par: `nonsphericity' valley param. par>=1 (method=3 only, default=1)
     par<=0: repeats w. extrapolation (narrow valley, all methods)
maxit: (max) number of steps (default=1000)
nerr: # of samples to calculate errors (default=0=do not calculate err.)
err: err=-1: uses stdev of data provided (assuming 1 if not provided)
     err=-2: as err=-1 if stdev of data available, err=0 otherwise (default)
     err=-3: stdev of dy_i is proportional to 3rd column,
             rescaled by stdev calculated from residual sum of squares
             (for data with known weights~1/dy_i^2 but unknown precision)
     err>=0: uses stdev calc. from resid. sum sq., multiplied by y^err:
     err=0: the same abs.err. (of stdev above) of all data
     err=1: the same rel.err. of all data
     err=0.5: error is proportional to y^0.5 (suited for histogram data)
*/

#include "ground.h"
#include "minimize.h"
#include "rndgen.h"
#include "qmin.h"

// number of parameters:
#define N 1

double A[N],A0[N];

/* Chebyshev polynomials T_i(x) in [-1,1]; x must be double */
#define T0 1
#define T1 x
#define T2 (x*x*2-1)
#define T3 x*((x*x)*4-3)
#define T4 (8*(x*x)*((x*x)-1)+1)
#define T5 x*(5-(x*x)*(20-16*(x*x)))

# if 0
/* normal polynomials */
#define T2 (x*x)
#define T3 x*(x*x)
#define T4 (x*x)*x*x)
#define T5 x*((x*x)*x*x))
#endif

// FUNCTION TO FIT
double func(double *P,double x)
{
  return exp(-x*P[0]);
}

#undef T0
#undef T1
#undef T2
#undef T3
#undef T4
#undef T5

#include "fitcode1.c"

    prt("Enter data (?? for list of data):");
    getdata
      getvec(A,,N)
      get(method)
      get(maxit) get(eps) get(par) get(D) 
      get(nerr) get(err)
    checkdata enddata

#include "fitcode2.c"

