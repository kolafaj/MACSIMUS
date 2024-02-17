/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   support for the conjugate gradient method to solve the constraint eqns
   note: directiter added 6/95
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sparse symmetric matrix n*n is coded by:
M_t M[(# of nonzero elements on or below diagonal) + 1]
where elements M[a,b] are in the order of
ab = 00 11 10 22 20 21 33 30 31 32 44 40 41 42 43 55 50 51 ... 0n
and some of them (not the diagonal ones) may be missing.
End of list is recognized by the additional term with b=n
*/

#include "conjgrad.h"

void Mvec(int n, double *Mv, void *voidM, double *v);
M_t *defineM (specinfo_t *spec, int *Mlen);

/* added: */
int directiter(int nc,double *x,M_t *M0,double *b,double omega,double eps);

#if 0
/*** debugging ***/

#ifdef PAR

void printM(int n, M_t *M) /***************************************** printM */
{
int a,b;
M_t *M0=M;
static char z[256];

debugf("*** matrix M");
loop (a,0,n) {
  z[0]=0;
  sprintf(z+strlen(z),"[%i] %i:%6.3f",a,a,M->Mab); /* a=b */
  while ( (b=(++M)->b) < a ) sprintf(z+strlen(z),"  %i:%6.3f",b,M->Mab);
  debugf(z); }
dput(M-M0);
} /* printM */

#else

void printM(int n, M_t *M) /***************************************** printM */
{
int a,b;
M_t *M0=M;

printf("*** matrix M\n");
loop (a,0,n) {
  printf("[%i] %i:%6.3f",a,a,M->Mab); /* a=b */
  while ( (b=(++M)->b) < a ) printf("  %i:%6.3f",b,M->Mab);
  printf("\n"); }
put(M-M0);
} /* printM */
#endif
#endif
