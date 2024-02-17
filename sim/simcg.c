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

#include "ground.h"
#include "simglob.h"
#include "simcg.h"

void Mvec(int n, double *Mv, void *voidM, double *v) /**************** Mvec */
/***
  Mv=M*v is performed where Mv and v are vectors[n] and matrix M[n,n]
  is sparse and coded by M_t
***/
{
  int a,b;
  M_t *M=(M_t*)voidM;

  loop (a,0,n) {
    Mv[a]=M->Mab*v[a]; /* a==b */
    while ( (b=(++M)->b) < a ) {
      Mv[a] += M->Mab*v[b];
      Mv[b] += M->Mab*v[a]; } }
} /* Mvec */

M_t *defineM (specinfo_t *spec, int *Mlen) /************************ defineM */
/***
    array containing definition of sparse matrix M is allocated and built,
    address of its 1st element is returned
***/
{
  int n,a,b,pass,i,j,s;
  M_t *M=NULL; /* initialized to suppress compile warning */

  if (!spec->nc) return NULL;

  *Mlen=0;
  loop (pass,0,2) {
    n=0;
    loop (a,0,spec->nc) {
      if (pass) {
        M[n].b=a;
        /* added in V3.1l: */
        M[n].imass = spec->si[spec->si[a].pair[0]].imass
                   + spec->si[spec->si[a].pair[1]].imass; }
      n++;
      loop (b,0,a) loop (i,0,2) loop (j,0,2)
        if ((s=spec->si[a].pair[i])==spec->si[b].pair[j]) {
          if (pass) {
            M[n].b=b;
            M[n].imass = spec->si[s].imass*(i==j ? 1 : -1); }
          n++; }  } /*a*/

    if (pass==0) rallocarrayzero(M,n+1); }

  M[n].b=spec->nc;
  *Mlen=n;

  return M;
} /* defineM */


/* ============================================================== directiter */
/*
  Solves the set of linear equations  M0 x = b  by the diagonal-based
  direct iteration method.  `omega' is the relaxation parameter.
*/
int directiter(int nc,double *x,M_t *M0,double *b,double omega,double eps)
{
  int it=0,a;
  double *slope,*Mx;
  double d,dd,xx,ee=eps*eps;
  M_t *M;

  allocarray(slope,nc);
  allocarray(Mx,nc);
  M=M0;
  loop (a,0,nc) {
    slope[a]=omega/M->Mab;
    while ( (++M)->b < a ); }

  do { it++;
    if (it>100*nc*nc) Error("directiter: too many iterations");
    dd=xx=0;
    Mvec(nc,Mx,M0,x);
    loop (a,0,nc) {
      x[a] += d=(b[a]-Mx[a])*slope[a];
      xx += x[a]*x[a];
      dd += d*d; }
  } while(dd/xx>ee);

  free(Mx); free(slope);

  return it;
} /* directiter */
