/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 The set of linear equations
   Mx=b
 is solved, where  M  is a symmetric positive-definite matrix.

 Efficient if M is sparse (and  Mvec()  makes use of it) and also if
 x  is very close to the solution.

 METHOD:
 Quadratic form
   f(x) = xMx/2 - bx
 is minimized by the conjugate gradient method.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "ground.h"
#include "conjgrad.h"

int conjgrad(
  int     nc,   /* order of vectors and matrix  M  */
  double *x,    /* x[nc]: initial estimate on input - solution returned */
  Mvec_f *Mvec, /* pointer to user-supplied procedure
                   void Mvec(int nc, double *Mv, void *Mdef, double *v)
                   returning product  M.v  in  Mv[nc] */
  void   *Mdef, /* pointer to sparse-coded matrix  M  passed as parameter
                   to  Mvec  */
  double *b,    /* b[nc]: rhs vector */
  double  eps   /* accuracy:
                   if  0<eps<=1  then relative accuracy  |delta x|/|x| < eps
                   if  eps>=2  then # of iterations */
  )             /* RETURN VALUE: number of *Mvec calls */


{
  double *h,*Mh,*g;
  double gg,gMh,lambda,gamma,t,xMh,bh,MhMh,hh,xx,err=3e33;
  int i,it=1;
  int ifit=eps>=2;
  static int maxitq=4;
  int maxit = ifit ? (int)(eps*maxitq/4) : nc*maxitq;

  allocarray(h,nc*3);
  Mh=(g=h+nc)+nc;

  /* initial direction h */
  (*Mvec)(nc,Mh,Mdef,x);
  loop (i,0,nc) h[i]=g[i]=b[i]-Mh[i];

  /* iterations */
  do { it++;

    (*Mvec)(nc,Mh,Mdef,h);

    gMh=xMh=gg=bh=MhMh=0;
    loop (i,0,nc) {
      xMh+=x[i]*Mh[i]; MhMh+=Mh[i]*Mh[i]; gMh+=g[i]*Mh[i]; 
      gg+=g[i]*g[i]; bh+=b[i]*h[i]; }
    if (gg==0) goto freeret;

    t=(bh-xMh)/gMh; /* f(x+t*h) is min for this t */
    lambda=gg/gMh;
    gamma=lambda*MhMh/gMh-1;

    if (ifit) {
      loop (i,0,nc) {
        x[i]+=t*h[i];
        g[i]-=lambda*Mh[i]; h[i]=g[i]+gamma*h[i]; } }
    else {
      xx=hh=0; /* to observe accuracy */
      loop (i,0,nc) {
        x[i]+=t*h[i];
        xx+=x[i]*x[i]; hh+=h[i]*h[i];
        g[i]-=lambda*Mh[i]; h[i]=g[i]/*new!*/+gamma*h[i]; }
      err=t*t*hh/xx; }
    
  } while (err>eps*eps && it<maxit);

  if (!ifit && it==maxit) {
    if (maxitq<125)
      WARNING(("conjgrad: too many iterations (%d), limit increased",it))
    else
      ERROR(("conjgrad: really too many iterations (%d)",it)) 
    maxitq=maxitq*5/4; }

 freeret:
  free(h);
  return it;
} /* conjgrad */
