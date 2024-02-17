#define BUG // probably gcc bug - void code added
/* blend while diagonalizing the dipole moment crashes at command 
   "loop (j,0,i) {"
   happens only with gcc -O2
   adding an empty command helps
*/   

#if 0 /* test */
/* cc -O2 -o jacobi jacobi.c -lm ; jacobi
*/

#undef BUG

#include "include.h"
int sig;
#define prt printf

double Jacobi(int n,double **A,double **R,double eps);

int main()
{
  int i,j,n=3;
  double **A,**R;

  alloc(A,sizeof(A[0])*n);
  alloc(R,sizeof(A[0])*n);

  loop (i,0,n) {
    alloc(A[i],sizeof(A[0][0])*n);
    alloc(R[i],sizeof(R[0][0])*n); }

  A[0][0]=2;
  A[0][1]=A[1][0]=A[1][1]=1;
  A[2][0]=A[0][2]=.5;
  A[1][2]=A[2][1]=.25;
  A[2][2]=2;
//  loop (i,0,n)loop (j,0,n) A[i][j]=0;
  loop (i,0,n) put3(A[i][0],A[i][1],A[i][2]);

  Jacobi(n,A,R,-1e-9);

  loop (i,0,n) printf("%.10f ",A[i][i]);
  _n
  loop (i,0,n) put3(R[i][0],R[i][1],R[i][2]);
  return 0;
}
#endif

/*
  threshold Jacobi diagonalization of symmetric matrix A
  [A. Ralston: The first course in numerical analysis, Chap. 10.4]
  input: 
    n = rank
    A[n][n] (as **A) = symmetric matrix
    eps = error; eps>0 : silent, eps<0 : verbose
  output:
    the error is returned
    A[n][n] (diagonalized)
    if (R!=NULL) then eigenvectors R[n][n] are returned
  iterations are interruptable by ^C
  A (and R if used) must be allocated in advance
  NEW 5/2004 negative n suppresses interruption because of sig
  NEW 5/2012 automatic threshold
  note that the convergence becomes quadratic when almost iterated
*/

double Jacobi(int n,double **A,double **R,double eps)
{
  int i,j,p,q,nit=0;
  double m,x,y,lambda,mu,nu,cost,sint,Q,Qlim,QQ;
  double *tp,*tq,*rp,*rq;
  int ifprt=eps<0;
  int activesig=n>0;

  n=abs(n);
  
  eps=fabs(eps);
  allocarray(tp,n);
  allocarray(tq,n);
  if (R) {
    allocarray(rp,n);
    allocarray(rq,n); }

  /* max out-diag and on-diag elements; init R */
  m=x=-1;
  loop (i,0,n) {
    if (fabs(A[i][i])>x) x=fabs(A[i][i]);
    if (R) R[i][i]=1;
    loop (j,0,i) {
#ifdef BUG
      if (m==-999) put(m) // never happens
#endif      
      if (R) R[i][j]=R[j][i]=0;
      if (fabs(A[i][j])>m) m=fabs(A[i][j]); } }
  if (m+x==0) return -1;

  Q=6./n;
  Qlim=1./n;
  QQ=Q;

  for (;;) {
    m*=exp(-Q); /* changing the threshold */
    nit=0;
    loop (q,0,n) loop (p,0,q) if (fabs(A[p][q])>m) {
      nit++;

      /* p<q */

      lambda=-A[p][q];
      mu=(A[p][p]-A[q][q])/2;
      nu=sqrt(Sqr(lambda)+Sqr(mu));
      cost=sqrt((nu+fabs(mu))/(2*nu));
      sint=(mu<0?-1:1)*lambda/((2*nu)*cost);

      copyarray(tp,A[p],n);
      copyarray(tq,A[q],n);
      if (R) {
	copyarray(rp,R[p],n);
	copyarray(rq,R[q],n); }

      loop (i,0,p) {
	A[p][i]=A[i][p]=tp[i]*cost-tq[i]*sint;
	A[q][i]=A[i][q]=tp[i]*sint+tq[i]*cost; }

      loop (i,p+1,q) {
	A[p][i]=A[i][p]=tp[i]*cost-tq[i]*sint;
	A[q][i]=A[i][q]=tp[i]*sint+tq[i]*cost; }

      loop (i,q+1,n) {
	A[p][i]=A[i][p]=tp[i]*cost-tq[i]*sint;
	A[q][i]=A[i][q]=tp[i]*sint+tq[i]*cost; }

      A[p][p]=tp[p]*Sqr(cost)+tq[q]*Sqr(sint)-tp[q]*(cost*sint*2);
      A[q][q]=tp[p]*Sqr(sint)+tq[q]*Sqr(cost)+tp[q]*(cost*sint*2);
      A[p][q]=A[q][p]=(tp[p]-tq[q])*(cost*sint)+tp[q]*(Sqr(cost)-Sqr(sint));

      if (R) loop (i,0,n) {
	R[p][i]=rp[i]*cost-rq[i]*sint;
	R[q][i]=rp[i]*sint+rq[i]*cost; } }

    x=-1;
    loop (i,0,n) if (fabs(A[i][i])>x) x=fabs(A[i][i]);

    y=sqrt(n/(nit+0.5));
    Min(y,1.25) Max(y,0.8)
    Q*=y;
    QQ=QQ*0.875+Q*0.125; /* not needed here for a symmetric matrix */
    //    prt("%d %g %g %g %g  nit y=Q_quotient m/maxL=err Q=quotient smoothed",nit,y,m/x,Q,QQ);

    if (x<0 || m/x<eps || (Q<Qlim && m/x<sqrt(eps)) || (sig && activesig)) break;
    if (ifprt) fprintf(stderr,"Jacobi iter=%d, err=%g, Q=%g\n",nit,m/x,Q); }

  free(tq); free(tp);
  if (R) { free(rp); free(rq); }

  if (sig && activesig) {
    prt("! Jacobi interrupted, err=%g",m/x);
    sig=0; }

  return m/x;
}
