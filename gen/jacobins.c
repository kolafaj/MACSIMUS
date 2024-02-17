/* #includable module, see also jacobi.c
   2019 test removed, GJlineq removed from here (see old+misc/test-jacobins.c)
   2018 unsuccessful attempt for superrelaxation: see below
   2013 Version 1

  The Jacobi method [A. Ralston: The first course in numerical analysis,
  Chap. 10.4] is generalized to a nonsymmetrix matrix with real
  eigenvalues.  The Schur form is produces.  Let's call this method
  Jacobi-Schur.

  INPUT:
    n = rank
    A[n][n] (as **A) = matrix
    stops iterations if error<eps OR if error<sqrt(eps) AND "smart convergence criterion"
    eps = error; eps>0 : silent, eps<0 : verbose
    key=-1: in addition, stops if the first eigenvalue becomes negative
    key>0: stops  key  consecutive iterations do not converge
  OUTPUT:
    the error is returned
    A[n][n] in the Schur form (upper diagnal matrix) is produced
    if (R!=NULL) then JacobiEigenvectors must be followed, then
     eigenvectors R[n][n] are returned:
     R[0][] is the eigenvector of eigenvalue A[0][0]
     R[1][] is the eigenvector of eigenvalue A[1][1]
     ...
  The iterations are interruptable by ^C
  A (and R if used) must be allocated in advance
  The iterations are linear (unlike for original Jacobi method for
  symmetric matrices)
*/

static int JacobiVerbose;

double Jacobins(int n,double **A,double **R,double eps,int key)
/* transforms A into the Schur form by the Jacobi method */
{
  int i,j,p,q,nit=0,nstall=0;
  double m,x,y,cost,sint,Q,Qlim,QQ,lasteps=999;
  double *tip,*tiq,*tpi,*tqi,*rp=NULL,*rq=NULL; /* initialized to suppress compiler warning */
  int ifprt=eps<0;
  int activesig=n>0;
  double tt,app_aqq,mindiag;
  char line[256];

  n=abs(n);

  JacobiVerbose=eps<0;

  eps=fabs(eps);
  allocarray(tip,n);
  allocarray(tpi,n);
  allocarray(tiq,n);
  allocarray(tqi,n);
  if (R) {
    allocarray(rp,n);
    allocarray(rq,n); }

  /* max out-diag and on-diag elements; init R */
  m=x=-1;
  loop (i,0,n) {
    if (fabs(A[i][i])>x) x=fabs(A[i][i]);
    if (R) R[i][i]=1;
    loop (j,0,i) {
      if (R) R[i][j]=R[j][i]=0;
      if (fabs(A[i][j])>m) m=fabs(A[i][j]); } }
  if (m+x==0) return -1;

  Q=6./n;
  Qlim=1./n;

  QQ=Q;

  if (ifprt) {
    fprintf(stderr,"#nit y=Q_quot err=m/maxL  Q*n smoothed  min_Aii\n");
    prt_("#nit y=Qquot err=m/maxL   Q*n   =smoothed  min_Aii\n"); }

  for (;;) {
    m*=exp(-Q); /* changing the threshold */
    nit=0;
    loop (q,0,n) loop (p,0,q) if (fabs(A[q][p])>m) {
      nit++;

      /* here p<q ; A[q][p]->0 */
      app_aqq=A[p][p]-A[q][q];
      tt=Sqr(app_aqq)+4*A[p][q]*A[q][p];
      if (tt<0) continue;
      tt=atan2(app_aqq+sqrt(Sqr(app_aqq)+4*A[p][q]*A[q][p]),2*A[p][q]);

      /* TEST: (column D should go to 0)
         prt("%14.6g %14.6g  %14.6g %14.6g  %g TT",A[p][p],A[p][q],A[q][q],A[q][p],tt);
         fgrep TT  TEST.prt | plot -:0:"D*exp(n*a)"
      */
      /*
The following attempt of superrelaxation does not work:
#define SUPER 1.05
#define TH 0.1
      if (fabs(tt)<TH) tt*=SUPER;
      if (tt>PI-TH) { tt-=PI; tt*=SUPER; tt+=PI; }
      if (tt<-PI+TH) { tt+=PI; tt*=SUPER; tt-=PI; }
SUPER Niter
0.5 75400
0.8 39334
0.9 31929
0.95 28745
0.99 27023
1   26964
1.01 27228
1.05 28922
1.1 31993
1.2 38950
1.5 74491
       */

      sint=sin(tt);
      cost=cos(tt);

      loop (i,0,n) {
        tpi[i]=A[p][i];
        tqi[i]=A[q][i];
        tip[i]=A[i][p];
        tiq[i]=A[i][q]; }

      if (R) {
	copyarray(rp,R[p],n);
	copyarray(rq,R[q],n); }

      loop (i,0,p) {
	A[p][i]=cost*tpi[i]-sint*tqi[i];
	A[i][p]=cost*tip[i]-sint*tiq[i];
	A[q][i]=cost*tqi[i]+sint*tpi[i];
	A[i][q]=cost*tiq[i]+sint*tip[i]; }


      loop (i,p+1,q) {
	A[p][i]=cost*tpi[i]-sint*tqi[i];
	A[i][p]=cost*tip[i]-sint*tiq[i];
	A[q][i]=cost*tqi[i]+sint*tpi[i];
	A[i][q]=cost*tiq[i]+sint*tip[i]; }

      loop (i,q+1,n) {
	A[p][i]=cost*tpi[i]-sint*tqi[i];
	A[i][p]=cost*tip[i]-sint*tiq[i];
	A[q][i]=cost*tqi[i]+sint*tpi[i];
	A[i][q]=cost*tiq[i]+sint*tip[i]; }

      A[p][p]=tpi[p]*Sqr(cost)+tqi[q]*Sqr(sint)-(tpi[q]+tqi[p])*(cost*sint);
      A[q][q]=tpi[p]*Sqr(sint)+tqi[q]*Sqr(cost)+(tpi[q]+tqi[p])*(cost*sint);
      A[p][q]=tpi[q]*Sqr(cost)-tqi[p]*Sqr(sint)+(tpi[p]-tqi[q])*(cost*sint);
      A[q][p]=tqi[p]*Sqr(cost)-tpi[q]*Sqr(sint)+(tpi[p]-tqi[q])*(cost*sint);

      if (R) loop (i,0,n) {
	R[p][i]=rp[i]*cost-rq[i]*sint;
	R[q][i]=rp[i]*sint+rq[i]*cost; } }

    x=-1;
    mindiag=3e33;
    loop (i,0,n) {
      if (fabs(A[i][i])>x) x=fabs(A[i][i]);
      Min(mindiag,A[i][i]) }

    y=sqrt(n/(nit+0.5));
    Min(y,1.25) Max(y,0.8)
    Q*=y;
    QQ=QQ*0.875+Q*0.125;
    if (ifprt) {
      sprintf(line,"%4d %7.5f %.3e %8.5f %8.5f %g",nit,y,m/x,Q*n,QQ*n,mindiag);
      fprintf(stderr,"%s\n",line);
      prt("%s =J",line); }

    if (x<0 ||
        m/x<eps ||
        (QQ<Qlim && m/x<sqrt(eps)) ||
        (sig && activesig) ||
        (key==-1 && mindiag<0) ||
        (key>0 && nstall>key) ) break;

    if ((m/x)/lasteps>1-eps) nstall++;
    else nstall=0;
    lasteps=m/x;
  }

  free(tqi); free(tpi);
  free(tiq); free(tip);
  if (R) { free(rp); free(rq); }

  if (sig && activesig) {
    prt("! Jacobi interrupted, err=%g",m/x);
    sig=0; }

  return m/x;
}

void JacobiEigenvectors(int n,double **A0,double **A,double **R)
/*
 input:
   A0 = original matrix
   A  = Schur form (diagonal = eigenvectors)
   R  = orthogonal vectors of the Schur form
 returned:
   R[i][] = eigenvector coresponding to eigenvalue A[i][i]
*/
{
  double **C;
  double *P;
  double *X;
  int i,j,k,a,b;

  if (!R) return;

  alloc2Darray(C,n,n);
  allocarray(P,n);
  allocarray(X,n);

  /* R[0][] unchanged */

  loop (i,1,n) {
    loop (k,0,n) {

      P[k]=0; /* inefficient */
      loop (a,0,n) loop (b,0,n)
        P[k]+=R[k][a]*(A0[a][b]-A[i][i]*(a==b))*R[i][b];

      loop (j,0,i) {
        C[k][j]=0;
        loop (a,0,n) C[k][j]+=R[k][a]*R[j][a];
        C[k][j] *= A[i][i]-A[j][j]; } }

    GJlineq(i,C,P,X);
    loop  (j,0,i)
      loop (a,0,n) R[i][a]+=X[j]*R[j][a];
    if (JacobiVerbose && i>100 && (i%20==0)) fprintf(stderr,"JacobiEigenvectors: %d of %d finished at %s",i,n,myctime(mytime()));
  }

  free(X);
  free(P);
  free2Darray(C);
}

double JacobiTest(int N,double lambda,double **A, double *R)
/*
  returns |A.R - lambda R|/|R|, i.e., rel. error of R as eigenvector of A
  (NB: this error may be large for lambda=0)
*/
{
  int i,j;
  double *RR,s1=0,s2=0;

  allocarrayzero(RR,N);
  loop (i,0,N) loop (j,0,N) RR[i]+=A[i][j]*R[j];

  loop (i,0,N) s1+=Sqr(RR[i]);
  loop (i,0,N) RR[i]-=lambda*R[i];
  loop (i,0,N) s2+=Sqr(RR[i]);

  free(RR);

  return sqrt(s2/s1);
}

/*
   Gauss-Jordan elimination with partial pivoting
   matrix  a  is replaced by the inverse
   see also invmat.c
*/
double invmat(int N, double **a) /*********************************** invmat */
{
  double *b,*c;
  int *d;
  double det=1.0,y,w;
  int i,j,k;

  allocarray(b,N);
  allocarray(c,N);
  allocarray(d,N);

  loop (j,0,N) d[j]=j;

  loop (i,0,N) {
    k=i; y=a[i][i];
    loop (j,i+1,N) {
      w=a[i][j];
      if (fabs(w) > fabs(y)) { k=j; y=w; } }

    if (y==0) goto ret;
    det*=y;
    y=1.0/y;

    loop (j,0,N) {
      c[j]=a[j][k]; a[j][k]=a[j][i]; a[j][i]= -c[j]*y;
      b[j]=a[i][j]*=y; }
    a[i][i]=y;
    j=d[i]; d[i]=d[k]; d[k]=j;
    loop (k,0,N) if (k!=i) loop (j,0,N) if (j!=i) a[k][j]-=b[j]*c[k]; }

  loop (i,0,N) while ((k=d[i])!=i) {
    loop (j,0,N) {
      w=a[i][j]; a[i][j]=a[k][j]; a[k][j]=w; }
    j=d[i]; d[i]=d[k]; d[k]=j;
    det= -det; }

 ret:
  free(d);
  free(c);
  free(b);

  return det;
}

/* matrix . vector, R := A.V
  V,R=N-vectors, A=N*N matrix
*/
void matdotvec(int N,double *R,double **A,double *V)
{
  int i,j;
  loop (i,0,N) {
    R[i]=0;
    loop (j,0,N) R[i]+=A[i][j]*V[j]; }
}
