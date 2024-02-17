/***
    Solves a system of  n  linear equations   AX=B  by the Gauss-Jordan
    elimination with a partial selection of the pivotal element.
    The function GJlineq() returns -1 in the case of singularity.

    Directly #includable after defining REAL (e.g., using prec.h)
***/

int GJlineq(int n, REAL **a, REAL *b, REAL *x)
{
  REAL aa,t;
  int i,j=-1,k;

  /*** Gauss-Jordan elimination ***/
  loop (k,0,n-1) {
    aa=0.0;
    loop (i,k,n) {
      t=a[i][k];
      if (fabs(t)>fabs(aa)) {
        aa=t; j=i; } }
    
    if (aa==0.0 || j<0) return -1;
    
    if (j!=k) {
      loop (i,k,n) {
        t=a[j][i]; a[j][i]=a[k][i]; a[k][i]=t; }
      t=b[j]; b[j]=b[k]; b[k]=t; }
    
    loop (i,k+1,n) {
      t=a[i][k]/aa;
      loop (j,k+1,n) a[i][j]-=t*a[k][j];
      b[i]-=t*b[k]; } }

  /*** back substitution ***/
  for (i=n-1; i>=0; i--) {
    if (a[i][i]==0.0) return -1;
    x[i]=b[i]/a[i][i];
    loop (j,0,i) b[j]-=x[i]*a[j][i]; }

  return 0;
}
