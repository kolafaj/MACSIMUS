/* 
   Gauss-Jordan elimination with partial pivoting
   matrix a (see matrix.c) is replaced by the inverse
   determinant is returned (0 indicates singular matrix or underflow)
   to be directly #included 
   rewritten from PL/1 version
   see also: ludecomp.c (more efficient, but more difficult to use)
*/
double GJinvmat(int il,int ih, double **a)
{
  double *b=newvector(il,ih);
  double *c=newvector(il,ih);
  int *d=newintvector(il,ih);
  double det=1.0,y,w;
  int i,j,k;

  loop (j,il,ih) d[j]=j;

  loop (i,il,ih) {
    k=i; y=a[i][i];
    loop (j,i+1,ih) {
      w=a[i][j];
      if (fabs(w) > fabs(y)) { k=j; y=w; } }

    if (y==0) goto ret;
    det*=y;
    y=1.0/y;

    loop (j,il,ih) {
      c[j]=a[j][k]; a[j][k]=a[j][i]; a[j][i]= -c[j]*y;
      b[j]=a[i][j]*=y; }
    a[i][i]=y;
    j=d[i]; d[i]=d[k]; d[k]=j;
    loop (k,il,ih) if (k!=i) loop (j,il,ih) if (j!=i) a[k][j]-=b[j]*c[k]; }

  loop (i,il,ih) while ((k=d[i])!=i) {
    loop (j,il,ih) {
      w=a[i][j]; a[i][j]=a[k][j]; a[k][j]=w; }
    j=d[i]; d[i]=d[k]; d[k]=j;
    det= -det; }

 ret:
  freeintvector(il,ih,d);
  freevector(il,ih,c);
  freevector(il,ih,b);

  return det;
}
