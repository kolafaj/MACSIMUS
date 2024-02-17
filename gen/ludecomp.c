/* LU decomposition:
   - matrix inversion
   - determinant
   - set of lin.eqs. (esp. good for several right-hand sides)
   module for direct #include
   module required: matrix.c
   HIGHPREC support 6/2004
   procedures:
     matcopy matprod
     LUdecomposition LUbacksubst
     LUlineq
     LUinvmat
     LUdet
   see also: lineq.c invmat.c (Gauss-Jordan w. partial pivot)
*/

#ifndef REAL
#define REAL double
#endif

void matcopy(int il,int ih,REAL **A,REAL **B)
/*
  B:=A
*/
{
  int i,j;

  loop (i,il,ih)
    loop (j,il,ih) B[i][j]=A[i][j];
}

void matprod(int il,int ih, REAL **A,REAL **B, REAL **C)
/* 
   C:=A.B (matrix product)
*/
{
  int i,j,k;
  REAL sum;

  loop (i,il,ih)
    loop (j,il,ih) {
      sum=0;
      loop (k,il,ih) sum+=A[i][k]*B[k][j];
      C[i][j]=sum; }
}

REAL LUdecomposition(int il,int ih, int *indx, REAL **a)
/* matrix a is replaced by its LU decomposition, 
   array index keeps track of rearranged rows */
{
  int i,imax=il-1,j,k;
  REAL big,x,s,d=1;
  REAL *v=newvector(il,ih);

  loop (i,il,ih) {
    big=0;
    loop (j,il,ih) if ((x=fabs(a[i][j]))>big) big=x;
    if (big==0) Error("zero matrix");
    v[i]=1/big; }

  loop (j,il,ih) {
    loop (i,il,j) {
      s=a[i][j];
      loop (k,il,i) s-=a[i][k]*a[k][j];
      a[i][j]=s; }
    big=0;
    loop (i,j,ih) {
      s=a[i][j];
      loop (k,il,j) s-=a[i][k]*a[k][j];
      a[i][j]=s;
      if ((x=v[i]*fabs(s))>=big) { big=x; imax=i; } }

    if (imax<il) Error("garbage matrix");

    if (j!=imax) {
      loop (k,il,ih) {
	x=a[imax][k]; a[imax][k]=a[j][k]; a[j][k]=x; }
      d= -d;
      v[imax]=v[j]; }

    indx[j]=imax;

    if (a[j][j]==0) Error("singular matrix"); /* or a[j][j]=small number */
    x=1/a[j][j];
    loop (i,j+1,ih) a[i][j]*=x; }

  freevector(il,ih,v);

  return d;
}

void LUbacksubst(int il,int ih, REAL **a, int *indx, REAL *b)
/* b is replaced by the solution x of a*x=b.  
   LUdecomposition must be called first; then, LUbacksubst may be called
   several times
*/
{
  int i,ii=il-1,ip,j;
  REAL s;

  loop (i,il,ih) {
    ip=indx[i];
    s=b[ip]; b[ip]=b[i];
    if (ii>=il) loop (j,ii,i) s-=a[i][j]*b[j];
    else if (s!=0) ii=i;
    b[i]=s; }

  for (i=ih-1; i>=il; i--) {
    s=b[i];
    loop (j,i+1,ih) s-=a[i][j]*b[j];
    b[i]=s/a[i][i]; }
}

REAL LUinvmat(int il,int ih, REAL **a, REAL **inva) /******** LUinvmat */
/* inva := inversion of a; a is replaced by its LU decomposition 
   determinant is returned */
{
  int i,j,*indx=newintvector(il,ih);
  REAL *col=newvector(il,ih);
  REAL d=LUdecomposition(il,ih,indx,a);

  loop (j,il,ih) {
    d*=a[j][j];
    loop (i,il,ih) col[i]=0;
    col[j]=1;
    LUbacksubst(il,ih,a,indx,col);
    loop (i,il,ih) inva[i][j]=col[i]; }

  freevector(il,ih,col);
  freeintvector(il,ih,indx);

  return d;
}

int LUlineq(int il,int ih, REAL **a, REAL *b, REAL *x) /****** LUlineq */
{
  /* x := solution of a*x=b; a is replaced by the LU decomposition */
  /* returns -1 if singular */
  int *indx=newintvector(il,ih);
  int sing;

  sing=LUdecomposition(il,ih,indx,a)==0;

  if (!sing) {
    /* x:=b */
    copy(x+il,b+il,(ih-il)*sizeof(REAL));

    LUbacksubst(il,ih,a,indx,x); }

  freeintvector(il,ih,indx);

  return -sing;
}

REAL LUdet(int il,int ih, REAL **a)
/* det(a) is returned, a is replaced by the LU decomposition */
{
  int i,*indx=newintvector(il,ih);
  REAL d;

  d=LUdecomposition(il,ih,indx,a);
  loop (i,il,ih) d*=a[i][i];
  freeintvector(il,ih,indx);

  return d;
}

