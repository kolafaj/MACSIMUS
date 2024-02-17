/*** routines to handle vectors and matrices ***/
/* HIGHPREC support 6/2004 */

#ifndef REAL
#define REAL double
#endif

REAL *newvector(int il,int ih)
{
  REAL *vec;

  if (ih<=il) Error("newvector");
  allocarray(vec,ih-il);

  return vec-il;
}

void freevector(int il,int ih,REAL *vec)
{
  if (ih<=il) Error("freector");
  vec+=il;
  free(vec);
}

int *newintvector(int il,int ih)
{
  int *vec;

  if (ih<=il) Error("newvector");
  allocarray(vec,ih-il);

  return vec-il;
}

void freeintvector(int il,int ih,int *vec)
{
  if (ih<=il) Error("freector");
  vec+=il;
  free(vec);
}

REAL **newmatrix(int il1,int ih1,int il2,int ih2)
{
  REAL **mat;
  int i;

  if (ih1<=il1 || ih2<=il2) Error("newmatrix");
  allocarray(mat,ih1-il1);
  mat-=il1;
  loop (i,il1,ih1) {
    allocarray(mat[i],ih2-il2);
    mat[i]-=il2; }

  return mat;
}

void freematrix(int il1,int ih1,int il2,int ih2,REAL **mat)
{
  int i;

  if (ih1<=il1 || ih2<=il2) Error("freematrix");
  loop (i,il1,ih1) { mat[i]+=il2; free(mat[i]); }
  mat+=il1;
  free(mat);
}
