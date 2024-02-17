#define DIM 2

typedef real vector[DIM];

#define VO(A,O) { A[0] O; A[1] O; }
#define VV(A,B) { A[0] B[0]; A[1] B[1]; }
#define VVO(A,B,O) { A[0] B[0] O; A[1] B[1] O; }
#define VVV(A,B,C) { A[0] B[0] C[0]; A[1] B[1] C[1]; }
#define VVVO(A,B,C,O) { A[0] B[0] C[0] O; A[1] B[1] C[1] O; }
#define VVVV(A,B,C,D) \
  { A[0] B[0] C[0] D[0]; A[1] B[1] C[1] D[1]; }
#define SQR(A) (A[0]*A[0]+A[1]*A[1])
#define SQRD(A,B) (Sqr(A[0]-B[0])+Sqr(A[1]-B[1]))
#define SCAL(A,B) (A[0]*B[0]+A[1]*B[1])
#define SUM(A) (A[0]+A[1])
#define PROD(A) (A[0]*A[1])
#define VECT(A,B) { A[0]=-B[1]; A[1]=B[0]; }

/* #define VSET(A,X,Y) { A[0]=X; A[1]=Y; } */

#define VARG(A) (double)A[0],(double)A[1]

#ifdef putv
#undef putv
#endif

#define putv(_X) prt("%11s=(%13.6g %13.6g)", #_X,(double)_X[0],(double)_X[1]);
