typedef void Mvec_f(
  int nc,       /* order of vectors and matrix  M[nc][nc]  */
  double *Mv,   /* Mv[nc]: returning product  M.v  in  Mv[nc] */
  void *Mdef,   /* pointer to sparse-coded matrix  M[nc][nc]  */
  double *v     /* v[nc]: returning product  M.v  in  Mv[nc] */
  );

int conjgrad(
  int     nc,   /* order of vectors and matrix  M[nc][nc]  */
  double *x,    /* x[nc]: initial estimate on input - solution returned */
  Mvec_f *Mvec, /* pointer to user-supplied procedure
                   void Mvec(int nc, double *Mv, void *Mdef, double *v)
                   returning product  M.v  in  Mv[nc] */
  void   *Mdef, /* pointer to sparse-coded matrix  M[nc][nc]
                   passed as parameter to  Mvec  */
  double *b,    /* b[nc]: rhs vector */
  double  eps   /* accuracy:
                   if  0<eps<=1  then relative accuracy  |delta x|/|x| < eps
                   if  eps>=2  then # of iterations */
  );            /* RETURN VALUE: number of *Mvec calls */
