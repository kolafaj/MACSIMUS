/* Fool-proof test: gradient */
{
  int i=0;
  double d=1e-5,Up,Um,fp,fm,aux,m=1,f,grad;
  static double H=1.008,C=12.011,N=14.0067,O=15.9994;
  int quit=0;
  ToIntPtr A,B,V;

  sdsalloc(A,cfg[0]->size); sdscopy(A,cfg[0]);
  sdsalloc(B,cfg[0]->size);
  sdsalloczero(V,cfg[0]->size);

  measure=1;
  for (;;) {
    prt_("enter i=index (i=-1 for xi)  d=step for num.deriv.  m=mass in g/mol  quit=1");
    getdata
      get(i) get(d) get(aux) get(m)
      get(H) get(C) get(N) get(O)
      get(f) get(grad) get(quit)
    checkdata enddata

    if (quit) break;

    sdszero(B);

    ((double*)(A->rp))[i]+=d;
    rhs(B,A,V);
    fp=((double*)(B->rp))[i]; Up=En.pot;
    _n
    put(B->rp[0][0])

    ((double*)(A->rp))[i]-=2*d;
    rhs(B,A,V);
    fm=((double*)(B->rp))[i]; Um=En.pot;

    ((double*)(A->rp))[i]+=d;

    f=(fp+fm)/2*m/Munit;
    grad=(Um-Up)/2/d;
    prt("f=%.9g  -grad=%.9g  dif=%g  U=%g", f,grad, f-grad, (Up+Um)/2);
  }

  free(V);
  free(B);
  free(A);
}
