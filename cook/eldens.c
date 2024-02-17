void eldensities(ToIntPtr A)
{
  int n,m;
  vector dr;
  real *r1,*r2;
  real rr;
  double E;
  sitesite_t *ss;

  if (eldens.ns) {
    if (eldens.ns!=No.s) ERROR(("eldensities: internal")) }
  else {
    eldens.ns=No.N;
    if (No.N!=No.s) ERROR(("eldensities: internal"))
    allocarray(eldens.rho,No.s);
    alloc2Darray(eldens.rhorho,No.s,No.s);
    loop (n,0,No.N) eldens.rhorho[n][n]=0; }

  if (!A) return;

  eldens.rref=A->rp;

  /* simplest code - no cache optimization */
  loop (n,0,No.N) {
    r1=(real*)(A->rp+n); /* atoms only */
    loop (m,0,n) {
      ss=&sstab[spec[molec[n].sp]->si->st][spec[molec[m].sp]->si->st];
      r2=(real*)(A->rp+m); /* atoms only */
#ifdef CUT
      NI(0) if ((rr =dr[0]*dr[0])<box.cq) {
        NI(1) if ((rr+=dr[1]*dr[1])< box.cq) {
          NI2(2)     rr+=dr[2]*dr[2];  if (rr<box.cq)
        eldens.rhorho[n][m]=eldens.rhorho[m][n]=Sqr(ss->a.xi)*exp(-(2*ss->a.q)*(sqrt(rr)/ss->a.r0-1));
      } }
#else
      beginCellImage(0,box.cq)
        eldens.rhorho[n][m]=eldens.rhorho[m][n]=Sqr(ss->a.xi)*exp(-(2*ss->a.q)*(sqrt(rr)/ss->a.r0-1));
      endCellImage
#endif
    } }

  loop (n,0,No.N) {
    eldens.rho[n]=0;
    loop (m,0,No.N) eldens.rho[n]+=eldens.rhorho[m][n];
    E+=sqrt(eldens.rho[n]); }

  En.pot -= E;
}
