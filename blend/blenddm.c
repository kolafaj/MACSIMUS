/* Dipole moment of a molecule. 
   see also totalcharge() in blendmed
   #included by blendmin.c 
   the center of dipole moment is atomic charge centroid
   (irrelevant for neutral molecule)

   also prtdist()
*/

double dipolemoment(species_t *spec,double charge,double *Qdiag)
                                                           /*** dipolemoment */
{
  int i;
  int ns=spec->ns;
  vector dm,cZ;
  double sZ=0,Z;
  int a,b;
  vector Q[3];
  double *passQ[3];

  VVO(cZ,=dm,=0)
  site=spec->site;

  loop (i,0,ns) {
    switch (Xopt.C) {
      case -1: Z=atom[site[i].type].Z; break;
      case -2: Z=atom[site[i].type].mass; break;
      case -3: Z=1; break;
      default: Z=Xopt.C==i; }
    sZ+=Z;
    VV(cZ,+=Z*site[i].r)
    VV(dm,+=site[i].charge*site[i].r) 
#ifdef POLAR
    VV(dm,+=mu[i])
#endif
    }

  if (sZ) {

    VO(cZ,/=sZ) 

    /* quadrupole moment
      Qab' = 3/2*SUM_i q_i Da Db
      Qab  = Qab' - 1/3 Tr Qab (remove trace)
    */
    loop (a,0,3) loop (b,0,3) Q[a][b]=0;

    loop (i,0,ns) {
      vector dr;
      double rr;
#ifdef POLAR
      double mudr;
#endif

      VVV(dr,=site[i].r,-cZ)
      rr=SQR(dr);
#ifdef POLAR
      mudr=SCAL(mu[i],dr);
#endif

      loop (a,0,3) loop (b,0,3)
        Q[a][b] += 3*site[i].charge*(
#ifdef POLAR
          dr[a]*mu[i][b]+dr[b]*mu[i][a]+
#endif
          dr[a]*dr[b]);

      loop (a,0,3) 
        Q[a][a] -= 
#ifdef POLAR
          2*mudr+
#endif
          site[i].charge*rr; }

    loop (a,0,3) {
      passQ[a]=Q[a];
      if (option('v')&4) prt("!%9s (%12.8f %12.8f %12.8f)",
                             a==0?"Q [eAA2] =":"",Q[a][0],Q[a][1],Q[a][2]); }
    Jacobi(3,passQ,NULL,spec->Xopt.Jeps);
    /* factor 1/2 in version 1.9e */
    loop (a,0,3) Qdiag[a]=Q[a][a]/2; }

  else {
    WARNING(("dipole and quadrupole: undefined center\n\
*** (Z undefined in par-file or wrong atom in option -C)"))
    return 0; }

  VV(dm,-=charge*cZ)

  return sqrt(SQR(dm));
}

#ifdef POLAR
/* this should become a module... */
#include "matrix.c"
#include "invmat.c"

void totalpolarizability(species_t *spec) /************* totalpolarizability */
{
  double alpha=0;
  int i,ns=spec->ns;

  if (polar&4) {
    /* J. Applequist etal, JACS 94, 2952 (1972) */
    int n,nps=0;
    int j,k,l;
    int ii,jj;
    double **A;
    vector a[3];
    double *passa[3];

    loop (i,0,ns) if (site[i].polar) {
      if (site[i].polar&POL_AXI) return; /* sorry, not supported */
      else if (site[i].polar&POL_ISO
               && ((isotropicparm_t*)site[i].pol)->alpha!=0) {
	nps++; 
	site[i].polar|=POL_NZISO; } }

    if (!nps) return;
    n=nps*3;

    A=newmatrix(0,n,0,n);

    ii=0;
    loop (i,0,ns) if (site[i].polar&POL_NZISO) {
      isotropicparm_t *pol=(isotropicparm_t*)site[i].pol;
      double *ri=site[i].r;

      loop (k,ii,ii+3)
        loop (l,ii,ii+3)
          A[k][l]=(k==l)/pol->alpha;

      jj=ii+3;
      loop (j,i+1,ns) if (site[j].polar&POL_NZISO) {
        double *rj=site[j].r;
	vector dr;
        double rr,f,rr3;

	VVV(dr,=rj,-ri)
        rr=SQR(dr); 
	f=-3/(rr*rr*sqrt(rr));
	rr3=rr/3;
        
        loop (k,0,3)
          loop (l,0,3)
            A[ii+k][jj+l]=A[jj+l][ii+k]=f*(dr[k]*dr[l]-rr3*(l==k)); 
	jj+=3; } 
      ii+=3; }

    if (!GJinvmat(0,n,A)) return;

    loop (k,0,3)
      loop (l,0,3) a[l][k]=0;

    loop (i,0,nps)
      loop (j,0,nps)
        loop (k,0,3)
          loop (l,0,3) a[l][k]+=A[i*3+k][j*3+l];

    alpha=(a[0][0]+a[1][1]+a[2][2])/3;

    prt("! total polarizability tensor [AA^3]:");
    loop (k,0,3) {
      prt_("! ");
      prt(sites_f,VARG(a[k]));
      passa[k]=a[k]; }

    Jacobi(3,passa,NULL,spec->Xopt.Jeps);
    prt_("! diagonalized: ");
    prt(sites_f,a[0][0],a[1][1],a[2][2]);

    freematrix(0,n,0,n,A); }

  else { /* ! polar&4; NOW JUST SUM, which is not correct - should include 1-4 and more */
    prt ("! WARNING: just sum of polarizabilities - should include self-field of 1-4 etc.");
    loop (i,0,ns) if (site[i].polar) {
      if (site[i].polar&POL_AXI) return; /* sorry, not supported */
      else if (site[i].polar&POL_ISO) {
        isotropicparm_t *pol=(isotropicparm_t*)site[i].pol;
        alpha+=pol->alpha; } } }

  prt_("! total polarizability =");
  prt_(charge_f,alpha);
  prt(" AA^3 = %g a.u.", alpha*6.74833396);

}
#endif

