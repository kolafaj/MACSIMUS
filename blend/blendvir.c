void virial(species_t *spec1,species_t *spec2,species_t *spec) /**** virial */
/* 
  spec must be merged spec1+spec2 
  molecules spec1 and spec2 must have been be minimized
*/
{
  int i,ns=spec->ns,hns,iv,nv;
  vector c1,c2;
  vector *r0,*r,*f;
  double sum,sumq,E,xrg=1/spec->Xopt.core,RT,Ezero;
  double maxexpE=0,corecheck;
#define EPS 1e-4 /* energy imprecise (with weight R^4) for r>1/EPS => 1/EPS is the distance limit */
  
  RT=0.0019872065*spec->Xopt.T; /* R=0.0019872065 kcal/mol/K */
  prt("\n*** Second virial coefficient by random shooting integration ***");
  if (spec->Xopt.V<0) {
    nv=-spec->Xopt.V; 
    spec->Xopt.V=-nv*2; /* if odd... */
    prt("inversion (r -> -r) included, good for both molecules with a dipole"); }
  else {
    prt("no inversion, good for symmetric molecules, bad for both molecules dipolar");
    nv=spec->Xopt.V; }
  
  put3(spec1->ns,spec2->ns,ns)
  if (ns!=spec1->ns+spec2->ns) ERROR(("number of sites in the merged molecule != sum of both molecules"))
  
  prt_("reference point = ");
  switch (spec->Xopt.C) {
    case -1: prt("atomic charge centroid"); break;
    case -2: prt("center-off-mass"); break;
    case -3: prt("geometric center"); break;
    default: if (spec->Xopt.C>=0) prt("atom %d",spec->Xopt.C); 
             else ERROR(("wrong option -C#d",spec->Xopt.C)) }
  
  site=spec->site;

  allocarray(r0,ns);
  allocarray(r,ns);
  allocarray(f,ns);

  hns=spec1->ns;
  loop (i,0,hns) VVV(r0[i],=r[i],=spec1->site[i].r)
  loop (i,hns,ns) VVV(r0[i],=r[i],=spec2->site[i-hns].r)

  Ezero=spec1->Emin+spec2->Emin;
  put3(spec1->Emin,spec2->Emin,Ezero)
  if (fabs(Ezero)>1e-12) prt("WARNING: molecules do not have zero energy - will be corrected");

  testout=0;
  
  VVO(c1,=c2,=0)
  
  if (spec->Xopt.C>=0) {
    if (spec->Xopt.C<hns) VV(c1,=r0[spec->Xopt.C])
    else ERROR(("B2: atom -C%d out of range for the 1st molecule",spec->Xopt.C))
    if (spec->Xopt.C<ns-hns) VV(c2,=r0[hns+spec->Xopt.C])
    else ERROR(("B2: atom -C%d out of range for the 2nd molecule",spec->Xopt.C)) }
  else {    
    double m1=0,m2=0,Z=1;
    
    loop (i,0,ns) {
      switch (spec->Xopt.C) {
        case -1: Z=atom[site[i].type].Z; break;
        case -2: Z=atom[site[i].type].mass; }
      if (i<hns) {
        m1+=Z; 
        VV(c1,+=Z*r0[i]) }
      else {
        m2+=Z; 
        VV(c2,+=Z*r0[i]) } }
      
    if (m1==0 || m2==0) ERROR(("B2: bad center (-C)"))
    VO(c1,/=m1)
    VO(c2,/=m2) }

  putv(c1)
  putv(c2)
    
  loop (i,0,hns)  VVV(r[i],=r0[i],-=c1)
  loop (i,hns,ns) VVV(r[i],=r0[i],-=c2)

  /* minimized molecules (may overlap) are in r0 */
//  loop (i,0,ns) prt("%f %f %f",VARG(r0[i])); exit(0);
  
  sum=sumq=0;
  corecheck=spec->Xopt.core*1.003;
  
  loop (iv,0,nv) {
    double x=rnd()*xrg,xx=x,R;
    int lim=x<EPS;
    int iax=irnd(3),iay=(iax+1)%3;
    double phi=rndcos()*PI;
    double cf=cos(phi),sf=sin(phi),ct=rndcos(),st=sqrt(1-ct*ct);
     
    if (lim) xx=EPS; /* rounding errors*r^4 too big for r>1/EPS */
     
    /* rotate molecule 1 in place */
    loop (i,0,hns) {
      double x=r[i][iax],y=r[i][iay];
      r[i][iax]=x*cf-y*sf;
      r[i][iay]=x*sf+y*cf; }
        
    /* random position molecule 2: random phi and cos(theta) */
    phi=rndcos()*PI;
    cf=cos(phi),sf=sin(phi);
    R=1/xx;
    loop (i,hns,ns) {
      r[i][0]=r0[i][0]+cf*st*R;
      r[i][1]=r0[i][1]+sf*st*R;
      r[i][2]=r0[i][2]+ct*R; }

    E=(Upot(r,f,spec,0)-Ezero)/RT;
    if (E<37) E=exp(-E); else E=0;
  
    if (R<corecheck) Max(maxexpE,E)
    E=(E-1)*Pow4(R);
    if (lim && spec->Xopt.V>0) E/=2; /* rough guess */
    sum+=E;

    if (spec->Xopt.V<0) {
      /* inversion: good for dipole-dipole */
      loop (i,0,hns) VV(r[i],=-r[i])
      st=E;
      E=(Upot(r,f,spec,0)-Ezero)/RT;
      if (E<37) E=exp(-E); else E=0;

      if (R<corecheck) Max(maxexpE,E)
      E=(E-1)*Pow4(R);
      if (lim && spec->Xopt.V>0) E/=2; /* rough guess */
      sum+=E;
      sumq+=Sqr(E+st); }
    else {
      sumq+=E*E; } }

  if (spec->Xopt.V<0) sum/=2,sumq/=4;
  sum*=xrg/nv;
  sumq*=Sqr(xrg)/nv;
  sumq=sqrt((sumq-Sqr(sum))/(nv-1));
  put2(maxexpE,sumq)
  if (maxexpE>sumq*0.1) 
    WARNING(("possibly too big core (-R), close to core exp(-E/RT)=%g",maxexpE));
  sum=2*PI*(Cub(spec->Xopt.core)/3-sum);
  sumq*=2*PI;
  
  prt("B2(%g) = %g +- %g AA^3 (core=%g, -C%d)",spec->Xopt.T,sum,sumq,spec->Xopt.core,spec->Xopt.C);
  prt("B2(%g) = %g +- %g cm3/mol (core=%g, -C%d)",spec->Xopt.T,sum*0.60221418,sumq*0.60221418,spec->Xopt.core,spec->Xopt.C);

  free(f);
  free(r);
  free(r0);
  testout=1;
}

#undef EPS
