/*
   To be #included from simmeas.c
   Measuring CoM, velocity, momentum
   CAVEAT: to some extent duplicates removedrifts
*/

void measuredrifts(void) /************************************ measuredrifts */
/* for whole configuration */
{
  int n,i,k,ns;
  molecule_t *mn;
  vector *r,*v;
  siteinfo_t *si;
  real mi,mass;
  ToIntPtr centered=cfg[0];

#ifdef SLAB
  if (slab.sp>=nspec) return;
#endif

  if (lag.CM || lag.AM) {
#ifndef FREEBC
  /* move molecules periodically so that their CMs are max Lh from box.center */
    sdsalloc(centered,cfg[0]->size);
    sdscopy(centered,cfg[0]);
    loop (n,FROM,No.N) {
      double m=0;
      vector c={0,0,0};

      mn=molec+n;
#ifdef SLAB
      if (slab.sp<0) { if (mn->sp>abs(slab.sp)) continue; }
      else { if (mn->sp!=slab.sp) continue; }
#endif
      si=spec[mn->sp]->si;
      ns=mn->ns;
      r=rof(mn,centered->rp);
      loop (i,0,ns) {
        m+=si[i].mass;
        VV(c,+=si[i].mass*r[i]) }
      VO(c,/=m)
        loop (k,0,3) {
        if (c[k]-box.center[k]> box.Lh[k]) loop (i,0,ns) r[i][k]-=box.L[k];
        if (c[k]-box.center[k]<-box.Lh[k]) loop (i,0,ns) r[i][k]+=box.L[k]; }
    }
#endif

    /* center of mass */
#ifdef SLAB
    /* copy of CoM with selection of species */
    VO(En.CM,=0)
    mass=0;

    loop (n,FROM,No.N) {
      mn=molec+n;
      if (slab.sp<0) { if (mn->sp>abs(slab.sp)) continue; }
      else { if (mn->sp!=slab.sp) continue; }
      si=spec[mn->sp]->si;
      r=rof(mn,centered->rp);

      loop (i,0,mn->ns) {
        mi=si[i].mass;
        mass+=mi;
        VV(En.CM,+=mi*r[i]) } }

    if (mass) VO(En.CM,/=mass)
    // NOT VV(CM,-=box.center): to be applied elsewhere
#else
    CoM(En.CM,centered);
#endif

    if (lag.CM) {
      StaSet(0,lag.CM,2,lag.n);
      StaAdd("CMx",En.CM[0]);
      StaAdd("CMy",En.CM[1]);
      StaAdd("CMz",En.CM[2]);
      StaAdd("CM^2",SQR(En.CM)); } }

  /* linear momentum (velocity drift) */
  if (lag.LM) {
    mass=0;

    VO(En.LM,=0)
    loop (n,FROM,No.N) {
      mn=molec+n;
#ifdef SLAB
      if (slab.sp<0) { if (mn->sp>abs(slab.sp)) continue; }
      else { if (mn->sp!=slab.sp) continue; }
#endif
      si=spec[mn->sp]->si;
      ns=mn->ns;
      v=rof(mn,cfg[1]->rp);
      loop (i,0,ns) {
        VV(En.LM,+=si[i].mass*v[i])
        mass+=si[i].mass; } }

    //    VO(En.LM,*=(massunit*lengthunit/timeunit))
    // cannot record in SI .cp because close to out of float range
    VO(En.LM,/=h) /* there is r(t)-r(t-h) in cfg[1]->rp */
    StaSet(0,lag.LM,2,lag.n);
    StaAdd("LMx [kg m/s]",En.LM[0]*(massunit*lengthunit/timeunit));
    StaAdd("LMy [kg m/s]",En.LM[1]*(massunit*lengthunit/timeunit));
    StaAdd("LMz [kg m/s]",En.LM[2]*(massunit*lengthunit/timeunit));

    /* translational temperature */
    if (mass) StaAdd("Ttr cluster (TLM)",En.TLM=SQR(En.LM)/(3*mass)); }

  /* angular momentum with respect to the center of mass */
  if (lag.AM) {
    vector I[3];
    vector rv,omega;
    real det,rr;
    int j,k;

    VO(En.AM,=0)
    memset(I,0,sizeof(I));

    loop (n,0,No.N) {
      mn=molec+n;
#ifdef SLAB
      if (slab.sp<0) { if (mn->sp>abs(slab.sp)) continue; }
      else { if (mn->sp!=slab.sp) continue; }
#endif
      si=spec[mn->sp]->si;
      ns=mn->ns;
      r=rof(mn,centered->rp);
      v=rof(mn,cfg[1]->rp);

      loop (i,0,ns) {
        mi=si[i].mass;
        VVV(omega,=r[i],-En.CM) /* omega=temporary vector here = r wrt center */
        VECT(rv,omega,v[i])
        VV(En.AM,+=mi*rv)
        rr=SQR(omega);
        loop (j,0,DIM) loop (k,j,DIM)
          /* matrix is symetric: j<=k */
          //          I[k][j]= // debug only
          I[j][k]+=mi*((j==k)*rr-omega[j]*omega[k]); } } /* n */

    //    putv(I[0])    putv(I[1])    putv(I[2])

    VO(En.AM,/=h)
    //    VO(En.AM,*=(massunit*lengthunit*lengthunit/timeunit)/h)
    // cannot record in SI in .cp because out of float range
    StaSet(0,lag.AM,2,lag.n);
    StaAdd("AMx [kg m2/s]",En.AM[0]*(massunit*lengthunit*lengthunit/timeunit));
    StaAdd("AMy [kg m2/s]",En.AM[1]*(massunit*lengthunit*lengthunit/timeunit));
    StaAdd("AMz [kg m2/s]",En.AM[2]*(massunit*lengthunit*lengthunit/timeunit));

    /* matrix inversion */
    det =
       I[0][0]*I[1][1]*I[2][2] + 2*I[0][1]*I[1][2]*I[0][2]
     - I[1][1]*Sqr(I[0][2]) - I[2][2]*Sqr(I[0][1]) - I[0][0]*Sqr(I[1][2]);
    /* omega/det is the angular velocity */
    omega[0] =
      ( En.AM[0]*(I[1][1]*I[2][2]-Sqr(I[1][2]))
      + En.AM[1]*(I[1][2]*I[0][2]-I[0][1]*I[2][2])
      + En.AM[2]*(I[0][1]*I[1][2]-I[1][1]*I[0][2]) ) ;
    omega[1] =
      ( En.AM[1]*(I[2][2]*I[0][0]-Sqr(I[0][2]))
      + En.AM[2]*(I[0][2]*I[0][1]-I[1][2]*I[0][0])
      + En.AM[0]*(I[1][2]*I[0][2]-I[2][2]*I[0][1]) ) ;
    omega[2] =
      ( En.AM[2]*(I[0][0]*I[1][1]-Sqr(I[0][1]))
      + En.AM[0]*(I[0][1]*I[1][2]-I[0][2]*I[1][1])
      + En.AM[1]*(I[0][2]*I[0][1]-I[0][0]*I[1][2]) ) ;

    //    prt("%g %g %g  %g %g %g @LMomega",VARG(En.LM),omega[0]/det,omega[1]/det,omega[2]/det);

    //    prt("%g %g %g  %g %g %g  %g %g %g @I",VARG(I[0]),VARG(I[1]),VARG(I[2]));

    /* rotational temperature */
    StaAdd("Trot cluster (TAM)",En.TAM=SCAL(En.AM,omega)/(3*det)); }

  if (cfg[0]!=centered) free(centered);
}

void measure1drift(int n) /*********************************** measure1drift */
/* for 1 molecule, directly printed */
{
  molecule_t *mn=molec+n;
  int i;
  vector *r=rof(mn,cfg[0]->rp),*v=rof(mn,cfg[1]->rp);
  siteinfo_t *si=spec[mn->sp]->si;
  real mi,mass;
  vector CM,LM,AM;
  vector I[3];
  vector rv,omega;
  double det,rr;
  int j,k;

  int nx,ix;
  int ixmin,imin,nxmin,nsx=0;
  double rrmin=9e99;

  underline("evaporated molecule");

  loop (nx,0,No.N) if (n!=nx) {
    molecule_t *mnx=molec+nx;
    vector *rx=rof(mnx,cfg[0]->rp),dr;

    loop (i,0,mn->ns) {
      loop (ix,0,mnx->ns) {
        VVV(dr,=rx[ix],-r[i])
#ifndef FREEBC
        loop (k,0,DIM) {
          while (dr[k]>box.Lh[k]) dr[k]-=box.L[k];
          while (dr[k]<-box.Lh[k]) dr[k]+=box.L[k]; }
#endif
        rr=SQR(dr);
        if (rr<rrmin) {
          rrmin=rr;
          nsx=mnx->ns;
          imin=i;
          ixmin=i;
          nxmin=nx; } } } }
    prt("molecule %d has evaporated and will be removed:\n\
minimum distance = %g (%d.%d/%d from %d.%d/%d) MINDIST",
        n,
        sqrt(rrmin),nxmin,ixmin,nsx,n,imin,mn->ns);
  
  VO(CM,=0)
  mass=0;

  loop (i,0,mn->ns) {
    mi=si[i].mass;
    mass+=mi;
    VV(CM,+=mi*r[i]) }

  if (!mass) ERROR(("the evaporated molecule has mass=0"))
  VO(CM,/=mass) 

  prt("mass = %.9g p.u. = %.9g g/mol = %.9g kg",mass,mass*Munit,mass*massunit);
  prt("center of mass = [ %9g %9g %9g ] AA",VARG(CM));

  /* linear momentum (velocity drift) */
  VO(LM,=0)
  loop (i,0,mn->ns) VV(LM,+=si[i].mass*v[i])

  VO(LM,/=h) /* there is r(t)-r(t-h) in cfg[1]->rp */
  prt("linear momentum = [ %9g %9g %9g ] kg m/s",
      LM[0]*(massunit*lengthunit/timeunit),
      LM[1]*(massunit*lengthunit/timeunit),
      LM[2]*(massunit*lengthunit/timeunit));

  prt("translational temperature = %.9g K",SQR(LM)/(3*mass));
  
  /* angular momentum with respect to the center of mass */
  VO(AM,=0)
  memset(I,0,sizeof(I));

  loop (i,0,mn->ns) {
    mi=si[i].mass;
    VVV(omega,=r[i],-CM) /* omega=temporary vector here = r wrt center */
    VECT(rv,omega,v[i])
    VV(AM,+=mi*rv)
    rr=SQR(omega);
    loop (j,0,DIM) loop (k,j,DIM)
        /* matrix is symetric: j<=k */
        //          I[k][j]= // debug only
      I[j][k]+=mi*((j==k)*rr-omega[j]*omega[k]); }

  VO(AM,/=h)
  prt("angular momentum = [ %9g %9g %9g ] kg m2/s",
      AM[0]*(massunit*lengthunit*lengthunit/timeunit),
      AM[1]*(massunit*lengthunit*lengthunit/timeunit),
      AM[2]*(massunit*lengthunit*lengthunit/timeunit));

  /* matrix inversion */
  det =
      I[0][0]*I[1][1]*I[2][2] + 2*I[0][1]*I[1][2]*I[0][2]
    - I[1][1]*Sqr(I[0][2]) - I[2][2]*Sqr(I[0][1]) - I[0][0]*Sqr(I[1][2]);
  /* omega/det is the angular velocity */
  omega[0] =
    ( AM[0]*(I[1][1]*I[2][2]-Sqr(I[1][2]))
    + AM[1]*(I[1][2]*I[0][2]-I[0][1]*I[2][2])
    + AM[2]*(I[0][1]*I[1][2]-I[1][1]*I[0][2]) ) ;
  omega[1] =
    ( AM[1]*(I[2][2]*I[0][0]-Sqr(I[0][2]))
    + AM[2]*(I[0][2]*I[0][1]-I[1][2]*I[0][0])
    + AM[0]*(I[1][2]*I[0][2]-I[2][2]*I[0][1]) ) ;
  omega[2] =
    ( AM[2]*(I[0][0]*I[1][1]-Sqr(I[0][1]))
    + AM[0]*(I[0][1]*I[1][2]-I[0][2]*I[1][1])
    + AM[1]*(I[0][2]*I[0][1]-I[0][0]*I[1][2]) ) ;

  /* rotational temperature */
  prt("rotational temperature = %.9g K\n\
(for more info, see remove1mol below)",SCAL(AM,omega)/(3*det));

}
