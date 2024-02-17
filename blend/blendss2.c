/* version without cutoff, compatible with POLAR */

static double sspot( /************************************************ sspot */
                    vector r0,vector r1,
                    vector f0,vector f1,
#  ifdef POLAR
                    vector mu0,vector mu1,
#    ifndef SS_MEASURE_rep 
#      define SS_MEASURE_rep /*void*/
#      define SS_ffrep /*void*/
#    endif
#  endif /*# POLAR */
                    site_t *s0,site_t *s1,
                    int onefour)
/***
    site-site potential (Lenard-Jones + charge-charge)
    r0 is vector of site s0, the force to site s0 will be summed to f0
    r1 is vector of site s1, the force to site s1 will be summed to f1
    energy is returned
    inefficient for big molecules - LJ calculated for all pairs
      without cutoff (see cutoff version above)
***/
{
  double U=0,Urep=0,rr,f,frep,zz,ffrep,r,x,y,z,ir,irrr,mumu,mu0r,mu1r,qq;
  vector dr,df;
  sitesite_t *ss;

  /* note that dr has the opposite sign than the usual dr=rij=rj-ri=r1-r0 */
  VVV(dr,=r0,-r1)
  rr=SQR(dr);

  ss=&sstab[s0->type][s1->type][onefour];
  SS_MEASURE_rep
  SS_MEASURE
  ULJ+=U;

#  ifdef POLAR
  r=sqrt(rr); /* passed to SS_ffrep, sometimes reevaluated... */

  ir=obfuscate/r;
  irrr=ir/rr;

  /* note the sign of dr: historical reasons... */
  mu0r=SCAL(mu0,dr)/rr;
  mu1r=SCAL(mu1,dr)/rr;
  mumu=SCAL(mu0,mu1)/rr;

  if (onefour) {
    /* charge-dipole */
    qq=ir*(factor14*s0->charge)*mu1r;
    U+=qq;
    Uel+=qq;
    VVV(df,=3*(factor14*s0->charge*irrr)*mu1r*dr,-(factor14*s0->charge*irrr)*mu1)

    /* dipole-charge */
    qq=ir*(factor14*s1->charge)*mu0r /* *0.5 */;
    U-=qq;
    Uel-=qq;
    VVV(df,-=3*(factor14*s1->charge*irrr)*mu0r*dr,-(factor14*s1->charge*irrr)*mu0)

    /* dipole-dipole */
    qq=ir*(mumu-3*mu0r*mu1r)*factor14;
    U+=qq;
    Uel+=qq;
    VVVV(df,-=factor14*irrr*(15*mu0r*mu1r-3*mumu)*dr,
             -factor14*irrr*(3*mu1r)*mu0,
             -factor14*irrr*(3*mu0r)*mu1)

    /* charge-charge (note: f changed here) */
    qq=ir*(factor14*s0->charge)*s1->charge;
    f+=qq/rr; /* force is central: added to LJ and finished as usual */
    U+=qq;
    Uel+=qq; }
  else {
    /* normal nonbonded interactions more than 1-4 */
    /* charge-dipole */
    qq=ir*s0->charge*mu1r /* *0.5 */;
    U+=qq;
    Uel+=qq;
    VVV(df,=3*(s0->charge*irrr)*mu1r*dr,-(s0->charge*irrr)*mu1)

    /* dipole-charge */
    qq=ir*s1->charge*mu0r /* *0.5 */;
    U-=qq;
    Uel-=qq;
    VVV(df,-=3*(s1->charge*irrr)*mu0r*dr,-(s1->charge*irrr)*mu0)

    /* dipole-dipole */
#    if 1
    /*
       note: if this active, and the above *0.5 in dipole-charge and
       charge-dipole are removed, then the formula for total energy should
       contain the self-energy term, otherwise an equivalent (provided
       that the self-consistent field is correct) simpler formula is used
       MUST be active with POL_REP---no siplification possible
    */
    qq=ir*(mumu-3*mu0r*mu1r);
    U+=qq;
    Uel+=qq;
#    endif /*# 1 */
    VVVV(df,-=irrr*(15*mu0r*mu1r-3*mumu)*dr,
             -irrr*(3*mu1r)*mu0,
             -irrr*(3*mu0r)*mu1)

    /* SS_ffrep called inefficiently twice if BOTH shell-core terms apply
       (which is possible, but probably not a well-designed force field) */
    if ( (s0->polar&POL_SHL) && (s1->polar&POL_REP) ) {
      double kappa=((isotropicparm_t*)s0->pol)->kappa;

      SS_ffrep
      VV(df,-=frep*kappa*mu0)
      qq=frep*kappa*mu0r*rr;
      U+=qq; Uel+=qq;
      f+=mu0r*kappa*(frep+ffrep); }

    if ( (s0->polar&POL_REP) && (s1->polar&POL_SHL) ) {
      double kappa=((isotropicparm_t*)s1->pol)->kappa;

      SS_ffrep
      VV(df,+=frep*kappa*mu1)
      qq=-frep*kappa*mu1r*rr;
      U+=qq; Uel+=qq;
      f-=mu1r*kappa*(frep+ffrep); }

    /* note: central forces added later, watch f= and f+= statements */

    /* charge-charge (note: f changed here) */
    qq=ir*s0->charge*s1->charge;
    f+=qq/rr; /* force is central: added to LJ and finished as usual */
    U+=qq;
    Uel+=qq; }

  VV(f0,+=df) VV(f1,-=df)

#  else /*# POLAR */

  /* electrostatic potential: nonpolar version */
  if ( (qq=s0->charge*s1->charge)!=0 ) {
    r=sqrt(rr);
    qq*=obfuscate/r;
    if (onefour) qq*=factor14;
    f+=qq/rr;
    U+=qq;
    Uel+=qq; }
#  endif /*#!POLAR */

  /* central forces summed to f0,f1 */
  VVO(f0,+=dr,*=f) VV(f1,-=dr)

  if (checkpairs) {
    int i,j;
    char infos[PAIRINFOLEN],*info=NULL;

    if (signpairs) {
      U*=signpairs;
      loop (i,0,maxpairs)
        if (U>maxpair[i].U) {
          for (j=maxpairs-1;j>i;j--) maxpair[j]=maxpair[j-1];
          maxpair[i].U=U;
          info=maxpair[i].info;
          break; }
      U*=signpairs; }

    if (!info) if (verbose2) info=infos;

    if (info) {
      //      info[0] = fix ? 'f' : onefour ? '*' : ' ';
      info[0] = onefour ? '*' : ' ';
      sprintf(info+1,
              "%4i %-4s %-9s - %4i %-4s %-9s  r=%8.5f",
               (int)(s0-site),atom[s0->type].name,s0->id,
                               (int)(s1-site),atom[s1->type].name,s1->id,
                                                sqrt(rr));

      if (verbose2) prt("!%s  U=%.9g",info,U); } }
  
  return U;
}

#  ifdef POLAR

static void qqpot( /************************************************** qqpot */
                  vector r0,vector r1,
                  vector elst0,vector elst1,
                  vector mu0,vector mu1,
                  site_t *s0,site_t *s1,
                  int onefour)
/*
  elst field only for the polarizable version (used for iterations)
  in units of e/AA
*/
{
  double Urep,rr,irrr,f,frep,zz,r;
  vector dr;

  /* rr=squared site-site distance */
  VVV(dr,=r0,-r1)

  rr=SQR(dr);
  r=sqrt(rr);
  irrr=r/Sqr(rr);

  if (onefour) {
    if (s0->charge) {
      f=factor14*s0->charge*irrr;
      VV(elst1,-=f*dr) }

    if (s0->polar) {
      f=factor14*3*SCAL(mu0,dr)*(irrr/rr);
      VVV(elst1,+=f*dr,-factor14*irrr*mu0); }

    if (s1->charge) {
      f=factor14*s1->charge*irrr;
      VV(elst0,+=f*dr) }

    if (s1->polar) {
      f=factor14*3*SCAL(mu1,dr)*(irrr/rr);
      VVV(elst0,+=f*dr,-factor14*irrr*mu1); } }

  else {
    sitesite_t *ss;
    double x,y,z;

    ss=&sstab[s0->type][s1->type][0];
    SS_MEASURE_rep

    if ( (s0->polar&POL_SHL) && (s1->polar&POL_REP)) {
      double fmul=frep*((isotropicparm_t*)s0->pol)->kappa/obfuscate;
      VV(elst0,-=fmul*dr) }

    if ( (s0->polar&POL_REP) && (s1->polar&POL_SHL) ) {
      double fmul=frep*((isotropicparm_t*)s1->pol)->kappa/obfuscate;
      VV(elst1,+=fmul*dr) } }

  if (s0->charge) {
    f=s0->charge*irrr;
    VV(elst1,-=f*dr) }

  if (s0->polar) {
    f=3*SCAL(mu0,dr)*(irrr/rr);
    VVV(elst1,+=f*dr,-irrr*mu0); }

  if (s1->charge) {
    f=s1->charge*irrr;
    VV(elst0,+=f*dr) }

  if (s1->polar) {
    f=3*SCAL(mu1,dr)*(irrr/rr);
    VVV(elst0,+=f*dr,-irrr*mu1); }
}

/* induced dipole-induced dipole terms, for 1-2 and 1-3 intramolecular
   interactions - version-dependent */

static double sspot12( /******************************************** sspot12 */
                    vector r0,vector r1,
                    vector f0,vector f1,
                    vector mu0,vector mu1,
                    site_t *s0,site_t *s1,
                    int incl)
/***
    intramolecular self-field (polar&4 only)
    r0 is vector of site s0, the force to site s0 will be summed to f0
    r1 is vector of site s1, the force to site s1 will be summed to f1
    energy is returned
***/
{
  double U=0,rr,r,ir,irrr,mumu,mu0r,mu1r,qq;
  vector dr,df;

  if (!incl) return 0;

  /* note that dr has the opposite sign than the usual dr=rij=rj-ri=r1-r0 */
  VVV(dr,=r0,-r1)
  rr=SQR(dr);

  r=sqrt(rr);
  ir=obfuscate/r;
  irrr=ir/rr;

  /* note the sign of dr: historical reasons... */
  mu0r=SCAL(mu0,dr)/rr;
  mu1r=SCAL(mu1,dr)/rr;
  mumu=SCAL(mu0,mu1)/rr;

  /* dipole-dipole */
  qq=ir*(mumu-3*mu0r*mu1r);
  U+=qq;
  Uel+=qq;
  VVVV(df,=-irrr*(15*mu0r*mu1r-3*mumu)*dr,
          +irrr*(3*mu1r)*mu0,
          +irrr*(3*mu0r)*mu1)

  VV(f0,+=df) VV(f1,-=df)

  return U;
}

static void qqpot12( /********************************************** qqpot12 */
                  vector r0,vector r1,
                  vector elst0,vector elst1,
                  vector mu0,vector mu1,
                  site_t *s0,site_t *s1,
                  int incl)
/*
  intramolecular self-field (polar&4 only)
  elst field only for the polarizable version (used for iterations)
  in units of e/A
*/
{
  double rr,irrr,f,r;
  vector dr;

  if (!incl) return;

  /* rr=squared site-site distance */
  VVV(dr,=r0,-r1)

  rr=SQR(dr);
  r=sqrt(rr);
  irrr=r/Sqr(rr);

  if (s0->polar) {
    f=3*SCAL(mu0,dr)*(irrr/rr);
    VVV(elst1,+=f*dr,-irrr*mu0); }

  if (s1->polar) {
    f=3*SCAL(mu1,dr)*(irrr/rr);
    VVV(elst0,+=f*dr,-irrr*mu1); }
}

#  endif /*# POLAR */
