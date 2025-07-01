/*
   #included in constrd.c, which is ugly

   POLAR only: #included from constrd.c
   ADDED: FQ4 support (POLAR&32)
   REMOVED in 2011 (see old+misc/scf+lambda.c):
      LAMBDA (convergence)
      optimization "if (iomegap)"
   REMOVED in 9/2017: divergence detection for epsp<0
*/

int selffield(ToIntPtr B,ToIntPtr A, /**************************** selffield */
              double epsp,double omegap,int run)
/*
  One iteration of the self-consistent induced dipoles.
  Forces are passed to the reference atom from the Drude charge.
  Virial of force and the pressure tensor are calculated.

  Calling scheme:
    scf.nit=0;
    do {
      zeroEn();
      forces(FORCES,POSITIONS);
    } while (!selffield(FORCES,POSITIONS,epsp,omegap,0|1));

  Parameters:
    B = forces (output)
    A = coordinates (1st part) not changed, Drude positions (2nd part) updated
    epsp = end-of-iteration criterion, see the Return value
    omegap = mixing iteration parameter
    run = affects how statistics is reported:
      run=1: called from simulation
      run=0: called from measurement of P by virtual volume change (or similar)
    
  Return value = 1 if:
    epsp>0:    if sufficiently accurate (max. dipole change < epsp)
    epsp<=-1:  fixed number of iterations |epsp| reached
    -1<epsp<0: at least 1+scf.maxit/20 iterations &&
               one-step error has increased (from the previous iteration) &&
               one-step error < |epsp| 
    epsp=0:    at least 2+scf.maxit/5 iterations && 
               one-step error has increased (from the previous iteration)
*/
{
  molecule_t *mn;
  siteinfo_t *si;
  int ns,n,i,nmax=-1,imax=-1;
  vector *fc,*fpol,*rpol,dr,dd;
#if POLAR&(8|32)
  vector *rc;
#endif /*# POLAR&(8|32) */
  int ierr=0,iret=1;
  double rr,f,maxrrold=0,maxrr=0,fdr;
  double iomegap=1-omegap;
  double err=0; /* standard deviation of induced dipole moments (in 1st iteration) */
  double maxerr=0; /* max error of induced dipole moments (in 1st iteration) */
#ifndef FREEBC
  double byerd;
#endif /*# FREEBC */
  static int option_v=-3333;
  static int nbad;
  static int scfmaxit0;

  if (!scfmaxit0) scfmaxit0=scf.maxit;
  scf.nit++;
  Max(scf.cycmaxit,scf.nit)

  if (epsp<0 && fabs(epsp)>scf.maxit) {
    WARNING(("Negative eps in selffield (%g) in abs.val. > scf.maxit=%d\n\
*** scf.maxit increased to %d",epsp,scf.maxit,(int)(1+fabs(epsp))))
    scf.maxit=1+fabs(epsp); }

  En.self=0;
#if POLAR&32
  En.fqself=0;
#endif /*# POLAR&32 */

  if (option_v==-3333) option_v=option('v'); /* |=8 if bad convergence */

  loop (n,FROM,No.N) {
    mn=molec+n;
    ns=mn->ns;
    si=spec[mn->sp]->si;
#if POLAR&(8|32)
    rc=rof(mn,A->rp);
#endif /*# POLAR&(8|32) */
    rpol=polarrof(mn,A->rp);
    fc=rof(mn,B->rp);
    fpol=polarrof(mn,B->rp);

#if POLAR&32
    {
      fq4_t *fq4;

      looplist (fq4,spec[mn->sp]->fq4) {
        int H1=fq4->indx[0];
        int H2=fq4->indx[1];
        int L1=fq4->indx[2];
        int L2=fq4->indx[3];
        double dqH1=fq4->H1H1*fpol[H1][0]+fq4->H1H2*fpol[H2][0]+fq4->HL*(fpol[L1][0]+fpol[L2][0]);
        double dqH2=fq4->H1H1*fpol[H2][0]+fq4->H1H2*fpol[H1][0]+fq4->HL*(fpol[L1][0]+fpol[L2][0]);
        double dqL1=fq4->L1L1*fpol[L1][0]+fq4->L1L2*fpol[L2][0]+fq4->LH*(fpol[H1][0]+fpol[H2][0]);
        double dqL2=fq4->L1L1*fpol[L2][0]+fq4->L1L2*fpol[L1][0]+fq4->LH*(fpol[H1][0]+fpol[H2][0]);

        /* note that later En.self is divided by 2 */
        En.fqself += (fq4->AHH*(Sqr(dqH1)+Sqr(dqH2))
                     +fq4->ALL*(Sqr(dqL1)+Sqr(dqL2))
                     +fq4->AHL2*(dqH1+dqH2)*(dqL1+dqL2));

        VV(dd, =rpol[H1][0]*rc[H1])
        VV(dd,+=rpol[H2][0]*rc[H2])
        VV(dd,+=rpol[L1][0]*rc[L1])
        VV(dd,+=rpol[L2][0]*rc[L2])

        /* relaxation */
        rpol[H1][0]=iomegap*(rpol[H1][0])+omegap*(si[H1].charge+dqH1);
        rpol[H2][0]=iomegap*(rpol[H2][0])+omegap*(si[H2].charge+dqH2);
        rpol[L1][0]=iomegap*(rpol[L1][0])+omegap*(si[L1].charge+dqL1);
        rpol[L2][0]=iomegap*(rpol[L2][0])+omegap*(si[L2].charge+dqL2);

        VV(dd,-=rpol[H1][0]*rc[H1])
        VV(dd,-=rpol[H2][0]*rc[H2])
        VV(dd,-=rpol[L1][0]*rc[L1])
        VV(dd,-=rpol[L2][0]*rc[L2])

        f=SQR(dd)/Sqr(omegap);
        Max(maxerr,f)
        err+=f; ierr++;
      }
    }
#endif /*# POLAR&32 */

    loop (i,0,ns) if (si[i].qtype&QTYPE_DRUDE) {
      VVV(dd,=dr,=rpol[i])
      rr=SQR(dr);
      Max(maxrrold,rr)
      if (rr>=Sqr(el.minqq))
        ERROR(("DRUDE old cfg: %d[%d] aux charge dist = %g",n,i,sqrt(rr)))

      /* NOTE (9/2001): the virial of force equals -2*self energy and
         is treated separately because of the Ewald corrections */

#if defined(COULOMB) && COULOMB<0
      /* Ewald correction for dipole self-term */
      f=si[i].charge*si[i].chargepol;
      En.el+=f*exacterud_sqrt_1(rr,&byerd);
      f*=byerd;
#  define ENPVIR En.Pvir
      MADDFORCES(fpol[i],fc[i],f,dr)
#  undef ENPVIR
#endif /*# defined(COULOMB) && COULOMB<0 */

#if POLAR&8
      /* axial polarizability */
      if ( (tozz=si[i].tozz)>=0 ) {
#  if POLAR&2
        if (si[i].Esat)
          /* saturation */
          ERROR(("axial saturated polarizability not supported"))
        else
#  else /*# POLAR&2 */
        {
          /* axial+no saturation */
          vector drtozz,Ezz,Exx,rpolnew;
          double zz,rr,rf,a;

          VVV(drtozz,=rc[tozz],-rc[i])
          rr=SQR(drtozz); rf=SCAL(fpol[i],drtozz);
          zz=rf/rr;
          VV(Ezz,=zz*drtozz)
          VVV(Exx,=fpol[i],-Ezz)
          VVV(rpolnew,=si[i].alpha_qq*Exx,+si[i].alphazz_qq*Ezz)

          /* relaxation */
          VVVV(dd,-=rpol[i],=iomegap*rpol[i],+omegap*rpolnew)

          /* forces passed to both atoms */
          /* WARNING: not for negative -p */
          VV(fc[i],+=fpol[i])
          VVV(drtozz,=fpol[i],-zz*drtozz)
          a=rf/rr*(si[i].alpha_qq-si[i].alphazz_qq);
          VO(drtozz,*=a)
          VV(fc[i],+=drtozz)
          VV(fc[tozz],-=drtozz) }
#    if PRESSURETENSOR
#      error not implemented...
#    endif /*# PRESSURETENSOR */
#  endif /*#!POLAR&2 */
        /* forces passed to both atoms -- eventually move here after axial
           saturated is implemented */
      }
      else
#endif /*# POLAR&8 */
        {
#if POLAR&2
          /* new pos of auxiliary charges = new dipole moments */
          if (si[i].Esat) {
            /* saturation */
            double asat=si[i].alpha_qq
              / (0.5+sqrt(0.25+si[i].alpha_qq*SQR(fpol[i])/si[i].Esat));

            /* relaxation */
            VVVV(dd,-=rpol[i],=iomegap*rpol[i],+(omegap*asat)*fpol[i]) }
        else
#endif /*# POLAR&2 */
          {
            /* standard: no saturation */
            VVVV(dd,-=rpol[i],=iomegap*rpol[i],+(omegap*si[i].alpha_qq)*fpol[i])
          }
        /* forces passed to the central atom */
        VV(fc[i],+=fpol[i]) }

      /* twice the dipole self-energy */
      f=SCAL(rpol[i],fpol[i]);
#if POLAR&2
      if (si[i].Esat)
        En.self += si[i].Esat*log1(f/si[i].Esat);
#  if PRESSURETENSOR&PT_VIR
#    error POLAR&2 & PRESSURETENSOR&PT_VIR - no support
#  endif /*# PRESSURETENSOR&PT_VIR */
      else
#endif /*# POLAR&2 */
        En.self += f; /* = -virial (will be added to En.vir later) */

#if PRESSURETENSOR&PT_VIR
      /* see the manual, Sec. "Pressure tensor for Drude oscillators" */
      VVV(En.Pvir,-=rpol[i],*fpol[i])
#  if PRESSURETENSOR&PT_OFF
#    if 1
      En.Pvir[3]-=rpol[i][1]*fpol[i][2];
      En.Pvir[4]-=rpol[i][2]*fpol[i][0];
      En.Pvir[5]-=rpol[i][0]*fpol[i][1];
#    else /*# 1 */
      /* debug: should be the same */
      En.Pvir[3]-=fpol[i][1]*rpol[i][2];
      En.Pvir[4]-=fpol[i][2]*rpol[i][0];
      En.Pvir[5]-=fpol[i][0]*rpol[i][1];
#    endif /*#!1 */
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_VIR */

#if 0
      /* old version: equivalent if SCF exactly solved */
      En.self += si[i].alpha_qq*SQR(fpol[i]);
#endif /*# 0 */

#if 0 /* MOVED and updated for axial... */
      /* forces passed to the central atom */
      VV(fc[i],+=fpol[i])
#endif /*# 0 */

      /* check that aux charges are close enough to the center */

      rr=SQR(rpol[i]);
      if (rr>=Sqr(el.minqq))
        ERROR(("DRUDE new cfg: %d[%d] aux charge dist = %g",n,i,sqrt(rr)))

      if (rr>maxrr) { maxrr=rr; imax=i; nmax=n; }

      /* error of the iteration */
      f=SQR(dd)/Sqr(omegap)*si[i].qqpol;
      //      put3(f,epsp,omegap)
      Max(maxerr,f)
      err+=f; ierr++; } }

  maxerr=sqrt(maxerr);
  err=sqrt(err/ierr);

  if (option('a') && option_v&8) {
    static int first=1;
    static double lastmaxerr=1,lasterr=1;
    double ratio=0,maxratio=0;

    if (first) {
      prt("omegap A      maxr  -> maxr  mol[site]  err_max  err_std  rate_max  rate_std");
      first=0; }
    if (lastmaxerr) maxratio=maxerr/lastmaxerr; else lastmaxerr=0;
    if (lasterr) ratio=err/lasterr; else lasterr=0;
    prt("%4.2f %.3g->%.3g %d[%d] %.5g %.5g %.7f %.8f SCF",
        omegap, /* option('@') removed */
        sqrt(maxrrold),sqrt(maxrr),nmax,imax,
        maxerr,err,
        maxratio,ratio);
    lastmaxerr=maxerr,lasterr=err; }

  En.self/=2;
#if POLAR&32
  En.fqself/=2;
  En.self += En.fqself; /* moved here from constrd.c */
#endif /*# POLAR&32 */

  if (scf.nit==1) {
    /* here we record the errors of the predictor (in the 1st iteration) */
    StaSet(0,2,2,lag.n);
    if (run) {
      En.sumpolstderr+=En.polstderr=err;
      StaAdd("polar 1st iter stderr",err);
      En.sumpolmaxerr+=En.polmaxerr=maxerr;
      StaAdd("polar 1st iter maxerr",maxerr); }
    nbad=0; }
  else /* scf.nit>1 */
    StaAdd(run?"polar iter rate":"polar iter rate aux",scf.rate=maxerr/scf.maxerr);

  if (epsp<=-1) {
    /*** fixed number of iterations |epsp| ***/
    scf.maxdr=sqrt(maxrr);
    iret=scf.nit>=fabs(epsp);
    goto ret; }

  if (epsp>-1 && epsp<=0) {
    /*** try to reach maximum precision ***/
    if (option('v')&4) put3(scf.nit,maxerr,maxerr/scf.maxerr)
    if (scf.nit<=1+scf.maxit/20) goto not_yet;
    if (epsp==0 && scf.nit<=2+scf.maxit/5) goto not_yet;
    if (epsp<0 && maxerr>fabs(epsp)) goto not_yet;
    if (maxerr==0 || maxerr/scf.maxerr>=1) goto ret;
   not_yet:
    scf.maxerr=maxerr;
    if (scf.nit>scf.maxit) {
      scf.maxit=scf.maxit*8/5+1;
      WARNING(("selffield max. precision algorithm (epsp=%g):\n\
*** %d iterations > limit => limit increased to %d",epsp,scf.nit,scf.maxit))
      if (scf.maxit>scfmaxit0*4) ERROR(("too bad...")) }
    iret=0;
    goto ret; }

  /*** here eps>0: divergence detected ***/
  if (scf.nit>3) {
    if (maxerr>scf.maxerr) nbad++;
    else nbad=0; }

  if (nbad>2 || scf.nit>scf.maxit) {
    static double dmax[10];
    int nmax[10],imax[10],im;
    int maxnbad=5;

    if (epsp<-5) maxnbad=-epsp/2;
    prt("scf info: epsp=%g maxnbad=%d\n",epsp,maxnbad);

    loop (im,0,10) dmax[im]=0;
    loop (n,FROM,No.N) {
      mn=molec+n;
      ns=mn->ns;
      si=spec[mn->sp]->si;
      rpol=polarrof(mn,A->rp);
#define ND 10
      loop (i,0,ns) {
        double d=SQR(rpol[i]);

        if (d>dmax[ND-1]) loop (im,0,ND) if (d>dmax[im]) {
          dmax[im]=d;
          imax[im]=i;
          nmax[im]=n;
          break; } } }

    prt("mol  site  |rpol| (%d worst)",ND);
    loop (im,0,ND) prt("%3d %3d  %.9g",nmax[im],imax[im],sqrt(dmax[im]));
#undef ND

    WARNING(("selffield: %d iterations, %d consecutive bad\n\
*** (maxerr=%g, scf.maxit=%d)",scf.nit,nbad,maxerr,scf.maxit))

    if (scf.maxit<scfmaxit0*6 && nbad<maxnbad) {
      scf.maxit += 2+scf.maxit/3;
      prt("scf.maxit increased to %d",scf.maxit); }
    else ERROR(("selffield: too bad - stop"))

    if (!(option_v&8)) {
      static int pass;

      if (pass++) {
        prt("Bad for the second time - verbose output (-v8) set");
        option_v|=8; } } }

  if (epsp>0 ? maxerr<epsp : scf.nit>=fabs(epsp)) {
    scf.maxdr=sqrt(maxrr);
    goto ret; }

  iret=0;

 ret:

  scf.maxerr=maxerr; // used also in nmc,nm
  //  scf.stderr=err;

  //  fprintf(stderr,"t=%g scf.nit=%d maxerr=%g iret=%d\n",t,scf.nit,maxerr,iret);
  
  if (iret) {
    StaAdd("polar last iter stderr",En.pollaststderr=err);
    StaAdd("polar last iter maxerr",En.pollastmaxerr=maxerr); }

  return iret;
}

void mechpolar(ToIntPtr B,ToIntPtr A) /************************** mechpolar */
/*
  DEPRECATED (Drude subsystem not thermostated)
  mechanical model of polarizable dipoles: just rhs
  NOTE: this was created as a simplified copy of selffield(),
        to be rewritten so that selffield() calls mechpolar()...
*/
{
  molecule_t *mn;
  siteinfo_t *si;
  int ns,n,i;
  vector *fc,*fpol,dr;
  double rr,maxrrold=0;
#if defined(COULOMB) && COULOMB<0
  double f;
#endif /*# defined(COULOMB) && COULOMB<0 */
#ifndef FREEBC
  double byerd;
#endif /*# FREEBC */

#if POLAR&30
  ERROR(("mechpolar: POLAR=%d not supported",POLAR))
#endif /*# POLAR&30 */

  zeroEn();

  forces(B,A);

  loop (n,FROM,No.N) {
    mn=molec+n;
    ns=mn->ns;
    si=spec[mn->sp]->si;
    // ?   rpol=polarrof(mn,A->rp);
    fc=rof(mn,B->rp);
    fpol=polarrof(mn,B->rp);
    loop (i,0,ns) {
      // ?   VVV(dd,=dr,=rpol[i])
      rr=SQR(dr);
      Max(maxrrold,rr)
      if (rr>=Sqr(el.minqq))
        ERROR(("DRUDE (mechpolar): %d[%d] aux charge dist = %g",n,i,sqrt(rr)))

      /* NOTE (9/2001): the virial of force equals -2*self energy and
         is treated separately because of the Ewald corrections */

#if defined(COULOMB) && COULOMB<0
      /* Ewald correction for dipole self-term */
      f=si[i].charge*si[i].chargepol;
      En.el+=f*exacterud_sqrt(rr,&byerd);
      f*=byerd;
      VVO(fc[i],-=dr,*=f) VV(fpol[i],+=dr) /* note sign: dr=pol-center */
#endif /*# defined(COULOMB) && COULOMB<0 */

      /* forces passed to the central atom */
      VV(fc[i],+=fpol[i]) } }
}

int scforces(ToIntPtr B, ToIntPtr A) /***************************** scforces */
/*
  accurate self-consistent forces, accuracy control: el.epsx (see the manual)
  number of iterations is returned
*/
{
  /* iterate self-field */
  scf.nit=0;
  do {
    zeroEn();
    forces(B,A);
  } while (!selffield(B,A,scf.epsx,scf.omegax,0));

  En.pot += En.self + En.el;

  return scf.nit;
}

int scfautoset(int icyc,int noint) /***************************** scfautoset */
{
  if (scf.domega) {
    static int lastcycmaxit;
    
    prt("%d %.6f %.6f %.6f %d SCF AUTOSET: icyc omega maxerr stderr nit",
        icyc,scf.omega,En.sumpolmaxerr/noint,En.sumpolstderr/noint,scf.nit);
    
    if (icyc && scf.cycmaxit>lastcycmaxit) {
      prt("! %d SCF iterations (%d in the last cycle) => divergence detected",
          scf.cycmaxit,lastcycmaxit);
      return 1; }
    
    En.sumpolmaxerr=En.sumpolstderr=0;
    scf.omega+=scf.domega;
    lastcycmaxit=scf.cycmaxit;
    scf.cycmaxit=0; }

  return 0;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 POLAR ERRORS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

void testSCF(void) /************************************************ testSCF */
/*
  For testing quality of SCF integration: full iterations
  FREEBC also supported
  See also measureP() where some test are also performed (except FREEBC)
  WARNING: 1st evalation of forces unnecessarily calculated twice
*/
{
  En_t En0;
  FILE *frun=NULL,*fex=NULL,*ferr=NULL;
  int n,i;

  if (!scf.test) return;

  En0=En;

  box.V=PROD(box.L);
  box.cq=Sqr(box.cutoff);

  if (!locA) sdsalloc(locA,cfg[0]->size);
  if (!locB) sdsalloc(locB,cfg[0]->size);

  sdscopy(locA,cfg[0]);

  if (!(constrd.mode&RESCALE_CM)) depend_r(locA,1); /* meaningfull for Rowlinson only */

  scforces(locB,locA);
  depend_f(locA,locB); /* warning: param order ! */
  StaAdd("scf.nit",scf.nit);
  StaAdd("scf.maxdr",scf.maxdr);

  En=En0;
  En.Polmaxerr=En.Polstderr=0;
  En.Fstderr=En.Fmaxerr=0;

  if (option('v')&32) {
    char *fmode=init<2 || (init_append&2) ? "at" : "wt";

    frun=fopen(Fn("run.pol"),fmode);
    fex=fopen(Fn("ex.pol"),fmode);
    ferr=fopen(Fn("err.pol"),fmode);
    fprintf(frun,"\
#  t=%.4f  Running induced dipole moments in p.u. (%.11g D)\n\
#  dip_x              dip_y              dip_z               |dip|\n",t,Debye);
    fprintf(fex,"\
#  t=%.4f  Iterated induced dipole moments in p.u. (%.11g D)\n\
#  dip_x              dip_y              dip_z               |dip|\n",t,Debye);
    fprintf(ferr,"\
#  t=%.4f  Errors in induced dipole moments in p.u. (%.11g D)\n \
#  dip_x              dip_y              dip_z               |dip|\n",t,Debye); }

  loop (n,FROM,No.N) {
    molecule_t *mn=molec+n;
    int ns=mn->ns;
    siteinfo_t *si=spec[mn->sp]->si;
    vector *rpolRun=polarrof(mn,cfg[0]->rp);
    vector *rpolA=polarrof(mn,locA->rp);
    vector *fRun=rof(mn,cfg[2]->rp);
    vector *fB=rof(mn,locB->rp);
    vector dd;
    double f;

    loop (i,0,ns) if (si[i].chargepol) {

      VVV(dd,=rpolRun[i],-rpolA[i])
      f=SQR(dd)*si[i].qqpol;
      if (frun) fprintf(frun,"%18.12f %18.12f %18.12f  %18.12f\n",
                        rpolRun[i][0]*si[i].chargepol,
                        rpolRun[i][1]*si[i].chargepol,
                        rpolRun[i][2]*si[i].chargepol,
                        sqrt(SQR(rpolRun[i])*si[i].qqpol));
      if (fex) fprintf(fex,"%18.12f %18.12f %18.12f  %18.12f\n",
                       rpolA[i][0]*si[i].chargepol,
                       rpolA[i][1]*si[i].chargepol,
                       rpolA[i][2]*si[i].chargepol,
                       sqrt(SQR(rpolA[i])*si[i].qqpol));
      if (ferr) fprintf(ferr,"%18.12f %18.12f %18.12f  %18.12f\n",
                        dd[0]*si[i].chargepol,
                        dd[1]*si[i].chargepol,
                        dd[2]*si[i].chargepol,
                        sqrt(f));
      En.Polstderr+=f; Max(En.Polmaxerr,f)

      f=SQRD(fRun[i],fB[i]);
      En.Fstderr+=f; Max(En.Fmaxerr,f) } }

  if (frun) { fprintf(frun,"\n"); fclose(frun); }
  if (fex) { fprintf(fex,"\n"); fclose(fex); }
  if (ferr) { fprintf(ferr,"\n"); fclose(ferr); }

  // StaSet(0,lag.ierr,2,lag.in); ?
  StaSet(0,lag.err,2,lag.n);

  En.Polstderr=sqrt(En.Polstderr/No.pol);
  En.Polmaxerr=sqrt(En.Polmaxerr);
  StaAdd("Polar stderr",En.Polstderr);
  StaAdd("Polar maxerr",En.Polmaxerr);

  En.Fstderr=sqrt(En.Fstderr/No.pol);
  En.Fmaxerr=sqrt(En.Fmaxerr);
  StaAdd("Polar force stderr [N]",En.Fstderr*(massunit*lengthunit/Sqr(timeunit)));
  StaAdd("Polar force maxerr [N]",En.Fmaxerr*(massunit*lengthunit/Sqr(timeunit)));
} /* testSCF() */

void measureepspol(void) /************************************ measureepspol */
/*
   to determine eps^inf (optical permittivity)
   Ewald needed
*/
{
#  ifndef FREEBC
  //  static ToIntPtr locA=NULL,locB; // not local now
  vector M0,dM;
  int i,nit0;
  double E,aux;
  static int phase; /* x,-x,y,-y,z,-z */
  static int stride=0; /* not every cycle */

  /* save previous state */
  struct Eext_s Eext0;
  En_t En0;
  vector QM0;

  if (scf.Estride && stride++%scf.Estride) return;

  En0=En;
  Eext0=Eext;
  VV(QM0,=Q->M)

  box.cq=Sqr(box.cutoff);

  if (!locA) sdsalloc(locA,cfg[0]->size);
  if (!locB) sdsalloc(locB,cfg[0]->size);

  sdscopy(locA,cfg[0]);

  scforces(locB,locA);
  nit0=scf.nit;

  VVV(M0,=dM,=Q->M)

  E=Eext.Eepspol*(1-2*(phase&1)); /* in prog. units, w. sign */

  Eext.E[phase/2]+=E;
  Eext.isE=1;

  scforces(locB,locA);

  VVV(dM,=Q->M,-M0)
  aux=4*PI/PROD(box.L)/E;
  if (option('v')&4) prt("%g %g %g  %c E=%g p.u. it=%d,%d HIGHFREQ",
      dM[0]*aux,dM[1]*aux,dM[2]*aux,
      'x'+phase/2,
      E,nit0,scf.nit);

  StaSet(0,2,2,0);
  StaAdd("hf:scf.nit(0)",nit0);
  StaAdd("hf:scf.nit(E)",scf.nit);

  StaSet(0,2,2,lag.in);
  En.chi_hf=dM[phase/2]*aux;
  StaAdd("high-frequency chi",En.chi_hf);
  aux=4*PI/(3*T*PROD(box.L));
  StaAdd("total chi",En.chi_hf+aux*SQR(M0));

  phase=(phase+1)%6;
  Eext=Eext0;
  /* to fix total chi (original Eext) */
  if (Eext.isE && scf.Estride>1) {
    aux=sqrt(aux);
    putv(Eext.E)
    loop (i,0,3) if (Eext.E[i]) StaAdd(string("xchi_%c",'x'+i),aux*M0[i]); }
  En=En0;
  VV(Q->M,=QM0)
#  endif /*# FREEBC */
}
