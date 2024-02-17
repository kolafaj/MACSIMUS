/* restricted: no constraints, no dependants, no frame dump... */
void normalmodes(void) /**************************************** normalmodes */
{
  int N=No.s*3,n,i,j;
  double **A,**R=NULL,*rp,*mq,*fp,*fm;
  double *fpp=NULL,*fmm=NULL; /* initialized to suppress compiler warning */
  double rr,x,sumlnf;
  sort_t *sort;
  FILE *nmv;
  ToIntPtr CFG,FP,FM,FPP=NULL,FMM=NULL,M; /* initialized to suppress compiler warning */
  double starttime=mytime(),stoptime,mem;
  char *err;
#ifdef POLAR
  int nit,nitm;
  double summaxerr=0,summaxerrq=0,*p0,*pp;
#endif

  underline("fundamental frequencies");

#if defined(COULOMB) && COULOMB<0
  if ((el.rshift&2)==0) WARNING(("nm: el.rshift=3 recommended with normal mode calculations"))
#endif

  if (No.c) ERROR(("nm: No.c=%d: constraints not implemented\n\
*** try nm.method=1 or 3",No.c))
  if (No.depend[0]) ERROR(("nm: No.depend=%d: dependants not implemented\n\
*** try nm.method=1 or 3",No.depend))

#define D(X) ((double*)(X[0]))
  if (FROM) ERROR(("normalmodes and FROM=%d: not implemented"))
  err=string("nm: started at %s",myctime(mytime()));
  fprintf(stderr,"%s",err);
  prt_("%s",err);

  prt("normal mode analysis: %dx%d matrix of d2U/drdr, nm.dr=%g, No.c=%d\n\
%s order formula for numerical derivative of forces",N,N,nm.dr,No.c,
      nm.dr>0?"2nd":"4th");
  prt("nm.method=0: good for fully flexible models (no constraints)");

  if (equalize.cfg || equalize.mol)
    WARNING(("nm: masses of atoms have been equalized\n\
*** normal modes and frequences are affected"))

  mem=(sizeof(double))/1073741824.*(Sqr((double)N)*(1+!!nm.modes));
  prt("nm: memory needed by matrices (excl. O(N) terms) = %g GiB",mem);
  if (mem>nm.mem) ERROR(("%.4g GiB needed exceeds limit nm.mem = %.4g GiB",mem,nm.mem))

  ralloc2Darrayzero(A,N,N); /* this is the mark for release */
  sdsralloc(CFG,cfg[0]->size); sdscopy(CFG,cfg[0]); rp=D(rof(molec,CFG->rp));
  sdsralloczero(FM,cfg[0]->size); fm=D(rof(molec,FM->rp));
  sdsralloczero(FP,cfg[0]->size); fp=D(rof(molec,FP->rp));
  if (nm.dr<0) {
    sdsralloczero(FMM,cfg[0]->size); fmm=D(rof(molec,FMM->rp));
    sdsralloczero(FPP,cfg[0]->size); fpp=D(rof(molec,FPP->rp)); }
  sdsralloczero(M,cfg[0]->size);
  mq=D(rof(molec,M->rp));
#ifdef POLAR
  p0=D(polarrof(molec,cfg[0]->rp));
  pp=D(polarrof(molec,CFG->rp));
#endif

  if (nm.modes) ralloc2Darrayzero(R,N,N);

  j=0;
  loop (n,0,No.N) {
    molecule_t *mn=molec+n;
    siteinfo_t *si=spec[mn->sp]->si;
    int ns=mn->ns;

    loop (i,0,ns) {
      x=sqrt(si[i].mass);
      mq[j++]=x;
      mq[j++]=x;
      mq[j++]=x; } }

  if (j!=No.s*3) ERROR(("internal: j!=No.s*3"))

  measure=1;

  /* forces at minimum, Drude at minimum:
     - FP thrown out
     - polar part of cfg[0] used in the predictor
     - not needed for nonpolar  */
#ifdef POLAR
  nit=scforces(FP,cfg[0]);
  summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
  put2(En.pot,nit)
  nit=0;
#endif

  /* calculate matrix A = expansion of the potential at minimum */
  loop (i,0,N) {
    rr=rp[i];

    rp[i]=rr+nm.dr;
#ifdef POLAR
    /* original dipoles */
    loop (j,0,No.s*3) pp[j]=p0[j];
    nitm+=scforces(FP,CFG);
    summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
#else
    zeroEn();
    forces(FP,CFG);
#endif

    rp[i]=rr-nm.dr;
#ifdef POLAR
    /* predictor */
    loop (j,0,No.s*3) pp[j]=2*p0[j]-pp[j];
    nit+=scforces(FM,CFG);
    summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
#else
    zeroEn();
    forces(FM,CFG);
#endif

    rp[i]=rr;

    if (nm.dr>0)
      loop (j,0,N) {
        /* NB: symmetrization (of numerical errors) here */
        x=(fp[j]-fm[j])/(nm.dr*mq[i]*mq[j]);
        A[i][j]-=x;
        A[j][i]-=x; }

    else {
      /* 4th order formula */
      rp[i]=rr+nm.dr*2;
#ifdef POLAR
      loop (j,0,No.s*3) pp[j]=3*p0[j]-2*pp[j];
      nitm+=scforces(FPP,CFG);
      summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
#else
      zeroEn();
      forces(FPP,CFG);
#endif
      rp[i]=rr-2*nm.dr;
#ifdef POLAR
      /* predictor */
      loop (j,0,No.s*3) pp[j]=2*p0[j]-pp[j];
      nit+=scforces(FMM,CFG);
      summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
#else
      zeroEn();
      forces(FMM,CFG);
#endif
      rp[i]=rr;

      loop (j,0,N) {
        /* NB: symmetrization (of numerical errors) here */
        x=(-fpp[j]+8*(fp[j]-fm[j])+fmm[j])/(6*nm.dr*mq[i]*mq[j]);
        A[i][j]-=x;
        A[j][i]-=x; } }

    if (i && (i%100==0)) fprintf(stderr,"nm: A[%d] of %d finished at %s",i,N,myctime(mytime()));
  }

  /* ... A is four times the correct (1/2)d2U/dxdy/sqrt(mi mj) */
  fprintf(stderr,"nm: A finished at %s",myctime(mytime()));

#ifdef POLAR
    prt("SCF iteration summary, scf.epsx=%g:\n\
averaged # of iterations for r=r-dr: %g  for r=r+dr(predicted): %g\n\
<maxerr> = %g  stdev(maxerr) = %g",
        scf.epsx,(double)nitm/N,(double)nit/N,
        summaxerr/(2*N+1),
        sqrt((summaxerrq/(2*N+1)-Sqr(summaxerr/(2*N+1)))*(2*N+1)/(2*N)));
    if (scf.epsx==0) prt("\
Hint: scf.epsx=2*<maxerr>  or so may be more efficient than scf.epsx=0,\n\
      yet accurate enough");
#endif

  sig=0;
  rr=Jacobi(N,A,R,nm.eps);
  /* ... and now:
     A=diag(omega_i^2),
     R[i][j]/sqrt(mj) = i-th normal mode eigenvector */
  prt("Jacobi diagonalization error = %g",rr);
  fprintf(stderr,"nm: Jacobi finished at %s",myctime(mytime()));

  rallocarray(sort,N);
  loop (i,0,N) {
    sort[i].l=A[i][i];
    sort[i].i=i; }
  qsort(sort,N,sizeof(struct sort_s),sortnmf);

  nmv=fopen(Fn("nmf"),"wt");

  fprintf(nmv,"\
# normal mode frequencies\n\
# matrix n=%d   d2U/drdr=-df/dr, dr=2*%g   Jacobi diag.err.=%g\n\
#  i        [THz]          [1/cm]\n",N,nm.dr,rr);
  rr=1e24/Sqr(4*PI); /* 1e24 = energyunit/lengthunit^2/massunit */
  sumlnf=0;
  j=0;
  loop (i,0,N) {
    x=sort[i].l*rr;
    if (x<0) x=-sqrt(-x); else x=sqrt(x);
    if (i>=nm.zero) sumlnf+=log(x);
    else if (fabs(x)>1e10) j++;
    fprintf(nmv,"%4d %15.8f %15.8f%s\n",
            i,x*1e-12,x/29979245800,i<nm.zero?" excluded from Tnu":"");
  }

  if (j) WARNING(("dr.zero=%d, but %d frequencies > 0.01 THz found\n\
*** check the frequences in %s\n\
*** make sure a potential minimum has been reached (T=0, thermostat, tau.T)\n\
*** check nm.dr, nm.eps\n\
*** check accuracy of Ewald (el.epsr, el.grid) and SCF (scf.epsx)",nm.zero,j,lastFn))

  j=No.s*3-nm.zero;

  fprintf(nmv,"# geometric average of vibrational temperatures Tnu = %.13g K\n",
          PLANCK_BOLTZMANN*exp(sumlnf/j));

  fclose(nmv);

  prt("nm.zero=%d frequences omitted, %d nonzero fundamental frequencies:\n\
sum ln(nu/s^-1) = %.14g\n\
geometric average of vibrational temperatures Tnu = %.13g K\n",
      nm.zero,j,sumlnf,
      PLANCK_BOLTZMANN*exp(sumlnf/j));

  /* write vibrations */
  if (R) {
    int from=0,to=N,k,indx,ix,j,i;
    FILE *plb;
    float r[3];
    double m,a,ampl;
    char *sys=string("rm %s.nm[0-9][0-9][0-9][0-9].plb",simils.simname);

    /* remove old vibration plb's */
    prt("%s: retc=%d",sys,system(sys));

    if (nm.modes>0) to=nm.modes;
    else from=N+nm.modes;
    Max(from,0)
    Min(to,N)
    if (from<to) {
      prt("normal modes: writing frames from %d to %d (incl)",from,to-1);
      loop (j,from,to) {
        indx=sort[j].i;
        plb=fopen(Fn(string("nm%04d.plb",j)),"wb");
        if (!plb) { ERROR(("cannot open %s",lastFn)) break; }
        ampl=nm.ampl;
        /* normal modes:
           ampl = max amplitude of moves over all atoms */
        if (ampl<0) ampl=-ampl*sqrt(No.s);
        m=a=0;
        loop (i,0,N) {
          R[indx][i]/=mq[i];
          a+=Sqr(R[indx][i]);
          if (i%3==2) { Max(m,a) a=0; } }
        m=nm.ampl/sqrt(m);
        r[0]=No.s;
        r[1]=-3;
        fwrite(r,4,2,plb);
        loop (ix,0,abs(nm.frames)) {
          double x=ix/(double)(abs(nm.frames)-1);

          if (nm.frames>0) x=m*cos(x*PI);
          else x=m*(1-2*x);

#ifdef FREEBC
          VO(r,=0)
#else
          VV(r,=box.L)
#endif
          fwrite(r,4,3,plb);
          loop (i,0,No.s) {
            loop (k,0,3) r[k]=CFG->rp[i][k]+R[indx][i*3+k]*x;
            fwrite(r,4,3,plb); } }
        fclose(plb); }
    prt("normal modes: %s.nm%04d.plb--%s.nm%04d.plb written",simils.simname,from,simils.simname,to-1);
    prt("recommended command to show normal modes:\n\
  show -h6 -g500x400 -d30 %s %s.nm%%04d.plb -I\\$wwir",simils.simname,simils.simname);
    }
  } /* **R */

  release(A);
  stoptime=mytime();
  err=string("nmc: duration = %.6g s, finished at %s\n",
             stoptime-starttime,myctime(stoptime));
  fprintf(stderr,"%s",err);
  prt_("%s",err);
}
