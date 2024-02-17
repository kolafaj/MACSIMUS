/* 
   6/2018 general version
   5/2012 version without B4 (rigid models only)
*/
void normalmodesc(void) /**************************************** normalmodes */
{
  int N,Ni,n,i,j,k,kk,inum,sp,nc,ii,jj,ir,a,b,ic,nitm,nit=0;
  int aoffset;
  int *S2M,*M2S; /* renumbering: see below */
  double **A,**B,**B1,**B2,**B3,**B4=NULL,**P,**Ai,**Ai0,**R=NULL;
  double *rp,*fp,*fm,*fpp=NULL,*fmm=NULL; /* initialized to suppress compiler warning */
  double **MM,**gm=NULL,***H=NULL; /* initialized to suppress compiler warning */
  double *Lambda,*G,*imass,*aux=NULL; /* initialized to suppress compiler warning */
  double rr,x,sumlnf,finerr=-1;
  double det,norm=0,scal;
  double octaveim[2]={0,0};
  int noctaveim=0;
  sort_t *sort;
  FILE *nmv;
  ToIntPtr CFG,FP,FM,FPP=NULL,FMM=NULL,F0; /* initialized to suppress compiler warning */
  molecule_t *mn;
  siteinfo_t *si;
  vector *r,*ra;
  M_t *M;
  double starttime=mytime(),stoptime,mem;
  char *err;
#ifdef POLAR
  double summaxerr=0,summaxerrq=0,*p0,*pp;
#endif

  if (No.c==0) WARNING(("no constraint found, it is better to use nm.method=0"))

  underline("fundamental frequencies");

#if defined(COULOMB) && COULOMB<0
  if ((el.rshift&2)==0) WARNING(("nm: el.rshift=3 recommended with normal mode calculations"))
#endif

  /* matrix M for constrained dynamics */
  if (!Ms) ERROR(("no constraints (rigid bonds) found\n\
*** (no matrix Ms has been allocated)\n\
*** try nm.method=0"))

#define D(X) ((double*)(X[0]))
  if (FROM) ERROR(("normalmodes and FROM=%d: not implemented"))
  err=string("nmc: started at %s",myctime(mytime()));
  fprintf(stderr,"%s",err);
  prt_("%s",err);

  if (equalize.mol || equalize.cfg)
    WARNING(("nmc: masses of atoms have been equalized\n\
*** normal modes and frequences are affected"))

  N=(No.s-No.depend[0])*DIM; /* vector size, 3* number of sites */
  Ni=N-No.c;
  mem=(sizeof(double))/1073741824.*(Sqr((double)N)*(6+!!nm.modes)+Sqr((double)No.c)+2*Sqr((double)Ni));
  prt("nmc: memory needed by matrices (excl. O(N) terms) = %g GiB",mem);
  if (mem>nm.mem) ERROR(("%.4g GiB needed exceeds limit nm.mem = %.4g GiB",mem,nm.mem))

  ralloc2Darrayzero(A,N,N); /* this is the mark for release */
  ralloc2Darrayzero(Ai,Ni,Ni);
  ralloc2Darrayzero(Ai0,Ni,Ni);
  rallocarray(ra,No.maxc);
  rallocarrayzero(imass,No.s-No.depend[0]);
  ralloc2Darrayzero(B,N,N);
  ralloc2Darrayzero(B1,N,N);
  ralloc2Darrayzero(B2,N,N);
  ralloc2Darrayzero(B3,N,N);
  ralloc2Darrayzero(P,N,N); /* 1st index filled to N, note different order vs. gam */
  sdsralloc(CFG,cfg[0]->size); sdscopy(CFG,cfg[0]); rp=D(rof(molec,CFG->rp));
  sdsralloczero(FM,cfg[0]->size); fm=D(rof(molec,FM->rp));
  sdsralloczero(FP,cfg[0]->size); fp=D(rof(molec,FP->rp));
  if (nm.dr<0) {
    sdsralloczero(FMM,cfg[0]->size); fmm=D(rof(molec,FMM->rp)); 
    sdsralloczero(FPP,cfg[0]->size); fpp=D(rof(molec,FPP->rp)); }
  sdsralloczero(F0,cfg[0]->size);

#ifdef POLAR
  p0=D(polarrof(molec,cfg[0]->rp));
  pp=D(polarrof(molec,CFG->rp));
#endif

  ralloc2Darrayzero(MM,No.c,No.c);
  rallocarrayzero(G,No.c);
  rallocarrayzero(Lambda,No.c); /* = g of paper, -g in constrdyn (?sign) */
  if (nm.method&2) {
    ralloc2Darrayzero(B4,N,N);
    ralloc2Darrayzero(gm,N,No.c);
    rallocarray(aux,No.c);
    rallocarray(H,N);
    loop (i,0,N) ralloc2Darrayzero(H[i],No.c,No.c); }

  prt("normal mode analysis: %dx%d matrix of d2U/drdr, nm.dr=%g, No.c=%d\n\
%s order formula for numerical derivative of forces",N,N,nm.dr,No.c,
      nm.dr>0?"2nd":"4th");
  prt("nm.method=%d: %s", nm.method,
      nm.method&2 ? nm.method&1 ? "inefficient fool-proof general code, use nm.method=2" : "general" : "good for rigid bodies");

  /* S2M[i] is the consecutive number of a mass site (not dependant),
     given site number i in the configuration (numbered from 0)
     The corresponding indices of matrix A are:
       S2M[i]*3, S2M[i]*3+1, S2M[i]*3+2 */
  rallocarray(S2M,No.s);
  loop (i,0,No.s) S2M[i]=-1;
  /* M2S[j] is the site number of M2S[j]-th mass site */
  rallocarray(M2S,No.s-No.depend[0]);

  j=inum=0;
  loop (n,FROM,No.N) {
    mn=molec+n;
    sp=mn->sp;
    r=rof(mn,cfg[0]->rp);
    si=spec[sp]->si;
    nc=mn->nc;
    loop (i,0,mn->ns) {
      if (si[i].mass) {
        if (inum>=No.s-No.depend[0]) ERROR(("internal: number of dependants"))
        M2S[inum]=j; /* inum = consecutive mass site number */
        S2M[j]=inum;
        imass[inum]=si[i].imass; /* numbered by mass sites */
        inum++; }
      j++; } }

  if (j!=No.s || inum!=No.s-No.depend[0])
    ERROR(("internal: sites: %d != %d, mass sites: %d != %d",
           j,No.s,inum,No.s-No.depend[0]))

  measure=1;

  /* forces at minimum, Drude at minimum:
     - F0 used in calc. constraints,
     - polar part of cfg[0] used in the predictor */
  depend_r(cfg[0],1);
#ifdef POLAR
  nit=scforces(F0,cfg[0]);
  summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
#else
  zeroEn();
  forces(F0,cfg[0]);
#endif
  depend_f(cfg[0],F0);

  put2(En.pot,nit)

  forDUMP=CFG;

  nit=nitm=0;
  /* calculate matrix A = expansion of the potential at minimum */
  loop (jj,0,N) {
    /* 2nd order central formula */
    j=M2S[jj/3]*3+(jj%3); /* j=jj if no dependants */
    //    if (jj!=j) ERROR(("temp.dep. jj=%d j=%d",jj,j))
    rr=rp[j];

    rp[j]=rr+nm.dr;
    depend_r(CFG,1);
#ifdef POLAR
    /* original dipoles */
    loop (i,0,No.s*3) pp[i]=p0[i];
    nitm+=scforces(FP,CFG);
    summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
#else
    zeroEn();
    forces(FP,CFG);
#endif
    depend_f(CFG,FP);

    rp[j]=rr-nm.dr;
    depend_r(CFG,1);
#ifdef POLAR
    /* linear predictor */
    loop (i,0,No.s*3) pp[i]=2*p0[i]-pp[i];
    nit+=scforces(FM,CFG);
    summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
#else
    zeroEn();
    forces(FM,CFG);
#endif
    depend_f(CFG,FM);

    rp[j]=rr;

    if (nm.dr>0) 
      loop (ii,0,N) {
        i=M2S[ii/3]*3+(ii%3);
        A[ii][jj]=-(fp[i]-fm[i])/(2*nm.dr); }

    else {
      /* 4th order formula */
      rp[j]=rr+nm.dr*2;
      depend_r(CFG,1);
#ifdef POLAR
      loop (i,0,No.s*3) pp[i]=3*p0[i]-2*pp[i];
      nitm+=scforces(FPP,CFG);
      summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
#else
      zeroEn();
      forces(FPP,CFG);
#endif
      depend_f(CFG,FPP);

      rp[j]=rr-2*nm.dr;
      depend_r(CFG,1);
#ifdef POLAR
      /* predictor */
      loop (i,0,No.s*3) pp[i]=2*p0[i]-pp[i];
      nit+=scforces(FMM,CFG);
      summaxerr+=scf.maxerr; summaxerrq+=Sqr(scf.maxerr);
#else
      zeroEn();
      forces(FMM,CFG);
#endif
      depend_f(CFG,FMM);
      rp[j]=rr; 

      loop (ii,0,N) {
        i=M2S[ii/3]*3+(ii%3);
        A[ii][jj]=-(-fpp[i]+8*(fp[i]-fm[i])+fmm[i])/(12*nm.dr); } }

    if (jj && (jj%100==0)) {
      err=string("nmc: A[%d] of %d finished at %s",jj,N,myctime(mytime()));
      fprintf(stderr,"%s",err);
      prt_("%s",err); } 
  } /* A: n */

  err=string("nmc: A finished at %s",myctime(mytime()));
  fprintf(stderr,"%s",err);
  prt_("%s",err);

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
  if (nm.dr<0) prt("4th order formula: r=r+dr includes r=r+2dr, r=r-dr includes r=r-2dr-dr");
#endif

  inum=aoffset=0;
  loop (n,FROM,No.N) {
    mn=molec+n;
    sp=mn->sp;
    si=spec[sp]->si;
    r=rof(mn,CFG->rp);
    ir=r-CFG->rp;

    if (ir!=inum) ERROR(("%d %d internal",ir,inum))

    if ( (nc=mn->nc) ) {

      /* because of B4, we need all ra[] in advance */
      loop (a,0,nc) {
        i=si[a].pair[0]; j=si[a].pair[1];
        VVV(ra[a],=r[i],-r[j]) }

      M=Ms[option('c')&1 ? n : sp];
      if (!M) ERROR(("internal: no matrix M"))

      loop (a,0,nc) {
        i=si[a].pair[0]; j=si[a].pair[1];

        /* ii,jj=numbers of mass sites of bond a of molecule n */
        ii=S2M[i+inum]; jj=S2M[j+inum];

        if (ii<0 || jj<0) ERROR(("constraint to a dependant not supported"))
        loop (k,0,3) {
          P[a+aoffset][ii*3+k]=ra[a][k];
          P[a+aoffset][jj*3+k]=-ra[a][k]; }

        loop (k,0,3) G[aoffset+a]+=(F0->rp[inum+i][k]*si[i].imass-F0->rp[inum+j][k]*si[j].imass)*ra[a][k];

        /* copying sparse matrix M to the square form MM */
        //        MM[a+aoffset][a+aoffset] = (si[i].imass+si[j].imass)*SQR(ra[a]); /* a==b */
        MM[a+aoffset][a+aoffset] = M->imass*SQR(ra[a]); /* a==b */
        while ( (b=(++M)->b) < a )
          MM[a+aoffset][b+aoffset] =
          MM[b+aoffset][a+aoffset] = M->imass*SCAL(ra[a],ra[b]);
        } /* a */

      if (nm.method&2) {
        /* general constraints */
        if (nm.method&1) {
          /* inefficient fool-proof version of H = grad_j M_ab */
          int K,J;

          M=Ms[option('c')&1 ? n : sp];
          loop (a,0,nc) {
            i=si[a].pair[0]; j=si[a].pair[1];
            ii=S2M[i+inum]; jj=S2M[j+inum];
            loop (b,0,nc) {
              int ib=si[b].pair[0], jb=si[b].pair[1];
              int iib=S2M[ib+inum], jjb=S2M[jb+inum];

              loop (J,0,N/3) 
                loop (K,0,N/3) {
                  loop (k,0,3) {
                    H[3*J+k][a+aoffset][b+aoffset] += imass[K]*
                      ((K==ii)-(K==jj))*((K==iib)-(K==jjb))
                      *(
                        ((J==ii)-(J==jj))*ra[b][k]
                        +
                        ((J==iib)-(J==jjb))*ra[a][k]
                       );
                } } } } }
        else {
          /* H = grad_j M_ab */
          M=Ms[option('c')&1 ? n : sp];
          loop (a,0,nc) {
            i=si[a].pair[0]; j=si[a].pair[1];
            ii=S2M[i+inum]; jj=S2M[j+inum];

            loop (k,0,3) {
              H[ii*3+k][a+aoffset][a+aoffset] += 2*M->imass*ra[a][k];
              H[jj*3+k][a+aoffset][a+aoffset] -= 2*M->imass*ra[a][k]; }

            while ( (b=(++M)->b) < a ) {
              int ib=si[b].pair[0], jb=si[b].pair[1];
              int iib=S2M[ib+inum], jjb=S2M[jb+inum];

              /* inefficient design: 
                 I do not known which site is common -
                 this info could have been prepared in M-> 
                 ic = common site of bonds a,b
                 xa = 2nd site of bond a
                 xb = 2nd site of bond a
                 ic,xa,xb are in global site numbering, [0,N/3)
                 sga = sign of ra[a] wrt site ic
                 sgb = sign of ra[b] wrt site ic */

              int ic=0,xa=0,xb=0,sga=0,sgb=0; /* except sga initialized to suppress compiler warning */

              if (ii==iib) ic=ii,xa=jj,xb=jjb,sga=1,sgb=1;
              if (ii==jjb) ic=ii,xa=jj,xb=iib,sga=-1,sgb=1;
              if (jj==iib) ic=jj,xa=ii,xb=jjb,sga=1,sgb=-1;
              if (jj==jjb) ic=jj,xa=ii,xb=iib,sga=-1,sgb=-1;

              if (!sga)
                ERROR(("internal: bad M: bonds a=%d-%d b=%d-%d do not have a common site",ii,jj, iib,jjb))

              loop (k,0,3) {
                H[ic*3+k][a+aoffset][b+aoffset] += M->imass*(sga*ra[a][k]+sgb*ra[b][k]);
                H[xa*3+k][a+aoffset][b+aoffset] -= M->imass*sgb*ra[b][k];
                H[xb*3+k][a+aoffset][b+aoffset] -= M->imass*sga*ra[a][k]; } }
          } /* a */
        } 
      } /* method&2 */

    aoffset+=nc; } /* nc */
  inum+=mn->ns; } /* n */

  if (aoffset!=No.c) ERROR(("internal: No.c"))
  if (inum!=No.s) ERROR(("internal: No.s"))

  /* Gram-Schmidt orthonormalization of the space of constraints
     and the rest (random vectors added) */
  inum=-N; /* counts repeated random vectors */
  loop (ic,0,N) {
    /* ... vectors P[0] .. P[ic-1] are already orthonormal */
    do {
      inum++;
      if (ic>=No.c) {
        /* random vector */
        norm=0;
        loop (i,0,N) {
          P[ic][i]=rndcos();
          norm+=Sqr(P[ic][i]); } }
      loop (i,0,ic) {
        scal=0; loop (j,0,N) scal+=P[ic][j]*P[i][j];
        loop (j,0,N) P[ic][j]-=scal*P[i][j]; }
      scal=0; loop (j,0,N) scal+=Sqr(P[ic][j]);
    } while (ic>=No.c && scal*N<norm*(N-ic));
    /* repeated if not orthogonal enough */

    /* normalize */
    scal=sqrt(scal); loop (j,0,N) P[ic][j]/=scal; }


  /* determining the error of the basis */
  norm=0;
  loop (i,0,N) {
    loop (j,0,N) {
      scal=-(i==j);
      loop (k,0,N) scal+=P[i][k]*P[j][k];
      Max(norm,fabs(scal)) } }

  prt("%d additional vectors in Gram-Schmidt for better numeric stability",inum);
  prt("max error of the basis orthogonal to the constraints=%g",norm);
  /* matrix P[N][N]
     line 0:      1st constraint, normalized (mass sites only)
     line 1:      2nd constraint, orthonormalized
     ...
     line No.c-1: last constraint, orthonormalized
     line No.c:   1st perpendicular vector, normalized         \
     ...                                                        } Ni vectors
     line N-1:    last perpendicular vector, orthonormalized   /

     S2M[No.s] : range=No.s-No.depend[0]=N/3
     M2S[N/3] : range=No.s
  */

  /* inverting the constraint dynamics matrix */
  det=invmat(No.c,MM);

  prt("nmc: invert constraint dynamics matrix (No.c=%d): det=%g%s",
      No.c,det,det==0?" (underflow)":"");

  /* Lagrange multipliers at minimum: the same as -g[] from constraint dynamics */
  loop (i,0,No.c)
    loop (j,0,No.c) Lambda[i]-=MM[i][j]*G[j];

  //  matdotvec(No.c,Lambda,MM,G); loop (i,0,No.c) Lambda[i]*=-1;

  if (nm.method&2)
    loop (j,0,N) {
      if ((nm.method&1)==0) /* symmetrization of H */
        loop (a,0,No.c)
          loop (b,0,a) H[j][a][b]=H[j][b][a]=x=H[j][a][b]+H[j][b][a];
      matdotvec(No.c,aux,H[j],Lambda);
      matdotvec(No.c,gm[j],MM,aux); }

  //  if (nm.method&2)  loop (j,0,N) loop (a,0,No.c) if (gm[j][a]) prt("gm[%d][%d]=%g",j,a,gm[j][a]);

  err=string("nmc: M,P,Lambda(,gm) finished at %s",myctime(mytime()));
  fprintf(stderr,"%s",err);
  prt_("%s",err);

  /* constraint part of A */
  inum=aoffset=0;
  r=rof(molec,rp);
  loop (n,FROM,No.N) {
    mn=molec+n;
    sp=mn->sp;
    si=spec[sp]->si;

    ir=rof(mn,rp)-CFG->rp;
    if (ir!=inum) ERROR(("%d %d internal",ir,inum))

    if ( (nc=mn->nc) ) {

      M=Ms[option('c')&1 ? n : sp];
      if (!M) ERROR((""))

      loop (a,0,nc) {
        int ia=si[a].pair[0];
        int ja=si[a].pair[1];

        VVV(ra[a],=r[inum+ia],-r[inum+ja]) }

      loop (a,0,nc) {
        int ia=si[a].pair[0];
        int ja=si[a].pair[1];
        int iA=S2M[inum+ia]*3;
        int jA=S2M[inum+ja]*3;

        loop (k,0,3) {
          B3[iA+k][iA+k] -= si[ia].imass*Lambda[a+aoffset];
          B3[iA+k][jA+k] += si[ia].imass*Lambda[a+aoffset];
          B3[jA+k][iA+k] += si[ja].imass*Lambda[a+aoffset];
          B3[jA+k][jA+k] -= si[ja].imass*Lambda[a+aoffset]; }

        loop (b,0,nc) {
          int ib=si[b].pair[0];
          int jb=si[b].pair[1];
          int iB=S2M[inum+ib]*3;
          int jB=S2M[inum+jb]*3;

          /* iM=M^-1[a][b] */
          double iM=MM[a+aoffset][b+aoffset];

          loop (kk,0,3) {
            double dfm=iM*(F0->rp[inum+ib][kk]*si[ib].imass
                          -F0->rp[inum+jb][kk]*si[jb].imass);

            loop (k,0,3) {
              B1[iA+k][iB+kk] += si[ia].imass*ra[a][k]*dfm;
              B1[iA+k][jB+kk] -= si[ia].imass*ra[a][k]*dfm;
              B1[jA+k][iB+kk] -= si[ja].imass*ra[a][k]*dfm;
              B1[jA+k][jB+kk] += si[ja].imass*ra[a][k]*dfm; } }

          loop (j,0,N) {
            /* to optimize... */
            double raAiB=ra[b][0]*A[iB][j]+ra[b][1]*A[iB+1][j]+ra[b][2]*A[iB+2][j];
            double raAjB=ra[b][0]*A[jB][j]+ra[b][1]*A[jB+1][j]+ra[b][2]*A[jB+2][j];
            double rax=iM*(raAiB*si[ib].imass-raAjB*si[jb].imass);

            loop (k,0,3) {
              B2[iA+k][j] -= ra[a][k]*rax*si[ia].imass;
              B2[jA+k][j] += ra[a][k]*rax*si[ja].imass; } }

        } /* b */

        /* B4 */
        if (nm.method&2)
          loop (j,0,N) {
            loop (k,0,3) {
              B4[iA+k][j] += ra[a][k]*gm[j][a]*si[ia].imass;
              B4[jA+k][j] -= ra[a][k]*gm[j][a]*si[ja].imass; } }

      } /* a */

      aoffset+=nc; }
    inum+=mn->ns; }

  if (aoffset!=No.c) ERROR((""))
  if (inum!=No.s) ERROR((""))


  /* A[i][j]*imass[i/3] is called B_0 in the manual */
  if (nm.eps<0) {
    double norm0=0,norm1=0,norm2=0,norm3=0,norm4=0;

    loop (i,0,N)
      loop (j,0,N) {
        norm0+=Sqr(A[i][j]*imass[i/3]);
        norm1+=Sqr(B1[i][j]);
        norm2+=Sqr(B2[i][j]);
        norm3+=Sqr(B3[i][j]);
        B[i][j] = A[i][j]*imass[i/3]+B1[i][j]+B2[i][j]+B3[i][j];
        if (nm.method&2) {
          norm4+=Sqr(B4[i][j]);
          B[i][j] += B4[i][j]; } }
    prt("norms (see the manual, not projected):\n\
  |B0|=%g, |B1|=%g, |B2|=%g, |B3|=%g, |B4|=%g",
        sqrt(norm0),sqrt(norm1),sqrt(norm2),sqrt(norm3),sqrt(norm4)); }
  else
    loop (i,0,N)
      loop (j,0,N) {
        B[i][j] = A[i][j]*imass[i/3]+B1[i][j]+B2[i][j]+B3[i][j];
        if (nm.method&2) B[i][j] += B4[i][j]; }

#if 0
    finerr=Jacobins(N,B,R,nm.eps,nm.key);
    prt("Jacobi diagonalization error = %g",finerr);
    err=string("nmc: Jacobi(Schur) finished at %s",myctime(mytime()));
    fprintf(stderr,"%s",err);
    prt_("%s",err); 
    loop (i,0,N) prt("%d %g",i,B[i][i]);
    exit(0);
#endif

  //  loop (i,0,N) loop (j,0,N) fprintf(stderr,"%8.3g%c",B4[i][j],j==N-1?'\n':' ');

  /*
    filter out (nonphysical) motions along the constraints
    Note: in old+misc/nmc-finaldebug.c you can find a version with
          these motions available (via #if)
  */

  /* projection first */
  loop (i,0,N) loop (j,0,N) B1[i][j]=B2[i][j]=0;

  loop (i,0,N)
    loop (j,0,N)
      loop (k,0,N) B1[i][j]+=P[i][k]*B[k][j];

  loop (i,0,N)
    loop (j,0,N)
      loop (k,0,N) B2[i][j]+=B1[i][k]*P[j][k];

  /* select the projected matrix (above some parts calculated unnecessarily) */
  loop (i,0,Ni)
    loop (j,0,Ni) Ai0[i][j]=Ai[i][j]=B2[No.c+i][No.c+j];

  if (nm.modes) {
    if (nm.method&8) {
      WARNING(("nmc: eigenvectors not supported (yet) with octave, nm.modes ignored"))
      nm.modes=0; }
    else
      ralloc2Darrayzero(R,Ni,Ni); }

  err=string("nmc: B0..B4 finished at %s",myctime(mytime()));
  fprintf(stderr,"%s",err);
  prt_("%s",err);

  sig=0;

  if (nm.method&8) { /* octave */
    char *octave_in=Fn("octave.in");
    char *octave_out=Fn("octave.out");
    FILE *f=fopen(octave_in,"wt");
    char line[256],*tok;
    int i,j,sg;
    double lambda[2]={0,0}; /* initialized to suppress compiler warning */

    fprintf(f,"\
format long\n\
N=%d\n\
A=[",Ni);
    loop (i,0,Ni)
      loop (j,0,Ni)
        fprintf(f,"%.15g%s",Ai[i][j],
                j<Ni-1 ? ",\\\n" : i<Ni-1 ? "\n" : "];\n");
    fprintf(f,"\nE=eig(A)\n");
    fclose(f);
    tok=string("octave < %s > %s",octave_in,octave_out);
    fprintf(stderr,"executing %s ...\n",tok);
    i=system(tok);

    err=string("nmc: octave:eig(A) errcode=%d finished at %s",i,myctime(mytime()));
    fprintf(stderr,"%s",err);
    prt_("%s",err);

    f=fopen(octave_out,"rt");
    if (!f) ERROR(("nmc: no %s",octave_out))

    while (fgets(line,256,f)) {
      if (!memcmp("E =",line,3)) break; }

    if (!fgets(line,256,f)) ERROR(("nmc: read %s: unexpected EOF",octave_out));
    loop (i,0,Ni) {
      if (!fgets(line,256,f)) ERROR(("nmc: read %s: unexpected EOF",octave_out))
      tok=strtok(line," \t\n");
      lambda[0]=lambda[1]=0;
      loop (j,0,2) {
        sg=1;
        if (!strcmp("+",tok)) tok=strtok(NULL," \t\n");
        if (!strcmp("-",tok)) tok=strtok(NULL," \t\n"),sg=-1;
        if (!tok) break;
        lambda[j]=sg*atof(tok);
        tok=strtok(NULL," \t\n");
        if (!tok) break; }

      if (lambda[1]) {
        Ai[i][i]=sqrt(Sqr(lambda[0])+Sqr(lambda[1]));
        if (fabs(lambda[1])>fabs(octaveim[1])) {
          noctaveim++;
          octaveim[0]=lambda[0];
          octaveim[1]=lambda[1]; } }
      else
        Ai[i][i]=lambda[0]; }

    finerr=fabs(lambda[1]);

    if (noctaveim) WARNING(("%d complex eigenvalues detected\n\
*** the eigenvalue with max. |Im| is %g + %g i\n\
*** check octave output file %s\n",
                            noctaveim,octaveim[0],octaveim[1],lastFn))
    fclose(f); }

  else { /* extended Jacobi */

    if (nm.dr<0) { 
      int i,j;
      double sum=0;

      loop (j,0,Ni) 
        loop (i,0,Ni) sum+=fabs(Ai[i][j]-Ai[j][i]);
      prt("nmc: matrix antisymmetry index=%g",sum/Sqr(Ni));
    }
                         
    finerr=Jacobins(Ni,Ai,R,nm.eps,nm.key);
    prt("Jacobi diagonalization error = %g",finerr);
    err=string("nmc: Jacobi(Schur) finished at %s",myctime(mytime()));
    fprintf(stderr,"%s",err);
    prt_("%s",err); }

  //    loop (i,0,Ni) prt("%d %g",i,Ai[i][i]);

  if (R) {
    JacobiEigenvectors(Ni,Ai0,Ai,R);
    err=string("nmc: JacobiEigenvectors finished at %s",myctime(mytime()));
    fprintf(stderr,"%s",err);
    prt_("%s",err); }

  if (R) {
    double *RR,s=0,x;

    allocarray(RR,Ni);

    loop (i,1,Ni) if (Ai[i][i]>1.001*Ai[i-1][i-1]) {
      x=JacobiTest(Ni,Ai[i][i],Ai0,R[i]);
      Max(s,x); }
    prt("Eigenvector test of nondegenerate eigenvalues (ratio>1.001) = %g (should be small)",s);
  }

  rallocarray(sort,Ni);
  loop (i,0,Ni) {
    sort[i].l=Ai[i][i];
    sort[i].i=i; }
  qsort(sort,Ni,sizeof(struct sort_s),sortnmf);

  nmv=fopen(Fn("nmf"),"wt");

  fprintf(nmv,"\
# normal mode frequencies\n\
# matrix n=%d   d2U/drdr=-df/dr, dr=2*%g   err=%g\n\
#  i        [THz]          [1/cm]\n",Ni,nm.dr,finerr);
  rr=1e24/Sqr(2*PI); /* 1e24 = energyunit/lengthunit^2/massunit */
  sumlnf=0;
  j=0;
  loop (i,0,Ni) {
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

  j=Ni-nm.zero;

  fprintf(nmv,"# geometric average of vibrational temperatures Tnu = %.12g K\n",
          PLANCK_BOLTZMANN*exp(sumlnf/j));

  fclose(nmv);

  prt("nm.zero=%d frequences omitted, %d nonzero fundamental frequencies:\n\
sum ln(nu/s^-1) = %.13g\n\
geometric average of vibrational temperatures Tnu = %.12g K\n",
      nm.zero,j,sumlnf,
      PLANCK_BOLTZMANN*exp(sumlnf/j));

  /* write vibrations */
  if (R) {
    int from=0,to=Ni,k,indx,ix,j,i;
    FILE *plb;
    float r[3];
    double m,a,ampl;
    double *RP;
    char *sys=string("rm %s.nm[0-9][0-9][0-9][0-9].plb",simils.simname);

    rallocarrayzero(RP,N);

    /* remove old vibration plb's */
    prt("%s: retc=%d",sys,system(sys));

    if (nm.modes>0) to=nm.modes;
    else from=N+nm.modes;
    Max(from,0)
    Min(to,Ni)
    if (from<to) {
      prt("normal modes: writing frames from %d to %d (incl)",from,to-1);
      loop (j,from,to) {
        indx=sort[j].i;
        loop (i,0,N) {
          RP[i]=0;
          loop (k,0,Ni) RP[i]+=R[indx][k]*P[No.c+k][i]; }
        plb=fopen(Fn(string("nm%04d.plb",j)),"wb");
        if (!plb) { ERROR(("cannot open %s",lastFn)) break; }
        ampl=nm.ampl;
        /* normal modes:
           ampl = max amplitude of moves over all atoms */
        if (ampl<0) ampl=-ampl*sqrt(No.s);
        m=a=0;
        loop (i,0,N) {
          a+=Sqr(RP[i]);
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
          loop (ii,0,N) {
            i=M2S[ii/3]*3+(ii%3);
            fm[i]=rp[i]+RP[ii]*x; }
          depend_r(FM,1);
          loop (i,0,No.s) {
            loop (k,0,3) r[k]=FM->rp[i][k];
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
