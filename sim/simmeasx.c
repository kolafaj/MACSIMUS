/*
  MORE MEASUREMENTS
  * Structure factor (only from playback, option -m1)
    The structure factor is weighted by `atomic masses': edit the
    ble-file to use `scattering amplitudes' b_i instead!
  * Mean square displacement (to determine diffusivity), see also plb2diff.c
  * Mean square charge diffusivity (to determine conductivity)
  * Cluster analysis (modules #included), #define CLUSTERS required
  * Cross section, #define XSECTION required

  See old+misc/sfdx for attempts to calculate viscosity from velocities
  Basic measurements are in simmeas.c
*/

#include "ground.h"
#include "sds.h"
#include <time.h>
#include "simglob.h"
#include "ewald.h"
#include "simils.h"
#include "simmeasx.h"
#include "cpmark.h"
#include "units.h"
#include "maxjump.h"
#include "statics.h"
#include "norm.h" /* int iscube(void) */
#include "asksig.h"
#include "simdef.h"

#define CMSHIFTERR 3e-7

#ifdef CLUSTERS
#  include "intermac.h"
#  include "cluster.c" /* updated and fixed 9/2009 */
extern cl_t cl;
#endif /*# CLUSTERS */

static vector oldL,lastL; /* ? - probably the same */

/*** structure factor ***/

#ifndef FREEBC
static int kk,nsf;
static double sumw=0,sumwq=0;
vector SFsumL;
#endif

void initSF(void) /************************************************** initSF */
{
#ifndef FREEBC
  int nn;

  if (el.sf<=0) return;

  prt("Initializing Ewald for structure factor: alpha=%g kappa=%g",el.alpha,el.kappa);

  if (FROM) WARNING(("obsolete option -j and SF"))

  if (el.kappa<=0) {
    if (el.sf) WARNING(("el.sf=%d ignored because el.kappa=%g",el.sf,el.kappa))
    el.sf=0; }

  nsf=0;
  VO(SFsumL,=0)

  /* ugly patch to force Ewald to use weights instead of charges */
  sumw=sumwq=0;
  loop (nn,0,No.N) {
    molecule_t *mn=molec+nn;
    siteinfo_t *si=spec[mn->sp]->si;
    int ns=mn->ns;
    int i;

    loop (i,0,ns) {
      sumw+=si[i].charge=si[i].sfweight;
      sumwq+=Sqr(si[i].charge); } }
  prt("SF: sum weights=%g  sum weights^2=%g",sumw,sumwq);

  /* the structure factor - initialize */
  if (el.sf==1) {
    if (!iscube()) ERROR(("initSF: el.sf=%d\n\
*** sphericalized structure factor can be calculated for a cube only",el.sf))
    el.sfsize=kk=(int)(Sqr(el.kappa*box.L[2])+1.000001);
    allocarrayzero(el.sfr,kk); }
  else if (el.sf==3) {
    /* upper estimate - WASTING MEMORY! */
    kk=(int)(el.kappa*box.L[0]+1.000001)
      *(int)(el.kappa*box.L[1]+1.000001)
      *(int)(el.kappa*box.L[2]+1.000001);
    allocarrayzero(el.sf3d,el.sfsize=kk*4); }

  if (tau.P) {
    VV(oldL,=box.L)
    /* WARNING: this will fail probably if not cube ... */
    if (iscube())
      WARNING(("initSF: el.sf=%d tau.P=%g\n\
*** variable-size cubic box will be scaled to a unit cube",el.sf,tau.P))
    else
      ERROR(("initSF: el.sf=%d tau.P=%g\n\
*** variable non-cubic box not supported",el.sf,tau.P))
    VO(box.L,=1) VO(box.Lh,=0.5)
    Ewald(0,NULL,NULL);
    VV(box.L,=oldL) VV(box.Lh,=0.5*box.L) }
  else
    /* fixed box size */
    Ewald(0,NULL,NULL);
#endif /*# FREEBC */
}

void calculateSF(void) /**************************************** calculateSF */
{
#ifndef FREEBC
  time_t stoptime;

  if (el.sf<=0) return;

  measure=1;

  if (tau.P) {
    /* WARNING: box scaled to 1-cube, structure factor = ??? */
    VV(oldL,=box.L)
    if (!iscube())
      ERROR(("calculateSF: variable non-cubic box not supported"))
    VO(box.L,=1) VO(box.Lh,=0.5)

    /* the structure factor */
    Ewald(1,NULL,cfg[0]->rp);
    Ewald(2,NULL,cfg[0]->rp);

    VV(box.L,=oldL) VV(box.Lh,=0.5*box.L) }
  else {
    /* the structure factor - fixed box */
    Ewald(1,NULL,cfg[0]->rp);
    Ewald(2,NULL,cfg[0]->rp); }

  nsf++;
  VV(SFsumL,+=box.L); /* for SF: cube only */

  if (option('t')) {
    time(&stoptime);
    fprintf(stderr,"%i SF done at %lu = %s",
                    nsf,           stoptime,ctime(&stoptime)); }
#endif /*# FREEBC */
}

void printfSF(void) /************************************************ printfSF */
{
#ifndef FREEBC
  int i;

  if (el.sf<=0) return;

  if (el.sf==1) {
    FILE *SF=fopen(Fn("sfr"),"wt");
    double W=0;

    VO(SFsumL,/=nsf)
    fprintf(SF,"# sphericalized (radial) structure factor\n");
    fprintf(SF,"# %d frames read, <L>=%.8f\n",nsf,SFsumL[2]);
    fprintf(SF,"# k=1/wavelength (NB: in older version, there was factor 2 pi)\n");
    fprintf(SF,"# k/AA^-1    S(k)      smooth[3]   smooth[5]   smooth[7]\n");
    if (option('v')&4) header(" k^2    S(k)     ");
    loop (i,1,kk)
      if (el.sfr[i].nk) Max(W,el.sfr[i].q/el.sfr[i].nk)
    loop (i,1,kk)
      if (el.sfr[i].nk) {
 /*.....      double sfi=el.sfr[i].q/(el.sfr[i].nk*sumwq);*/
        double sfi=1+No.s*(el.sfr[i].q/el.sfr[i].nk-sumwq)/Sqr(sumw);
        static double W1[3]={0.25,0.5,0.25};
        static double W2[5]={0.0625,0.25,0.375,0.25,0.0625};
        static double W3[7]={0.015625,0.09375,0.234375,0.3125,0.234375,0.09375,0.015625};
        double sm1,sm2,sm3;
        double qq,nnk;
        int is;

#  define DOIT(M) \
        qq=nnk=0; \
        loopto (is,-M,M) {\
          int j=i+is; \
          if (j>=1 && j<kk) qq+=el.sfr[j].q*W##M[is+M],nnk+=el.sfr[j].nk*W##M[is+M]; }\
          sm##M=1+No.s*(qq/nnk-sumwq)/Sqr(sumw);

        DOIT(1)
        DOIT(2)
        DOIT(3)

        /* SF: cube only ... */
        fprintf(SF,"%8.5f %11.7f %11.7f %11.7f %11.7f %d\n",
                sqrt((double)i)/SFsumL[2],sfi,sm1,sm2,sm3,el.sfr[i].nk);
        if (option('v')&4) {
          prt_("%3i %9.6f ",i,sfi);
          graph(el.sfr[i].q/el.sfr[i].nk/W,65); } }
      else
        if (option('v')&4) prt("%3i",i);
    header("");
    fclose(SF); }

  if (el.sf==3) {
    FILE *SF=fopen(Fn("sf3d"),"wt");

    VO(SFsumL,/=nsf)
    fprintf(SF,"# 3D structure factor\n");
    fprintf(SF,"# %d frames read, <L>=%.8f %.8f %.8f\n",nsf,VARG(SFsumL));
    fprintf(SF,"# k=(unit vector of direction)/wavelength\n");
    fprintf(SF,"#kx/AA^-1 ky/AA^-1 kz/AA^-1  k/AA^-1      S(k)      nk\n");
    loop (i,0,kk*4)
      if (SQR(el.sf3d[i].k)) {
        double sfi=1+No.s*(el.sf3d[i].q/nsf-sumwq)/Sqr(sumw);

        fprintf(SF,"%8.5f %8.5f %8.5f  %8.5f  %11.7f  %d\n",
                el.sf3d[i].k[0]/SFsumL[0],
                el.sf3d[i].k[1]/SFsumL[1],
                el.sf3d[i].k[2]/SFsumL[2],
                sqrt(Sqr(el.sf3d[i].k[0]/SFsumL[0])+Sqr(el.sf3d[i].k[1]/SFsumL[1])+Sqr(el.sf3d[i].k[2]/SFsumL[2])),
                sfi,nsf); }
    header("");
    fclose(SF); }
#endif /*# FREEBC */
}

static double **charge; /* copy of charges; note that si[].charge is rewritten! */
static ToIntPtr lastcfg,firstcfg;
static vector lastCM;
static vector *rcenter;
static FILE *qcp,*mcp;
static float *rcp;
static int nL,ndiff,nd;
static double sumV,sumL,sumL2,dt,firstt;

static struct maxjump_s
  maxjump= {{0,0,0},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
  max1jump={{0,0,0},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

void initdiff(double dtplb) /************************************** initdiff */
{
  int n,i;
  vector *rrow=rof(molec,cfg[0]->rp); /* r: access by No.s sites */

  if (dtplb==0) dtplb=1; /* info only */

  if (diff.mode==0) return;

  if (tau.P) {
    loop (i,0,No.s) VV(rrow[i],/=box.L)
    sumL+=box.L[0];
    sumL2+=Sqr(box.L[0]);
    sumV+=box.L[0]*box.L[1]*box.L[2];
    nL++; }

  sdsalloc(lastcfg,cfg[0]->size);
  sdsalloc(firstcfg,cfg[0]->size);

  VV(lastL,=box.L)
  if (tau.P) VO(lastL,=1)
  VV(oldL,=box.L) /* ??? duplicated */

  nd=nL=ndiff=0;
  dt=dtplb;
  firstt=t;
  sumV=sumL=sumL2=0;

  if (diff.mode&4) alloczero(rcenter,No.s*sizeof(vector));
  alloczero(charge,sizeof(charge[0])*nspec); // WHY THIS? IN SF???

  loop (n,0,nspec) {
    siteinfo_t *si=spec[n]->si;
    int ns=spec[n]->ns;

    alloc(charge[n],sizeof(*si)*ns);

    loop (i,0,ns) {
#ifdef POLAR
      /* it is recommended to set -a0 unless Gaussian charges */
      charge[n][i]=si[i].charge+si[i].chargepol;
#else /*# POLAR */
      charge[n][i]=si[i].charge;
#endif /*#!POLAR */
    } }

  if (diff.mode&1) mcp=fopen(Fn("m.cp"),"wb");
  if (diff.mode&2) qcp=fopen(Fn("q.cp"),"wb");

  allocarrayzero(rcp,nspec+1);
  rcp[0]=CPmark;
  *(int4*)(rcp+1)=nspec+1;
  loop (i,2,nspec) sprintf((char*)(rcp+i),"%d",i);
  if (nspec>1) copy(rcp+nspec,"mdif",4); /* WARNING: no header for nspec=1 */
  if (mcp) fwrite(rcp,sizeof(*rcp),nspec+1,mcp);
  if (nspec>1) copy(rcp+nspec,"cond",4); /* WARNING: no header for nspec=1 */
  if (qcp) fwrite(rcp,sizeof(*rcp),nspec+1,qcp);

  header(tau.P ? " t-t0[ps] t[ps]       msd           mscd     L [AA]  "
               : " t-t0[ps] t[ps]       msd [AA^2]    mscd     L [AA]  ");
}

void calculatediff(int n,int no) /**************************** calculatediff */
{
  int i,j,sp,ns;
  vector *rrow=rof(molec,cfg[0]->rp); /* r: access by No.s sites */
  vector dr,*firstr,*r=NULL; /* initialized to suppress compile warning */
#ifndef FREEBC
  vector *lastr;
  vector CM;
#endif
  struct msd_s {
    double m,q;
    int n;
  } *msd,*msp;
  vector msdm,msdq;
  double totm,totq,CMshift;
  siteinfo_t *si;
  struct maxjump_s maxthisjump= {{0,0,0},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

  if (tau.P) {
    loop (i,0,No.s) VV(rrow[i],/=box.L)
    if (!iscube()) ERROR(("not cube (2nd check)"))

    sumL+=box.L[0];
    sumL2+=Sqr(box.L[0]);
    sumV+=PROD(box.L);
    nL++; }

  /* remove nearest image jumps in the molecule moves */
#ifndef FREEBC
  if (!ndiff) {
    sdscopy(firstcfg,cfg[0])
    CoM(lastCM,firstcfg);
    fprintf(stderr,"%g %g %g\n",VARG(lastCM));
  }
  else {
    if (tau.P) { VV(oldL,=box.L) VO(box.L,=1) VO(box.Lh,=0.5) }

    loop (j,FROM,No.N) {
      firstr=rof(molec+j,firstcfg->rp);
      r=     rof(molec+j,cfg[0]->rp);
      lastr= rof(molec+j,lastcfg->rp);
      ns=molec[j].ns;

      loop (i,0,ns) {
        int k;
        double D;

        loop (k,0,3) {
          D=r[i][k]/box.L[k]-lastr[i][k]/lastL[k];
          while (D>0.5)  { r[i][k]-=box.L[k]; D-=1; }
          while (D<-0.5) { r[i][k]+=box.L[k]; D+=1; }
          /* now r (cfg[0]) moves continuously over periodic boxes */

          if (fabs(D)>fabs(maxthisjump.xi[k])) {
            maxthisjump.xi[k]=D;
            maxthisjump.frame[k]=n;
            maxthisjump.n[k]=j;
            maxthisjump.i[k]=i;
            maxthisjump.no=no; }

          if (fabs(D)>fabs(max1jump.xi[k])) {
            max1jump.xi[k]=D;
            max1jump.frame[k]=n;
            max1jump.n[k]=j;
            max1jump.i[k]=i;
            max1jump.no=no; }

          D=r[i][k]/box.L[k]-firstr[i][k]/lastL[k];
          if (fabs(D)>fabs(maxjump.xi[k])) {
            maxjump.xi[k]=D;
            maxjump.frame[k]=n;
            maxjump.n[k]=j;
            maxjump.i[k]=i;
            maxjump.no=no; } } } }

    CoM(CM,cfg[0]);
    CMshift=sqrt(SQRD(lastCM,CM));
    if (CMshift>CMSHIFTERR) {
      double maxD=0;
      int k,kk=-1;

      putv(max1jump.xi)
      putv(max1jump.frame)
      putv(max1jump.n)
      putv(max1jump.i)
      put(max1jump.no)

      putv(maxthisjump.xi)
      putv(maxthisjump.frame)
      putv(maxthisjump.n)
      putv(maxthisjump.i)
      put(maxthisjump.no)

      WARNING(("center of mass shifted by %g\n\
*** first(n=%d) = %g %g %g\n\
*** this(no=%d) = %g %g %g\n\
*** I will try to fix it assuming that one molecule\n\
*** moved by more than L/2 in one direction.\n",
                                        CMshift,
            n,    VARG(lastCM),
            no,   VARG(CM)))
        // sleep(1);

      loop (k,0,3) if (fabs(maxthisjump.xi[k])>fabs(maxD)) {
        j=maxthisjump.n[k];
        maxD=maxthisjump.xi[k];
        kk=k; }
      if (k<0) ERROR(("internal"))
      k=kk;
      put3(j,k,maxD)
      ns=molec[j].ns;

      loop (i,0,ns) r[i][k]-=sign(maxD)*box.L[k];
      maxD-=sign(maxD);
      maxjump.xi[k]=maxD;
      max1jump.xi[k]=maxD;

      CoM(CM,cfg[0]);
      CMshift=sqrt(SQRD(lastCM,CM));
      put(CMshift)
      if (CMshift>CMSHIFTERR)
        ERROR(("sorry, the fixup failed"))
    }

    VV(lastCM,=CM)

    if (tau.P) {
      VV(box.L,=oldL)
      VV(box.Lh,=0.5*box.L) } }
#endif /*# FREEBC */

  ndiff++;
  fprintf(stderr,"\r%d %d:%d ",ndiff,n,no);

  allocarrayzero(msd,nspec);
  VVO(msdm,=msdq,=0)

  loop (j,FROM,No.N) {
    vector vm,vq;
    double summ=0;

    VVO(vm,=vq,=0)
    firstr=rof(molec+j,firstcfg->rp);
    r=rof(molec+j,cfg[0]->rp);
    ns=molec[j].ns;
    sp=molec[j].sp;
    msp=msd+sp;
    msp->n++;
    si=spec[sp]->si;

    loop (i,0,ns) {
      VVV(dr,=r[i],-firstr[i])

      VV(vq,+=charge[sp][i]*dr)
      VV(vm,+=si[i].mass*dr)
      summ += si[i].mass; }

    msp->m+=SQR(vm)/Sqr(summ)/6;
    msp->q+=SQR(vq)/6;
    VV(msdm,+=vm)
    VV(msdq,+=vq) } /* molecules */

  if (diff.mode&4) loop (i,0,No.s) VV(rcenter[i],+=rrow[i])
  nd++;

  putv(msdm)

  totm=SQR(msdm)/6;
  totq=SQR(msdq)/6;

  prt("%7.3f %8.3f %12.5g %12.5g (%9.5f,%9.5f,%9.5f)",
      (n-1)*dt,firstt,totm,totq,VARG(box.L));
  loop (i,0,nspec) rcp[i]=msd[i].m/(msd[i].n+(msd[i].n==0));
  rcp[nspec]=totm;
  if (mcp) fwrite(rcp,sizeof(*rcp),nspec+1,mcp);

  loop (i,0,nspec) rcp[i]=msd[i].q;
  rcp[nspec]=totq;
  if (qcp) fwrite(rcp,sizeof(*rcp),nspec+1,qcp);

  sdscopy(lastcfg,cfg[0]);
  VV(lastL,=box.L)
  if (tau.P) VO(lastL,=1)

  free(msd);
}

void printdiff(void) /******************************************** printdiff */
{
  FILE *aux=fopen(Fn("aux"),"wb"); /* separate files in V3.0 */
  FILE *fcenter=NULL;
  float fr[3];
  float q;
  int i,k;

  header("");
  if (mcp) fclose(mcp);
  if (qcp) fclose(qcp);

  if (diff.mode&4) fcenter=fopen("center.plb","wb");

  loop (k,0,3) if (maxjump.i[k]>=0) {
    prt("%c-axis:\n\
max jump over periodic boxes = %g*L (frames 1->%d, mol.s=%d.%d)",
        k+'x',maxjump.xi[k],maxjump.frame[k],maxjump.n[k],maxjump.i[k]);
    prt("max jump between consecutive frames: %g*L (frames %d->%d, mol.s=%d.%d)",
        max1jump.xi[k],max1jump.frame[k]-1,max1jump.frame[k],max1jump.n[k],max1jump.i[k]);
    if (fabs(max1jump.xi[k])>0.4)
      WARNING(("\
max jump between consecutive frames is too close to L/2\n\
*** diffusions etc. may be wrong - check whether msd and column difm are tiny\n\
*** if wrong then simulate again with shorter dt.plb")) }

  fwrite(&maxjump,sizeof(maxjump),1,aux);
  fwrite(&max1jump,sizeof(max1jump),1,aux);
  fclose(aux);

  if (tau.P) prt("multiplication factor for diffusions and conductivities: <L^2>=%f\n\
  <L>^2=%f <V>^2/3=%f",sumL2/nL,Sqr(sumL/nL),pow(sumV/nL,2./3));

  if (fcenter) {
    fr[0]=No.s;
    if (tau.P) {
      fr[1]=sqrt(sumL2/nL);
      q=fr[1]/nd; }
    else {
      fr[1]=-3;
      q=1./nd; }
    fwrite(fr,4,2,fcenter);
    VV(fr,=box.L)
      fwrite(fr,4,3,fcenter);
    loop (i,0,No.s) {
      VV(fr,=q*rcenter[i])
        fwrite(fr,4,3,fcenter); }
    fclose(fcenter); }

  free(lastcfg);
  free(firstcfg);
}


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef XSECTION
#  include "sphint.c"

typedef vector matrix[DIM];

static void rndmat(matrix M) /*************************************** rndmat */
/* random orientational matrix, orthonormal */
{
  double rr;

  rndsphere(M[0]);
  do rr=rndball(M[2]); while (fabs(SCAL(M[0],M[2]))<0.0625);
  VECT(M[1],M[2],M[0])
  rr=sqrt(SQR(M[1]));
  VO(M[1],/=rr);
  VECT(M[2],M[0],M[1])
}

static void mpl(vector a,matrix M,vector b) /*************************** mpl */
/* a:=M.b */
{
  int i;

  loop (i,0,DIM) a[i]=SCAL(M[i],b);
}

#  define XSECTION_M 1
/* molecule-based version */
#  include "xsection.c"

#  undef XSECTION_M
#  define XSECTION_M 0
/* whole configuration version */
#  include "xsection.c"

#  ifdef CLUSTERS
#    undef XSECTION_M
#    define XSECTION_M 2
/* cluster version */
#    include "xsection.c"
#  endif /*# CLUSTERS */

#endif /*# XSECTION */

#ifdef BJERRUM
#  include "bjerrum.c"
#endif

#include "drifts.c"

#ifdef SPCTCF
# include "spctcf.c"
#endif
