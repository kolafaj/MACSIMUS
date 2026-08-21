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
#ifdef CLUSTERS
#  include "intermac.h"
#  include "cluster.c" /* updated and fixed 9/2009 */
extern cl_t cl;
#endif /*# CLUSTERS */


/*** structure factor ***/

#ifndef FREEBC
static int kk,nsf;
static double sumw=0,sumwq=0;
vector SFsumL;
static vector lastL;
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
    VV(lastL,=box.L)
    /* WARNING: this will fail probably if not cube ... */
    if (iscube())
      WARNING(("initSF: el.sf=%d tau.P=%g\n\
*** variable-size cubic box will be scaled to a unit cube",el.sf,tau.P))
    else
      ERROR(("initSF: el.sf=%d tau.P=%g\n\
*** variable non-cubic box not supported",el.sf,tau.P))
    VO(box.L,=1) VO(box.Lh,=0.5)
    Ewald(0,NULL,NULL);
    VV(box.L,=lastL) VV(box.Lh,=0.5*box.L) }
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
    VV(lastL,=box.L)
    if (!iscube())
      ERROR(("calculateSF: variable non-cubic box not supported"))
    VO(box.L,=1) VO(box.Lh,=0.5)

    /* the structure factor */
    Ewald(1,NULL,cfg[0]->rp);
    Ewald(2,NULL,cfg[0]->rp);

    VV(box.L,=lastL) VV(box.Lh,=0.5*box.L) }
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
static vector CMcfglast;
static vector *rcenter;
static FILE *qcp,*mcp;
static float *rcp;
static int nL,ndiff,nd;
static double sumV,sumL,sumL2,dt,firstt;

/* NB: "jump" here means the displacement extended by following
       the trajectory over periodic boundaries */
static struct maxjump_s
  /* max jump over the whole trajectory: */
maxjump= {{0,0,0},{-1,-1,-1},{-1,-1,-1},-1},
  /* max of all 1-frame jumps: */
  max1jump={{0,0,0},{-1,-1,-1},{-1,-1,-1},-1};

void initMSD(double dtplb) /**************************************** initMSD */
/* preparation of data structures
   - the first cfg is read on the 1st call of calculateMSD() */
{
  int n,i;

  if (dtplb==0) dtplb=1; /* info only */

  if (MSD.mode==0) return;

  if (box.follow>2e-6) prt("\
WARNING: box.follow=%g (plb2diff -F%g) is likely too long.\n\
         Please use as short value as possible.",box.follow,box.follow);

  sdsalloc(lastcfg,cfg[0]->size);
  sdsalloc(firstcfg,cfg[0]->size);

  nd=nL=ndiff=0; /* ndiff==0 tells calculateMSD() to store */
  dt=dtplb;
  sumV=sumL=sumL2=0;

  if (MSD.mode&4) alloczero(rcenter,No.s*sizeof(vector));
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

  if (MSD.mode&1) mcp=fopen(Fn("m.cp"),"wb");
  if (MSD.mode&2) qcp=fopen(Fn("q.cp"),"wb");

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

void calculateMSD(int frame,int to) /************************** calculateMSD */
{
  int i,n,k,sp,ns;
  vector dr,*firstr,*r=NULL; /* initialized to suppress compile warning */
#ifndef FREEBC
  vector *lastr;
  vector CMcfg;
  static vector firstL;
#endif
  struct msd_s {
    double m,q;
    int n;
  } *msd,*msp;
  vector msdm,msdq,CMshift;
  double totm,totq;
  siteinfo_t *si;
  struct maxjump_s
    maxthisjump= {{0,0,0},{-1,-1,-1},{-1,-1,-1},-1},
    nextthisjump= {{0,0,0},{-1,-1,-1},{-1,-1,-1},-1};

  if (tau.P) {
    if (!iscube()) ERROR(("not cube (2nd check)"))

    sumL+=box.L[0];
    sumL2+=Sqr(box.L[0]);
    sumV+=PROD(box.L);
    nL++; }

  if (!ndiff) sdscopy(firstcfg,cfg[0])

#ifndef FREEBC
  /* remove nearest image jumps in the molecule moves */
  if (!ndiff) {
    CoM(CMcfglast,firstcfg);
    VV(firstL,=box.L)
    firstt=t;
    VV(CMcfglast,/=box.L) }
  else {
    double D;

    loop (n,FROM,No.N) {
      vector CMmol={0,0,0},CMlast={0,0,0},CMfirst={0,0,0}; /* inefficient */
      molecule_t *mn=molec+n;
      siteinfo_t *si=spec[mn->sp]->si;

      r=     rof(mn,cfg[0]->rp);
      firstr=rof(mn,firstcfg->rp);
      lastr= rof(mn,lastcfg->rp);
      ns=molec[n].ns;

      /* here we rescale to unit box */
      loop (i,0,ns) {
        VVV(CMmol,+=si[i].mass*r[i],/box.L)
        /* inefficient: unnecessarily recalculated, but more compact code */
        VVV(CMlast,+=si[i].mass*lastr[i],/lastL)
        VVV(CMfirst,+=si[i].mass*firstr[i],/firstL) }
      VO(CMmol,/=spec[mn->sp]->mass)
      VO(CMlast,/=spec[mn->sp]->mass)
      VO(CMfirst,/=spec[mn->sp]->mass)

      loop (k,0,3) {
        D=CMmol[k]-CMlast[k];
        while (D>0.5)  { loop (i,0,ns) r[i][k]-=box.L[k]; D-=1; CMmol[k]-=1; }
        while (D<-0.5) { loop (i,0,ns) r[i][k]+=box.L[k]; D+=1; CMmol[k]+=1; }
        /* now r (cfg[0]) moves continuously over periodic boxes */

        /* two max jumps after 1 frame */
        if (fabs(D)>fabs(maxthisjump.r[k])) {
          nextthisjump.r[k]=maxthisjump.r[k];
          maxthisjump.r[k]=D;
          nextthisjump.frame[k]=maxthisjump.frame[k];
          maxthisjump.frame[k]=frame;
          nextthisjump.n[k]=maxthisjump.n[k];
          maxthisjump.n[k]=n;
          nextthisjump.to=maxthisjump.to;
          maxthisjump.to=to; }

        /* record the absolute maximum */
        if (fabs(D)>fabs(max1jump.r[k])) {
          max1jump.r[k]=D;
          max1jump.frame[k]=frame;
          max1jump.n[k]=n;
          max1jump.to=to; }

        /* displacement from 1st configuration */
        D=CMmol[k]-CMfirst[k];
        if (fabs(D)>fabs(maxjump.r[k])) {
          maxjump.r[k]=D;
          maxjump.frame[k]=frame;
          maxjump.n[k]=n;
          maxjump.to=to; } } } /*n*/

    CoM(CMcfg,cfg[0]);
    VV(CMcfg,/=box.L)
    VVV(CMshift,=CMcfg,-CMcfglast)

    putv(CMshift)

    loop (k,0,3) {

      if (fabs(CMshift[k])>box.follow) {
        /* center of mass has shifted too much */

        WARNING(("\
Center of mass in coordinate %c (k=%d) has shifted by %g*L[k]\n\
*** in frame %d (t=%g) in the block starting from frame %d.\n\
*** CM(prev)= %13.9f  CM(this)= %13.9f  CMshift = %13.9f\n\
*** Double-check box.follow=%g (plb2diff -F%g)!\n\
*** I'll try to fix it assuming that a single molecule has moved by > L/2.",
                 'x'+k,k,CMshift[k],
                 to, frame, t,
               CMcfglast[k],CMcfg[k],CMshift[k],box.follow,box.follow))

        prt("max1jump.r=%g frame=%d n=%d no=%d (max jump so far)",
            max1jump.r[k], max1jump.frame[k],
            max1jump.n[k], max1jump.to);
        prt("maxthisjump.r=%g frame=%d n=%d no=%d (max jump from prev frame)",
            maxthisjump.r[k], maxthisjump.frame[k],
            maxthisjump.n[k], maxthisjump.to);
        prt("nextthisjump.r=%g frame=%d n=%d no=%d (2nd max jump from prev frame)",
            nextthisjump.r[k], nextthisjump.frame[k],
            nextthisjump.n[k], nextthisjump.to);
        D=fabs(maxthisjump.r[k])-fabs(nextthisjump.r[k]);
        if (D<0.5-box.jump)
          ERROR(("A molecule to follow is likely not unique (D=%g),\n\
*** check box.jump (plb2diff -J).",D))
        if (D<2*(0.5-box.jump))
          WARNING(("A molecule to follow may not be unique (D=%g),\n\
*** check box.jump (plb2diff -J).",D))
        else
          prt("The next next longest jump is short enough (D=%g).",D);

        /* fix using the longest jump */
        n=maxthisjump.n[k];
        D=maxthisjump.r[k];
        prt("The longest jump found for mol=%d jump[%d]=%g L",n,k,D);
        ns=molec[n].ns;

        r=rof(molec+n,cfg[0]->rp);
        loop (i,0,ns) r[i][k]-=sign(D)*box.L[k];
        D-=sign(D);
        maxjump.r[k]=D;
        max1jump.r[k]=D;

        CoM(CMcfg,cfg[0]);
        VV(CMcfg,/=box.L)
        VVV(CMshift,=CMcfg,-CMcfglast)
        prt("Corrected CoM shift = %g",CMshift[k]);
        if (CMshift[k]<box.follow)
          prt("*** The fixup has been successful. ***\n");
      else
        ERROR(("The fixup failed.")) }
    } /* k */

    VV(CMcfglast,=CMcfg) }

#endif /*# FREEBC */

  ndiff++;
  fprintf(stderr,"\r%d %d:%d ",ndiff,frame,to);

  allocarrayzero(msd,nspec);
  VVO(msdm,=msdq,=0)

  loop (n,FROM,No.N) {
    vector vm,vq;
    double summ=0;

    VVO(vm,=vq,=0)
    firstr=rof(molec+n,firstcfg->rp);
    r=rof(molec+n,cfg[0]->rp);
    ns=molec[n].ns;
    sp=molec[n].sp;
    msp=msd+sp;
    msp->n++;
    si=spec[sp]->si;

    loop (i,0,ns) {
      loop (k,0,3) dr[k]=r[i][k]/box.L[k]-firstr[i][k]/firstL[k];
      VV(vq,+=charge[sp][i]*dr)
      VV(vm,+=si[i].mass*dr)
      summ += si[i].mass; }

    msp->m+=SQR(vm)/Sqr(summ)/6;
    msp->q+=SQR(vq)/6;
    VV(msdm,+=vm)
    VV(msdq,+=vq) } /* molecules */

  //  putv(msdm)

  if (MSD.mode&4) loop (i,0,No.s) VVV(rcenter[i],+=rof(molec,cfg[0]->rp)[i],/box.L)
  nd++;

  totm=SQR(msdm)/6;
  totq=SQR(msdq)/6;

  prt("%7.3f %8.3f %12.5g %12.5g (%9.5f,%9.5f,%9.5f)",
      (frame-1)*dt,firstt,totm,totq,VARG(box.L));
  loop (i,0,nspec) rcp[i]=msd[i].m/(msd[i].n+(msd[i].n==0));
  rcp[nspec]=totm;
  if (mcp) fwrite(rcp,sizeof(*rcp),nspec+1,mcp);

  loop (i,0,nspec) rcp[i]=msd[i].q;
  rcp[nspec]=totq;
  if (qcp) fwrite(rcp,sizeof(*rcp),nspec+1,qcp);

  sdscopy(lastcfg,cfg[0]);
  VV(lastL,=box.L)

  free(msd);
}

void printMSD(void) /********************************************** printMSD */
{
  FILE *aux=fopen(Fn("aux"),"wb"); /* separate files in V3.0 */
  FILE *fcenter=NULL;
  float fr[3];
  float q;
  int i,k;

  header("");
  if (mcp) fclose(mcp);
  if (qcp) fclose(qcp);

  if (MSD.mode&4) fcenter=fopen("center.plb","wb");

  loop (k,0,3) if (maxjump.n[k]>=0) {
    prt("%c-axis:\n\
max jump over periodic boxes = %g*L (frames 1->%d, mol=%d)",
        k+'x',maxjump.r[k],maxjump.frame[k],maxjump.n[k]);
    prt("max jump between consecutive frames: %g*L (frames %d->%d, mol=%d)",
        max1jump.r[k],max1jump.frame[k]-1,max1jump.frame[k],max1jump.n[k]);
    if (fabs(max1jump.r[k])>box.jump)
      WARNING(("\
The max raw jump between consecutive frames in %c = %g*L ~ L/2\n\
*** Diffusions/conductivity may be wrong!\n\
*** Check MSD and column difm; if not tiny, decrease dt.plb and simulate again!\n\
*** Current setup: box.follow=%g (plb2diff -F), box.jump=%g (plb2diff -J)",
               'x'+k,max1jump.r[k],
               box.follow,box.jump)) }

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
