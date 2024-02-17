/*
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                  initialize/load/save configuration                   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "ground.h"
#include "sds.h"
#include "simglob.h"
#include "varfile.h"
#include "units.h"
#include "siminit.h"
#include "norm.h"
#include "forces.h"
#include "interpot.h"
#ifdef SLAB
#  include "simpot.h"
#endif  /*# SLAB */
#include "simdef.h"
#include "simils.h"
#include "options.h"
#ifdef POLAR
#  include "constrd.h"
#endif /*# POLAR */

#ifdef Reverse
#  error compile-time option Reverse no longer supported
#endif /*# Reverse */

#include <unistd.h> /* sleep */

extern volatile int sig; /* defined in main.c */

typedef float Epairreal; /* to spare space */

char *lastFn;

#define FNLEN 256

char *Fn(char *ext) /**************************************************** Fn */
/***
 - pointer to a static string containing SIMNAME.EXT is returned
 - two static strings are switched => you can use Fn() twice in one expression
 - if EXT="plb","vlb",NUMBER, then PLBNAME.EXT is returned (if PLBNAME not NULL)
***/
{
  static char n1[FNLEN],n2[FNLEN],*n;
  char *name=simils.simname;

  n = n==n1 ? n2 : n1; /* trick to be able to use the result of Fn twice */

  if (simils.plbname) {
    if (!strcmp("plb",ext) || !strcmp("vlb",ext) || isdigit(ext[0]))
      name=simils.plbname; }
  if (strlen(name)+strlen(ext)>FNLEN-4)
    ERROR(("Fn: file name %s.%s is too long",name,ext))
  sprintf(n,"%s.%s",name,ext);

  return lastFn=n;
} /* Fn */

const char *getext(const char *fn) /********************************* getext */
  /* returns extension incl. ., NULL if no extension */
{
  const char *c,*dot=NULL;

  for (c=fn; *c; c++ ) {
    if (*c=='.') dot=c;
    if (*c=='/') dot=NULL; }

  return dot;
}

char *stripext(const char *fn) /*********************************** stripext */
  /* make a permanent copy of fn without 4 characters of .EXT */
{
  int l=strlen(fn)-3;
  char *c=malloc(l);

  memcpy(c,fn,l);
  c[l-1]=0;
  return c;
}

void backup(char *ext) /********************************************* backup */
{
  char *fn,*bak;
  char *bakext;

  alloc(bakext,strlen(ext)+2);
  sprintf(bakext,"%s~",ext);

  fn=Fn(ext); bak=Fn(bakext);
  if (!remove(bak)) prt("%s removed",bak);
  if (!rename(fn,bak)) prt("%s renamed to %s",fn,bak);
  free(bakext);
}

void waitfordiskspace(int blocks) /************************ waitfordiskspace */
/*
  waits until at least blocks+option('w') kB is available on the file system
  active only if option('w')>1
*/
{
  if (option('w')>1)
    for (;;) {
      int delay;
      static int lastdelay,sum;
      FILE *f=NULL;
      char line[128],*tok="0";
      blocks+=option('w');

      line[0]=0;
      if (system(string("df . > %s",Fn("df")))) delay=2;
      else if (!(f=fopen(lastFn,"rt"))) delay=2;
      else if (!fgets(line,128,f)) delay=2;
      else if (!fgets(line,128,f)) delay=2;
      else if (!(tok=strtok(line,"\n\t "))) delay=2;
      else if (!(tok=strtok(NULL,"\n\t "))) delay=2;
      else if (!(tok=strtok(NULL,"\n\t "))) delay=2;
      else if (!(tok=strtok(NULL,"\n\t "))) delay=2;
      else if (atoi(tok)<blocks/2) delay=2;
      else if (atoi(tok)<blocks) delay=1;
      else delay=0;

      if (f) fclose(f);
      remove(Fn("df"));
      if (!delay) {
        if (lastdelay)
          WARNING(("disk has been full, now enough space has been detected\n\
*** I have had to wait for disk space for %d s, now I am continuing",sum))
        sum=lastdelay=0;
        return; }

      lastdelay=delay;
      delay=60*(delay+lastdelay)+irnd(60);
      sum+=delay;
      sleep(delay); }
}

int MC=0;
int mirror;
struct simils_s simils={ {1,1,1} };

static void rotate(vector *c,int ns,double drot) /******************* rotate */
/*** rotation by random angle in (-drot,drot) around random axix ***/
{
  int j=irnd(DIM), k=(j+1)%DIM, i;
  double z=rndcos()*drot, cf=cos(z), sf=sin(z);

  loop (i,0,ns) {
    z=c[i][j];
    c[i][j]=cf*z-sf*c[i][k];
    c[i][k]=sf*z+cf*c[i][k]; }
}


void rndorientation(vector *c,int ns) /********************** rndorientation */
/*** random orientation; if mirror then also random mirror reflection ***/
{
  vector o[3],r;
  double a;
  int i,j;

  rndsphere(o[0]);
  do {
    rndsphere(o[1]);
    a=SCAL(o[0],o[1]); }
  while (fabs(SCAL(o[0],o[1])>0.875));
  VV(o[1],-=a*o[0])
  a=sqrt(SQR(o[1]));
  VO(o[1],/=a)

  if (mirror && rnd()>0.5) VECT(o[2],o[1],o[0])
  else VECT(o[2],o[0],o[1])

  loop (j,0,ns) {
    VO(r,=0)
    loop (i,0,3) r[i]+=SCAL(o[i],c[j]);
    VV(c[j],=r) }
}

#ifdef FREEBC
static double sticky;

static double stickypot(int ns1,vector *r1,int ns2,vector *r2)
{
  int i,j;
  double rr,minrr=9e9;

  loop (i,0,ns1) loop (j,0,ns2) {
    rr=SQRD(r1[i],r2[j]);
    Min(minrr,rr) }

  return minrr/sticky;
}
#endif /*# FREEBC */

#if 0 /* debug */
static void checkpairE(Epairreal **Epair,double initqf)
{
  int m,n;
  double EE;
  vector *rp=a[0]->rp,*rpf=a[2]->rp;

  loop (n,0,No.N)
    loop (m,0,n) {
      En.el=0;
      EE= (*pot[molec[n].sp][molec[m].sp]) (
             rof(cfg+n,rpf),rof(cfg+m,rpf), /*dummy*/
             cfg+n, cfg+m, rp ) + En.el*initqf;
#  ifdef FREEBC
      if (sticky!=0 && m<FROM) EE+=stickypot(molec[n].ns,rof(cfg+n,rp),molec[m].ns,rof(cfg+m,rp))*T;
#  endif /*# FREEBC */
      if (fabs(Epair[n][m]-EE)>1e-9)
        ERROR(("(%d,%d) Epair=%.5f real=%.5f dif=%g",
               n,m,(double)Epair[n][m],EE,Epair[n][m]-EE)) }
}
#endif /*# 0 */

static void doMC(Epairreal **Epair,int maxns,double Emax,double initqf,int nplb)
/* WARNING: MC may require a LOT of memory - pair table */
{
#ifdef GOLD
  ERROR(("MC for GOLD not implemented"))
#else /*# GOLD */
  static struct {
    double d;
    double ar;
    int all,acc; } move[3]={ {0}, {.2}, {.5} };
  vector rdis,*old,*rp=cfg[0]->rp,*rpf=cfg[2]->rp;
  int m,n=0,i,maxn,nn,nMC=0;
  double Enew,Eold,EE,maxE;
  Epairreal *Epairnew;

  move[0].d=fmin(fmin(box.L[0],box.L[1]),box.L[2]);
  prt("\nMonte Carlo          shoot              displace            rotate");
  prt("___#_______maxE__  _all_acc.r__length  _all_acc.r__length  _all_acc.r.__angle");
  alloc(old,sizeof(vector)*maxns);
  alloc(Epairnew,No.N*sizeof(Epairnew[0]));

  maxn=No.N+1+(No.N-FROM)/4;

  do {
    if (nplb) if (nMC%nplb==0) writeplayback();
    nMC++;

    maxE=-3e33;
    loop (i,0,3) move[i].acc=move[i].all=0;

/*.....checkpairE(Epair,initqf);*/

    loop (nn,FROM,maxn) { /* sweep */

      if (nn<No.N)
        /* normal move, molecules selected sequentially */
        n=nn;
      else {
        /* find the molecule with highest energy (violates Metropolis!) */
        double maxU=-3e33,sum;

        loop (i,FROM,No.N) {
          sum=0;
          loop (m,0,i) sum+=Epair[i][m];
          loop (m,i+1,No.N) sum+=Epair[m][i];
          if (sum>maxU) maxU=sum,n=i; } }

      {
        molecule_t *mn=molec+n;
        int ns=mn->ns;
        vector *r=rof(mn,rp);
        enum MCmove_e { SHOOT,DISPLACE,ROTATE } MCmove; /* this order! */

        if (ns>maxns) ERROR(("%d %d",ns,maxns))

        copy(old,r,sizeof(vector)*ns);

        MCmove=(enum MCmove_e)(irnd(2+(ns>1)));
        move[MCmove].all++;

        switch (MCmove) {
          case SHOOT:
            if (n>=No.rotatefrom) rndorientation(r,ns);
#  ifdef FREEBC
            VV(rdis,=rndgauss()*box.Lh)
#  else /*# FREEBC */
            VV(rdis,=rnd()*box.L)
#    ifdef SLAB
            if (wall.n) rdis[2]=(wall.z[1]-wall.z[0])*rdis[2]+wall.z[0]*box.L[2];
#    endif /*# SLAB */
#  endif /*#!FREEBC */
            goto ANYDISPLACE;
          case DISPLACE:
            rndball(rdis);
            VO(rdis,*=move[MCmove].d);
         ANYDISPLACE:
            loop (i,0,ns) VV(r[i],+=rdis)
            break;
          case ROTATE:
        rotate(r,ns,move[MCmove].d); }

        normalize(n);

        Enew=Eold=0;

        //#ifdef FREEBC
        //        /* central force */
        //        if (center.K) {
        //          double r0q=Sqr(center.r0);
        //
        //          if (center.K<0) WARNING(("CHECK THIS -- sticky?"))
        //          if (center.r0) loop (m,0,ns) {
        //            double rr=SQR(r[m]);
        //
        //            if (rr>r0q) Enew+=sqr(sqrt(rr)-center.r0);
        //            rr=SQR(old[m]);
        //            if (rr>r0q) Eold+=sqr(sqrt(rr)-center.r0); }
        //          else loop (m,0,ns) {
        //            Enew+=SQR(r[m]);
        //            Eold+=SQR(old[m]); }
        //
        //          Enew*=center.K;
        //          Eold*=center.K; }
        //#endif

        loop (m,0,No.N) if (m!=n) {
          En.el=0;
          EE=(*pot[mn->sp][molec[m].sp]) (
                rof(mn,rpf),rof(molec+m,rpf), /*dummy*/
                mn, molec+m, rp ) + En.el*initqf;
#  ifdef FREEBC
          if (sticky!=0 && m<FROM) EE+=stickypot(mn->ns,rof(mn,rp),molec[m].ns,rof(molec+m,rp))*T;
#  endif /*# FREEBC */
          Epairnew[m]=EE;
          if (m<n) Eold+=Epair[n][m]; else Eold+=Epair[m][n];
          Enew+=EE; }

        if (Enew-Eold<-T*log(rnd()+1e-33)) {
          /* accepted */
          move[MCmove].acc++;
          loop (m,0,n) Epair[n][m]=Epairnew[m];
          loop (m,n+1,No.N) Epair[m][n]=Epairnew[m]; }
        else {
          /* rejected */
          Enew=Eold;
          copy(r,old,sizeof(vector)*ns); }
        Max(maxE,Enew); } } /* nn */

    loop (i,0,3) if (move[i].all) {
      move[i].ar=(double)move[i].acc/move[i].all;
      if (i) {
        /* to reach acceptance ratio exp(-1)=0.368 */
        double q=move[i].ar;
        if (q>0.60653066) q=2;
        else if (q<0.13533528) q=0.5;
        else q=-1/log(q);
        move[i].d*=pow(q,1-0.9/sqrt((double)move[i].all)); } }
    loop (i,0,DIM) Min(move[1].d,box.Lh[i])
    Min(move[2].d,PI)

    prt_("%4d%12.4g",nMC,maxE);
    loop (i,0,3) prt_(" %5d %.3f %7.4f",move[i].all,move[i].ar,move[i].d);
    _n

  } while (!( (maxE<Emax && MC<0) || sig || --MC==0 ));

  sig=0;
  free(Epairnew);
  free(old);
#endif /*#!GOLD */
}


#define declare_cryst(NS) struct { \
  int ns;      /* number of sites in a basic cubic cell */ \
  vector r[NS];/* [ns] site coordinates */ }

/* generic type: */
typedef declare_cryst(1) cryst_t;

/* basic crystal structures:
  note that the argument of declare_cryst must be the same as member ns
  for SC,BCC,FCC: reference points for a molecule
  coordinates relative - to be multiplied by cell size */
static declare_cryst(1) SC =  { 1,{{0.5,0.5,0.5}}};
static declare_cryst(2) BCC = { 2,{{0.25,0.25,0.25},
                                   {0.75,0.75,0.75}}};
static declare_cryst(4) FCC = { 4,{{0.25,0.25,0.25},
                                   {0.25,0.75,0.75},
                                   {0.75,0.75,0.25},
                                   {0.75,0.25,0.75}}};

#define NPINS 5
static cryst_t *cryst[NPINS]={
  NULL, /* n.a. */
  (cryst_t*)& SC,
  (cryst_t*)& BCC,
  (cryst_t*)& FCC,
  NULL /* 4=external cfg */
};

static double cellsize;

static int Nvert(int nL[3],int npins,char *info) /******************** Nvert */
{
  int ns=cryst[npins]->ns; /* sites in elementary cubic lattice */
  double x1=box.L[0]+box.L[1]+box.L[2],x0=x1*1e-10,x;
  int k,nvert=0;
  int nL0[3];
  double avert;

  if (npins<1 || npins>=NPINS) ERROR(("initcryst: %d is bad pins",npins))

  cellsize=0;

  /* determine elem. box size by interval halving */
  do {
    x=sqrt(x0*x1);
    avert=ns;
    loop (k,0,3) avert*=(nL0[k]=(int)(box.L[k]/x));
    if (avert<No.N) x1=x;
    else { nvert=avert; cellsize=x0=x; VV(nL,=nL0) }
  } while (fabs(x1/x0-1)>1e-13);

  if (option('v')&1) prt("  %d:%-3s %d=%dx%dx%d vertices %s",
                            npins,
                               npins==1?"sc":
                               npins==2?"bcc":
                               npins==3?"fcc":"from cfg",
                                    nvert,nL[0],nL[1],nL[2],info);
  if (cellsize==0 || nvert==0)
    ERROR(("initcryst: npins=%d cannot calculate # of sites",npins))

  return nvert;
}

void initcryst(int npins,double Emax) /*************************** initcryst */
/***
  Initial cfg is a crystal with possible holes.
  Molecules/holes are distributed randomly.
  Orientations are random unless -q-XXX
  npins=0 auto selection (min vacancies)
  npins=1 sc
  npins=2 bcc
  npins=3 fcc
  npins>3 config with npins molecules in ASCII file cfg<ns>.<npins>
  TWODIM not supported.
***/
{
  int i,n,nx,sp= -1, insertkey,nscfg=0,k;
  vector *c=NULL,*cc=NULL;
  vector *rp=cfg[0]->rp,*rpv=cfg[1]->rp,*rpf=cfg[2]->rp,*r;

  int nL[DIM],ir[DIM],ivu,nscell,inscell;
  double Lcell,initqf=1.0;

  char **vertexused; /* [site in cell] [molecule] */

  underline("initial configuration");

  prt("pins=%d: initial configuration = %s",
      npins, npins<4 ? "crystal" : "repeating cell from cfg");

  VV(box.Lh,=0.5*box.L)
  box.cq=Sqr(box.cutoff);
  measure=1;

  if (npins<=0) {
    /* automatic selection of the best crystal structure */
    int ipins,Nv,Nmax=1;

    prt("cfginit: auto selection");

    loopto (ipins,1,3)
      if ( (Nv=Nvert(nL,ipins,"tried"))*(ipins>1) < Nmax ) {
        Nmax=Nv;
        npins=ipins; } }

  if (npins>3) {
    char fn[16];
    FILE *wat;
    char line[256],*ch;
    int nmol;

    nscfg=molec[No.N-1].ns;
    sprintf(fn,"cfg%d.%d",nscfg,npins);
    if ( !(wat=fopen(fn,"rt")) ) ERROR(("cfg file %s not found",fn))
    if (!fgets(line,256,wat)) ERROR(("%s: no data",fn))

    ch=line; if (*ch=='#') ch++; if (*ch==' ') ch++;
    if (sscanf(ch,"ns=%d nmol=%d L=%lf\n",&nscell,&nmol,&Lcell)!=3)
      ERROR(("%s: format 1st line",fn))
    put3(nscell,nmol,Lcell)
    if (nscell!=nmol*nscfg)
      ERROR(("%s: nscell=%d ns=%d nmol=%d inconsistent",fn,nscell,nscfg,nmol))
    alloc(cryst[4],sizeof(cryst_t)-sizeof(vector) + nscell*sizeof(vector));
    cryst[4]->ns=nmol;
    loop (i,0,nscell) {
      do {
        if (!fgets(line,256,wat)) ERROR(("%s: not enough data",fn))
        } while (line[0]=='#');

      if (5!=sscanf(line,"%d%d" realfmt realfmt realfmt,
                    &nmol,&nmol,
                    cryst[4]->r[i], cryst[4]->r[i]+1, cryst[4]->r[i]+2))
      ERROR(("%s line %d: format",fn,i)) }

    fclose(wat);
    npins=4; }

  n=Nvert(nL,npins,"selected");

  if (npins==4 && fabs(Lcell/cellsize-1)>1e-13) {
    VO(box.L,*=Lcell/cellsize)
    prt("box rescaled %.14g times to match elementary cell size",Lcell/cellsize);
    cellsize=Lcell; }

  nscell=cryst[npins]->ns;
  ralloc2Darrayzero(vertexused,nscell,PROD(nL));

  prt("cell size=%f (or more),  %d sites/cell, %d molecules on %d vertices",
     cellsize,nscell,No.N,n);

  n=0;
  loop (nx,0,No.N) {
    molecule_t *mn=molec+n;
    int ns=mn->ns;

    if (spec[mn->sp]->config)
      ERROR(("config in a ble-file: crystal initialization not supported"))

   again:
    do {
      loop (k,0,DIM) ir[k]=irnd(nL[k]);
      inscell=irnd(nscell);
    } while (vertexused[inscell][ivu=(ir[0]*nL[1]+ir[1])*nL[2]+ir[2]]);

    vertexused[inscell][ivu]=1;

    if (mn->sp != sp) {
      /* new species */
      if (cc) free(cc);
      cc=initr(sp=mn->sp); /* newly allocated */
      c=cc;
      Scorrect(mn,c,rpv/*dummy*/,1e-7); }

    r=rof(mn,rp);

    if (npins==4 && n>=FROM) {
      if (ns!=nscfg) ERROR(("species %d with %d sites is not in cfg (%d sites)\n\
*** (wrong option -j ?)",sp,ns,nscfg))
      loop (i,0,ns) VV(r[i],=cryst[npins]->r[ns*inscell+i]) }
    else {
      int i;

      if (n>=No.rotatefrom) rndorientation(c,ns);

      if (npins==4 && n<FROM)
        loop (i,0,ns) VV(r[i],=c[i])
      else
        loop (i,0,ns) loop (k,0,3)
          r[i][k]=c[i][k]+(box.L[k]/nL[k])*cryst[npins]->r[inscell][k]; }

    loop (i,0,ns) loop (k,0,DIM) r[i][k]+=(box.L[k]/nL[k])*ir[k];

    /* molecule has been inserted: what about overlaps? */
    insertkey=1;
    {
      int m,ton=min(FROM,n);
      double Esum=0;

      loop (m,0,ton) {
        En.el=0;
        Esum += (*pot[mn->sp][molec[m].sp]) (
              rof(mn,rpf),rof(molec+m,rpf), /*dummy*/
              mn, molec+m, rp ) + En.el*initqf;
        if (Esum > Emax) {
          /* cannot insert - overlap */
          insertkey=0; break; } } }

    if (!insertkey && n<FROM) {
      static int count,maxcount=10;

      /* overlap detected between first FROM molecules - trying again */
      vertexused[inscell][ivu]=0;

      if (++count>maxcount) {
        WARNING(("cannot insert first %d molecules within %d trials - still trying",FROM,maxcount))
        maxcount*=2; }
      goto again; }

    n += insertkey;
  } /* nx,n */

  if (n!=No.N) {
    int dif=No.N-n;

    WARNING(("%d molecules removed (old=%d new=%d)\n\
*** The re-initialization is not compatible with equalize.cfg.\
",dif,No.N,n))
    if ( (spec[nspec-1]->N-=dif) < 0 ) Error("");
    prt(">>>>> to be included to %s :",Fn("def"));
    loop (i,0,nspec) prt_(" N[%d]=%d",i,spec[i]->N);
    _n _n
    No.N=n;
    initNo();
    dif *= spec[nspec-1]->ns*sizeof(vector);
    loop (i,0,9) if (cfg[i]) cfg[i]->size-=dif; }

  Maxwell(FROM,No.N,2);

  if (cc) free(cc);

  ivu=0;
  loop (n,0,nL[0]*nL[1]*nL[2])
    loop (inscell,0,nscell) ivu += !vertexused[inscell][n];
  prt("%d vacancies in the lattice",ivu);

  if (cryst[4]) free(cryst[4]);

  release(vertexused);

  sdszero(cfg[2]);
  normalize(-1);

#ifdef SLAB
  if (wall.n) {
    rescalecfg(cfg[0],RESCALE_Z|RESCALE_CM,wall.z[1]-wall.z[0],NULL);
    loop (n,0,No.N) {
      molecule_t *mn=molec+n;
      vector *r=rof(mn,cfg[0]->rp);
      int ns=mn->ns,i;
      loop (i,0,ns) r[i][2]+=wall.z[0]*box.L[2]; }
    prt("configuration rescaled to [%g,%g]*Lz",wall.z[0],wall.z[1]); }
#endif /*# SLAB */

  removedrifts(1);
  initqf=constrainterror(cfg[0],cfg[1]);
  prt("(v)constraint error=%g (%g)",initqf,vconstrainterror);
  cfg[0]->dep=0;
  prt("<<< cfg initialized >>>");
} /* initcryst */

#ifndef FREEBC
static double xrnd(int center) /*************************************** xrnd */
/*
   returns a random number:
     center=0:  uniform in [0,1)
     center=1:  in [0.1,0.9) enhanced at center
     center=-1: in [0,1) enhanced outside center
*/
{
  if (center) {
    double x=rndcos();
    x*=1+Pow6(x);
    if (center>0) return x/5+0.5;
    else return x/4+(x<0); }
  else
    return rnd();
}
#endif /*# FREEBC */

void initcfg(double pins,double Emax,int nplb,int slab_geom) /****** initcfg */
/***
  pins = insertion probablility:
    pins<0: increase density if prob(insert)>sqrt(|pins|)
            decrease density if prob(insert)<|pins|
    pins>0: decrease density if prob(insert)<|pins|
            accepted if inserted
    pins=0: try to insert forever
    note: also sticky increased if density decreases (if -o-# and FREEBC)
  Emax = energy limit for one molecule insertion
  nplb = for playback while MC initializing
  slab_geom = slab.geom (see the manual, SLAB only)

  An attempt is made to fill the simulation box with randomly
  distributed molecules.  If it does not succeed, a new attempt is
  made with a lower density.  If it does but the ratio of successful
  insertions is higher than pins (i.e. the system is too dilute), a
  new attempt is made with a higher density; this does not apply if
  pins<=0.
  (An insertion of a molecule is considered unsuccessful if the
  interaction energy of the molecule being inserted with a certain set
  of other molecules (taken in a loop) exceeds Emax.
  An attempt to fill the box is considered unsuccessful if the
  probability of successful insertions falls below pins^2.)

  NOTE the long-range k-space part of elst interaction energy is not
  taken into account while calculating energy.  The configuration
  produced by this procedure is still dilute (even if pins is low) and
  is far from equilibrium Hint: for pins=0 all probability checks are
  disabled

  WARNINGS: central force (center.K, center.r0) IS NOT included in this
            initializer, but it is included in MC (for FREEBC only!)
            WALL potential IS included in the initializer, but NOT in MC
***/
{
  double swell=0.5;
  double shrink=pow(swell,-0.6180339887498948/*golden ratio*/);
  double initqf=1.0; /* OLD: sqr(option('q')*0.01) */
  int n,m,nins,sp= -1,maxns=0,nmod;
  vector *c=NULL,*cc=NULL;
  vector *rp=cfg[0]->rp,*rpv=cfg[1]->rp,*rpf=cfg[2]->rp;
  Epairreal **Epair=NULL,EE;
  int cfgindex=0;
  int countagain=0;

  // REMOVED in V2.9f
  // #ifdef FREEBC
  //   /* sticky potential, e.g., water around protein */
  //   if (center.r0<0) {
  //     sticky=Sqr(center.r0);
  //     if (pins<0)
  //       WARNING(("pins<0 ignored with the sticky potential (center.r0<0)"))
  //     pins=fabs(pins);
  //     prt("sticky potential for MC, range=%d",-center.r0);
  //     swell=0.8;
  //     center.r0=0; }
  // #endif

  box.cq=Sqr(box.cutoff);
  measure=1;

  if (MC) {
    alloc(Epair,No.N*sizeof(Epair[0]));
    loop (n,1,No.N) alloc(Epair[n],n*sizeof(Epair[0][0])); }

 TRYAGAIN:

  if (countagain++>100) ERROR(("more than 100 unsuccessful box size trials"))

  nins=0;
  VV(box.Lh,=0.5*box.L)
  nmod=No.N/8+1;

  loop (n,simils.frommol,No.N) {
    molecule_t *mn=molec+n;
    int ns=mn->ns;
    double Esum;
    int config=spec[mn->sp]->config;
#ifdef ANCHOR
    int insanchor=-DIM; /* needed by ANCHOR only */
#endif /*# ANCHOR */

    if (n) if (simils.frommol || (n%nmod==0)) {
      nmod=(No.N-n)/4+1;
      prt("inserting %d prob=%f",n,(n-simils.frommol)/(nins+1e-9)); }

    if (config) {
      if (mn->sp != sp) cfgindex=0;
      sp=mn->sp;
      c=initr(sp)+cfgindex*ns; /* reference */
      cfgindex++; }
    else if (mn->sp != sp) {
      /* new species */
      if (n>=FROM) Max(maxns,ns)
      if (cc) free(cc);
      cc=initr(sp=mn->sp); /* newly allocated */
      c=cc; }
    /* Scorrect(mn,c,rpv,1e-7*(init==3)); now init=4=slab */
    Scorrect(mn,c,rpv/*dummy*/,1e-7);

  INSERT: {
      int i;
      vector *r=rof(mn,rp);

      if (!config && n>=No.rotatefrom) rndorientation(c,ns);
      copy(r,c,ns*sizeof(vector));
      if (fabs(pins)*nins++>n) {
        prt("L = %f %f %f, nins=%d: low insert prob for n = %d", VARG(box.L),nins,n);
        VO(box.L,/=swell)
        if (box.L[0]+box.L[1]+box.L[2]>1e6) ERROR(("box has swelled too much"))
#ifdef FREEBC
        if (sticky!=0) {
          sticky/=Sqr(swell); /* typical distance *=swell */
          prt("new sticky range=%f",sqrt(sticky)); }
#endif /*# FREEBC */
        goto TRYAGAIN; }

      if (!config) loop (i,0,DIM) {
        int j;
#ifdef GOLD
        double dz=0;
#endif /*# GOLD */

#ifdef FREEBC /* centered Gauss */
        double inL=rndgauss()*box.L[i]*(n!=0);
#else /*# FREEBC */
#  ifdef SLAB
        double inL = xrnd((init==4 && i>=(slab_geom%3))*(slab_geom/3?-1:1))*box.L[i];
        if (i==2 && wall.n) inL=(wall.z[1]-wall.z[0]-0.02)*inL+(wall.z[0]+0.01)*box.L[2];
#  else /*# SLAB */
        /* old meaning (only simple slab) */
        double inL = xrnd(init==4 && i==2)*box.L[i];
#  endif /*#!SLAB */
#endif /*#!FREEBC */

#ifdef ANCHOR
        double d0=0.1*insanchor/(DIM*ns);
        /*
          anchored sites/CM are first placed at their position
          if not successful, in the neighborhood
        */
        { double d=0;

          if (mn->anchor & ANCHOR_c)
            inL=mn->r0[i]; /* CM - init. r is already centered */
          else
            if (mn->anchor & ANCHOR_s) inL=mn->r0[i]-r[mn->iaxis[0]][i];
          if (insanchor++>=DIM*ns) {
            d=rndgauss()*d0;
            if (i==0 && (insanchor%10==0))
              prt("%d-th attempt to insert mol[%d] in range %g",
                  insanchor,n,d0);
            if (d>=box.Lh[i]) d=box.Lh[i];
            inL+=d; }
        }
#endif /*# ANCHOR */

#ifdef GOLD
        if (i==2) {
          inL=fabs(inL)*wall.scalez;
          if (n<FROM) {
            double m=9e9;
            loop (j,0,ns) Min(m,r[j][2])
              inL=wall.minz-m+1e-5; } }
      tryagain:
#endif /*# GOLD */
        loop (j,0,ns) r[j][i] += inL;
#ifdef GOLD
        if (i==2) {
          double m=9e9;
          loop (j,0,ns) Min(m,r[j][2])
            if (m<wall.minz) {
              inL=-inL+fabs(gaussh(L))+dz; dz+=1; goto tryagain; } }
#endif /*# GOLD */
        }

#ifdef SLAB
      if (wall.is) {
        Esum=mol_wall(rof(molec+n,rpf) /*dummy*/,mn,rp);
        if (Esum > Emax) goto INSERT; }
      else
        Esum=0;
#else /*# SLAB */
      Esum=0;
#endif  /*#!SLAB */

#ifndef FREEBC
      if (center.on) {
        int i,k;

        if (center.on&2) loop (i,0,ns) loop (k,0,DIM) {
          double dr=r[i][k]-box.Lh[k];

          if (dr<0) dr+=center.r0[k]; else dr-=center.r0[k];

          Esum+=center.K[k]*Sqr(dr); }

#  ifdef SLAB
        if (center.on&1) loop (i,0,ns) loop (k,0,NCENTER) {
          if (slab.nn[k]<=0) break;
          if (n<slab.nn[k]) {
            double dz=r[i][2]-(box.Lh[2]+slab.z[k]);

            if (fabs(dz)>slab.z0[k]
                && (slab.ns[k]==0
                    || (slab.ns[k]>0 && i<slab.ns[k])
                    || (slab.ns[k]<0 && i>=slab.ns[k]))) {

              if (dz<0) dz+=slab.z0[k]; else dz-=slab.z0[k];

              if (slab.z1[k]) {
                /* UNFINISHED */
                static int pass;
                if (!pass++) WARNING(("slab bias forces (slab.z1) are not (yet) included in the initializer")) }
              else
                Esum +=slab.Kz[k]*Sqr(dz); } } }
#  endif /*# SLAB */
      }


      if (Esum > Emax) goto INSERT;
#endif /*# FREEBC */

      if (!config) loop (m,0,n) {
        En.el=0;
        EE= (*pot[mn->sp][molec[m].sp]) (
            rof(mn,rpf),rof(molec+m,rpf), /*dummy*/
            mn, molec+m, rp ) + En.el*initqf;
#ifdef FREEBC
        if (sticky!=0 && m<FROM)
          EE+=stickypot(mn->ns,rof(mn,rp),molec[m].ns,rof(molec+m,rp))*T;
#endif /*# FREEBC */
        Esum+=EE;

        if (MC) {
          Epair[n][m]=EE;
          if ( (n<FROM || m<FROM) && EE>Emax ) goto INSERT; }
        else
          if (Esum > Emax) goto INSERT; }
      } /* INSERT: */

    } /* n */

  prt("L = %f %f %f  prob(insert) = %.3f", VARG(box.L),No.N/(double)nins);
  if (pins<0 && nins*sqrt(-pins)<No.N) {
    VO(box.L,/=shrink)
    goto TRYAGAIN; }

  if (MC) {
    doMC(Epair,maxns,Emax,initqf,nplb);
    loop (n,1,No.N) free(Epair[n]);
    free(Epair); }

  if (!simils.frommol) Maxwell(FROM,No.N,2);

  if (cc) free(cc);

  sdszero(cfg[2]);

  normalize(-1);
  removedrifts(1);
  swell=constrainterror(cfg[0],cfg[1]);
  prt("(v)constraint error=%g (%g)",swell,vconstrainterror);
  cfg[0]->dep=0;
  prt("<<< cfg initialized >>>");
} /* initcfg */

static int OpenCfg(int n,char *mode) /* ---------------------------- OpenCfg */
/* n>=0: extension is n (in decimal)
   n<0: extension is cfg */
{
  char suff[16]="cfg"; /* [8] would be enough, but [16] suppresses verbose warnings */
  char *aux=simils.simname;

  if (n>=0) {
    if (*mode=='r' && simils.cfgname) simils.simname=simils.cfgname;
    if (!simils.simname) ERROR(("OpenCfg: no NAME to read NAME.1, NAME.2,..."))
    sprintf(suff,"%d",n); }

  if (mode[0]=='w') backup(suff);
  VarOpen(Fn(suff),mode);

  if (n>=0) simils.simname=aux;

  return n<0 ? option('m') : option('r');
}

/* see also loadcfg0.c for a very old version of loadcfg */

static int CFGKEY = 1 /* New version */
#ifdef POLAR
                  + 2
#endif /*# POLAR */
                  ;

/* new key-based record; must be less than 7 doubles */
struct rec_s {
  int key;
  int intval;
  vector vecval;
};

/* removed: SCALEDMOL +4 "old" new +8 +16*/

void static swapxyz(real *x,int swapload) /************************* swapxyz */
{
  real sw;

  switch (swapload) {
    case 1: sw=x[0]; x[0]=x[1]; x[1]=sw; break; /* x <-> y */
    case 2: sw=x[2]; x[2]=x[1]; x[1]=sw; break; /* y <-> z */
    case 3: sw=x[2]; x[2]=x[0]; x[0]=sw; break; /* x <-> z */
    case 4: sw=x[2]; x[2]=x[1]; x[1]=x[0]; x[0]=sw; break; /* z -> y -> x -> z */
    case 5: sw=x[2]; x[2]=x[0]; x[0]=x[1]; x[1]=sw; break; /* z -> x -> y -> z */
    default: ERROR(("bad swapload=%d (option -[-%d0000)",swapload,swapload)) }
}

static int swapload,replicate=1;
struct load_s load= { {1,1,1} };

int initreplicate(void) /************************************* initreplicate */
{
  int i;

  if (!simils.Ns) ERROR(("initreplicate: number of sites not initialized"))

  VV(simils.pc,=load.n)
  swapload=load.tr;

  loop (i,0,DIM) if (simils.pc[i]<=0) {
    ERROR(("load.n[%d]=%d: bad cell replicate factor",i,simils.pc[i]))
    simils.pc[i]=1; }

  replicate = PROD(simils.pc);

  if (simils.Ns!=No.s) {
    prt("Number of sites loaded = %d, number of sites specified = %d",simils.Ns,No.s);
    if (simils.N)
      prt("Number of molecules loaded = %d, number of molecules specified = %d",simils.N,No.N); }

  if (!simils.N) {
    simils.N=(long)No.N*simils.Ns/No.s;
    prt("Unknown number of molecules loaded, calculated %d will be used",simils.N); }

  if ((No.N==simils.N) != (No.s==simils.Ns))
    ERROR(("molecules/sites numbers inconsistent (composition changed?)\n\
*** No.N=%d simils.N=%d  No.s=%d simils.Ns=%d", No.N,simils.N,No.s,simils.Ns))

  if (simils.Ns>No.s) {
    if ((load.N&2)==0) ERROR(("The configuration loaded is longer than requested:\n\
*** loaded simils.Ns=%d, requested=No.s=%d\n\
*** Specify load.N&2 to suppress this message and omit the extra molecules.",simils.Ns,No.s))

    WARNING(("The configuration loaded is longer than requested,\n\
*** the extra molecules/sites will be omitted\n\
*** this is safe only for the last species omitted")) }

  if (simils.Ns<No.s && replicate==1) {
    if (nspec>1) ERROR(("The configuration loaded is shorter than requested.\n\
*** Sorry, I cannot fix this for more than 1 species."))

    if (load.N&2) ERROR(("The configuration loaded is shorter than requested.\n\
*** Specify load.N&1 to insert missing molecules.")) }

  if (replicate!=1) prt("Load with periodic copies: replicate=%d*%d*%d=%d",simils.pc[0],simils.pc[1],simils.pc[2],replicate);
  if (replicate<1) ERROR(("internal: replicate<1"))

  simils.changecfg=simils.Ns!=No.s || replicate>1;

  return replicate;
}

void replicatecfg(void) /************************************** replicatecfg */
/* make periodic copies */
{
  int i,j,ii,sp,off=0,chunk;

  if (swapload) {
    static char *info[]={"","x <-> y", "y <-> z","x <-> z","z -> y -> x -> z","z -> x -> y -> z"};

    swapxyz(box.L,swapload);
    VV(box.Lh,=0.5*box.L)
    loop (i,0,simils.to)
      loop (j,0,simils.Ns) swapxyz(simils.cfg[i]->rp[j],swapload);
    prt("box rotated: swapload=%d (%s)",swapload,info[swapload]);
    prt("box changed to: L=[%.14g %.14g %.14g]",box.L[0],box.L[1],box.L[2]); }

  if (!simils.changecfg) return;

  if (replicate>1) {
    /* replicate cell */
    prt("replicating cell %d times:\n\
  %d molecules loaded -> %d molecules requested\n\
  %d sites loaded -> %d sites requested",
        replicate, simils.N,No.N, simils.Ns,No.s);
    if (No.N!=simils.N*replicate || No.s!=simils.Ns*replicate)
      ERROR(("cannot replicate cell: not a multiple: total %d molecules expected\n\
*** (or mixture composition changed)",simils.N*replicate))
    En.tot*=replicate;

    loop (i,0,simils.to) {
      off=0;
      loop (sp,0,nspec) {
        chunk=spec[sp]->ns*spec[sp]->N/replicate;
        if (spec[sp]->ns*spec[sp]->N!=chunk*replicate)
          ERROR(("species %d: target %d sites is not %dx replicated chunk=%d\n\
*** (or mixture composition changed)",sp,spec[sp]->ns*spec[sp]->N,replicate,chunk))
        loop (j,0,replicate) {
          vector dr;

          if (i)
            VO(dr,=0)
          else {
            dr[0]=(j/(simils.pc[1]*simils.pc[2]))*box.L[0];
            dr[1]=(j/simils.pc[2])%simils.pc[1]*box.L[1];
            dr[2]=j%simils.pc[2]*box.L[2]; }
          loop (ii,0,chunk)
            VVV(cfg[i]->rp[off*replicate+j*chunk+ii],=simils.cfg[i]->rp[off+ii],+dr) }
        off+=chunk; }
      if (off*replicate!=No.s) ERROR(("internal: %d %d",off,No.N)) }

    prt("%d x %d x %d periodic copies made",simils.pc[0],simils.pc[1],simils.pc[2]);
    VV(box.L,*=simils.pc)
    VV(box.Lh,=0.5*box.L)
    replicate=1; }

  else if (simils.Ns>No.s) {
    /* truncate too long configuration */
    int size=cfg[0]->size;

    put2(simils.size,size)

    if (simils.size<=size) ERROR(("internal: sizes"))

    loop (i,0,simils.to) {
#ifdef POLAR
      int hdr=(char*)(cfg[i]->rp)-(char*)(cfg[i]);
      int newblk=(size-hdr)/2;
      int oldblk=(simils.size-hdr)/2;

      copy(cfg[i],simils.cfg[i],size); /* more than needed */
      copy((char*)(cfg[i])+hdr+newblk,(char*)simils.cfg[i]+hdr+oldblk,newblk);
#else /*# POLAR */
      copy(cfg[i],simils.cfg[i],size);
#endif /*#!POLAR */
      cfg[i]->size=size; } }

  else if (simils.Ns<No.s) {
    /* insert extra molecules (postpone) */
    int size=cfg[0]->size;
    int hdrsize=sizeof(ToInt)-sizeof(vector);
    int offset=simils.size-hdrsize;
    int n;

    Fn("cfg"); /* set lastFn */

    put2(simils.size,size)
    if (simils.size>=size) ERROR(("internal: sizes"))

    loop (n,0,No.N) if (molec[n].ir==offset) simils.frommol=n;
    if (!simils.frommol)
      ERROR(("The configuration in %s is truncated in the middle of a molecule.",lastFn))
    else
      WARNING(("The configuration in %s is shorter than requested,\n\
*** the missing molecules from %d will be inserted at random places\n\
*** with energy limit Emax.",lastFn,simils.frommol))

    loop (i,0,simils.to) {
#ifdef POLAR
      int hdr=(char*)(cfg[i]->rp)-(char*)(cfg[i]);
      int newblk=(size-hdr)/2;
      int oldblk=(simils.size-hdr)/2;

      copy(cfg[i],simils.cfg[i],simils.size); /* more than needed */
      copy((char*)(cfg[i])+hdr+newblk,simils.cfg[i]+hdr+oldblk,newblk);
#else /*# POLAR */
      copy(cfg[i],simils.cfg[i],simils.size);
#endif /*#!POLAR */
      cfg[i]->size=size; /* restored: have been overwritten */ } }
  else
    ERROR(("internal"))

  loop (i,0,simils.to) free(simils.cfg[i]);
  simils.changecfg=0;
}

void loadcfg(int n,double *sigvdW) /******************************* loadcfg */
/***
   Configuration and some additional information (No.N, L,...) are loaded
     n<0: from SIMNAME.cfg
     n>=0: from PLBNAME.n (changed in V3.3s)
   Should be called in two passes:
     1st: sigvdW!=NULL, info is read, the file remains opened
     2nd: sigvdW==NULL, the configuration is read, the file is closed
   Consistency of data is checked.
   sigvdW = pointer to selected van der Waals diameter, see tau.i and tau.R
            bad habit: also denotes the pass
   NB: endian change and float/double conversion no longer supported,
       real=double now assumed
***/
{
  static int cfgkey,i,j,cfgol[32],m,persite;
  int calcns=0,locnspec;
  static real cfgh;
#ifdef POLAR
  static double cfgtaudip;
#endif /*# POLAR */
  static struct varfile_s usedVarFile;

  if (sigvdW) { /***** 1st pass *****/
    if (usedVarFile.file) ERROR(("loadcfg: file already opened in 1st pass"))

    if (n<0 || option('v')&4) prt_("loadcfg: 1st pass (read header), cfgkey=");

    simils.to=OpenCfg(n,"r");
    VarRead(&cfgkey,sizeof(cfgkey));
    if (n<0 || option('v')&4) prt("%d",cfgkey);
    if (!(cfgkey&1)) ERROR(("Sorry, I cannot load the old (<V2.7) format of %s.\n\
*** Use  cfgconv ",lastFn,simils.simname))
    VarRead(cfgol,sizeof(cfgol));
    if (n<0) loop (i,0,32) if (optionlist[i]!=cfgol[i]) {
      prts_("option"); prtoption(i,cfgol)
      prts_(" read from file, "); prtoption(i,optionlist)
      prts(" in this run");
      if (init==0 && strchr("cejkmpquwxy[]^",i+'`'))
        WARNING(("option -%c has changed and init=0",i+'`')) }

    VarRead(&locnspec,sizeof(locnspec));

    if (locnspec==nspec) {
      int dif=0;

      if (n<0 || option('v')&4) prt("%s: locnspec=%d read (OK), simils.specsize=%d",lastFn,locnspec,simils.specsize);

      if (simils.specsize<locnspec) {
        if (simils.spec) free(simils.spec);
        allocarray(simils.spec,simils.specsize=locnspec); }
      VarRead(simils.spec,sizeof(struct spec_s)*nspec);

      loop (i,0,min(locnspec,nspec)) {
        if (simils.spec[i].N!=spec[i]->N || simils.spec[i].ns!=spec[i]->ns) dif++;
        calcns+=simils.spec[i].N*simils.spec[i].ns; }

      if (dif || option('v')&4) {
        prt("%s: inconsistent number of species or sites:",lastFn);
        header("i   loaded: N   ns  specified: N   ns");
        loop (i,0,max(locnspec,nspec)) {
          prt("%d       %3d    %4d    %3d    %4d",
              i,i<locnspec?simils.spec[i].N:0,i<locnspec?simils.spec[i].ns:0,
              i<nspec?spec[i]->N:0,i<nspec?spec[i]->ns:0); }
        header(""); } }
    else {
      prt("%s: WARNING: locnspec=%d read disagrees with system, ble-file number applies",lastFn,locnspec);
      if (load.N&4) WARNING(("This disagreement for load.N=4 is unexpected and may cause a crash"));
      locnspec=nspec; }

    if (init<2 && option('e')!=cfgol['e' & 31]) {
      if (option('e')) ERROR(("CP record length changed (wrong option -e)"))
      else prt("-e disabled"); }

    if (cfgkey&8) ERROR(("%s: cfgkey & 8 no longer supported, use cfgconv",lastFn))

#ifdef POLAR
    cfgtaudip=tau.dip;
#endif /*# POLAR */

    if (VarFile.size!=sizeof(struct rec_s))
      ERROR(("VarFile.size=%d (%d expected) in %s\n\
*** possible causes: real!=double, different nspec (check %s.ble)",
             (int)VarFile.size,(int)sizeof(struct rec_s),lastFn,simils.sysname))

    while (VarFile.size==sizeof(struct rec_s)) {
      struct rec_s rec;

      VarRead(&rec,sizeof(rec));
      prt_("loadcfg key=%d: ",rec.key);
      switch (rec.key) {
        case 1: simils.N=rec.intval;
          if (simils.N!=No.N)
            prt("%s contains %d molecules but %d specified (will be solved later)",
                lastFn,simils.N,No.N);
          t=rec.vecval[0];
          cfgh=rec.vecval[1];
          En.tot=rec.vecval[2];
          break;

        case 2:
          simils.Ns=rec.intval;
          VV(box.L,=rec.vecval)
          prt(" Ns=%d L=%g %g %g loaded",simils.Ns,VARG(box.L));
          break;

        case 3:
          /* reserved */
          break;

        case 4:
          *sigvdW=rec.vecval[0];
          if (*sigvdW) prt_(" sigvdW=%g",simils.Ns,*sigvdW);
          prt(" (thermostat=%d tau.T=%g tau.P=%g ignored) ",rec.intval,rec.vecval[1],rec.vecval[2]);
          break;

#if defined(SLAB) && SLAB & 2
        case 5:
          if (cleave.n!=rec.intval) prt(" cleave.n=%d (ignored)",rec.intval);
          cleave.z[0]=rec.vecval[0]; cleave.z[1]=rec.vecval[1];
          loop (i,0,cleave.n) prt("cleave[%d]=%.9f loaded",i,cleave.z[i]);
          break;
#endif /*# defined(SLAB) && SLAB & 2 */
#ifdef POLAR
        case 17:
          cfgtaudip=rec.vecval[0];
          prt(" cfgtaudip=%g",cfgtaudip);
          break;
#endif /*# POLAR */
        default:
          WARNING(("%s: record key=%d ignored",lastFn,rec.key)) }
    }

#ifndef FREEBC
    putv(box.L)
    loop (i,0,DIM) if (box.L[i]<=0) {
      WARNING(("L%c = 0 loaded: replaced by 1",i+'x'))
      box.L[i]=1; }
#endif /*# FREEBC */

    VV(box.Lh,=0.5*box.L)
    box.V=PROD(box.L);

    if (!iscube()) prt("NOT CUBE loaded: L=[%.14g %.14g %.14g]",box.L[0],box.L[1],box.L[2]);

    /* next record = configuration */
    simils.size=VarFile.size;
    simils.Neq=(VarFile.size-2*sizeof(int))/sizeof(real); // realsize
    put3(VarFile.size,simils.Neq,simils.Ns)

    /* bug fixed by Tomas Trnka 2/2011 */
    persite = (cfgkey & 2) ? 2 : 1; /* vectors per site (polar => atom + Drude particle) */

    if (persite==2) prt("%s: POLAR cfg",lastFn);

    if (simils.Ns) {
      if (calcns!=simils.Ns) {
        if (calcns)
          WARNING(("number of sites calculated from species (%d) differs from Ns=%d\n",calcns,simils.Ns))
        else
          prt("\
WARNING: no species/site information found in %s\n\
         (OK for cfgs obtained using plb2cfg, cfgconv or similar)",lastFn); }
      if ((simils.Neq-RPOFFSET)/DIM!=simils.Ns*persite)
        ERROR(("Neq=%d Ns=%d: number of sites inconsistent with record length %ld",
               simils.Neq,simils.Ns,(long)VarFile.size)) }
    else {
      simils.Ns=(simils.Neq-RPOFFSET)/DIM/persite;
      prt("NOTE: old cfg format read, number of sites %d calculated from record length", simils.Ns); }

    usedVarFile=VarFile; }

  else { /****** 2nd pass ******/
    struct { int size,padd; char c[1]; } *aux;

    prt("loadcfg: 2nd pass (read configuration)");

    VarFile=usedVarFile;
    if (!VarFile.file) ERROR(("loadcfg: file %s not opened in 2nd pass", lastFn))

    Min(simils.to,cfgol['m' & 31])

    if ( (cfgkey&2) != (CFGKEY&2) ) {
      prt("%sPOLAR conversion",cfgkey&2 ? "POLAR->NON" : "NONPOLAR->");
      if (No.N != simils.N
          || replicate>1) ERROR(("cannot combine with N change or replicate"))
      m=cfg[0]->size;
      j=VarFile.size;
      alloc(aux,j);
      loop (i,0,simils.to) {
        if (j!=VarFile.size) ERROR(("%s: records of different sizes",lastFn))
        VarRead(aux,j);
        copy(cfg[i],aux,min(m,j));
        cfg[i]->size=m; }
      free(aux); }

    else {
      /* no POLAR/nonPOLAR conversion; N may differ, though */

      loop (i,0,simils.to) {
        if (simils.changecfg) sdsalloc(simils.cfg[i],simils.size)
        else simils.cfg[i]=cfg[i]; }

      loop (i,0,simils.to) VarRead(simils.cfg[i],simils.cfg[i]->size); }

    if (h!=cfgh) {
      double q=1;

      if (fabs(cfgh/h-1)>1e-5) {
        prt("WARNING integration step changed %g -> %g",cfgh,h);
        if (!init) WARNING(("init==0 and h changed")) }
      loop (i,1,option('m')) {
        q*=h/cfgh;
        if (!simils.cfg[i]) {
          if (simils.changecfg) sdsalloc(simils.cfg[i],simils.size)
          else simils.cfg[i]=cfg[i];
          if (!simils.cfg[i]) ERROR(("internal: loadcfg: cfg[%d] not initialized",i)) }
        /* BUG (trapped above): rescalecfg uses No.N, but smaller needed if replicate */
        if (simils.changecfg && i)
          WARNING(("implementation limitation: cfg[%d] not rescaled\n\
*** because size of configuration changed",i))
        else
          rescalecfg(simils.cfg[i],RESCALE_XYZ|RESCALE_H,q,NULL); } }

    if (option('m')==2 && cfgol['m' & 31]>2) {
      /* special: more accurate Gear->Verlet(Shake) conversion */
      if (!simils.cfg[2]) {
        if (simils.changecfg) sdsalloc(simils.cfg[2],simils.size)
        else simils.cfg[2]=cfg[2];
        if (!simils.cfg[2]) ERROR(("internal: loadcfg: cfg[2] not initialized")) }
      VarRead(simils.cfg[2],simils.cfg[2]->size);
      if (h!=cfgh)
        rescalecfg(simils.cfg[2],RESCALE_XYZ|RESCALE_H,Sqr(h/cfgh),NULL);
      loop (j,0,simils.Neq) {
        (&simils.cfg[1]->logs)[j]-=(&simils.cfg[2]->logs)[j]; }
      prt("Gear->Verlet conversion corrected"); }

#ifdef POLAR
    if (option('p')/10%10==0) {
      /* Car-Parrinello-like */
      if (cfgol['p'&31]/10%10>0)
        /* SCF->CP conversion */
        Maxwell(FROM,No.N,-2);
      else if (h!=cfgh || tau.dip!=cfgtaudip) {
        double q=(h/cfgh);

        if (tau.dip!=0 && cfgtaudip!=0) q*=cfgtaudip/tau.dip;

        if (q!=1) {
          /* this correction assumes that cfg[0]=r(t), cfg[1]=r(t-h),
             which applies for positions of auxiliary (Drude) charges */
          loop (j,simils.Neq,simils.Neq*2-1)
            (&simils.cfg[1]->logs)[j]=(&simils.cfg[0]->logs)[j]
            +((&simils.cfg[1]->logs)[j]-(&simils.cfg[0]->logs)[j])*q;
          prt("Car-Parrinello-like, change of h or tau.dip: velocity of aux charges recalculated"); } } }
#endif /*# POLAR */

    VarClose();
    memset(&usedVarFile,0,sizeof(usedVarFile));
    /* cfg[i]->nr (re-)CALCULATED because of compatibility with OLD configs */
    /*.....loop (i,0,simils.to) cfg[i]->nr=(cfg[i]->size-2*sizeof(int))/sizeof(real)-1;*/

    if (n<0) {
      removedrifts(1);
      prt("<<< cfg loaded >>> L=[%g %g %g]",box.L[0],box.L[1],box.L[2]); }
    cfg[0]->dep=0;

    loop (j,0,DIM) {
      if (box.L[j]) cfg[0]->lambda[j]=log(box.L[j]);
      /* TO BE TURNED OFF : */
      //      loop (i,1,simils.to) cfg[i]->lambda[j]=0;
    }

  } /* pass==2 */
} /* loadcfg */

void replaceL(vector L) /************************************************** */
/*
  if permitted by load.L, replace box.L by L
*/
{
  int i;

  if (!SUM(simils.changeL)) return;

  loop (i,0,DIM)
    if (simils.changeL[i]) {
      if (abs(simils.changeL[i]-2)>1) ERROR(("simils.changeL[%d]=%d wrong key",i,simils.changeL[i]))
      if ( (simils.changeL[i]==1 && box.L[i]<L[i])
        || (simils.changeL[i]==2 && box.L[i]>L[i])
        || (simils.changeL[i]==3 && box.L[i]!=L[i]) ) {
        prt("L%c=%g (loaded or initialized) changed to %g",i+'x',box.L[i],L[i]);
        box.L[i]=L[i]; } }
}

void savecfg(int n,int4 timekey,double *sigvdW) /******************* savecfg */
/***
    Configuration and some additional information are saved to file.
    n=-1: save SIMNAME.cfg (full configuration)
    n=-2: as above but do not call removedrifts(1) and depend_r(cfg[0],0)
          (good after removing a molecule because the tables are broken)
    SIMNAME.cfg (n<0:  full configuration) or
    SIMNAME.n   (n>=0: see option -r.  n = rounded t/tcfg.
***/
{
  int i,to;
  int cfgkey=CFGKEY;
  char *info="";
  struct rec_s rec;

  if (n>=-1) {
    if ((drift&DRIFT_WHEN)==DRIFT_SAVE) removedrifts(option('v')&64);
    depend_r(cfg[0],0); }

  to=OpenCfg(n,"w");

  put(n)
  
  VarPut(&cfgkey,sizeof(cfgkey));
  VarPut(optionlist,sizeof(optionlist));
  VarPut(&nspec,sizeof(nspec));
  if (simils.specsize!=nspec) {
    if (simils.spec) free(simils.spec);
    allocarray(simils.spec,simils.specsize=nspec);  }
  loop (i,0,nspec) {
    simils.spec[i].N=spec[i]->N;
    simils.spec[i].ns=spec[i]->ns; }
  VarPut(simils.spec,sizeof(struct spec_s)*nspec);

  memset(&rec,0,sizeof(rec));

  rec.key=1;
  rec.intval=No.N;
  rec.vecval[0]=t; rec.vecval[1]=h; rec.vecval[2]=En.tot;
  VarPut(&rec,sizeof(rec));

  rec.key=2;
  rec.intval=No.s;
  VV(rec.vecval,=box.L)
  VarPut(&rec,sizeof(rec));

  rec.key=4;
  rec.intval=thermostat;
  rec.vecval[0]=*sigvdW;
  if (*sigvdW) prt("savecfg: sigvdW=%g",*sigvdW);
  rec.vecval[1]=tau.T;
  rec.vecval[2]=tau.P;
  VarPut(&rec,sizeof(rec));

#if defined(SLAB) && SLAB & 2
  rec.key=5;
  rec.intval=cleave.n;
  rec.vecval[0]=cleave.z[0]; rec.vecval[1]=cleave.z[1];
  loop (i,0,cleave.n) prt("cleave[%d]=%.9f saved",i,cleave.z[i]);
  VarPut(&rec,sizeof(rec));
#endif /*# defined(SLAB) && SLAB & 2 */
#ifdef POLAR
  rec.key=17;
  rec.vecval[0]=tau.dip;
  VarPut(&rec,sizeof(rec));
#endif /*# POLAR */

  loop (i,0,to) VarPut(cfg[i],cfg[i]->size);
  /* could be VarPutSds(cfg[i]), but for back-compatibility kept */

  KeyClose(timekey);
  prt("<<< %d saved %s>>>",n,info);
} /* savecfg */

#if defined(POLAR) && POLAR&32
void initfqcharges(void) /************************************ initfqcharges */
/* assign permanent charges to fluctuating charges
   (of cfg[0]: will propagate later)
*/
{
  int n,i;

  if (fabs(el.epsinf)<1e15)
    ERROR(("el.epsinf=%g for fluctuating charge\n\
*** use 3e33 (infinity) or re-consider the code",el.epsinf))

  loop (n,0,No.N) {
    molecule_t *mn=molec+n;
    siteinfo_t *si=spec[mn->sp]->si;

    loop (i,0,mn->ns) {
      vector *rpol=polarrof(mn,cfg[0]->rp);

      rpol[i][0]=si[i].charge; } }
}
#endif /*# defined(POLAR) && POLAR&32 */

/******************************** playback *********************************/

static struct playback_s {
  int nm;      /* # of molecules recorded */
  int ns;      /* # of sites */
  int nplb;    /* 1=plb, 2=also vlb, 3=also dlb */
  FILE *plb[3];/* [0]: playback file, [1]: velocities, [2]=Drude (abs.pos.) */
  char *ext[3];/* coresponding extension names */
} playback;

typedef float fvector[DIM];

void openplayback(int nm,int append) /************************* openplayback */
/*
  Opens a playback file called SIMNAME.plb
  nm>0: write first nm molecules only
  nm=-1: write whole configuration
  append=1: open at end for appending
  append=0: rewrite
  option('n') is a sum of flags:
    1 = also SIMNAME.vlb (velocities) opened
    2 = plb recorded at time t-h/2,
        to be used with 1 and leap-frog (option('m')==2)
    4 = plb is centered, 2nd float is -6, L is unchanged
    8 = also SIMNAME.dlb (drude charges) opened
  writeplayback() will actually write a frame (frames)
  closeplayback() will close the file(s)
  NB: not used for reading -- see readplayback()
*/
{
  float header[2];
  int i,iplb;

#ifndef POLAR
  if (option('n')&8) ERROR(("option -n8 (export Drude) illegal for nonpolar simulation"))
#endif /*# POLAR */
  if (simils.plbname) {
    prt_("%s changed to ",simils.plbname);
    simils.plbname=NULL; /* NB: simils.cfgname is kept */
    prt("%s%s%s for writing%s%s",Fn("plb"),
        option('n')&1?",vlb":"",
        option('n')&8?",dlb":"",
        option('n')&2?", plb shifted by -h/2":"",
        option('n')&4?", plb centered":""); }

  playback.ns=0;
  playback.nplb=0;
  playback.ext[playback.nplb++]="plb";
  if (option('n')&1) playback.ext[playback.nplb++]="vlb";
  if (option('n')&8) playback.ext[playback.nplb++]="dlb";

  if (nm==0) ERROR(("internal: nm==0"))
  nm=abs(nm);

  Min(nm,No.N)
  loop (i,0,nm) playback.ns+=molec[i].ns;
  playback.nm=nm;

  if (nm==No.N) {
    if (playback.ns!=No.s) ERROR(("internal")) }
  else
    prt("WARNING: only first %d molecules (of %d) recorded in the playback%c",
        nm,No.N,"  s"[playback.nplb]);

  loop (iplb,0,playback.nplb) {
    if (playback.plb[iplb]) ERROR(("internal: playback.%s already opened",playback.ext[iplb]))
    if (append) {
      playback.plb[iplb]=fopen(Fn(playback.ext[iplb]),"rb");
      if (playback.plb[iplb]) {
        if (2!=fread(header,sizeof(float),2,playback.plb[iplb]))
          ERROR(("%s is too short",lastFn))
        if (header[0]!=playback.ns)
          ERROR(("number of sites in %s (ns=%d) does not match ns=%d requested",
            lastFn,(int)header[0],playback.ns))
        fclose(playback.plb[iplb]); }
      else {
        ERROR(("playback file %s is about to be appended but it does not exist\n\
*** a new one will be created (if -n, also .vlb/.dlb EVEN IF ONE EXISTS)\n\
*** (the .cp and .plb/.vlb/.dlb files need not match)",lastFn))
        append=0; } }

    if (!append) backup(playback.ext[iplb]);
    playback.plb[iplb]=fopen(Fn(playback.ext[iplb]),append?"ab":"wb");

    if (!playback.plb[iplb]) ERROR(("open %s",lastFn))

    if (!append) {
      header[0]=playback.ns;
      header[1]=-3; /* forces variable L format of .plb */
      if (option('n')&4) header[1]=-6;
      fwrite(header,2,sizeof(float),playback.plb[iplb]);
      prt("playback file %s (%d sites) created",lastFn,playback.ns); }

    prt("%splayback opened t=%.3f L=[%8.5f %8.5f %8.5f] position=%ld",
        iplb?"vel.":"",t,VARG(box.L),ftell(playback.plb[iplb])); }
}

static double lastt=-1;

void writeplayback(void) /*********************************** writeplayback */
/*
   WARNING: with SHAKE, velocity is shifted by h/2
   use -n2 to shift position by h/2 backwards
*/
{
  vector d;
  fvector r;
  int iplb,i;

  loop (iplb,0,playback.nplb) {

    if (!playback.plb[iplb]) ERROR(("internal: playback.%s (iplb=%d) not opened",playback.ext[iplb],iplb))
    if (!playback.nm) ERROR(("internal"))

    // if (iplb==0) depend_r(cfg[0],0);
    depend_r(cfg[iplb],0); // also velocities (?)

    waitfordiskspace(No.s/40+2);

#ifdef FREEBC
    VO(r,=0)
#else /*# FREEBC */
    VV(r,=box.L)
#endif /*#!FREEBC */
    if (DIM!=fwrite(r,sizeof(float),DIM,playback.plb[iplb])) ERROR(("write playback"))

    loop (i,0,playback.ns) {

      if (playback.ext[iplb][0]=='v')
        /* velocity: Verlet(leap-frog)/SHAKE: v(t-h/2), Gear: v(t) */
        VV(r,=(1./h)*cfg[iplb]->rp[i])
      else {
        /* position r(t) */
        VV(d,=cfg[0]->rp[i])
#ifdef POLAR
        if (playback.ext[iplb][0]=='d')
          VV(d,+=cfg[0]->rp[i+No.s])
#endif /*# POLAR */
        if (option('n')&2)
          /* shift by -h/2: r(t-h/2) */
          VV(d,-=0.5*cfg[1]->rp[i]) /* r(t-h/2) */
        if (option('n')&4)
          /* center in box (for FREEBC not useful) */
          VV(d,-=box.Lh)
        /* convert to float */
        VV(r,=d) }

      if (DIM!=fwrite(r,sizeof(float),DIM,playback.plb[iplb]))
        ERROR(("write playback")) }

    if (fflush(playback.plb[iplb])) ERROR(("write playback"))
    if (!iplb && option('v')&4)
      prt("plb: t=%.3f r[0]=(%8.5f %8.5f %8.5f) L=[%8.5f %8.5f %8.5f]",
          t,VARG(cfg[iplb]->rp[0]),VARG(box.L)); }

  lastt=t;
}

void closeplayback(void) /************************************ closeplayback */
{
  int iplb;

  loop (iplb,0,playback.nplb) {
    char *fn=Fn(playback.ext[iplb]);

    if (playback.plb[iplb]) {
      long ft=ftell(playback.plb[iplb]);
      long frame=ft/(DIM*sizeof(float)*(1+playback.ns));

      underline(iplb?"closing velocity playback file and crash-restart hints"
                    :"closing playback file and crash-restart hints");

      prt("closing %s at position %ld\nlast frame=%ld written at t=%.14g",
          fn,ft,frame,lastt);
      prt("(this has been determined by ftell from really written data)\n\
# SIMPLE RESTART in case of crash during the next sweep (data set ended by ;):");
      prt("  truncate -s %ld %s",ft,fn);
      prt("# OR more slowly with a full backup and plb file check:\n\
# mv %s C~%s\n\
# plbcut C~%s %s 0 1:%ld",fn,fn,fn,fn,frame);
      prt("# + optional partial backup (keep the cut-off part only, remove full backup):\n\
# plbcut C~%s C%ld-%s 0 %ld:2147483647\n\
  rm C~%s",fn,frame+1,fn,frame+1,fn);
      prt("# then, the following test in bash should pass:");
      /*      prt("  if [ `ls -lG %s | cut -d \" \" -f4` = %ld ]\n    then echo OK\n    else echo WRONG SIZE\n  fi",fn,ft);*/
      prt("  if [ `stat -c%%s %s` = %ld ] ; then echo OK ; else echo WRONG SIZE; fi",fn,ft);
      prt("# remove lock-file\n\
  rm %s\n\
#and restart simulation with init=0",Fn("loc"));
        if (fclose(playback.plb[iplb])) ERROR(("close %s\n\
*** hint: try restart.sh",lastFn)) }
    else
      ERROR(("%s: nothing to close",fn))
    playback.plb[iplb]=NULL; }
}

void writeasc(void) /********************************************** writeasc */
/*
  option -l[-]EGFAVC (decimal digits, can be combined)
     - :  writes the configuration to file SIMNAME.atm instead of SIMNAME.asc
     C :  writes the configuration to file SIMNAME.asc
     V :  writes the velocities to file SIMNAME.vel
     A :  writes the accelerations to file SIMNAME.acc
     F :  writes the forces to file SIMNAME.for
     G :  writes the gradient of energy (numerical derivative)
          to file SIMNAME.gra (to compare with SIMNAME.for) SLOW!
     E :  writes the sum of the numerical gradient and force
          (should be numerical zero)
     P :  writes the Drude amplitides (POLAR only)
   digit=1 denotes the maximum precision in the g-format
   digit=2..9: number of decimal digits in the f-format
   Example: -l111 will output SIMNAME.{asc,vel,acc} with max. precision
   option -n2 does not apply
*/
{
  FILE *asc;
  char fmt[48]; /*e.g., "%4d %2d %14.10f %14.10f %14.10f  %14.10f\n" */
  char *FMT; /* =fmt or a bit more */
  int nnn;
  int optionl=abs(option('l')),order;
  int i,n,ns,ll,llf,pos;
  molecule_t *mn;
  vector *r,*R;
  char *fn="???";
  char *rinfo="???";
  enum key_e { ASC,VEL,ACC,FOR,GRA,ERR,DRU,ATM } key;
  double q;

  if (!optionl) return;

  for (key=ASC,pos=1; key<=DRU; key++,pos*=10) {

    ll=llf=(optionl/pos)%10;

    if (!ll) continue;

    if (llf==1) ll=12;

    order=0;

    if (key!=ASC) order=option('m')+1,ll+=2;
    if (key==VEL) order=1;

    if (key==FOR || key==ACC || key==DRU) {
#ifdef POLAR
      scforces(cfg[order],cfg[0]);
#else /*# POLAR */
      zeroEn();
      forces(cfg[order],cfg[0]);
      En.pot+=En.el;
#endif /*#!POLAR */
    }

    if (key==DRU) order=0; /* we have to calculate SCF, but use Drude */

    if (key==ASC && option('l')<0) key=ATM;

    q=1;
    switch (key) {
      case ASC: fn="asc"; rinfo="r"; break;
      case VEL: fn="vel"; rinfo="veloc"; q=1/h; break;
      case ACC: fn="acc"; rinfo="accel"; break;
      case FOR: fn="for"; rinfo="force"; break;
      case GRA: fn="gra"; rinfo="grad";  break;
      case ERR: fn="err"; rinfo="f+grad";  break;
#ifdef POLAR
      case DRU: fn="dru"; rinfo="Drude"; break;
#else /*# POLAR */
      case DRU: ERROR(("no POLAR and dump of Drude amplitudes requested")) break;
#endif /*#!POLAR */
      case ATM: fn="atm"; rinfo=""; break;
      default: ERROR(("internal")) }

    fn=Fn(fn);

    asc=fopen(fn,"wt");
    if (!asc) ERROR(("cannot write to %s",fn))
    if (key==ATM)
      fprintf(asc,"%d\n",No.s);
    else {
      fprintf(asc,"# ns=%d nmol=%d L=%.15f %.15f %.15f\n",
              No.s,No.N,VARG(box.L));
      fprintf(asc,"# logs=%.16g lambda=%.16g %.16g %.16g\n",
              q*cfg[order]->logs,VARG(q*cfg[order]->lambda));
      if (key==FOR || key==ACC) fprintf(asc,"# Epot=%.15g K\n",En.pot);
      if (key==GRA) fprintf(asc,"# eps=%g (step for numerical derivative)\n",eps);

      fprintf(asc,"# i spec%*sx    %*sy    %*sz         |%s|\n",
              ll,rinfo, ll,rinfo, ll,rinfo, rinfo); }
    ll+=4;

#ifndef TWODIM
    if (llf==1) strcpy(fmt,"%d %d %.16g %.16g %.16g %.16g\n");
    else sprintf(fmt,"%%%dd %%2d %%%d.%df %%%d.%df %%%d.%df  %%%d.%df\n",
                 (int)(log10((double)No.s)+0.9999999),
                 ll,llf, ll,llf, ll,llf, ll,llf);
#else /*# TWODIM */
    if (llf==1) strcpy(fmt,"%d %d %.16g %.16g %.16g\n");
    else sprintf(fmt,"%%%dd %%2d %%%d.%df %%%d.%df  %%%d.%df\n",
                 (int)(log10((double)No.s)+0.9999999),
                 ll,llf, ll,llf, ll,llf);
#endif /*#!TWODIM */

    FMT=fmt;
    if (key==ATM) {
      FMT=strchr(fmt+2,'%');
      FMT=strchr(FMT+1,'%')-5;
      memcpy(FMT,"%-4s ",4);
      strcpy(strend(fmt)-7,"\n"); }

    fprintf(stderr,"writing %s with format %s",fn,FMT);

    nnn=0;

    if (key==ATM)
      fprintf(asc,FMT," box",VARG(box.L));

    loop (n,0,No.N) {
      siteinfo_t *si;

      mn=molec+n;
      ns=mn->ns;
      r=rof(mn,cfg[order]->rp);
      R=rof(mn,cfg[0]->rp);
#ifdef POLAR
      if (key==DRU) r=polarrof(mn,cfg[order]->rp);
      if (key==DRU) R=polarrof(mn,cfg[0]->rp);
#endif /*# POLAR */
      si=spec[mn->sp]->si;

      loop (i,0,ns) {
        vector grad;

        if (key==GRA || key==ERR) {
          double Eplus;
          vector fplus;
          int k;

          loop (k,0,DIM) {
            double R0=R[i][k];

            R[i][k]+=eps;
#ifdef POLAR
            scforces(cfg[order],cfg[0]);
#else /*# POLAR */
            zeroEn();
            forces(cfg[order],cfg[0]);
#endif /*#!POLAR */
            Eplus=En.pot;
            VV(fplus,=r[i])

            R[i][k]=R0-eps;
#ifdef POLAR
            scforces(cfg[order],cfg[0]);
#else /*# POLAR */
            zeroEn();
            forces(cfg[order],cfg[0]);
#endif     /*#!POLAR */
            R[i][k]=R0;

            grad[k]=(Eplus-En.pot)/(2*eps); }
          if (key==ERR) VVV(grad,+=0.5*fplus,+0.5*r[i])
          VV(r[i],=grad) }

        if (key==ACC) VO(r[i],*=si[i].imass)

        if (key==ATM) fprintf(asc,FMT,sitedef[si[i].st].name,VARG(r[i]));
        else fprintf(asc,fmt,nnn,mn->sp,VARG(q*r[i]),sqrt(SQR(q*r[i])));
        nnn++; } /* i */
    } /* n */

    fclose(asc);
    prt("%s writen",fn);

    if (key==ATM && option('l')<0) key=ASC; }
}

void readplayback(int frame,int all) /************************* readplayback */
/*
  Read SIMNAME.plb (to cfg[0]->rp)
  if option -n then read also SIMNAME.vlb (velocities, to cfg[1]->rp)
  all=0:   read L and ns only, close the file(s) even if frame>=0
  all=1:   read also cfg and velocities, if frame>0 then close the file(s)
  frame=0: open file(s) if necessary, read the next frame, and leave the
           file(s) opened
  frame>0: read the specified frame from file(s) and leave the file(s) opened
  frame<0: read the specified frame, and close the file(s);
           if not -n then assign velocities from Maxwell-Boltzmann at given T
  CHANGES in V3.0a:
           frame>0 and <0 swapped
           checking all frames for frame<0 (former frame>0) removed
  NOTES:
    cannot mix calls to frame>1 and frame=0!
    Maxwell velocities not assigned for frame>=0 even if no .vlb
    (??? if -[ then now reading .plb moved to loadcfg())
    -n2 not considered
*/
{
  static FILE *plb,*vlb;
  fvector r,v;
  int i,k,repl;
  static int NS, varL=0;
  int closeplb=frame<0,maxwellplb=closeplb;

  frame=abs(frame);

  if (!plb) {
    /* open and read header */
    plb=fopen(Fn("plb"),"rb");
    if (!plb) ERROR(("open %s",lastFn))
    if (option('n')&1 && all) {
      vlb=fopen(Fn("vlb"),"rb");
      maxwellplb=0;
      if (!vlb) ERROR(("open %s (option -n)",lastFn)) }

    if (fread(r,sizeof(float),2,plb)!=2)
      ERROR(("%s: read header",Fn("plb")))
    if (vlb) {
      if (fread(v,sizeof(float),2,vlb)!=2) ERROR(("%s: read header",Fn("vlb")))
      if (v[0]!=r[0]) ERROR(("%s:%s inconsistent ns",Fn("vlb"),Fn("plb")))
      if (v[1]!=r[1]) ERROR(("%s:%s inconsistent header",Fn("vlb"),Fn("plb"))) }

    NS=r[0];
    if (r[1]>=0) {
      VO(box.L,=r[1]) VV(box.Lh,=0.5*box.L)
      prt("L=%g (cube) read from playback file",box.L[0]); }
    else if (r[1]==-3)
      varL=1;
    else
      ERROR(("wrong header L=%g in the playback file",r[1]))

    if (simils.Ns!=NS)
      WARNING(("readplayback: simils.Ns=%d != NS=%d (from %s)\n\
*** Fixing and continuing with fingers crossed..",simils.Ns,NS,Fn("plb")))
    simils.Ns=NS;
        
    if (!all && simils.Ns!=No.s)
      prt("Reading %s: %d sites found but %d specified (will be solved later)",
          lastFn,simils.Ns,No.s); }

  if (frame) {
    long int pos=8+(long)(NS+varL)*(frame-1)*12;

    fseek(plb,pos,SEEK_SET);
    if (vlb) fseek(vlb,pos,SEEK_SET); }

  /* reading box size */
  if (varL) {
    if (fread(r,sizeof(float),3,plb)!=3) {
      if (frame) ERROR(("%s: reading frame %d, L[3]",Fn("plb"),frame))
      else ERROR(("%s: reading next frame, L[3]",Fn("plb"))) }

    if (vlb) {
      if (fread(r,sizeof(float),3,plb)!=3) {
        if (frame) ERROR(("%s: reading frame %d, L[3]",Fn("vlb"),frame))
        else ERROR(("%s: reading next frame, L[3]",Fn("vlb"))) }

      if (fabs(v[0]-r[0])+fabs(v[1]-r[1])+fabs(v[2]-r[2])>0)
        ERROR(("%s: L=%g %g %g\n\
*** %s: L=%g %g %g inconsistent",Fn("vlb"),VARG(r),Fn("plb"),VARG(v))) } }

  repl=0;
  prt_("%s: box read = [",Fn("plb"));
  loop (k,0,DIM) {
    if (r[k]<=0) {
      box.L[k]/=simils.pc[k];
      repl++;
      prt_("R%g ",box.L[k]); }
    else {
      box.L[k]=r[k];
      prt_("%g ",r[k]); }
    box.Lh[k]=box.L[k]/2; }
  prt(repl?"] (R=replaced by L[])":"]");

  if (all) {
    simils.to=1;
    if (vlb) simils.to=2;
    if (NS!=simils.Ns) ERROR(("internal: simils.Ns=%d != NS(from plb)=%d",simils.Ns,NS))
#ifdef POLAR    
    simils.size=sizeof(ToInt)-sizeof(vector)+simils.Ns*2*sizeof(vector);
#else
    simils.size=sizeof(ToInt)-sizeof(vector)+simils.Ns*sizeof(vector);
#endif
    loop (i,0,simils.to) {
      if (simils.changecfg) sdsalloc(simils.cfg[i],simils.size)
      else simils.cfg[i]=cfg[i]; }

    if (simils.cfg[0]->size != simils.size)
      ERROR(("readplayback: Expected record size=%d != actual record size=%d\n\
*** (as calculated from %s).  Sorry, I cannot fix this.\n\
*** NB: Whereas load.N=4 allows replacing the number of molecules given in\n\
*** the def-file by the cfg-file, the plb-file must match the cfg-file.",simils.cfg[0]->size,simils.size,Fn("plb")))
    if (NS!=No.s) WARNING(("readplayback: previously calculated No.s=%d != actual simils.Ns=NS=%d",simils.Ns,NS))

    loop (i,0,NS) {
      if (fread(r,sizeof(float),3,plb)!=3) {
        if (frame) ERROR(("%s: reading frame %d, site %d",Fn("plb"),frame,i))
        else ERROR(("%s: reading next frame, site %d",Fn("plb"),i)) }
      VV(simils.cfg[0]->rp[i],=r)
      if (vlb) {
        if (frame) ERROR(("%s: reading frame %d, site %d",Fn("vlb"),frame,i))
        else ERROR(("%s: reading next frame, site %d",Fn("vlb"),i))
        VV(simils.cfg[1]->rp[i],=v) } }
      cfg[0]->dep=0; }

  if (closeplb) {
    fclose(plb); plb=NULL;
    if (vlb) { fclose(vlb); vlb=NULL; } }

  if (maxwellplb) {
    //       normalize(-1); ?? - cannot be used if simils.a
    // what if option('m')<2 ?
    Maxwell(FROM,-NS,2); }
}

void readasc(int order,int all) /*********************************** readasc */
/* Reads configuration (order=0), velocities (order=1),
   or accelerations (order=2) in ascii
   Configuration should have been created with option -l
   all=0 means that L only is read
   # info-lines
   # i spec ...
   0 sp x y [z]
   1 sp x y [z]
   ...
   WARNING: limited checks - the asc-files should have been
     saved under the same circumstances
*/
{
  FILE *asc;
  vector *rp=cfg[order]->rp;
  int n,j,sp;
  char line[128];
  char *fn,*c;

  if (simils.changecfg) ERROR(("cfg change not implemented with readasc"))

  if (order<0 || order>2) ERROR(("readasc: bad order=%d",order))

  fn=Fn("asc\0vel\0acc"+order*4);
  asc=fopen(fn,"rt");
  if (!asc) ERROR(("%s not found",fn));
  prt("reading %s",fn);

  for (;;) {
    if (!fgets(line,128,asc)) ERROR(("%s too short",fn));
    if (line[0]!='#') ERROR(("%s: bad format",fn))
    if (!memcmp(line,"# i spec",7)) break;
    prts(line);
    if ( (c=strstr(line,"ns=")) )
      if ((n=atoi(c+3))!=No.s) ERROR(("%s: ns=%d, should be %d",fn,n,No.s))
    if ( (c=strstr(line,"nmol=")) )
      if ((n=atoi(c+5))!=No.N) ERROR(("%s: nmol=%d, should be %d",fn,n,No.N))
    if ( (c=strstr(line,"Epot=")) )
      prt("Epot=%g read",En.pot=atof(c+5));
    if ( (c=strstr(line,"L=")) ) {
      sscanf(c+2,"%lf%lf%lf",box.L,box.L+1,box.L+2);
      prt("L=%g %g %g read",VARG(box.L));
      VV(box.Lh,=0.5*box.L) } }

  if (all) loop (n,0,No.s) {
    if (!fgets(line,128,asc)) ERROR(("readasc: bad line %d",n))
    if (sscanf(line,"%d%d"realfmt realfmt realfmt,&j,&sp,
               &rp[n][0],&rp[n][1]
#ifndef TWODIM
               ,&rp[n][2]
#endif /*# TWODIM */
               )!=DIM+2) Error("asc file: format");
    if (n!=j) ERROR(("%s: n=%d expected, %d read",lastFn,n,j)); }
  fclose(asc);
}

void Maxwell(int from,int to,double prob) /************************ Maxwell */
/*
  assign random Maxwell-Boltzmann velocities (temperature=T) to all
    atoms in molecules from<=n<to
    (prob=1: all, prob=0.5: one half of atoms assigned, prob=0:none...)
  to<0 refers to |to| sites, they are counted from molecule 0,
    but assigned only if molecule>=from
  prob>1 && init==3: always assigned + calls Scorrect
  POLAR only: prob<0: assigns aux masses only (Car-Parrinello)
              with temperature multiplied by -prob
  from<0 is the same as from=FROM but silent
*/
{
  vector *rp=cfg[0]->rp,*rpv=cfg[1]->rp;
  int n;
  int maxns=0,ins=0,verbose=prob>=1;
  double sigma, eps=0;

  if (from<0) from=FROM,verbose=0;

#ifdef POLAR
  int auxonly=0;
  double aux=1;

  if (prob<0) { auxonly++; aux=-prob; prob=1; }
#endif /*# POLAR */

  if (prob>1 && init==3) eps=1e-7;

  if (from<0) from=0;
  if (to<0) maxns=-to, to=No.N;

  Min(to,No.N)

  loop (n,0,to) {
    molecule_t *mn=molec+n;
    int ns=mn->ns;
    vector *r=rof(mn,rpv);
#ifdef POLAR
    vector *rpol=polarrof(mn,rp);
    vector *rvpol=polarrof(mn,rpv);
#endif /*# POLAR */
    int i;
    siteinfo_t *si=spec[mn->sp]->si;

    loop (i,0,ns) {
      ins++;
      if (n>=from && (!maxns || ins<=maxns)) if (rnd()<prob) {
#ifdef POLAR
        if (option('p')/10%10==0) {
          /* Car-Parrinello-like - init velocities */
          sigma=sqrt(aux*T/Sqr(tau.dip)*si[i].alpha_qq)*h;
          VVO(rvpol[i],=rpol[i],+rndgauss()*sigma) }
        if (!auxonly)
#endif /*# POLAR */
          {
            sigma=sqrt(T*si[i].imass)*h;
            VO(r[i],=rndgauss()*sigma)
          }
      } }
    Scorrect(mn,rof(mn,rp),rpv,eps); }

  if (verbose) prt("molecules %d<=n<%d velocities assigned T=%g",from,to,T);
}

void MaxwellCM(double prob) /************************************* MaxwellCM */
/*
  assign random Maxwell-Boltzmann velocities (temperature=T) to all
    molecules from FROM (center-of-mass)
*/
{
  int n;

  loop (n,FROM,No.N) {
    molecule_t *mn=molec+n;
    int ns=mn->ns;
    vector *r=rof(mn,cfg[1]->rp);
    int i;
    siteinfo_t *si=spec[mn->sp]->si;

    if (rnd()<prob) {
      vector cmx;
      double im=1./spec[mn->sp]->mass;
      double sigma=sqrt(T*im)*h;

      VO(cmx,=0)
      loop (i,0,ns) VV(cmx,-=si[i].mass*r[i])
      VVO(cmx,=im*cmx,+rndgauss()*sigma)
      loop (i,0,ns) VV(r[i],+=cmx) } }

  removedrifts(1); /* added in V3.6k */
}

void makefixa(void) /********************************************** makefixa */
{
  if (fixsites0) {
    sdsralloczero(fixa,cfg[0]->size);
    sdscopy(fixa,cfg[0]); }
}

void loadfixa() /************************************************** loadfixa */
{
  if (fixsites0) {
    fixsites_t *f;
    FILE *fxc=fopen(Fn("fxc"),"rb");

    if (fxc) {
      fclose(fxc);
      VarOpen(Fn("fxc"),"rb");
      prt("loading %s (atom positions to be fixed)",lastFn);
      if (fixa->size!=VarFile.size) ERROR(("%s: bad number of sites",lastFn))
      VarRead(fixa,fixa->size);
      VarClose();
      prt("%s loaded",lastFn); }
    else
      prt("WARNING: %s missing\n\
         the positions to be fixed will be taken from configuration\n\
         or from %s (if present there)\n", Fn("fxc"),Fn("fix"));

    looplist (f,fixsites0) {
      int i;
      char *src=fxc?Fn("fxc"):"configuration";

      loop (i,f->from,f->to) {
        if (f->isr) {
          if (option('v')&4)
            prt("site %d fixed by spring at %g %g %g (specified in %s)",
              i,VARG(f->r),Fn("fix"));
          VV(fixa->rp[i],=f->r) }
        else {
          if (option('v')&4)
            prt("site %d fixed by spring at %g %g %g (taken from %s)",
              i, VARG(fixa->rp[i]), src); } } } }

}

void savefixa(void) /********************************************** savefixa */
{
  if (fixsites0) {
    VarOpen(Fn("fxc"),"wb");
    VarPut(fixa,fixa->size);
    VarClose();
    prt("%s saved",lastFn); }
}

void dumpall(ToIntPtr B,ToIntPtr A) /******************************* dumpall */
{
  static int pass=0;
  FILE *dump=fopen(Fn(string("%d.dump",pass)),"wt");
  int n,i;

  fprintf(dump,"%d -3 NS -DIM\n%17.12f %17.12f %17.12f L[]\n",No.s,VARG(box.L));

  loop (n,0,No.N) {
    molecule_t *mn=molec+n;
    siteinfo_t *si=spec[mn->sp]->si;
    vector *r=rof(mn,A->rp);
    vector *f=rof(mn,B->rp);
#ifndef POLAR
    loop (i,0,mn->ns)
      fprintf(dump,"%17.12f %17.12f %17.12f  %17.12f  %16.9f %16.9f %16.9f %d %d r[] q f[] mol site\n",
              VARG(r[i]),                    si[i].charge, VARG(f[i]),     n, i);
#elif POLAR&32
    vector *rp=polarrof(mn,A->rp);
    vector *fp=polarrof(mn,B->rp);
    loop (i,0,mn->ns) if (si[i].qtype&FQ)
      fprintf(dump,"%17.12f %17.12f %17.12f  %17.12f %17.12f  %16.9f %16.9f %16.9f %d %d r[]  FQ Phi  f[] mol site\n",
              VARG(r[i]),                    rp[i][0],fp[i][0],VARG(f[i]),          n, i);
    else
      fprintf(dump,"%17.12f %17.12f %17.12f  %17.12f  %17.12f %17.12f %17.12f  %17.12f  %16.9f %16.9f %16.9f   %16.9f %16.9f %16.9f %d %d r[] q rDrude[] qDrude f[] fDrude[] mol site\n",
              VARG(r[i]),                    si[i].charge,VARG(rp[i]),         si[i].chargepol,VARG(f[i]),     VARG(fp[i]),         n, i);
#else /*#!POLAR!POLAR&32 */
    /* normal + Drude */
    vector *rp=polarrof(mn,A->rp);
    vector *fp=polarrof(mn,B->rp);
    loop (i,0,mn->ns)
      fprintf(dump,"%17.12f %17.12f %17.12f  %17.12f  %17.12f %17.12f %17.12f  %17.12f  %16.9f %16.9f %16.9f   %16.9f %16.9f %16.9f %d %d r[] q rDrude[] qDrude f[] fDrude[] mol site\n",
              VARG(r[i]),                    si[i].charge,VARG(rp[i]),         si[i].chargepol,VARG(f[i]),     VARG(fp[i]),             n, i);
#endif /*#!POLAR!POLAR&32 */
  }
  fclose(dump);

  if (++pass>=10) ERROR(("limit of 10 dumps reached"))
}

void zerocfg(int zero) /******************************************** zerocfg */
{
  int i,k;
  vector *r;

  if (!zero) return;

  loop (k,0,9) if (cfg[k]) {
    r=rof(molec,cfg[k]->rp);
    if (zero&1) loop (i,0,No.s) r[i][0]=0;
    if (zero&2) loop (i,0,No.s) r[i][1]=0;
    if (zero&4) loop (i,0,No.s) r[i][2]=0; }

  //  normalize(-1);
}

double Lfromfile(int k,double t) /******************************** Lfromfile */
{
  static struct history_s {
    double t;
    vector L;
  } one,hist[3];
  static int nl=0;
  static FILE *rhof=NULL;
  char line[128];
  double rho;
  int i;

  if (!rhof) {
    rhof=fopen(Fn("box"),"rt");
    if (!rhof) ERROR(("tau.rho<0 specified and no %s found",lastFn)) }

  while (t>=hist[2].t || nl<3) {
    if (!fgets(line,128,rhof)) break;
    if (tau.rho==-1) {
      if (sscanf(line,"%lf%lf",&one.t,&rho)<2) continue;
      one.L[0]=one.L[1]=one.L[2]=cbrt(No.mass*rhounit/rho); }
    else
      if (sscanf(line,"%lf%lf%lf%lf",&one.t,one.L,one.L+1,one.L+2)<4) continue;
    hist[0]=hist[1]; hist[1]=hist[2]; hist[2]=one;
    nl++; }

  if (nl<3) ERROR(("%s: at least 3 data lines expected",lastFn))
  i=t>hist[1].t;

  return hist[i].L[k]+(t-hist[i].t)/(hist[i+1].t-hist[i].t)*(hist[i+1].L[k]-hist[i].L[k]);
}

void remove1mol(int n) /***************************************** remove1mol */
/***
    remove molecule number n
    - move and shorten all allocated cfg[i] (but does not re-allocate them)
    - change spec[i]->N, No.N
    ! does NOT change No.f, molec[n], etc. - to be used to save cfg only
    ! does NOT change fixed/tethered/ANCHORed sites
***/
{
  int i,j,ns;
  molecule_t *mn;

  if (n<0 || n>=No.N) {
    WARNING(("Molecule n=%d to be removed is out of range and ignored",n))
    /* no stop simulation signal: */
    sig=0;
    return; }

  /* this cannot be in savecfg because internal tables (as rof) are broken */
  removedrifts(1);
  depend_r(cfg[0],0);

  mn=molec+n;
  ns=mn->ns;

  if (spec[mn->sp]->N)
    spec[mn->sp]->N--;
  else {
    ERROR(("request to remove a molecule but spec[%d]->N=0",mn->sp))
    return; }

  loop (i,0,option('r')) if (cfg[i]) {
    vector *r=rof(mn,cfg[i]->rp);
    //    prt("remove1mol: size=%d ns=%d %p %p %p %ld",cfg[i]->size,ns,cfg[i]->rp,r,r+ns,((cfg[i]->rp+No.s)-(r+ns))*sizeof(vector));
    prt("remove1mol: a[%d]: molecule=%d spec=%d ns=%d",i,n,mn->sp,ns);
    loop (j,0,ns) prt("%.16g %.16g %.16g %d %d REMOVE1MOL%d", VARG(r[j]),j,n,i);
    memmove(r,r+ns,((cfg[i]->rp+No.s)-(r+ns))*sizeof(vector));
    cfg[i]->size-=sizeof(vector)*ns;
#ifdef POLAR
    // not tested !!!
    memmove(cfg[i]->rp+No.s-ns,cfg[i]->rp+No.s,(r-cfg[i]->rp)*sizeof(vector));
    memmove(r+No.s-ns,r+No.s+ns,((cfg[i]->rp+No.s)-(r+ns))*sizeof(vector));
    cfg[i]->size-=sizeof(vector)*ns;
#endif /*# POLAR */
  }
  No.s-=ns;
  No.N-=1;
  /* WARNING: internal tables are broken, only savecfg is allowed */

  prt("NOTE: molecule %d removed, to be restarted with changed N[] or load.N=4",n);

  return;
}
