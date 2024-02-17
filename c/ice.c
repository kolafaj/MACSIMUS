/* \make ice
 */
#include "ground.h"
#include "rndgen.h"
#include "statics.h"
typedef float real;
#include "vector3d.h"
#include <signal.h>
#include "datetime.c"

#define VERSION "V2.9"

volatile int sig;
void asksig(int signr) /* =========================================== asksig */
{
  if (sig++) exit(-1); /* ^C typed twice */

  signal(SIGINT,asksig);
  fprintf(stderr,
    "\nINTERRUPT(%d) ACCEPTED, wait to finish cycle or type ^C again to exit immediately\n",
    signr);
}

#include "icedata.c"

vector L;
int N,nomit,sphere,checkodd=1;

double ddist(float *a,float *b) /************************************* ddist */
{
  int i;
  double rr=0;
#if 0
  loop (i,0,3)
    rr+=Sqr(L[i]/2-fabs(L[i]/2-fabs(a[i]-b[i])));
#else
  /* needed for ice V */
  vector d;
  VVV(d,=a,-b)
  loop (i,0,3) {
    while (d[i]<-L[i]/2) d[i]+=L[i];
    while (d[i]>L[i]/2) d[i]-=L[i];
    rr+=d[i]*d[i]; }
#endif

  return rr;
}

/* oxygen positions */
struct cfg_s {
  vector r;
  int nnbr; /* should be 4; -1 means to omit this water */
  int nbr[4]; /* HHee */
} *cfg;

int ndef(int i) /****************************************************** ndef */
/* number of Bjerrum defects around molecule i */
{
  int n,d=0;

  n=cfg[i].nbr[0]; if (cfg[n].nbr[0]==i || cfg[n].nbr[1]==i) d++;
  n=cfg[i].nbr[1]; if (cfg[n].nbr[0]==i || cfg[n].nbr[1]==i) d++;
  n=cfg[i].nbr[2]; if (cfg[n].nbr[2]==i || cfg[n].nbr[3]==i) d++;
  n=cfg[i].nbr[3]; if (cfg[n].nbr[2]==i || cfg[n].nbr[3]==i) d++;

  return d;
}

int charge(int i) /************************************************** charge */
/* Bjerrum defects around molecule i w. sign */
{
  int n,d=0;

  n=cfg[i].nbr[0]; if (cfg[n].nbr[0]==i || cfg[n].nbr[1]==i) d++;
  n=cfg[i].nbr[1]; if (cfg[n].nbr[0]==i || cfg[n].nbr[1]==i) d++;
  n=cfg[i].nbr[2]; if (cfg[n].nbr[2]==i || cfg[n].nbr[3]==i) d--;
  n=cfg[i].nbr[3]; if (cfg[n].nbr[2]==i || cfg[n].nbr[3]==i) d--;

  return d;
}

double xM(int i) /******************************************************* xM */
/* x-component of the dipole moment of molecule i */
{
  int j,n;
  double dx,d=0;

  loop (j,0,2) {
    n=cfg[i].nbr[j];
    dx=cfg[n].r[0]-cfg[i].r[0];
    if (dx<-L[0]/2) dx+=L[0];
    if (dx> L[0]/2) dx-=L[0];
    d+=dx; }

  return d;
}

int ndef_lt(int i) /************************************************ ndef_lt */
/* number of Bjerrum defects around molecule i from molecules <i */
{
  int n,d=0;

  n=cfg[i].nbr[0]; if (n<i) if (cfg[n].nbr[0]==i || cfg[n].nbr[1]==i) d++;
  n=cfg[i].nbr[1]; if (n<i) if (cfg[n].nbr[0]==i || cfg[n].nbr[1]==i) d++;
  n=cfg[i].nbr[2]; if (n<i) if (cfg[n].nbr[2]==i || cfg[n].nbr[3]==i) d++;
  n=cfg[i].nbr[3]; if (n<i) if (cfg[n].nbr[2]==i || cfg[n].nbr[3]==i) d++;

  return d;
}

vector M; /* dipole moment */
int pass;

double dipmq(double OO) /********************************************* dipmq */
{
  int i,n,j,k;
  vector dr;
  double mu=2.3052; /* [D] dip.mom. of TIP4P/2005 water */
  double Mcomp=fabs(OO)*sqrt(4./3); /* dip.mom. of this tetrahedral model */

  VO(M,=0)

  loop (i,0,N) loop (j,0,2) {
    n=cfg[i].nbr[j];
    VVV(dr,=cfg[n].r,-cfg[i].r)
    loop (k,0,3) {
      if (dr[k]<-L[k]/2) dr[k]+=L[k];
      if (dr[k]> L[k]/2) dr[k]-=L[k]; }
    VV(M,+=dr) }

  VO(M,/=Mcomp)
  prt("  M = %.6f %.6f %.6f = %g mu_H2O (pass %d)",VARG(M),sqrt(SQR(M)),pass);
  VO(M,*=mu)
  prt("  M = %.6f %.6f %.6f = %g D (for mu_H2O=2.3052 D)",VARG(M),mu=sqrt(SQR(M)));

  pass++;
  return mu;
}

int indx(char *arg,int df) /******************************************* indx */
/* "8y" returns 1, "9" returns df */
{
  char c=strend(arg)[-1];
  if (strchr("xyz",c)) return c-'x';
  else return df;
}

struct cell_s *ice;
double Mmid,Merr;
int quiet=1;
char *args;
char *SITES="HHO";

void nextcfg(int i,int nr) /**************************************** nextcfg */
{
  int k;

  if (nr&1) {
    k=cfg[i].nbr[0],cfg[i].nbr[0]=cfg[i].nbr[2],cfg[i].nbr[2]=k;
    k=cfg[i].nbr[1],cfg[i].nbr[1]=cfg[i].nbr[3],cfg[i].nbr[3]=k; }
  else {
    k=cfg[i].nbr[0],cfg[i].nbr[0]=cfg[i].nbr[1],cfg[i].nbr[1]=cfg[i].nbr[2],cfg[i].nbr[2]=k; }
}

void randomize(int nrand) /*************************************** randomize */
/* shuffle molecules, nrand times each */
{
  int i,irand;

  loop (i,0,N)
    loop (irand,0,nrand) {
    int a=irnd(4),b=irnd(4),k;
      k=cfg[i].nbr[a], cfg[i].nbr[a]=cfg[i].nbr[b], cfg[i].nbr[b]=k; }
}

void annealing(void) /******************************************** annealing */
{
  /* simulated annealing if initial random lattice */
  int i,j,k;
  int a,b;
  int old[4];
  double T=1;
  int Uold,Unew;

  do {
    randomize(100);
    prt("original dipole moment of the box:");
    dipmq(ice->unit);

    /* Monte Carlo assignment of protons by simulated annealing */
    for (;sig==0;) {
      int nd=0;
      double exptab0[9],*exptab=exptab0+4;

      loop (i,0,N) nd+=ndef(i);
      nd/=2; /* every defect was included twice */
      if (!quiet) prt("%d %g",nd,T);
      if (nd==0) break;
      if (nd&1) ERROR(("%d=odd number of Bjerrum defects",nd))
      T*=1-2*T/N;

      loopto (k,-4,4) exptab[k]=exp(k/T);

      loop (j,0,8*N) {
        i=irnd(N);

        memcpy(old,cfg[i].nbr,sizeof(cfg[i].nbr));
        Uold=ndef(i);
        //a=irnd(4),b=irnd(4);
        a=irnd(2),b=2+irnd(2); // H,e
        k=cfg[i].nbr[a],cfg[i].nbr[a]=cfg[i].nbr[b],cfg[i].nbr[b]=k;
        Unew=ndef(i);
        if (rnd()>exptab[Uold-Unew])
          /* rejected: return back the old cfg */
          memcpy(cfg[i].nbr,old,sizeof(cfg[i].nbr)); } }

    prt("final dipole moment of the box:");
  } while (fabs(dipmq(ice->unit)-Mmid)>Merr && sig==0);
}

double MC(int ino,double T) /******************************************** MC */
{
  /* MC of given configuration; averaged energy returned */
  int i,j,k,nd;
  int a,b;
  int old[4];
  int Uold,Unew;
  double exptab0[9],*exptab=exptab0+4;

  loopto (k,-4,4) exptab[k]=exp(k/T);

  loop (j,0,ino*N) {
    i=irnd(N);

    memcpy(old,cfg[i].nbr,sizeof(cfg[i].nbr));
    Uold=ndef(i);
    //      do a=irnd(4),b=irnd(4); while (a==b);
    a=irnd(2),b=2+irnd(2); // H,e
    k=cfg[i].nbr[a],cfg[i].nbr[a]=cfg[i].nbr[b],cfg[i].nbr[b]=k;
    Unew=ndef(i);

    if (rnd()>exptab[Uold-Unew])
      /* rejected: return back the old cfg */
      memcpy(cfg[i].nbr,old,sizeof(cfg[i].nbr)); }

  nd=0;
  loop (i,0,N) nd+=ndef(i);
  //  if (nd&1) ERROR(("%d=odd number of Bjerrum defects",nd))
  nd/=2; /* every defect was included twice */

  return (double)nd/N;
}

double MCfield(int ino,double T,double E) /************************* MCfield */
/* MC in electric field */
{
  int i,j,k,nd;
  int a,b;
  int old[4];
  double Uold,Unew;

  loop (j,0,ino*N) {
    i=irnd(N);

    memcpy(old,cfg[i].nbr,sizeof(cfg[i].nbr));

    Uold=ndef(i)+E*xM(i);

    a=irnd(2),b=2+irnd(2); // H,e
    k=cfg[i].nbr[a],cfg[i].nbr[a]=cfg[i].nbr[b],cfg[i].nbr[b]=k;

    Unew=ndef(i)+E*xM(i);

    if (rnd()>exp(-(Unew-Uold)/T))
      /* rejected: return back the old cfg */
      memcpy(cfg[i].nbr,old,sizeof(cfg[i].nbr)); }

  nd=0;
  loop (i,0,N) nd+=ndef(i);
  //  if (nd&1) ERROR(("%d=odd number of Bjerrum defects",nd))
  nd/=2; /* every defect was included twice */

  return (double)nd/N;
}

double MCinside(int ino,double T) /******************************** MCinside */
/* simulate with the layers fixed */
{
  int i,j,k,nd;
  int a,b;
  int old[4];
  int Uold,Unew;
  double exptab0[9],*exptab=exptab0+4;

  loopto (k,-4,4) exptab[k]=exp(k/T);

  loop (j,0,ino*N) {
    i=nomit/2+irnd(N-nomit/2);

    memcpy(old,cfg[i].nbr,sizeof(cfg[i].nbr));
    Uold=ndef(i);
    //      do a=irnd(4),b=irnd(4); while (a==b);
    a=irnd(2),b=2+irnd(2); // H,e
    k=cfg[i].nbr[a],cfg[i].nbr[a]=cfg[i].nbr[b],cfg[i].nbr[b]=k;

    Unew=ndef(i);

    if (rnd()>exptab[Uold-Unew])
      /* rejected: return back the old cfg */
      memcpy(cfg[i].nbr,old,sizeof(cfg[i].nbr)); }

  nd=0;
  loop (i,0,N) nd+=ndef(i);
  //  if (nd&1) ERROR(("%d=odd number of Bjerrum defects",nd))
  nd/=2; /* every defect was included twice */

  return (double)nd/N;
}

void writecfg(char *fn) /****************************************** writecfg */
{
  FILE *plb=fopen(string("%s.plb",fn),"wb");
  FILE *mol=fopen(string("%s.mol",fn),"wt");
  FILE *def=fopen(string("%s.def",fn),"wt");
  FILE *fix=fopen(string("%s.fix",fn),"wt");
  vector hdr,rn,rn2;
  double wOH=0.96/ice->OO;
  double wOM=0.1559/ice->OO;
  double wOL=0.8/ice->OO;
  double cutoff;
  int i;
  int nsites=strlen(SITES),isM=!!strchr(SITES,'M'),isL=!!strchr(SITES,'L');
  int Nleft=N;

  if (sphere>=N) sphere=0; /* all molecules */

  if (sphere*nomit) Error("internal");

  if (sphere) {
    /* keep sphere molecules from the center */
    int j;

    loop (j,0,N-sphere) {
      double maxrr=0;
      int im=-1;

      loop (i,0,N) if (cfg[i].nnbr>=0) {
        double rr=Sqr(cfg[i].r[0]-L[0]/2)
                 +Sqr(cfg[i].r[1]-L[1]/2)
                 +Sqr(cfg[i].r[2]-L[2]/2)+1e-9*rnd();
        if (rr>maxrr) maxrr=rr,im=i; }
      if (im<0) Error("internal");
      cfg[im].nnbr=-1; }
    Nleft=sphere; }

  if (nomit) {
    /* omit z-layer (may be wrong for monoclinic ice V) */
    Nleft=N-nomit;
    double zmin=9e9,zmax=-9e9;
    int j;

    if (Nleft<=0) ERROR(("wrong option -o (nomit=%d)",Nleft ))

    loop (i,0,N) {
      Min(zmin,cfg[i].r[2])
      Max(zmax,cfg[i].r[2]) }

    loop (j,0,nomit) {
      double zz=9e9;
      int im=-1;

      loop (i,0,N) if (cfg[i].nnbr>=0) {
        double dz=min(cfg[i].r[2]-zmin,zmax-cfg[i].r[2]);
        if (zz>dz) zz=dz,im=i; }
      if (im<0) Error("internal");
      cfg[im].nnbr=-1; } }

  hdr[0]=Nleft*nsites;
  hdr[1]=-3;
  fwrite(hdr,4,2,plb);
  fwrite(L  ,4,3,plb);
  cutoff=fmin(L[0],L[1]);
  cutoff=fmin(cutoff,L[2])/2;
  fprintf(def,"! %s\n\
! Mcell = %.6f %.6f %.6f = %g D (for mu_H2O=2.3052 D)\n\
n=%d N[0]=n\n\
L[0]=%.7f\n\
L[1]=%.7f\n\
L[2]=%.7f\n\
cutoff=%.7f LJcutoff=cutoff\n",
          args,
          VARG(M),sqrt(SQR(M)),
          Nleft,VARG(L),cutoff);
  fclose(def);

  fprintf(mol,"ice %s\n\
\n\
parameter_set = dumb\n\
number_of_atoms = %d\n\
\n\
! Mcell = %.6f %.6f %.6f = %g D (for mu_H2O=2.3052 D)\n\
\n\
atoms\n\
! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n",
          args,Nleft*(int)strlen(SITES),VARG(M),sqrt(SQR(M)));

  loop (i,0,N) if (cfg[i].nnbr>=0) {
    int hb=0,hl=0,l,n;
    char *c,*cc;

    for (c=SITES; *c; c++) {
      n=i*nsites+(int)(c-SITES);
      switch (*c) {
        case 'O':
          VV(hdr,=cfg[i].r)
          fprintf(mol,"%3d R120/O%d- O 0.0 0 %d",n,i,2+isM+2*isL);
          for (cc=SITES; *cc; cc++) if (strchr("HML",*cc))
            fprintf(mol," %d",i*nsites+(int)(cc-SITES));
          fprintf(mol,"\n");
          fprintf(fix,"%4d %6.3f %6.3f %6.3f\n",n,hdr[0],hdr[1],hdr[2]);
          break;
        case 'H':
          if (hb>=2) Error("wrong SITES (too many H)");
          VV(rn,=cfg[cfg[i].nbr[hb]].r)
          loop (l,0,3) {
            if (rn[l]-cfg[i].r[l]>L[l]/2) rn[l]-=L[l];
            if (rn[l]-cfg[i].r[l]<-L[l]/2) rn[l]+=L[l]; }
          VVV(hdr,=(1-wOH)*cfg[i].r,+wOH*rn)
          hb++;
          fprintf(mol,"%3d W070/H%d-%d H 0.0 0 1",n,i,n%nsites);
          for (cc=SITES; *cc; cc++) if (*cc=='O')
            fprintf(mol," %d",i*nsites+(int)(cc-SITES));
          fprintf(mol,"\n");
          break;
        case 'L':
          if (hl>=2) Error("wrong SITES (too many L)");
          VV(rn,=cfg[cfg[i].nbr[2+hl]].r)
          loop (l,0,3) {
            if (rn[l]-cfg[i].r[l]>L[l]/2) rn[l]-=L[l];
            if (rn[l]-cfg[i].r[l]<-L[l]/2) rn[l]+=L[l]; }
          VVV(hdr,=(1-wOL)*cfg[i].r,+wOL*rn)
          hl++;
          fprintf(mol,"%3d B060/L%d-%d L 0.0 0 1",n,i,n%nsites);
          for (cc=SITES; *cc; cc++) if (*cc=='O')
            fprintf(mol," %d",i*nsites+(int)(cc-SITES));
          fprintf(mol,"\n");
          break;
        case 'M':
          VV(rn,=cfg[cfg[i].nbr[2]].r)
          loop (l,0,3) {
            if (rn[l]-cfg[i].r[l]>L[l]/2) rn[l]-=L[l];
            if (rn[l]-cfg[i].r[l]<-L[l]/2) rn[l]+=L[l]; }
          VVV(hdr,=(1-wOH)*cfg[i].r,+wOH*rn)
          VV(rn2,=cfg[cfg[i].nbr[3]].r)
          loop (l,0,3) {
            if (rn2[l]-cfg[i].r[l]>L[l]/2) rn2[l]-=L[l];
            if (rn2[l]-cfg[i].r[l]<-L[l]/2) rn2[l]+=L[l]; }
          VVVV(hdr,=(1+2*wOM)*cfg[i].r,-wOM*rn,-wOM*rn2)
          fprintf(mol,"%3d M110/M%d M 0.0 0 1",n,i);
          for (cc=SITES; *cc; cc++) if (*cc=='O')
            fprintf(mol," %d",i*nsites+(int)(cc-SITES));
          fprintf(mol,"\n");
          break;
        default:
          Error("bad SITES: only OHML accepted"); }
      VO(hdr,+=1) /* linked-cell dislikes negative coords */
      fwrite(hdr,4,3,plb); } }
  fclose(plb);
  fclose(mol);
  fclose(fix);

  prt("molcfg -%d %s-water %s",N,SITES,fn);
  prt("box = %.7f %.7f %.7f",VARG(L));

  {
    double epsr=95; /* 0 C: 96.5; 91.5; */
    double eps0=8.8541878e-12;
    double k=1.380649e-23;
    double V=L[0]*L[1]*L[2]*1e-30;
    double Debye=1./2.99792458e+29;

    double M=sqrt((epsr-1)*(2*epsr+1)/(3*epsr)*3*eps0*273*k*V);
    prt("M[epsr=%g] = %g D",epsr,M/Debye);
  }
}

struct cfg_s *cfg0;
long long int ndis=0,maxndis=0,info,branched;

int nbrsof(int n) /************************************************** nbrsof */
{
  int i;

  if (n>=N) {
    ndis++;
    if (ndis>=info) { prt("%lld FOUND",ndis); info+=1000+info/10; }
    if (ndis<maxndis)
      /* NEW in V2: numbered from 0 */
      writecfg(string("all-%lld",ndis-1));
    return 1; /* proton cfg found! */
  }

  cfg[n]=cfg0[n];

  loop (i,0,6) {
    nextcfg(n,i);
    if (ndef_lt(n)==0) nbrsof(n+1); }

  return 0;
}

int nx,ny,nz;

int shift(int i,int j,int k,int n) /********************************** shift */
{
  int jcell=j/n;

  switch (k) {
    case 2:
      jcell=jcell%nz;
      if (jcell==nz-1) jcell=-1;
      break;
    case 1:
      jcell=(jcell/nz)%ny;
      if (jcell==ny-1) jcell=-1;
      break;
    case 0:
      jcell=(jcell/(nz*ny))%nx;
      if (jcell==nx-1) jcell=-1;
      break;
  }
  if (jcell>1) ERROR(("shift %d %d %d %d  jcell=%d",i,j,k,n,jcell))

  return jcell;
}

#define NEST 20
int maxnest,nest,pertest;
unsigned cyc[NEST+1];
unsigned hist[NEST+1];

int findi,nestto,maxnestto,found;
unsigned cycto[NEST+1];

void propagateto(int i) /*************************************** propagateto */
/* one step of generating path to findi */
{
  int j;

  nestto++;
  cycto[nestto]=i;
  if (cycto[nestto]==findi) {
    //    int idr[3]={0,0,0};

    //    /* check periodic */
    //    loop (i,0,nestto) loop (j,0,3)
    //      idr[j]+=(int)((cfg[cyc[i+1]].r[j]-cfg[cyc[i]].r[j])/(L[j]/2));
    //
    //    if (SQR(idr)) goto ret;

    found++;
    goto ret; }
  else {
    loop (j,1,nestto) if (cycto[nestto]==cycto[j]) goto ret;
  }
  if (nestto<=maxnestto)
    loop (j,0,4) propagateto(cfg[i].nbr[j]);

 ret:
  nestto--;
}

int mindist(int i,int j) /****************************************** mindist */
/* min topological distance */
{
  if (i==j) return 0;
  findi=j;
  loop (maxnestto,0,nest) {
    nestto=0;
    found=0;
    loop (j,0,4) propagateto(cfg[i].nbr[j]);

    if (found) return maxnestto+1; }

  return -1;
}

void propagate(int i) /******************************************* propagate */
/* one step of generating cycles */
{
  int j,k,n;

  nest++;
  cyc[nest]=i;
  if (cyc[nest]==cyc[0]) {
    int idr[3]={0,0,0};

    if (pertest) {
      /* check periodic, fails for ice V */
      loop (i,0,nest) loop (j,0,3)
        idr[j]+=(int)((cfg[cyc[i+1]].r[j]-cfg[cyc[i]].r[j])/(L[j]/2));

      if (SQR(idr)) goto ret; }

    /* check branched */
    if (!branched) loop (j,0,nest)
      loop (k,j+1,nest) {
        n=k-j;
        if (n>nest/2) n=nest-n;
        if (mindist(cyc[j],cyc[k])<n) goto ret; }

#if 0
    if (nest==12) {
      loop (i,0,nest) printf("%d %d\n",i,cyc[i]);
      printf("\n");
    }
#endif

    hist[nest]++;
    goto ret; }
  else {
    loop (j,1,nest) if (cyc[nest]==cyc[j]) goto ret;
  }
  if (nest<=maxnest)
    loop (j,0,4) propagate(cfg[i].nbr[j]);

 ret:
  nest--;
}

int nnum,dden;

void numden(int inum,int iden) /************************************* numden */
{
  int n;

  for (;;) {
    if (inum>iden) n=inum,inum=iden,iden=n;
    iden-=inum;
    if (iden==0) {
      nnum/=inum; dden/=inum;
      return; }
  }
}

int Num(int inum,int iden) /******************************************** Num */
{
  nnum=inum; dden=iden;
  numden(inum,iden);
  return nnum;
}

int Den(int inum,int iden) /******************************************** Den */
{
  nnum=inum; dden=iden;
  numden(inum,iden);
  return dden;
}

char *wire;
void rewire(char *fn) /********************************************** rewire */
{
  int i;
  FILE *f=fopen(string("%s.wire",fn),"rt");
  /* file format:
! comment
d site site ! to disconnect
...
c site site ! to connect
...
! number of connect and disconnect must match and be in order
  */
  char line[128];

  loop (i,0,N) printf("%d -> %d %d %d %d original\n",i
                      ,cfg[i].nbr[0]
                      ,cfg[i].nbr[1]
                      ,cfg[i].nbr[2]
                      ,cfg[i].nbr[3]
                      );

  if (!f) Error("no ICENAME.wire file, stop");

  while (fgets(line,128,f)) {
    int s1,s2,err;

    if (line[0]=='!') continue;

    if (2!=sscanf(line+2,"%d%d",&s1,&s2)) {
      fprintf(stderr,"%swire line format error, skipping\n",line);
      continue; }

    if (s1==s2) {
      fprintf(stderr,"%swire site numbers must differ, skipping\n",line);
      continue; }

    if (line[0]=='d') {
      /* disconnect */
      err=1;
      loop (i,0,cfg[s1].nnbr) if (cfg[s1].nbr[i]==s2) cfg[s1].nbr[i]=-1,err=0;
      if (err) {
        fprintf(stderr,"%sno bond %d->%d\n",line,s1,s2);
        Error("wire"); }
      err=1;
      loop (i,0,cfg[s2].nnbr) if (cfg[s2].nbr[i]==s1) cfg[s2].nbr[i]=-1,err=0;
      if (err) {
        fprintf(stderr,"%sno bond %d->%d\n",line,s2,s1);
        Error("wire"); } }

    else if (line[0]=='c') {
      /* connect */
      err=1;
      loop (i,0,cfg[s1].nnbr) if (cfg[s1].nbr[i]==-1) cfg[s1].nbr[i]=s2,err=0;
      if (err) {
        fprintf(stderr,"%scannot wire %d->%d, use d[isconnect] first\n",line,s1,s2);
        Error("wire error"); }
      err=1;
      loop (i,0,cfg[s2].nnbr) if (cfg[s2].nbr[i]==-1) cfg[s2].nbr[i]=s1,err=0;
      if (err) {
        fprintf(stderr,"%scannot wire %d->%d, use d[isconnect] first\n",line,s2,s1);
        Error("wire error"); } }
  }

  loop (i,0,N) printf("%d -> %d %d %d %d rewired\n",i
                      ,cfg[i].nbr[0]
                      ,cfg[i].nbr[1]
                      ,cfg[i].nbr[2]
                      ,cfg[i].nbr[3]
                      );

  fclose(f);
}

int main(int narg,char **arg) /**************************************** main */
{
  int n,nn,ix,iy,iz,i,j,k,l,iarg,num=0;
  double V;
  int RNDSEED=0;
  double Mfrom=0, Mto=1e10, E=0;
  double T=-1;
  int no=1048576,noint=2,lag=11,eq=-128,usebeta=0,gennbrs=0;
  char *ICENAME="",*basename="NONAME";

  if (narg<2) {
    fprintf(stderr,"\
Generate ice structure with proton disorder, " VERSION ". Call by:\n\
  ice TYPE [NX[COORD] NY[COORD] NZ[COORD] [OPTIONS]\n\
TYPE = type of ice {2D,Ic,Ih,III,V,VI,VII} or clathrate {ClI,ClII,ClH}\n\
NX,NY,NZ = how many times to repeat the basic cell (missing=only info on ice)\n\
COORD = {x,y,z} denotes order of coordinates; default=xyz\n\
OPTIONS =\n\
  -sSITES generate plb-file with given structure and order of sites [df.=%s]\n\
          SITES=string of {H,H,O,M,L}; M=TIP4P-like site, L=lone pair\n\
  -v      verbose (annealing protocol to .prt, with -g connect data)\n\
  -mMIN   (with annealing) minimum dipole moment accepted\n\
  -MMAX   maximum dipole moment, in D with mu_H2O=mu_TIP4P/2005\n\
  -aN     all cfgs generated, max N are written (NEW in V2: numbered from 0)\n\
          BUG: incorrect for 1- and 2-cycles (e.g., ice 2D 2 2 1 -a1)\n\
  -oN     N<0, with annealing: |N/2| atoms will be omitted from both z-ends\n\
            Useful for ice V. Bad numbering: run mol2mol INNAME +a OUTNAME\n\
          N<0, with MC: field simulation\n\
          N>0, with annealing: make sphere with N water molecules\n\
  -zSEED  seed for rnd [default = time]\n\
  -TTEMP  temperature for Monte Carlo, 0=annealing to 0 [default=-1=nothing]\n\
  -bBETA  beta=1/T for Monte Carlo\n\
  -BNUM/DEN   beta=1/T for MC = NUM/DEN*ln(N) (N=no.of molecules,NUM=integer)\n\
  -nNO    number of MC cycles [1048576]\n\
  -iNOINT number of MC steps/molecule in one cycle [2]\n\
  -eEQUIL number of MC cycles for equilibration\n\
          negative: diffusion-based correlation time * |EQUIL| [-128]\n\
  -lLAG   lag for statistics [11]\n\
  -g      generate neighbors from distances (not VII, cf. ice-gentab)\n\
          (at least TYPE 3 3 3 -v -T0 -g to generate nbrs_t for icedata.c)\n\
  -cMX[p] cycle analysis upto MX-cycle [0], negative MX = include branched\n\
          p=exclude periodic cycles (this test fails for V)\n\
  -EE     elst field, only with MC + -o\n\
  -SSHIFT change zshift/x = cos(beta) to fit box (hint: use -S-.333.. for V)\n\
  -w      rewire H-bonds, provide file ICENAME.wire:\n\
          ! comment\n\
          d site1 site2 ! to remove bond site1-site2\n\
          ...\n\
          c site1 site2 ! to add bond site1-site2\n\
          ...\n\
Environment:\n\
  ICENAME: base name [default = derived from parameters]\n\
Algorithm:\n\
  Either MC (defect: E=1) or simulated annealing of initial undefined-oriented\n\
  water molecules until all Bjerrum defects annihilate\n\
Caveats:\n\
  * the output structure is topologically correct, but not exact\n\
  * water is tetrahedral (angle HOH=109.47 deg), HO=0.96, OM=0.18, OL=0.80\n\
  * ice V connect algorithm is inefficient for large lattices\n\
Example (basic cell repeated 5x in x, but written as y):\n\
  ice Ih 5y 3z 3x -sHOMH -T0\n\
  ice Ih 3 2 2 -T0.5\n\
See also:\n\
  lattice iceII iceI+III(old version) Ihinbox(with rotations)\n\
",SITES);
    exit(0); }

  initscroll(0);
  stringinit(0,8);
  signal(SIGINT,asksig);

  if (!strcmp(arg[1],"2D")) ice=&Ice2D;
  else if (!strcmp(arg[1],"Ih")) ice=&IceIh;
  else if (!strcmp(arg[1],"Ic")) ice=&IceIc;
  else if (!strcmp(arg[1],"II")) ERROR(("\
Ice II not supported because proton ordered.\n\
*** Use prepared iceIIhexcell.plb obtained using iceII.c. Example:\n\
  plbreplicate iceIIhexcell.plb aux.plb 2 1 4\n\
  wat2wat OHH HOMH < aux.plb > IceII214in.plb"))
  else if (!strcmp(arg[1],"III")) ice=&IceIII;
  else if (!strcmp(arg[1],"V")) checkodd=0,ice=&IceV;
  else if (!strcmp(arg[1],"VI")) ice=&IceVI;
  else if (!strcmp(arg[1],"VII")) ice=&IceVII;
  else if (!strcmp(arg[1],"ClI")) ice=&ClI;
  else if (!strcmp(arg[1],"ClII")) ice=&ClII;
  else if (!strcmp(arg[1],"ClH")) ice=&ClH;
  else ERROR(("%s: unsupported ice/clathrate",arg[1]))

  prt("ice %s basic cell:",arg[1]);
  prt("  %d molecules",ice->n);
  if (ice->unit>0) {
    prt("  a=[%.7f %.7f %.7f]*%.7g",VARG(ice->a),ice->unit);
    VO(ice->a,*=ice->unit) }
  V=PROD(ice->a)/ice->n;
  prt("  a=[%.7f %.7f %.7f] A",VARG(ice->a));
  prt("  volume/molecule=%g A^3  density=%g kg/m^3",V,29915.07113/V);

  if (narg<3) exit(0);

  nx=atoi(arg[2]); ix=indx(arg[2],0);
  ny=atoi(arg[3]); iy=indx(arg[3],1);
  nz=atoi(arg[4]); iz=indx(arg[4],2);

  if (ix+iy+iz!=3) Error("bad transformation; either omit modifiers xyz or use each just once");

  N=ice->n*nx*ny*nz;
  L[ix]=ice->a[0]*nx;
  L[iy]=ice->a[1]*ny;
  L[iz]=ice->a[2]*nz;

  prt("basic cell multiplied by %dx%dx%d = %g %g %g lattice\n\
%d vertices, transformation=%c%c%c",
         nx,ny,nz,         L[ix],L[iy],L[iz],       N,               ix+'x',iy+'x',iz+'x');

  {
    double minL=fmin(L[0],L[1]);
    double maxL=fmax(L[0],L[1]);
    Min(minL,L[2])
    Max(maxL,L[2])
    prt("max aspect=%5.3f  quality=%5.3f (minL^3/N/V, cube=1)",
        maxL/minL,Cub(minL)/N/V);
  }

  args=strdup(string("%s %s %s %s",arg[1],arg[2],arg[3],arg[4]));

  loop (iarg,5,narg)
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 's': SITES=arg[iarg]+2; break;
        case 'v': quiet=0; break;
        case 'z': RNDSEED=atoi(arg[iarg]+2); break;
        case 'm': Mfrom=atof(arg[iarg]+2); break;
        case 'M': Mto=atof(arg[iarg]+2); break;
        case 'a': T=-9; maxndis=atoi(arg[iarg]+2); break;
        case 'o':
          sphere=atoi(arg[iarg]+2);
          if (sphere<0) nomit=-sphere,sphere=0;
          break;
        case 'T': T=atof(arg[iarg]+2); usebeta=0; break;
        case 'b': T=1/atof(arg[iarg]+2); usebeta=-1; break;
        case 'B':
          num=atof(arg[iarg]+2);
          usebeta=atoi(strchr(arg[iarg],'/')+1);
          T=(double)usebeta/(num*log(N));
          break;
        case 'n': no=atoi(arg[iarg]+2); break;
        case 'i': noint=atoi(arg[iarg]+2); break;
        case 'l': lag=atoi(arg[iarg]+2); break;
        case 'e': eq=atoi(arg[iarg]+2); break;
        case 'g': gennbrs++; break;
        case 'c':
          maxnest=atoi(arg[iarg]+2);
          pertest=strchr(arg[iarg]+2,'p')!=NULL;
          break;
        case 'S': ice->shift=atof(arg[iarg]+2); break;
        case 'E': E=atof(arg[iarg]+2); break;
        case 'w': wire=arg[iarg]+2; break;
        default: ERROR(("%s: unknown option",arg[iarg]))
      }
    else ERROR(("%s: option should start with -\n\
*** (first 4 parameters, ICE NX NY NZ, are position without -)",arg[iarg]))


  ICENAME=getenv("ICENAME")
           ?getenv("ICENAME")
           :strdup(string("%s%d%d%d%s",arg[1],nx,ny,nz,SITES));

  if (T>0)
    basename=strdup(usebeta==0 ? string("%sT=%g",ICENAME,T) :
                    num ? string("%sB%d_%d",ICENAME,num,usebeta) :
                    string("%sb=%g",ICENAME,1/T));
  else
    basename=ICENAME;

  branched=maxnest<0;
  maxnest=abs(maxnest);
  if (maxnest>NEST || maxnest==1 || maxnest==2) ERROR(("maxnest=%d out of range",maxnest))

  prt("using %s as the basename",basename);

  if (ICENAME[0]) out=fopen(string("%s.prt",basename),"wt"); // else stdout

  prt_("ice" VERSION );
  loop (i,0,narg) prt_(" %s",arg[i]);
  _n

  if (ice==&Ice2D && nz>1) WARNING(("several independent sheets of 2D ice\n\
*** (2D ice has the unit cell in x,y)"))

  if (gennbrs && ice==&IceVII)
    WARNING(("-g with ice VII not allowed - will be wrong.\n\
*** Cf. ice-gentab.c and iceV1"))

  Mmid=(Mfrom+Mto)/2;
  Merr=(Mto-Mfrom)/1.999999999999;
  if (Merr<=0) ERROR(("bad or zero dipole moment range -m%g -M%g",Mfrom,Mto))

  if (RNDSEED) {
    rndinit(7,RNDSEED);
    prt("RNDSEED=%d 1st rnd=%11.9f",RNDSEED,rnd()); }
  else {
    rndinit(0,0);
    prt("RNDSEED=TIME 1st rnd=%11.9f",rnd()); }

  allocarrayzero(cfg,N);
  prt("N=%d",N);

  /* force using distances to generate neighbors */
  if (gennbrs) ice->nbrs=NULL;
  if (!ice->nbrs) prt("NOTE: H-bond structure from O-O distances (slow for large lattice)");

  nn=0;
  loop (i,0,nx)
    loop (j,0,ny)
      loop (k,0,nz)
        loop (n,0,ice->n) {
          if (ice->unit>0) {
            cfg[nn].r[ix]=ice->a[0]*(i+k*ice->shift+ice->r[n][0]);
            cfg[nn].r[iy]=ice->a[1]*(j+ice->r[n][1]);
            cfg[nn].r[iz]=ice->a[2]*(k+ice->r[n][2]); }
          else {
            cfg[nn].r[ix]=ice->a[0]*(i+k*ice->shift)+ice->r[n][0];
            cfg[nn].r[iy]=ice->a[1]*j+ice->r[n][1];
            cfg[nn].r[iz]=ice->a[2]*k+ice->r[n][2]; }
          if (ice->nbrs) {
            cfg[nn].nnbr=4;
            loop (l,0,4)
              cfg[nn].nbr[l]=(ice->nbrs[n][l].nbr +ice->n * (
               (k+nz+ice->nbrs[n][l].z)%nz
                 +nz*( (j+ny+ice->nbrs[n][l].y)%ny
                   +ny*( (i+nx+ice->nbrs[n][l].x)%nx ))))%N; }
          nn++; }

  if (nn!=N) ERROR(("nn(%d)!=N(%d)",nn,N))


  // loop (i,0,N) loop (j,0,N) prt("%g",ddist(cfg[i].r,cfg[j].r));

  if (!ice->nbrs) {
    /* neigbors from distances */
    double aa=Sqr(ice->OOHB);

    loop (i,0,N) loop (j,0,i)
      if (ddist(cfg[i].r,cfg[j].r)<aa) {
        if (cfg[i].nnbr<4) cfg[i].nbr[cfg[i].nnbr++]=j;
        else fprintf(stderr,"%d %d: WARNING more than 4 neighbors omitted\n",i,j);
        if (cfg[j].nnbr<4) cfg[j].nbr[cfg[j].nnbr++]=i;
        else fprintf(stderr,"%d %d: WARNING more than 4 neighbors omitted\n",i,j); }

    loop (i,0,N) if (cfg[i].nnbr!=4)
      fprintf(stderr,"%d: WARNING only %d neighbors\n",i,cfg[i].nnbr); }

  if (wire) rewire(basename);

  /* check neighbor table consistency */
  loop (i,0,N) loop (k,0,cfg[i].nnbr) {
    int n=cfg[i].nbr[k],j;
    if (n<0 || n>=N) ERROR(("%d %d neighbor index out of range",i,n))
    loop (j,0,cfg[n].nnbr) if (cfg[n].nbr[j]==i) goto OK;
      WARNING(("neigbor table inconsistency: %d->%d but not back",i,n))
    OK:; }

  if (!quiet)
    loop (i,0,N) {
      prt("%d:(%g,%g,%g) nbrs=%d %d %d %d",
          i,
          cfg[i].r[0],cfg[i].r[1],cfg[i].r[2],
          cfg[i].nbr[0],cfg[i].nbr[1],cfg[i].nbr[2],cfg[i].nbr[3]); }

#if 0
  nest=maxnest;
  loop (i,0,N) printf("|0..%d|=%d\n",i,mindist(0,i));
  exit(0);
#endif

  /* cycle analysis */
  if (maxnest) {
    /* central cell - makes sense for ice V */
    int ii=(nx/2*ny+ny/2)*nz+nz/2;

    loop (i,ii,ii+ice->n) {
      if (!quiet) fprintf(stderr,"site %d\n",i);
      nest=0;
      cyc[0]=i;
      loop (j,0,4) propagate(cfg[i].nbr[j]); }
    loopto (i,3,maxnest) if (hist[i])
      printf("%d %g = %d/%d\n",i,(double)hist[i]/(ice->n*2*i),Num(hist[i],ice->n*2*i),Den(hist[i],ice->n*2*i)); }

  /* find ALL disordered configurations, option -a */
  if (T==-9) {
    info=maxndis;
    allocarray(cfg0,N);
    copyarray(cfg0,cfg,N);
    nbrsof(0);
    prt("%lld proton disordered configurations found",ndis); }

  else if (T>0 && nomit==0) {
    /* equilibrium Monte Carlo simulation */
    int i;
    double U=0;
    FILE *cpa=fopen(string("%s.cpa",basename),"wt");

    if (E) ERROR(("elst field not supported here, needs -o"))

    StaSet(noint,lag,2,12);
    //    if (eq<0) eq=(201+60*exp(2./3/T))/noint; /* diffusion model V2 */
    if (eq<0) {
      if (ice==&Ice2D) eq=(201-eq*exp(1/T))/noint; /* 2D ice */
      else eq=(201-eq*exp(2./3/T))/noint; } /* diffusion model: V2.1 */


    prt("EQ %d cycles by %d*%d steps started at %s",eq,N,noint,datetime(0,0));

    loop (i,0,eq) {
      U=MC(noint,T);
      fprintf(cpa,"%d %g EQ\n",i,U);
      if (!quiet) prt("%d %g EQ",i,U); }
    if (eq) prt("MC equilibrated (%d cycles): U=%g T=%g noint=%d",eq,U,T,noint);
    else WARNING(("no equilibration (eq=0)"))

    prt("MC %d cycles by %d*%d steps started at %s",no,N,noint,datetime(0,0));
    fflush(out);
    fprintf(cpa,"# equilibration finished, will print also part of the run\n");
    loop (i,0,no) {
      U=MC(noint,T);
      if (cpa) {
        fprintf(cpa,"%d %g MC\n",i+eq,U);
        if (i>eq) { fclose(cpa); cpa=NULL; } }
      if (!quiet) prt("%d %g MC",i,U);
      StaAdd("U",U);
      if (sig) break; }
    if (cpa) fclose(cpa);
    prt("MC stopped at %s",datetime(0,0));
    StaPrintAll("");
    prt("%.10g %.10g %g %g %d %d %d  # 1/T <U> stderr T N no noint",
        1/T,StaMean("U"),StaStdErr("U"),T,N,StaN("U"),noint);
    StaSave(string("%s.sta",basename,1/T));
  }

  else if (T>0 && nomit!=0) {
    /* special: Monte Carlo simulation of defects */
    int i;
    double U=0;
    FILE *cpa=fopen(string("%s.cpa",basename),"wt");

    prt("DL %d cycles by %d*%d steps started at %s",eq,N,noint,datetime(0,0));

    loop (i,0,eq) {
      U=MCfield(noint,T,E);
      fprintf(cpa,"%d %g DL\n",i,U);
      if (!quiet) prt("%d %g DL",i,U); }

    if (eq) prt("MC Bjerrum defect layers generated (%d cycles): U=%g T=%g noint=%d",eq,U,T,noint);
    else  WARNING(("no equilibration (eq=0)"))

    i=nomit; nomit=0; writecfg(ICENAME[0]?ICENAME:"aux"); nomit=i;

    prt("MC %d cycles by %d*%d steps started at %s",no,N,noint,datetime(0,0));
    fflush(out);
    fprintf(cpa,"# Bjerrum defect layers generated\n");

    loop (i,0,no) {
      U=MCinside(noint,T);
      if (cpa) fprintf(cpa,"%d %g MC\n",i+eq,U);

      //loop (j,nomit/2,N-nomit/2) {
      loop (j,0,N) {
        int sg=charge(j);
        if (sg) {
          prt("%d %g %g %g D%d",i,VARG(cfg[j].r),sg); } }

      if (!quiet) prt("%d %g MC",i,U);
      if (sig) break; }
    if (cpa) fclose(cpa);
    prt("MC stopped at %s",datetime(0,0));
  }

  /* simulated annealing (no Bjerrum defect) */
  else if (T==0) {
    annealing();
    writecfg(ICENAME[0]?ICENAME:"aux");
  }
  else
    fprintf(stderr,"nothing to do, set option -T or -b or -B\n");

  if (!quiet)
    loop (i,0,ice->n) {
      prt_("{ ");
      loop (k,0,4) {
        prt_("{%d,%d,%d,%d}, ",
             cfg[i].nbr[k]%ice->n,
             shift(i,cfg[i].nbr[k],0,ice->n),
             shift(i,cfg[i].nbr[k],1,ice->n),
             shift(i,cfg[i].nbr[k],2,ice->n)); }
      prt("},"); }

  if (out) fclose(out);

  return 0;
}
