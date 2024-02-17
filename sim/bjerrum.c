#include "minimize.h"

//#define PER

/* max permitted O-O distance */
#define OO 5.5

/* H-bond limit */
#define HO 2.5

static void nearestdr(vector dr) /******************************** nearestdr */
/* must solve out-of-cell coordinates */
{
  int k;

  loop (k,0,3) {
    while (dr[k]>box.Lh[k]) dr[k]-=box.L[k];
    while (dr[k]<-box.Lh[k]) dr[k]+=box.L[k]; }
}

static struct bjerrum_t {
  int D1[2],D2[2]; /* [0]=D-defect, [1]=L-defect; D1,D2=both molecules */
  ToIntPtr A;      /* configuration */
  double sig;      /* Gauss' sigma */
  double isigq2;   /* 1/(2*sigma^2) */
  double cutoff;   /* Gauss renormalized at cutoff */
  double cq;       /* cutoff^2 */
  double shift;    /* Gauss at cutoff */
  double ishift;   /* 1/(1-shift) */
  int it;          /* # of minimization iterations counter */
} bjerrum={{-1,-1},{-1,-1},NULL,1,0};

static double renormgauss(double rr) /************************** renormgauss */
/* rr =|r1-r2|^2 */
{
  if (rr<bjerrum.cq) return (exp(-rr*bjerrum.isigq2)-bjerrum.shift)*bjerrum.ishift;
  else return 0;
}

static double maxm=0,minm=9e9;
static int iO,iH1,iH2;
static int (*hb)[2];

void bjerrum_read_oo(void) /******************************** bjerrum_read_oo */
{
  int ns,i;
  molecule_t *m1;
  siteinfo_t *si1;


  if (maxm==0) {
    FILE *oo;
    char line[64];
    int n=0,a,b;

    loop (i,0,2) bjerrum.D1[i]=-1,bjerrum.D2[i]=-1;

    allocarrayzero(hb,No.N*2);

    /* O,H sites */
    m1=molec;
    ns=m1->ns;
    si1=spec[m1->sp]->si;
    loop (i,0,ns) {
      if (si1[i].mass>maxm) {
        maxm=si1[i].mass; iO=i; }
      if (si1[i].mass>0 && si1[i].mass<minm) {
        minm=si1[i].mass*1.0001;
        iH1=iH2,iH2=i; } }
    prt("molecule 0 sites: O=%d H=%d,%d",iO,iH1,iH2);

    oo=fopen(Fn("oo"),"rt");
    if (!oo) ERROR(("bjerrum: neighbor list (H-bonded water-water pairs) %s missing",lastFn))

    while (fgets(line,64,oo)) {
      if (n>=2*No.N) ERROR(("%s: too many lines",lastFn))
      sscanf(line,"%d%d",&a,&b);
      if (a<0 || a>=No.N || b<0 || b>=No.N) ERROR(("%s: water molecule numbers out of range"))

      hb[n][0]=a;
      hb[n][1]=b;
      if (strchr(line,'D')) bjerrum.D1[0]=a,bjerrum.D2[0]=b;
      if (strchr(line,'L')) bjerrum.D1[1]=a,bjerrum.D2[1]=b;
      n++; }

    if (n<2*No.N) ERROR(("%s: too few lines (called because bj.mode&1)",lastFn))
    prt("%s read: D-defect: D1=%d D2=%d, L-defect: L1=%d L2=%d",
        lastFn,bjerrum.D1[0],bjerrum.D2[0],bjerrum.D1[1],bjerrum.D2[1]);

    fclose(oo); }
}

void bjerrum_topology(ToIntPtr A) /************************ bjerrum_topology */
/* bk.mode=4: OUT OF ORDER or only partially working */
{
  int n,n1,n2;
  molecule_t *m1;
  vector *r1,*r2,dr;
  int nHB1,nHB2,nL=0,nL1=-1,nL2=-1,nD=0,nD1=-1,nD2=-1;

  /* finding a defect - topology way */
  prt("bjerrum: t=%f",t);

  if (!hb) ERROR(("no hb: bjerrum_read_oo should be called first"))

  loop (n,0,2*No.N) {
    n1=hb[n][0];
    n2=hb[n][1];
    m1=molec+n1;
    r1=rof(m1,A->rp);
    r2=rof(molec+n2,A->rp);

    VVV(dr,=r1[iO],-r2[iO]) nearestdr(dr);

    if (SQR(dr)>Sqr(OO)) prt("O%d-O%d=%g exceeded the limit (%g)",n1,n2,sqrt(SQR(dr)),(double)OO);

    nHB1=0;
    VVV(dr,=r1[iO],-r2[iH1]) nearestdr(dr); if (SQR(dr)<Sqr(HO)) nHB1++;
    VVV(dr,=r1[iO],-r2[iH2]) nearestdr(dr); if (SQR(dr)<Sqr(HO)) nHB1++;

    nHB2=0;
    VVV(dr,=r2[iO],-r1[iH1]) nearestdr(dr); if (SQR(dr)<Sqr(HO)) nHB2++;
    VVV(dr,=r2[iO],-r1[iH2]) nearestdr(dr); if (SQR(dr)<Sqr(HO)) nHB2++;

    if (nHB1==0 && nHB2==0) {
      /* L-defect */
      nL++; nL1=n1; nL2=n2; }

    if (nHB1==1 && nHB2==1) {
      /* D-defect */
      nD++; nD1=n1; nD2=n2; } }

  prt("!nD=%d nD1=%d nD2=%d   nL=%d nL1=%d nL2=%d BJERRUM4",
        nD,   nD1,   nD2,     nL,   nL1,   nL2);
}

double SumGauss(vector R) /**************************************** SumGauss */
/*
   Sum of charges multiplied by Gaussian-like weights; w(center)=1.
   w = [Gaussian(distance_from_R,bj.sig)-cutoff]/[1-cutoff] (to be smooth).
   Charge is converted to e.
   D-charge is returned with opposite sign (to ease minimization).
   Only large cutoff/bj.sig are recommended.
   Since we use Minimize, for bj.sig>0 (D-defect), -sum is returned!
   WARNING: sign bugs in V3.3d and older, L-D swapped and inefficient
*/
{
  double sumq=0,rr;
  int n,i;
  molecule_t *m;
  siteinfo_t *si;
  vector *r,dr;
#ifdef POLAR
  double auxsum;
  vector *rpol;
#endif /*# POLAR */

  bjerrum.isigq2=0.5/Sqr(bjerrum.sig);
  bjerrum.cutoff=fmin(box.Lh[0],box.Lh[1]);
  bjerrum.cutoff=fmin(bjerrum.cutoff,box.Lh[2]);
#ifdef PER
  bjerrum.cutoff*=1.73;
#endif
  bjerrum.cq=Sqr(bjerrum.cutoff);
  bjerrum.shift=exp(-bjerrum.cq*bjerrum.isigq2);
  bjerrum.ishift=1/(1-bjerrum.shift);

  bjerrum.it++;

  loop (n,0,No.N) {
    m=molec+n;
    si=spec[m->sp]->si;
    r=rof(m,bjerrum.A->rp);
#ifdef POLAR
    rpol=polarrof(m,bjerrum.A->rp);
    loop (i,0,m->ns) {
      if (si[i].charge || si[i].chargepol) {
        /* nearest-image */
        VVV(dr,=R,-r[i])
        nearestdr(dr);
        rr=SQR(dr);
        //        sumq+=exp(-rr/sigq)*si[i].charge;
        auxsum=renormgauss(rr)*si[i].charge;
        VV(dr,+=rpol[i])
        rr=SQR(dr);
        //        sumq+=exp(-rr/sigq)*si[i].chargepol;
        sumq+=auxsum+renormgauss(rr)*si[i].chargepol; }
    }
#else /*# POLAR */
    loop (i,0,m->ns) if (si[i].charge) {
#  ifdef PER
      /* sum over 27 images: -1,0,+1 SLOW */
      /* DO NOT USE: the effective charge begins to decrease for large sigma's */
      int j[3];

      loopto (j[0],-1,1)
        loopto (j[1],-1,1)
          loopto (j[2],-1,1) {
            VVVVV(dr,=R,-r[i],+j,*box.L)
            rr=SQR(dr);
            //      sumq+=exp(-rr/sigq)*si[i].charge;
            sumq+=renormgauss(rr)*si[i].charge; }
#  else
      /* nearest-image (standard) */
      VVV(dr,=R,-r[i])
      nearestdr(dr);
      rr=SQR(dr);
      //      sumq+=exp(-rr/sigq)*si[i].charge;
      sumq+=renormgauss(rr)*si[i].charge;
#  endif /*#!0 */
    }
#endif /*#!POLAR */
  }

  /* WARNING: sign bugs in V3.3d and older, L-D swapped and inefficient */
  if (bjerrum.sig<0) return sumq/electron;
  else return -sumq/electron; /* D-defect charge changed sign */
}

void bjerrum_gauss(ToIntPtr A,double from,double to,double q,double eps)
                                                   /********** bjerrum_gauss */
{
  int i,k,idefect;
  molecule_t *m;
  vector *r,R0,R;
  double sig,sumq;
  double sig_to=cbrt(PROD(box.L))*to;

  bjerrum.A=A;

  if (sizeof(REAL)!=sizeof(double)) DISASTER(("Minimize: wrong REAL size"))

  loop (idefect,0,2) {
    /* idefect=0: D-defect, idefect=1: L-defect */

    if (bjerrum.D1[idefect]<0)
      VV(R0,=idefect*box.Lh)
    else {
      m=molec+bjerrum.D1[idefect];
      r=rof(m,A->rp);
      VO(R0,=0)
      loop (i,0,m->ns) VV(R0,+=r[i])
      VO(R0,/=m->ns) }

    if (bjerrum.D2[idefect]<0)
      VV(R,=idefect*box.Lh)
    else {
      m=molec+bjerrum.D2[idefect];
      r=rof(m,A->rp);
      VO(R,=0)
      loop (i,0,m->ns) VV(R,+=r[i])
      VO(R,/=m->ns) }

    /* CoM of two molecules comprising a defect */
    loop (k,0,DIM) {
      if (R0[k]>R[k]+box.Lh[k]) R0[k]-=box.L[k];
      if (R0[k]<R[k]-box.Lh[k]) R0[k]+=box.L[k];
      R0[k]=(R0[k]+R[k])/2;
      if (R0[k]<0) R0[k]+=box.L[k];
      if (R0[k]>box.L[k]) R0[k]-=box.L[k]; }

    prt("! Gauss %c-defect D1=%d D2=%d R=[%g %g %g]","DL"[idefect],bjerrum.D1[idefect],bjerrum.D2[idefect],VARG(R0));


#if 0
    /* DEBUG-3D data, D-defect only: */
    bjerrum.sig=from;
    {
      int i[3];
      double d=from*0.5,s;
#define RG 16

      loopto (i[0],-RG,RG) {
        R[0]=R0[0]+d*i[0];
        loopto (i[1],-RG,RG) {
          R[1]=R0[1]+d*i[1];
          loopto (i[2],-RG,RG) {
            if (SQR(i)<RG*(RG+1)) {
              R[2]=R0[2]+d*i[2];
              s=SumGauss(R);
              prt("%d %d %d  %.10f  %g %g %g  BJERRUMDEBUG",i[0],i[1],i[2],s,VARG(R)); } } } } }
    return;
#endif

    for (sig=from; sig<sig_to; sig*=q) {

      VV(R,=R0)

      bjerrum.sig=idefect?-sig:sig;
      /* WARNING: sign bugs in V3.3d and older, L-D swapped and inefficient */
      bjerrum.it=0;

#if 0
      {
        /* prescan over the whole box INEFFICIENTLY CALLED TWICE */
        vector D,RR;
        double d=sig/4,s,m=9e9;
        int i;

        loop (i,0,3) D[i]=box.L[i]/ceil(box.L[i]/d);

        for (R[0]=D[0]/2; R[0]<box.L[0]; R[0]+=D[0]) {
          for (R[1]=D[1]/2; R[1]<box.L[1]; R[1]+=D[1]) {
            for (R[1]=D[1]/2; R[1]<box.L[1]; R[1]+=D[1]) {
              s=SumGauss(R);
              if (s<m) { m=s; VV(RR,=R) }
              if (!idefect) s=-s;
              prt("%g %g %g  %.10f %g BJERRUMSCAN", VARG(R), s, sig); } } }
        VV(R,=RR)
      }
#endif

#if 1
      {
        /* prescan in the 2*sigma vicinity of the position read from .oo */
        int i[3];
        vector RR;
        double d=sig/4,s,m=9e9;
#define RG 8

        loopto (i[0],-RG,RG) {
          R[0]=R0[0]+d*i[0];
          loopto (i[1],-RG,RG) {
            R[1]=R0[1]+d*i[1];
            loopto (i[2],-RG,RG) {
              if (SQR(i)<RG*(RG+1)) {
                R[2]=R0[2]+d*i[2];
                s=SumGauss(R);
                if (s<m) { m=s; VV(RR,=R) } } } } }
        VV(R,=RR)
      }
#endif

      //      Minimize(MINIMIZE_MC,3,R,SumGauss,16,eps,bjerrum.sig*0.1,9);
      Minimize(MINIMIZE_NR,3,R,SumGauss,20,eps,0.001,2);

      sumq=SumGauss(R);

      /* WARNING: sign bugs in V3.3d and older, L-D swapped and inefficient */
      prt("%9.6f %9.6f [%10.6f %10.6f %10.6f] %d BJERRUM2",
          sig,idefect?sumq:-sumq,VARG(R),bjerrum.it); }

  prt("\t\tBJERRUM2"); }
}
