/*
This is the standard (not particle-mesh) version of Ewald summation.
This code uses ideas of J. Perram (making use of symmetry).
See also macsimus/c/ewald.c (fool-proof slow standalone version for debugging).

UPDATE 6/2022
  - background charges, Yeh-Berkowitz revisited

UPDATE 3/2011
  - fluctuating charges

UPDATE 2/2010
  - parallel support rehacked for pthreads
  - CLOCK, TIM_FILE removed

UPDATE 12/2008:
  - code cleaned and a bit optimized
  - obsolete DOS and HIT switches removed

UPDATE 9/2008:
  - obsolete switch SITEEWALD removed

UPDATE 2/2005:
  - L[3] (not cube) implemented
  - obsolete switches PAR, FASTRAM, POINTERS removed
  - k-vectors are allocated in a table and L must not change a lot

FIRST VERSION 1991

  pass -1: initialization
        0: initialization & protocol printed to out
        1: sums Q(k)
        2: energy (return value) and forces calculated
  frp  (in format described by cfg): to sum up forces
  rp   (in format described by cfg): coordinates of sites
  if pass==2 then energy is returned

WARNINGS: if alpha<0 then masses are used (and copied to the configuration!)
instead of charges.  This is for calculating the structure factor from
stored configurations and cannot be used for simulating!!!
Sequential version only!

Calling scheme:
^^^^^^^^^^^^^^^
vector *forces,*configuration;
Ewald(0,<any>,<any>);  (* initialization - may be repeated *)
repeat {
  Ewald(1,<any>,configuration->rp);
  E=Ewald(2,forces->rp,configuration->rp);
}

Formulas:
^^^^^^^^^
U =
    SUM_j<l qj ql Erfc(alpha rjl)/rjl [ r-space term: see module interpol.c ]
  + SUM_k!=0 exp(-pi^2/alpha^2*kk) / (2 V pi kk) |Q(k)|^2
  + 2 pi / (2 epsinf+1) *M^2/V
  - alpha/sqrt(pi)*SUM_j qj^2
(is divided by 4 pi eps0 in SI)
where
  kk denotes (kx/Lx)^2+(ky/Ly)^2+(ky/Ly)^2
  Q(k) = SUM_j qj exp[2 pi i (kx/Lx.xj+ky/Ly.yj+kz/Lz.zj)]
  M = SUM_j qj rj
  V = Lx Ly Lz

More:
^^^^^
For ready-to-code formulas, see Aguado, Madden: JCP 119, 7471 (2003).
Differences:
* h=2*pi*k/L (vector, by components)
* r-space parameter is called kappa (here: alpha, while kappa means k-space)
* constant term -alpha/sqrt(pi)*SUM_j qj^2 is missing in U

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "ground.h"
#include "sds.h"
#include "simglob.h"

#include "ewald.h"
#include "units.h"
#include "rhs.h"
#include "norm.h" /* iscube() only */

/* theses headers are in  macsimus/cook/ (to be cleaned?) */
#include "forces.h"

#include "cputime.h"
#include <time.h>

#include "simgear.h"
#ifdef POLAR
/* selffield() needed */
#  include "constrd.h"
#endif /*# POLAR */

#include "simdef.h" /* initrspace() only, more for GAUSSIANCHARGES */

/* global : */
Q_t *Q;

static double estimr,estimk;

/*** complex arithmetic ***/
#define CxE(A) { A.re=(c=A.re)*E->re-A.im*E->im; A.im=c*E->im+A.im*E->re; }
#define CxiE(A) { A.re=(c=A.re)*E->re+A.im*E->im; A.im=A.im*E->re-c*E->im; }
#define Cadd(A,B) { A.re+=B.re; A.im+=B.im; }
#define Csqr(A) (A.re*A.re+A.im*A.im)
#define Csqr2(A,B) (A.re*B.re+A.im*B.im)

#if PARALLEL==1 || PARALLEL==3
/* PARALLEL==3 no longer active */

/*** loops over charges ***/
#  define iqloop for (iqt=iqt0; iqt<iqt1; iqt++)
#  define iqfloop for (iqt=iqt0; iqt<iqt1; iqt++)
#  define iqeloop(E0) for (iqt=iqt0; E=E0+iqt,iqt<iqt1; iqt++)

#else /*# PARALLEL==1 || PARALLEL==3 */

/*** serial: loops over charges ***/
#  define iqloop for (iqt=0; iqt<Nq; iqt++)
#  define iqfloop for (iqt=0; iqt<Nq; iqt++)
#  define iqeloop(E0) for (iqt=0; E=E0+iqt,iqt<Nq; iqt++)
#endif /*#!PARALLEL==1 || PARALLEL==3 */

typedef struct {
  complex
    x,   /* q*exp[2*PI*i*x*kx/Lx] */
    xy,  /* q*exp[2*PI*i*(x*kx/Lx+y*ky/Ly)] */
    xY,  /* q*exp[2*PI*i*(x*kx/Lx-y*ky/Ly)] */
    xyz, /* q*exp[2*PI*i*(x*kx/Lx+y*ky/Ly+z*kz/Lz)] */
    xyZ, /* q*exp[2*PI*i*(x*kx/Lx+y*ky/Ly-z*kz/Lz)] */
    xYz, /* q*exp[2*PI*i*(x*kx/Lx-y*ky/Ly+z*kz/Lz)] */
    xYZ; /* q*exp[2*PI*i*(x*kx/Lx-y*ky/Ly-z*kz/Lz)] */
#ifdef GAUSSIANCHARGES
  int site; /* site type, QQTAB style needed */
  //  double sigmaq; /* Gauss sigma^2 - enable with optimization? */
#  if defined(POLAR) && POLAR&32
#    error "POLAR&32 + GAUSSIANCHARGES not implemented"
#  endif /*# defined(POLAR) && POLAR&32 */
#else /*# GAUSSIANCHARGES */
#  if defined(POLAR) && POLAR&32
  double q;     /* q; in pass=2 all q-factors above are 1 */
  double dummy; /* cache-line optimization */
#  else /*# defined(POLAR) && POLAR&32 */
  complex dummy; /* cache-line optimization */
#  endif /*#!defined(POLAR) && POLAR&32 */
#endif /*#!GAUSSIANCHARGES */
} qtab_t;

/*** some ugly globals with sums of charges etc. ***/
struct charges_s charges;

/* moved here because of PARALLEL */
static int Nq;               /* # of all charges */
static int kxmax;
static complex *ex,*ey,*ez;  /* array[Nq] : exp(2*PI*i*r/L) */
static qtab_t *qtab;

#if PARALLEL==1 || PARALLEL==3
#  ifdef GAUSSIANCHARGES
static struct Qk_s **Qkcell;
#  else /*# GAUSSIANCHARGES */
static complex **Qkcell;
#  endif /*#!GAUSSIANCHARGES */
struct pll_ewald_s *pll_ewald;
#endif /*# PARALLEL==1 || PARALLEL==3 */

#if PARALLEL==2
struct pll_ewald_s *pll_ewald;
#endif /*# PARALLEL==2 */

static real *ekk;    /* array[Nkk] : table of ...exp(...k*k) factors */
#ifdef GAUSSIANCHARGES
static real **sekk; /* array[Nkk][nsites]: as above, additional factor for Gaussian charge */
static real **ssekk; /* array[Nkk][nsites]: *sigma_i^2 (for pressure tensor) */
#endif /*# GAUSSIANCHARGES */
static struct ktab_s {
  int kymax;
  int *kzmax; /*[kymax(incl.)]*/
} *ktab /*[kxmax(incl.)]*/ ;
static int ktabsize; /* =kxmax+1 */
static vector *F;    /* array[Nq] */
#if defined(POLAR) && POLAR&32
static double *Phi; /* array[Nq] */
/* QTIMES is a POLAR&32 hack because qtab does not contain charges with FQ */
#  define QTIMES qtab[iqt].q*
#else /*# defined(POLAR) && POLAR&32 */
#  ifdef GAUSSIANCHARGES
// TO BE OPTIMIZED - moved to qtab
#    define QTIMES 0.5*
#  else /*# GAUSSIANCHARGES */
#    define QTIMES /* dummy */
#  endif /*#!GAUSSIANCHARGES */
#endif /*#!defined(POLAR) && POLAR&32 */
static double Ediag; /* systematic part of k-space cutoff error of energy */
static real KX,PIaq;


#if PARALLEL==1

struct par1_ew_s {
  struct sfr_s *sfr;
  struct sf3d_s *sf3d;
  int pass;
} par1_ew;

void *parEwald(void *arg) /**************************************** parEwald */
{
  complex *E,A,B,C,D;
#  ifdef GAUSSIANCHARGES
  complex AT,BT,CT,DT;
  complex AS,BS,CS,DS;
#  endif /*# GAUSSIANCHARGES */
  real FA,FB,FC,FD,e,c,kk;
  int iqt,kx,ky,kz,iQ=0,ikk=0,jkk,kymax,kzmax;
  int iqt0,iqt1;
  vector ek;
#  if PRESSURETENSOR&PT_VIR
  double ptf,E1,E1x;
#  endif /*# PRESSURETENSOR&PT_VIR */
  int pass=par1_ew.pass;
  int ith=(pthread_t*)arg-No.thread;
  clock_t clock0=0;

  if (partimes.kspace) clock0=clock();

  iqt0=Nq*ith/No.th;
  iqt1=Nq*(ith+1)/No.th;

#  define Energy pll_ewald[ith].E
#  define PVIR pll_ewald[ith].Pvir
#  define Qk Qkcell[ith]

  if (pass==1) {
#  include "ewpass1.c"
  }
  if (pass==2) {
    if (measure) {
#  include "ewpass2m.c"
    }
    else {
#  include "ewpass2.c"
    } }

#  undef Qk
#  undef Energy
#  undef PVIR
  /* will be used serially in Pvir contributions due to M */
#  define PVIR En.Pvir

  if (partimes.kspace) pll_ewald[ith].t+=clock()-clock0;

  return arg;
}
#endif /*# PARALLEL==1 */

#if COULOMB<0
static real machineprec;

static int almostint(real x) /************************************ almostint */
/* returns true if x is very close to an integer */
{
  x-=(int)x;
  return fabs(0.5-fabs(x-0.5))<machineprec;
}
#endif /*# COULOMB<0 */

double Ewald(int pass,vector *frp,vector *rp) /*********************** Ewald */
{
  static vector oldL={-3e33,-3e33,-3e33};
  static vector firstL;

  static int Nkk= -1; /* # of k-vectors irresp. of signs of components */
                      /* -1 forces initialization */

  real *ptr;
#ifdef POLAR
  vector rplusrpol;
  real *polptr;
#endif /*# POLAR */

  vector ek;
#ifdef GAUSSIANCHARGES
  struct Qk_s *Qk;
  int isites;
#else /*# GAUSSIANCHARGES */
  complex *Qk;
#endif /*#!GAUSSIANCHARGES */
  real e,c,kk,aux;
#if COULOMB<0
  real kxq,kyq,dK;
  static real Kk; /* K refers to z-coordinate integer vector */
  int nagain=0;
#endif /*# COULOMB<0 */
  vector PIL;
  double Energy=0,alsPI=el.alpha/sqrt(PI);
  int i,kx,ky,kz,iQ,ikk,kymax,kzmax;
  long size,qtabsize;
#if PARALLEL==1 || PARALLEL==3
  int ith,ipt;
#else /*# PARALLEL==1 || PARALLEL==3 */
  real FA,FB,FC,FD;
  complex *E,A,B,C,D;
#  ifdef GAUSSIANCHARGES
  complex AT,BT,CT,DT;
  complex AS,BS,CS,DS;
#  endif /*# GAUSSIANCHARGES */
  int jkk,iqt;
#endif /*#!PARALLEL==1 || PARALLEL==3 */

  molecule_t *mn;

  int iq,n,ns,sp,nmol;
  siteinfo_t *si;
  vector *r,*fr;
  static int oldkappa=-99;

  if (el.kappa<=0) return 0;

  VV(PIL,=2*PI/box.L)
  /* WARNING: the initialization logic is weird... */
  if (pass<=0) {
    /*
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %                        pass 0 - initialization                       %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    */

    if (Nkk==-1 || fabs(el.kappa-oldkappa)>1e-11) {

#if COULOMB<0
      /* Ewald only: Nq,Nkk,charges calculated */
      oldkappa=el.kappa;
      if (machineprec==0) {
        if (sizeof(machineprec)<8) machineprec=1e-6;
        else machineprec=1e-9; }
      dK=machineprec*2;

     AGAIN: if (nagain++>25) ERROR(("Ewald: cannot solve rounding problems"))

      if (pass>=0) prt("initializing Ewald: kappa=K/L=%.9g",el.kappa);
      KX=box.L[0]*el.kappa;
      Kk=Sqr(el.kappa);

      /*** if not the 1st initialization then free the previous arrays ***/
      if (Nkk != -1) {
        if (ekk) free(ekk);
#  ifdef GAUSSIANCHARGES
        if (sekk) free2Darray(sekk);
        if (ssekk) free2Darray(ssekk);
#  endif /*# GAUSSIANCHARGES */
        if (Q) free(Q);
        if (ktab) {
          loop (n,0,ktabsize) free(ktab[n].kzmax);
          free(ktab); } }
#endif /*# COULOMB<0 */

      /*** computing Nq and (sums of) charges ***/
      Nq=0;
#ifdef GAUSSIANCHARGES
      charges.sumaq=
#endif /*# GAUSSIANCHARGES */
#ifndef POLAR
        charges.molsumq=charges.molsumdq=
#endif /*# POLAR */
        charges.sum=charges.sumq=charges.max=charges.summq=0;

      /* this part was optimized in V3.3b:
         - sped up
         - cumulation of errors (esp. for almost point Drude dipoles) avoided */
      for (n=0;n<No.N;) {
        mn=molec+n;
        si=spec[mn->sp]->si;
        ns=mn->ns;

        sp=mn->sp;
        nmol=n;
        while (n<No.N && sp==molec[n].sp) n++; /* now n=next species */
        nmol=n-nmol; /* # of molecules of the same kind */
        loop (i,0,ns) {
#ifndef POLAR
          /* nonpolar */
          if ((e=si[i].charge)!=0) {
            Max(charges.max,fabs(e))
            charges.sum+=nmol*e;
            charges.sumq+=nmol*(e*e);
#  ifdef GAUSSIANCHARGES
            charges.sumaq+=nmol*(e*e)/sqrt(1/Sqr(el.alpha)+2*Sqr(sitedef[si[i].st].LJ[0].parm[SS_PARMS-1]));
#  endif /*# GAUSSIANCHARGES */
            if (si[i].mass) charges.summq+=nmol*Sqr(e*si[i].imass);
            Nq+=nmol; }
#elif POLAR&32
          /* FQ + Drude (Gaussian not supported) */
#  if POLAR&4
          if (si[i].qtype)
#  endif /*# POLAR&4 */
          if (si[i].qtype&FQ) {
            e=si[i].charge;
            Max(charges.max,fabs(e))
            charges.sum+=nmol*e;
            charges.sumq+=nmol*(e*e);
            if (si[i].mass) charges.summq+=nmol*Sqr(e*si[i].imass);
            if (si[i].charge) Nq+=nmol; }
          else {
            e=si[i].charge+si[i].chargepol;
            Max(charges.max,fabs(e))
            charges.sum+=nmol*e;
            charges.sumq+=nmol*(e*e);
            if (si[i].mass) charges.summq+=nmol*Sqr(e*si[i].imass);
            if (si[i].charge   ) Nq+=nmol;
            if (si[i].chargepol) Nq+=nmol; }
#else  /*#!POLAR!POLAR&32 */
          /* normal+Drude charges */
#  if POLAR&4
          if (si[i].qtype)
#  endif /*# POLAR&4 */
          {
            e=si[i].charge+si[i].chargepol;
            Max(charges.max,fabs(e))
            charges.sum+=nmol*e;
            charges.sumq+=nmol*(e*e);
#  ifdef GAUSSIANCHARGES
            /* gaussian sigma = sitedef[si[i].st].LJ[0].parm[SS_PARMS-1] */
            charges.sumaq+=nmol*(e*e)/sqrt(1/Sqr(el.alpha)+2*Sqr(sitedef[si[i].st].LJ[0].parm[SS_PARMS-1]));
#  endif /*# GAUSSIANCHARGES */
            if (si[i].mass) charges.summq+=nmol*Sqr(e*si[i].imass);
            if (si[i].charge   ) Nq+=nmol;
            if (si[i].chargepol) Nq+=nmol; }
#endif /*#!POLAR!POLAR&32 */
        }
      } /* n */

#ifdef ECC
      /* ECC support - EXPERIMENTAL, inefficient code */
      if (el.ecc && rp!=NULL) loop (n,0,No.N) {
        mn=molec+n;
        si=spec[mn->sp]->si;
        ns=mn->ns;

        double qsum=0;
        vector rqsum={0,0,0};
        r=rof(mn,rp);

        loop (i,0,ns) {
          if ((e=si[i].charge)!=0) {
            qsum+=e; /* for ECC support: interferes with SF, so the test */
            VV(rqsum,+=e*r[i]) } }

        if (fabs(qsum)<1e-7) charges.molsumdq+=SQR(rqsum);
        else charges.molsumq+=Sqr(qsum);
      } /* n */
#endif /*# ECC */

#if COULOMB<0
      /*** computing diagonal correction */
      Ediag=0;
      if (el.diag>=0) {
        double qkappax=Sqr(el.kappa)+Sqr(el.alpha);
        /* then relative accuracy of Ediag will be ~ exp(-PI^2) */
        double kappax=sqrt(qkappax);
        double kk;
        double q= -Sqr(PI/el.alpha);
        int Kx=kappax*box.L[0];
        int Ky=kappax*box.L[1];
        int Kz=kappax*box.L[2];

        loopto (kx,-Kx,Kx) loopto (ky,-Ky,Ky) loopto (kz,-Kz,Kz) {
          kk=Sqr(kx/box.L[0])+Sqr(ky/box.L[1])+Sqr(kz/box.L[2]);
          if (kk>Sqr(el.kappa) && kk<qkappax) Ediag+=exp(kk*q)/kk; }
        Ediag *= charges.sumq/(2*PI*box.L[0]*box.L[1]*box.L[2]); }

      /*** computing sizes of Q and ekk;
           readjusting K for rounding problems:
             because of optimization, an expression like
             L[2]*sqrt(kyq) may give once e.g. 2+1e-15 and in other part
             of the code 2-1e-15; then, (int) of such a number is not
             well defined: K is a bit increased if this might happen;
           (computing ekk moved to pass==1 with check of L change) ***/
      iQ=ikk=0;
      kxmax=(int)KX;
      loopto (kx,0,kxmax) {
        kxq=Kk-Sqr(kx/box.L[0]);
        aux=box.L[1]*sqrt(kxq);
        if (pass==0 && almostint(aux)) {
          el.kappa+=el.kappa*dK; dK*=1.4;
          prt("Ewald: possible rounding problem (x) => kappa increased to %.12g",el.kappa);
          goto AGAIN; }
        kymax=(int)aux;
        loopto (ky,0,kymax) {
          kyq=kxq-Sqr(ky/box.L[1]);
          aux=box.L[2]*sqrt(kyq);
          if (pass==0 && almostint(aux)) {
            el.kappa+=el.kappa*dK; dK*=1.4;
            prt("Ewald: possible rounding problem (y) => kappa increased to %.12g",el.kappa);
            goto AGAIN; }
          kzmax=(int)aux;
          //          put3(kx,ky,kzmax)
          loopto (kz,!(kx+ky),kzmax) {
            iQ += 1 << ((kx>0)+(ky>0)+(kz>0)-1); /* = 2^(..) */
            ikk++; } } } /* kz,ky,kx */

      /*** creating ktab (repeats the previous code to some extent) ***/
      ktabsize=kxmax+1;
      allocarray(ktab,ktabsize);
      loopto (kx,0,kxmax) {
        kxq=Kk-Sqr(kx/box.L[0]);
        kymax=(int)(box.L[1]*sqrt(kxq));
        ktab[kx].kymax=kymax;
        allocarray(ktab[kx].kzmax,kymax+1);
        loopto (ky,0,kymax) {
          kyq=kxq-Sqr(ky/box.L[1]);
          ktab[kx].kzmax[ky]=(int)(box.L[2]*sqrt(kyq)); } }

      VO(oldL,=-3e33) /* to pretend the old L is different and calculate ekk later */
      if (!pass) { /*** initialization protocol printed to out ***/
#  include "ewaldprt.c"
      } /* !pass */

      sdsalloc(Q,sizeof(Q_t)+(iQ-1)*sizeof(Q->Qk[0]));
      Q->N=iQ;
      VO(Q->M,=0) /* to have Q->M defined even if no charges
                     ... though it is stupid to allocate Q at all */
      Nkk=ikk;
      /* during pre-inititalization, Nkk may be zero, so prevent alloc error */
      allocarray(ekk,Nkk+!Nkk);
#  ifdef GAUSSIANCHARGES
      alloc2Darray(sekk,Nkk+!Nkk,nsites);
      alloc2Darray(ssekk,Nkk+!Nkk,nsites);
#  endif /*# GAUSSIANCHARGES */
#endif /*# COULOMB<0 */
    }

#if COULOMB<0
    return (double)Q->N/No.N;
#else /*# COULOMB<0 */
    return 0;
#endif /*#!COULOMB<0 */
  } /* initialization (pass<=0) */

  /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             passes 1 and 2                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */

  /*** computing ekk (moved here from initialization to enable L change) ***/
  /* selected k-vectors, table sizes, alpha, ... are unchanged! */
//  if (L[0]!=oldL[0] || L[1]!=oldL[1] || L[2]!=oldL[2]) {
  if (memcmp(box.L,oldL,sizeof(box.L))) {
    if (firstL[0]!=0)
      if (fabs(log(box.L[0]/firstL[0]))>el.diff
       || fabs(log(box.L[1]/firstL[1]))>el.diff
       || fabs(log(box.L[2]/firstL[2]))>el.diff) {
        if (charges.sumq) WARNING(("variable L with Ewald:\n*** box has changed by >%g%% (log-scale) and k-vectors remain unchanged\n*** L=[%g %g %g]",el.diff*100,VARG(box.L)))
        else prt("WARNING: box has changed by >%g%% (but Ewald inactive because zero charge)",el.diff*100);
        el.diff*=2; }

    if (oldL[0]>-1e33 && firstL[0]==0) {
      prt("WARNING: variable box L with Ewald:\n(only small changes in L are allowed without restarting the simulation)");
      VV(firstL,=oldL) }

    PIaq=-Sqr(PI/el.alpha);
    iQ=ikk=0;
    if (kxmax+1!=ktabsize)
      ERROR(("ktab size error: requested=%d  allocated with=%d  KX=%g",kxmax+1,ktabsize,KX))
    loopto (kx,0,kxmax) {
      kymax=ktab[kx].kymax;
      loopto (ky,0,kymax) {
        aux=Sqr(ky/box.L[1])+Sqr(kx/box.L[0]);
        kzmax=ktab[kx].kzmax[ky];
        loopto (kz,!(kx+ky),kzmax) {
          kk=Sqr(kz/box.L[2])+aux;
          ekk[ikk] = 4/(kk*(box.L[0]*box.L[1]*box.L[2]))*exp(PIaq*kk);
          /* NB: factor 4 above comes from: 2 from differentiation of Q^2
                                            2 from counting k and -k at once */
#ifdef GAUSSIANCHARGES
          loop (isites,0,nsites) {
            double sigmaq=Sqr(sitedef[isites].LJ[0].parm[SS_PARMS-1]); // optimization - prepare table?!

            sekk[ikk][isites] = exp(-Sqr(PI)*2*sigmaq*kk);
            ssekk[ikk][isites] = sekk[ikk][isites]*sigmaq; }
#endif /*# GAUSSIANCHARGES */
          iQ += 1 << ((kx>0)+(ky>0)+(kz>0)-1); /* = 2^(..) */
          ikk++; } } } /* kz,ky,kx */
    if (iQ!=Q->N || ikk!=Nkk) ERROR(("internal: Q->N=%d (%d expected)  Nkk=%d (%d expected)\n\
*** (may be caused by K rounding problems - try to change K a bit)",iQ,Q->N, ikk,Nkk))

    VV(oldL,=box.L) }

  /*** preparing tables of multiplying factors ex=exp[rx*..], ey, ez, etc.
       and the total dipole moment ***/

  if (Nq==0) return 0; /* nothing to do */

  qtabsize=Nq*sizeof(qtab_t);

  if (pass==1) {
    VO(Q->M,=0)
#ifdef POLAR
    charges.fqSumq=0;
#endif /*# POLAR */
    size=Nq*sizeof(complex);
    /* ex,ey,ez will be freed in pass 2 */
    alloc(qtab,qtabsize);
    alloc(ex,size); alloc(ey,size); alloc(ez,size); }
  else {
#if defined(POLAR) && POLAR&32
    allocarrayzero(Phi,Nq);
#endif /*# defined(POLAR) && POLAR&32 */
    allocarrayzero(F,Nq); }

  /*** configuration -> Ewald interface ***/

  iq=0;
  if (pass==1)
    loop (n,0,No.N) {
      mn=molec+n;
      si=spec[mn->sp]->si;
      ns=mn->ns;
      r=rof(mn,rp);
      loop (i,0,ns) {
#if defined(POLAR) && POLAR&32
        if (si[i].qtype&FQ) {
          aux=polarrof(mn,rp)[i][0]; /* fluctuating charge */
          charges.fqSumq+=Sqr(aux);
          ptr=r[i];
          VV(Q->M,+=aux*ptr)
          qtab[iq].x.re=aux; qtab[iq].x.im=0;
          ex[iq].re=cos(aux=PIL[0]*ptr[0]); ex[iq].im=sin(aux);
          ey[iq].re=cos(aux=PIL[1]*ptr[1]); ey[iq].im=sin(aux);
          ez[iq].re=cos(aux=PIL[2]*ptr[2]); ez[iq].im=sin(aux);
#  ifdef GAUSSIANCHARGES
#    error "not implemented"
#  endif /*# GAUSSIANCHARGES */
          iq++; }
        else
#endif /*# defined(POLAR) && POLAR&32 */
        if ((aux=si[i].charge)) {
          ptr=r[i];
          VV(Q->M,+=aux*ptr)
          qtab[iq].x.re=aux; qtab[iq].x.im=0;
          ex[iq].re=cos(aux=PIL[0]*ptr[0]); ex[iq].im=sin(aux);
          ey[iq].re=cos(aux=PIL[1]*ptr[1]); ey[iq].im=sin(aux);
          ez[iq].re=cos(aux=PIL[2]*ptr[2]); ez[iq].im=sin(aux);
#ifdef GAUSSIANCHARGES
          qtab[iq].site=si[i].st;
          //        qtab[iq].sigmaq=Sqr(sitedef[si[i].st].LJ[0].parm[SS_PARMS-1]); // needed?
#endif /*# GAUSSIANCHARGES */
          iq++; }
#ifdef POLAR
        if ((aux=si[i].chargepol)) {
          ptr=r[i]; /* for sure if chargepol but not charge */
          polptr=(real*)((char*)ptr+polar_off);
          VVV(rplusrpol,=ptr,+polptr)
          VV(Q->M,+=aux*rplusrpol)
          qtab[iq].x.re=aux; qtab[iq].x.im=0;
          ex[iq].re=cos(aux=PIL[0]*rplusrpol[0]); ex[iq].im=sin(aux);
          ey[iq].re=cos(aux=PIL[1]*rplusrpol[1]); ey[iq].im=sin(aux);
          ez[iq].re=cos(aux=PIL[2]*rplusrpol[2]); ez[iq].im=sin(aux);
#  ifdef GAUSSIANCHARGES
          qtab[iq].site=si[i].st;
          //        qtab[iq].sigmaq=Sqr(sitedef[si[i].st].LJ[0].parm[SS_PARMS-1]); // needed?
#  endif /*# GAUSSIANCHARGES */
          iq++; }
#endif /*# POLAR */
      } /* i */
      //      prt("%d %g %g %g QQQ", n,VARG(Q->M));
    } /* n */
  else /* pass 2 */
    loop (n,0,No.N) {
      mn=molec+n;
      si=spec[mn->sp]->si;
      ns=mn->ns;
      r=rof(mn,rp);
      loop (i,0,ns) {
#if defined(POLAR) && POLAR&32
        if (si[i].qtype&FQ) {
          /* fluctuating charge */
          qtab[iq].q=polarrof(mn,rp)[i][0];
          qtab[iq].x.re=1; qtab[iq].x.im=0; iq++; }
        else if ((aux=si[i].charge)) {
          /* normal charge (Drude may follow) */
          qtab[iq].q=aux;
          qtab[iq].x.re=1; qtab[iq].x.im=0; iq++; }
        if ((aux=si[i].chargepol)) {
          /* Drude */
          qtab[iq].q=aux;
          qtab[iq].x.re=1; qtab[iq].x.im=0; iq++; }
#else /*# defined(POLAR) && POLAR&32 */
        if ((aux=si[i].charge)) {
          /* normal charge (Drude may follow) */
          qtab[iq].x.re=aux; qtab[iq].x.im=0; iq++; }
#  ifdef POLAR
        if ((aux=si[i].chargepol)) {
          /* Drude */
          qtab[iq].x.re=aux; qtab[iq].x.im=0; iq++; }
#  endif /*# POLAR */
#endif /*#!defined(POLAR) && POLAR&32 */
      }
    }
  if (iq!=Nq) Error("Ewald:Nq");

/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %             main loop over k-vectors  (passes 1 and 2)              %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

  Qk=Q->Qk;
  iQ=ikk=0;
  kxmax=(int)KX;

#if PARALLEL==1 || PARALLEL==3

  if (pass==1) {
    if (sizeof(int)!=4) WARNING(("sizeof(int)=%d, but 4 is assumed - it will go wrong..",sizeof(int)))
    if (sizeof(struct pll_ewald_s)%CACHELINE) DISASTER(("sizeof(pll_ewald_s)=%d is not an integer multiple of the cache line=%d\n\
*** adjust padding in `struct pll_ewald_s' in sim/ewald.h and recompile!",sizeof(struct pll_ewald_s),CACHELINE))

    /* pll_ewald no longer allocated in forces.c */
    if (!pll_ewald) allocarray(pll_ewald,No.th); /* forever (No.th should not change) */

    arrayzero(pll_ewald,No.th);
    alloc2Darrayzero(Qkcell,No.th,Q->N); }

  par1_ew.sfr=el.sfr;
  par1_ew.sf3d=el.sf3d;
  par1_ew.pass=pass;

  PARALLELIZE(parEwald,No.th)

  if (pass==1) {
    //    loop (i,0,Q->N) { Qk[i].re=0; Qk[i].im=0; }

    arrayzero(Qk,Q->N); /* also for GAUSSIANCHARGES */
    loop (ith,0,No.th) {
      int i;

      loop (i,0,Q->N) {
#  ifdef GAUSSIANCHARGES
        Qk[i].Q.re+=Qkcell[ith][i].Q.re;
        Qk[i].Q.im+=Qkcell[ith][i].Q.im;
        Qk[i].T.re+=Qkcell[ith][i].T.re;
        Qk[i].T.im+=Qkcell[ith][i].T.im;
        Qk[i].S.re+=Qkcell[ith][i].S.re;
        Qk[i].S.im+=Qkcell[ith][i].S.im;
#  else /*# GAUSSIANCHARGES */
        Qk[i].re+=Qkcell[ith][i].re;
        Qk[i].im+=Qkcell[ith][i].im;
#  endif /*#!GAUSSIANCHARGES */
      } }
    /* warning: unnecessarilly(?) all Qkcell in pass 2 */
    loop (ith,0,No.th) {
      int i;

      loop (i,0,Q->N) Qkcell[ith][i]=Qk[i]; } }
  else {

    Energy=pll_ewald[0].E/(4*PI);
#  if PRESSURETENSOR&PT_VIR
    VO(En.Pvir,+=Energy)
    loop (kx,0,PT_DIM) En.Pvir[kx]-=2/(4*PI)*pll_ewald[0].Pvir[kx];
#  endif /*# PRESSURETENSOR&PT_VIR */
    /* calculated (unnecessarily) in all threads */

    loop (ith,0,No.th) {
      if (fabs(pll_ewald[0].E-pll_ewald[ith].E)>1e-6)
        ERROR(("different energies in threads (0 vs. %d): %g %g",
               ith,pll_ewald[0].E,pll_ewald[ith].E))
#  if PRESSURETENSOR
      loop (ipt,0,PT_DIM)
        if (fabs(pll_ewald[0].Pvir[ipt]-pll_ewald[ith].Pvir[ipt])>1e-6)
          ERROR(("different Pvir[%d] in threads 0 and %d",ipt,ith))
#  endif /*# PRESSURETENSOR */
    }
    free2Darray(Qkcell);
  }

  if (partimes.kspace) loop (ith,0,No.th)
    partimes.kspace[ith]+=pll_ewald[ith].t/(double)CLOCKS_PER_SEC;

#elif PARALLEL==2
  /* whole Ewald is a separate thread (parallel to pairforces)
     => must sum into separate memory and sum up later */
#  define Energy pll_ewald->E
#  define PVIR pll_ewald->Pvir
  /* ..pll_ewald->E,pll_ewald->Pvir summed after par2_func() has finished */

  /* to be defined even if not measured */
  pll_ewald->E=pll_ewald->Ereturned=0;

  if (pass==1) {
#  include "ewpass1.c"
    }
  if (pass==2)
    if (measure) {
#  if PRESSURETENSOR&PT_VIR
      double ptf,E1,E1x;

      Energy=0;
      memset(PVIR,0,sizeof(PVIR));
#  endif /*# PRESSURETENSOR&PT_VIR */
#  include "ewpass2m.c"
      Energy/=4*PI; }
    else {
#  include "ewpass2.c"
    }
//?? #  undef PVIR -- must be kept shadowed because of M-contrib
//?#  undef Energy  -- can keep using Energy==pll_ewald->E
//?  Energy=pll_ewald->E;

#else /*#!PARALLEL==1 || PARALLEL==3!PARALLEL==2 */
  /* serial version */

  if (pass==1) {
#  include "ewpass1.c"
  }
  if (pass==2) {
    if (measure) {
#  if PRESSURETENSOR&PT_VIR
      double PVIR[PT_DIM];
      double ptf,E1,E1x;

      memset(PVIR,0,sizeof(PVIR));
#  endif /*# PRESSURETENSOR&PT_VIR */
#  include "ewpass2m.c"
      Energy/=4*PI;
#  if PRESSURETENSOR&PT_VIR
      VO(En.Pvir,+=Energy)
      loop (kx,0,PT_DIM) En.Pvir[kx]-=2/(4*PI)*PVIR[kx];
#  endif /*# PRESSURETENSOR&PT_VIR */
    } else {
#  include "ewpass2.c"
    } }

  /* end of serial version */
  /* because the code below calculating the M^2-term constributions */
#  define PVIR En.Pvir
#endif /*#!PARALLEL==1 || PARALLEL==3!PARALLEL==2 */

  if (pass==2) {
    double Edip;
    int k;

    if (el.corr) {
      /*
         Extended "slab" correction, see doi:10.1063/1.479595
         [In-Chul Yeh and Max L. Berkowitz, JCP 111, 3155 (1999)].
         See below for addditional terms with charged background (el.bg).
         Use el.corr=4 for slab pseudo-2D b.c.
      */

      Edip=0; /* dipolar energy */
      VO(ek,=0) /* -> dipolar force */
      c = 4*PI/PROD(box.L);
      loop (k,0,3) if (el.corr & 1<<k) {
        if ((el.corr&8)==0) ek[k] = -c*Q->M[k]; /* standard Yeh-Berkowitz */
        Edip+=c/2*Sqr(Q->M[k]); } }
    else {
      /* isotropic dielectric b.c. */
      c = 4*PI/(2*el.epsinf+1)/PROD(box.L);
      VV(ek,=-c*Q->M)
      Edip = c/2*SQR(Q->M); }

#if PRESSURETENSOR&PT_VIR
    VO(PVIR,+=Edip)
#endif /*# PRESSURETENSOR&PT_VIR */
    iq=0;

    /* NB:
       - bug fixed in V2.7h: diag. of PTENS was wrong (but tr(PTENS) was OK)
       - if PARALLEL==2, PVIR is pll_ewald->Pvir, otherwise En.Pvir
       - in the following macros, aux will be a charge */

    if (!el.sf) {

      if (el.bg && el.corr) {
        /*
           Yeh-Berkowitz correction for CHARGED SYSTEMS compensated by
           electric background, see doi:10.1063/1.3216473, eq. (30)
           [V. Ballenegger, A. Arnold and J. J. Cerda, JCP 131, 094107 (2009)]
        */
        int k;
        vector ekx;

        c = 4*PI*charges.sum/PROD(box.L);

#if PRESSURETENSOR&PT_VIR
#  if PRESSURETENSOR&PT_OFF
#    define ADDFORCES(R) { \
       loop (k,0,3) { ekx[k]=ek[k]; if (el.corr & 1<<k) ekx[k]+=c*R[k]; } \
       VVV(ptr,+=F[iq],+aux*ekx) \
       if (measure) { \
         VVV(PVIR,+=aux*ekx,*R) \
         PVIR[3]+=aux*ekx[1]*R[2]; \
         PVIR[4]+=aux*ekx[2]*R[0]; \
         PVIR[5]+=aux*ekx[0]*R[1]; } }
#  else /*# PRESSURETENSOR&PT_OFF */
#    define ADDFORCES(R) { \
       loop (k,0,3) { ekx[k]=ek[k]; if (el.corr & 1<<k) ekx[k]+=c*R[k]; } \
       VVV(ptr,+=F[iq],+aux*ekx) \
         if (measure) { \
           VVV(PVIR,+=aux*ekx,*R) } }
#  endif /*#!PRESSURETENSOR&PT_OFF */
#else /*# PRESSURETENSOR&PT_VIR */
#  define ADDFORCES(R) { \
       loop (k,0,3) { ekx[k]=ek[k]; if (el.corr & 1<<k) ekx[k]+=c*R[k]; } \
       VVV(ptr,+=F[iq],+aux*ekx) }
#endif /*#!PRESSURETENSOR&PT_VIR */

#include "ewaldf.c"
#undef ADDFORCES

      } else {
        /* standard version for neutral systems */

#if PRESSURETENSOR&PT_VIR
#  if PRESSURETENSOR&PT_OFF
#    define ADDFORCES(R) \
       VVV(ptr,+=F[iq],+aux*ek) \
       if (measure) { \
         VVV(PVIR,+=aux*ek,*R) \
         PVIR[3]+=aux*ek[1]*R[2]; \
         PVIR[4]+=aux*ek[2]*R[0]; \
         PVIR[5]+=aux*ek[0]*R[1]; }
#  else /*# PRESSURETENSOR&PT_OFF */
#    define ADDFORCES(R) \
       VVV(ptr,+=F[iq],+aux*ek) \
         if (measure) { \
           VVV(PVIR,+=aux*ek,*R) }
#  endif /*#!PRESSURETENSOR&PT_OFF */
#else /*# PRESSURETENSOR&PT_VIR */
#  define ADDFORCES(R) \
       VVV(ptr,+=F[iq],+aux*ek)
#endif /*#!PRESSURETENSOR&PT_VIR */

#include "ewaldf.c"
      } }

#if defined(POLAR) && POLAR&32
    free(Phi);
#endif /*# defined(POLAR) && POLAR&32 */
    free(F); free(qtab);
    free(ez); free(ey); free(ex);

    if (el.kappa<=0) return 0; /* patch !!! */

    if (measure) {
#ifndef POLAR
      double sumq=charges.sumq;
#elif POLAR&32
      double sumq=charges.sumq+charges.fqSumq;
#else /*#!POLAR!POLAR&32 */
      double sumq=charges.sumq; /* using exacterud_sqrt_1(): Sumq removed */
#endif /*#!POLAR!POLAR&32 */

      return
#if PARALLEL==2
        /* cannot return value in parallel */
        pll_ewald->Ereturned=
#endif /*# PARALLEL==2 */
        Ediag*(el.diag==1) // diagonal correction
        + Energy
        + charges.bgE/PROD(box.L) // background correction term
#ifdef GAUSSIANCHARGES
        - charges.sumaq/sqrt(PI)
#else /*# GAUSSIANCHARGES */
        - alsPI*sumq /* NB: this term is HUGE for POLAR/Drude, therefore
                        if r-space (i.e., pairforces()) are called after
                        k-space (Ewald()), the contributions should not be
                        cumulated, but En.el should be reset to zero */
#endif /*#!GAUSSIANCHARGES */
        + Edip; }
    else
      return 0; } /* pass==2 */

  return (double)Q->N/No.N;

} /* Ewald */

#undef CxE
#undef CxiE
#undef Cadd
#undef Csqr
#undef iqeloop
#undef iqfloop
#undef iqloop

/*
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                  Ewald tests and parameter setting                    %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef LINKCELL
static double rkspeed=0; /* speed factor for optimization of the linkcell method
  CPU time = const + A*(K+2/3)^3 + B*cutoff^3
  rkspeed = B*L^3/Nq/A   (is dimensionless, typical value around 3 for linkcell)
  CPU time ~ const + (K+2/3)^3 + rkspeed*Nq*(cutoff/L)^3
  if rkspeed<>0 and el.test<0 then cutoff is determined to minimize time */
#endif /*# LINKCELL */

static double Ewaldparm(double cutoff,vector L) /***************** Ewaldparm */
{
  /* calculates alpha, kappa from cutoff, epsr, epsk */
  double x,kappa=el.kappa; /* ??? (for the case if K is single) */
  double V=PROD(L);
  int nit=10000;
/*.....  double minL=fmin(fmin(L[1],L[0]),L[2]);*/

  if (charges.sumq<=0) {
    /* just calculate charges CUMBERSOME, FIX!!!! */
    /* el.kappa=???: problem with CUTELST - alloc failed; is this needed at all? */
    el.kappa=0.0001;
    Ewald(-1,NULL,cfg[0]->rp);
    el.kappa=kappa; }

  if (charges.sumq>0) {
    /* alpha, random charges */
    el.alpha=log(2*charges.max*sqrt(charges.sumq/cutoff/V)/el.epsr);
    if (el.alpha<0) {
      ERROR(("Ewald auto set: cannot determine alpha (check epsr)\n\
*** (charges.max=%g charges.sumq=%g cutoff=%g V=%g el.epsr=%g)",charges.max,charges.sumq,cutoff,V,el.epsr));
      return -1; }
    el.alpha=sqrt(el.alpha)/cutoff;
    /* kappa, random charges */
    do {
      if (!nit--) {
        ERROR(("Ewald auto set: cannot determine kappa (check el.epsk)\n\
*** (charges.max=%g el.alpha=%g el.epsk=%g)",charges.max,el.alpha,el.epsk))
        return -1; }
      x=kappa;
      kappa=log(charges.max*el.alpha/(PI*el.epsk)*sqrt(8*charges.sumq/(kappa*V)));
      if (kappa<=0) ERROR(("Ewald auto set: cannot determine kappa (check epsk)"))
      kappa=el.alpha*sqrt(kappa)/PI;
    } while (fabs(x-kappa)>1e-9);
    el.kappa=kappa; }
#ifdef LINKCELL
  return V*Cub(kappa)+rkspeed*Nq*Cub(cutoff)/V;
#else /*# LINKCELL */
  return 0;
#endif /*#!LINKCELL */
}

// #ifndef NIBC ??? - charges probably needed anyway
int Ewaldtest(double *setL) /************************************* Ewaldtest */
/***
  Test module of electrostatics errors,
  calculate sums of charges (not cutoff electrostatics)
  Ewald: set parameters acording to el.test=
    0: quit Ewaldtest ( erfc and Ewald initialized only )
    1: calculate forces, energy, and errors (if reference is set)
    2: calculate forces and energy and store them as the reference for
       evaluating errors in next steps, sets el.test=1
   -1: calculate alpha and kappa from cutoff, epsr, epsk, then el.test=1
   -2: calculate alpha and kappa from cutoff, epsr, epsk, then el.test=2
   -3: calculate alpha and kappa from cutoff, epsr, epsk
  -10: calculate alpha and kappa from cutoff, epsr, epsk silently, then
       quit Ewaldtest (to be used for automatic alpha,kappa setting)

  1 is returned if the program is to be quitted
  NOTES:
  - Forces and energy include also all non-electrostatic interactions and
    these are used also to evaluate relative errors.  Using smaller (or big
    integer) eps is recommended, the old value is automatically restored on
    return (i.e., el.test=0).
  - Ewald is unnecessarily re-initialized even if the parameters do not change
  - if init=20 has been specified, forces read by readasc() are directly
    available as a reference (OBSOLETE?)
***/
{
  static double Epotref,Ediagref,dfref=1,maxfref,Evirref;
#if PRESSURETENSOR&PT_VIR
  static double Pvirref[PT_DIM];
#endif /*# PRESSURETENSOR&PT_VIR */
  double ff,df=1,maxf;
  double oldalpha=el.alpha,oldkappa=el.kappa,oldeps=eps,varf=1;
  double cutoff=box.cutoff;
  vector *d,*dref,dr;
  int ns=No.s,quit=0,etestref=0;

  int ii,no=1; /* just to do calculations several times to measure time */
  int i=-1,lasti=-1;
  vector ri={0,0,0};
  ToIntPtr A;
  time_t time0,time1;

  measure=1;

  df=0;
  loop (ii,0,DIM)
    Max(df,fabs(log(setL[ii]/box.L[ii])))
  if (df>1e-9)
    prt("WARNING Ewaldtest: reference box differs from the actual box by max %.2f%%\n\
(the reference box will be used for Ewald auto set if requested)",df*100);

  prt("Ewaldtest: el.test=%d, reference L=[%g %g %g]",el.test,setL[0],setL[1],setL[2]);

  if (init==12) {
    A=cfg[gear.order+2];
    sdscopy(A,cfg[2]);
    goto setref; }

 again:

  if (!(el.test==0 || el.test==-10)) {
#if COULOMB>=0
    WARNING(("el.test==0 || el.test==-10 only partly supported for COULOMB=%d",COULOMB))
#endif /*# COULOMB>=0 */

    _n
    putline('T',79);
    prt("el.test=%d: 0=end  1=calc  2=set ref  -3=set parm  -1,-2=parm+1,2  -10=parm+end",el.test);

    cutoff=box.cutoff;

    if (i>=0 && i<ns) {
      vector fi;
      VV(ri,=cfg[0]->rp[i])
      VV(fi,=cfg[gear.order+abs(el.test)]->rp[i])
      put(i)
      putv(ri)
      putv(fi) }

    getdata
      get(el.test)
      get(el.epsk) get(el.epsr) get(eps)
      get(el.grid) get(el.minqq)
      get(el.alpha) get(el.kappa) get(cutoff)
      get(quit)
#ifdef LINKCELL
      getvec(No.cell,,DIM) get(rkspeed)
#endif /*# LINKCELL */
#ifdef POLAR
      get(scf.eps) get(scf.maxit) get(scf.omega)
#endif /*# POLAR */
      get(no) get(varf)
      get(i) if (i!=lasti && i>=0 && i<ns) {
        vector fi;
        VV(ri,=cfg[0]->rp[i])
        VV(fi,=cfg[gear.order+abs(el.test)]->rp[i])
        putv(ri)
        putv(fi)
        lasti=i; }
      getvec(ri,,3)
      checkdata
    enddata
    box.cutoff=cutoff;

#ifdef POLAR
    if (scf.omega<0) {
      scf.omega=0.95;
      Max(scf.maxit,100)
      Min(scf.eps,1e-9)
      WARNING(("polar iterations not set, using scf.omega=%g scf.maxit=%d scf.eps=%g",
        scf.omega,scf.maxit,scf.eps)) }
#endif /*# POLAR */

    if (el.alpha<=0) el.alpha=oldalpha;
    if (el.kappa<=0) el.kappa=oldkappa;

    if (i>=0 && i<ns) VV(cfg[0]->rp[i],=ri)

    if (quit) {
#if PARALLEL
      printpartimes();
#endif /*# PARALLEL */
      return 1; }
  }

  if (el.test<0) {
#ifdef LINKCELL
    if (rkspeed) {
      double dc=box.cutoff*0.01,t,tt;
      int ic;

      if (dc<=0) dc=setL[2]*0.005;
      if (dc>0.085 && dc<0.7) dc=1;    /* typical if lengths in AA */
      else if (dc>0.03) dc=0.1;      /* typical if reduced units */

      ic=(int)(setL[2]*0.6/dc+1.5); /* max cutoff=0.6*L */
      t=3e33;
      do {
        ic--;
        tt=t;
        if ((t=Ewaldparm(ic*dc,setL))<0) goto again;
      } while (t<tt);
      box.cutoff=(ic+1)*dc; }
#endif /*# LINKCELL */

    if (Ewaldparm(cutoff,setL)<0) goto again;


#if COULOMB<0
    prt("\
Ewald auto set: cutoff=%g alpha=%g kappa=%g\n\
                Kx=%.3f Ky=%.3f Kz=%.3f",
        cutoff,el.alpha,el.kappa,
        VARG(el.kappa*setL));
#else /*# COULOMB<0 */
    prt("Short electrostatics: cutoff=%g alpha=%g (kappa=%g ignored)\n",
        cutoff,el.alpha,el.kappa);
#endif /*#!COULOMB<0 */
    if (el.test!=-10 && el.test<-2) goto again; }

  box.cq=Sqr(cutoff);
  initrspace();
  Ewald(0,NULL,cfg[0]->rp);

  if (el.test==0 || el.test==-10) {
    eps=oldeps;
    return 0; }

  if (abs(el.test)>2) {
    Error("el.test out of range, replaced by 1");
    el.test=1; }

  if (abs(el.test)==2) {
    if (cfg[gear.order+2]==NULL) sdsralloczero(cfg[gear.order+2],cfg[0]->size) }

  A=cfg[gear.order+abs(el.test)];

  time(&time0);

  loop (ii,0,no) {
#ifdef POLAR
    scf.nit=0;
    do { /* iterate self-field */
#endif /*# POLAR */
      zeroEn();
      forces(A,cfg[0]);
#ifdef POLAR
      put2(scf.nit,En.el)
    } while (!selffield(A,cfg[0],scf.eps,scf.omega,1));
#endif /*#!POLAR */
  }

#ifdef POLAR
  En.pot += En.self;
#endif /*# POLAR */
  En.pot += En.el; /* do not move to forces() */

  time(&time1);

#ifdef LINKCELL
  prt("No.cell=(%d,%d,%d) cutoff=%g finished in %li/%i s",
      No.cell[0],No.cell[1],No.cell[2],cutoff,time1-time0,no);
#else /*# LINKCELL */
  prt("cutoff=%5.1f finished in %li/%i s",
      cutoff,time1-time0,no);
#endif /*#!LINKCELL */

#ifdef POLAR
  prt("POLAR: %d iterations (scf.eps=%g scf.omega=%g)",scf.nit,scf.eps,scf.omega);
#endif /*# POLAR */

  if (abs(el.test)==1) { /* el.test == +-1 */
    df=maxf=0;
    d=A->rp;
    dref=cfg[gear.order+2]->rp;
    loop (ii,0,ns) {
      VVV(dr,=d[ii],-dref[ii])
      df += ff = SQR(dr);
      Max(maxf,ff) }
    prt("\nTTTTTT Ewald test TTTTTT\nEn.pot = %.7f",En.pot);

    if (etestref) {
      prt("\
En.pot-Epotref = %.3e K (with diagonal corr. %.3e)\n\
summarized standard deviation of forces = %.3e K/AA = %.3e rel.\n\
calc.deviation/estimate = %.4g (should be 1 if accurate reference)\n\
maximum  deviation of forces = %.3e K/AA = %.3e rel.",
          En.pot-Epotref,
          En.pot+Ediag-Epotref-Ediagref,
          sqrt(df),sqrt(df/dfref),
          sqrt(df/(1e-99+Sqr(estimr)+Sqr(estimk))),
          sqrt(maxf),sqrt(maxf/dfref));
      prt("deviation of Pvir = %g [Pa]",En.vir*Punit/(box.L[0]*box.L[1]*box.L[2]*DIM)-Evirref);
#if PRESSURETENSOR&PT_VIR
      prt_("deviation of virial tensor =");
      loop (ii,0,PT_DIM) prt_(" %g",En.Pvir[ii]-Pvirref[ii]);
      prt("\nPvir-diag(Pvir)/3 = %g BUG HERE (units? /V?)",Evirref-SUM(En.Pvir)/DIM);
#endif /*# PRESSURETENSOR&PT_VIR */
    }

    goto again;
  } /* el.test==+-1 */


 setref: /* el.test=+=2 or init=20 */
  etestref=1;
  el.test=el.test/abs(el.test);
  Epotref=En.pot; Ediagref=Ediag;
  dfref=maxfref=0;
  d=A->rp;
  loop (ii,0,ns) { dfref += ff = SQR(d[ii]); Max(maxfref,ff); }
  varf=sqrt(dfref/No.s);
  prt(
      "\nTTTTTT Ewald test TTTTTT reference set TTTTTT\n"
      "Epot = %.14g K\n\
sqrt(sum forces^2) = %.14g K/AA max = %.10g",
      Epotref,sqrt(dfref),maxfref);

  prt("Pvir (from En.vir) = %g (all pressures in [Pa])",
      Evirref=En.vir*Punit/(PROD(box.L)*DIM));
#if PRESSURETENSOR&PT_VIR
  prt_("virial tensor:");
  loop (ii,0,PT_DIM)
    prt_(" %s%.9g",ii==3?" ":"",Pvirref[ii]=En.Pvir[ii]);
  prt("\nPvir-diag(Pvir)/3 = %g",Evirref-SUM(En.Pvir)/DIM);
#endif /*# PRESSURETENSOR&PT_VIR */

  goto again;
} /* Ewaldtest */
