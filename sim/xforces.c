/* some extra forces */

/* WARNING: the empty userforces.c and xforces.c MUST NOT be placed in the
   same directory, otherwise userforces.c in MACSIMUS/cook/PROJECT will not
   be found
*/

#include "ground.h"
#include <time.h>
#include "simglob.h"
#include "norm.h"
#include "ewald.h"
#include "forces.h"
#include "xforces.h"
#include "simmeas.h"
#include "cputime.h"
#include "intermac.h" /* because of fixsites */
#include "statics.h" /* cleaving */
#include "units.h" /* cleaving, wall */

#ifdef SLAB
#  include "simdef.h" /* nsites */
#  include "interpot.h" /* sstab */
#  include "simpot.h"
#endif /*# SLAB */

void fixforces(ToIntPtr B, ToIntPtr A) /************************** fixforces */
{
  vector *rpB=B->rp,*rpA=A->rp,*rpf;
#ifdef GOLD
  double xxyy;
#endif /*# GOLD */
  real *r1,*r2,*f,rr;
  vector dr;
  int i;
  fixsites_t *fix;
  double U=0,K=option('k');

  if (!fixa || K==0) return;
  rpf=fixa->rp;

  looplist (fix,fixsites0)
    loop (i,fix->from,fix->to) {
      r1=rpA[i];
      r2=rpf[i];
      //    put(i) putv(r1) putv(r2 pf[i])

      beginCellImage(,box.cq)
        f=rpB[i];
        if (measure) U+=K*rr;
        VV(f,-=K*dr)
      } else { /* elseCellImage */
        static int count=0;
        double KK=K*box.cq/rr; /* smaller force outside cutoff sphere */

        f=rpB[i];
        /* U wrong for this reduced force */
        VV(f,-=KK*dr)
        if (count<10) {
          WARNING(("too far from the fix location (more than 10 warnings will be suppressed)\n\
*** fix=(%g %g %g), real=(%g %g %g)",
                   VARG(r2),VARG(r1)))
            count++; }
      endCellImage
      }
  U/=2;
  if (measure && option('v')&4) prt("fixE=%g", U);

  En.pot += En.fix=U;
}

void centerforces(ToIntPtr B, ToIntPtr A) /******************** centerforces */
/*
   additional forces, originally "to the center of the box"
*/
{
  if (center.cmn) {
    /* a harmonic force is adeed to the whole configuration (as slab) */
    int ns,i,k,n;
    vector CM,K2;
    double mi;
    molecule_t *mn;
    siteinfo_t *si;
    vector *r,*f;

    VO(CM,=0)

    if (center.cmn>No.N) center.cmn=No.N;

    loop (n,FROM,center.cmn) {
      mn=molec+n;
      ns=mn->ns;
      r=rof(mn,A->rp);
      f=rof(mn,B->rp);
      si=spec[mn->sp]->si;

      loop (i,0,ns) {
        mi=si[i].mass;
        VV(CM,+=mi*r[i]) } }

    VO(CM,/=center.cmmass)
#ifndef FREEBC
    VV(CM,-=box.Lh) /* in periodic b.c. centereded to the center of the slab */
#endif /*# FREEBC */

    loop (k,0,DIM) En.pot+=center.cmK[k]*Sqr(CM[k]);
    VV(K2,=2/center.cmmass*center.cmK)

    loop (n,FROM,center.cmn) {
      mn=molec+n;
      ns=mn->ns;
      r=rof(mn,A->rp);
      f=rof(mn,B->rp);
      si=spec[mn->sp]->si;

      loop (k,0,DIM) if (center.cmK[k])
        loop (i,0,ns)
          f[i][k]-=K2[k]*si[i].mass*CM[k]; }
    }

  /*
    harmonic force to the center (by particles);
    FREEBC: center=(0,0,0)
    periodic b.c.: center=(Lx,Ly,Lz)/2
  */

  if (center.on) {
    int ns,i,k,n;
    molecule_t *mn;
    vector *r,*f;

#ifdef SLAB
#  if SLAB & 2
    int ifrom,ito,icleave;
    double locWcleave=0;

    En.fcleave[0]=En.fcleave[1]=0;
    En.WPhi=0; /* PROBABLY NOT NEEDED... */
#  endif /*# SLAB & 2 */
#endif /*# SLAB */

    loop (n,FROM,No.N) {
      mn=molec+n;
      ns=mn->ns;
      r=rof(mn,A->rp);
      f=rof(mn,B->rp);

      /* forces to box center */
      if (center.on&2) loop (i,0,ns) loop (k,0,DIM) {
#ifdef FREEBC
        double dr=r[i][k];
#else /*# FREEBC */
        double dr=r[i][k]-box.Lh[k];
#endif /*#!FREEBC */

        if (fabs(dr)>center.r0[k]) {
          if (dr<0) dr+=center.r0[k]; else dr-=center.r0[k];
          En.pot+=center.K[k]*Sqr(dr);
          f[i][k]-=center.K2[k]*dr; } }

#ifdef SLAB
      /* z (slab) forces */
#  if SLAB & 2
      /* cumbersome code: */
      if (cleave.i<0) ifrom=0,ito=abs(cleave.i);
      else ifrom=cleave.i,ito=cleave.i;
      if (ito>=ns) ito=ns-1;

      loop (icleave,0,cleave.n) {
        loopto (i,ifrom,ito) {
          double zc=box.L[2]*cleave.z[icleave];
          double dz=r[i][2]-zc,x,ff,uu;

#    ifndef FREEBC
          /* WARNING: FREEBC version not checked */
          if (dz>box.Lh[2]) dz-=box.L[2];
          if (dz<-box.Lh[2]) dz+=box.L[2];
#    endif /*# FREEBC */
          x=dz/cleave.sigma;
          if (fabs(x)>=1) continue;
          if (fabs(x)<0.5) ff=x,uu=0.5-x*x;
          else if (x>0) ff=1-x,uu=ff*ff;
          else ff=-x-1,uu=ff*ff;
          ff*=cleave.K; /* now ff=force */
          En.fcleave[icleave]+=ff;
          f[i][2]+=ff;
          En.pot+=uu*=cleave.Ksigmah; /* now uu=cleaving potential Phi */
          En.vir+=dz*ff;
#    if PRESSURETENSOR&PT_VIR
          En.Pvir[2]+=dz*ff;
#    endif /*# PRESSURETENSOR&PT_VIR */
          En.WPhi+=(uu)/(cleave.n*box.L[0]*box.L[1]*cleave.sigma); /* PROBABLY NOT NEEDED... */
          locWcleave+=(uu+dz*ff)/(cleave.n*box.L[0]*box.L[1]*cleave.sigma); }
      }
#  endif /*# SLAB & 2 */

#  if SLAB & 4
      /* LJ wall, old simple and inconsistent code (use wall.n etc) */
      if (center.on&8) loop (i,0,ns) {
          double dz=r[i][2],dz3=dz*dz*dz;

          En.pot+=slab.wall.A/Cub(dz3)-slab.wall.B/dz3;
          f[i][2]+=(9*slab.wall.A/Cub(dz3)-3*slab.wall.B/dz3)/dz; }
#  endif /*# SLAB & 4 */

      if (center.on&1) loop (i,0,ns) loop (k,0,NCENTER) {
        if (slab.nn[k]<=0) break;
        if (n<slab.nn[k]) {
#  ifdef FREEBC
          double dz=r[i][2]-slab.z[k]; /* ??? */
#  else /*# FREEBC */
          double dz=r[i][2]-(box.Lh[2]+slab.z[k]);
#  endif /*#!FREEBC */

          if (fabs(dz)>slab.z0[k]
              && (slab.ns[k]==0
                  || (slab.ns[k]>0 && i<slab.ns[k])
                  || (slab.ns[k]<0 && i>=slab.ns[k]))) {

            if (dz<0) dz+=slab.z0[k]; else dz-=slab.z0[k];

//#    if SLAB & 1 !!! for SLAB always active since V2.8c, SLAB&1 now = zone melting
            if (slab.z1[k]) {
              /*                              ___     ___ _Kz2 */
              /* slab/outside bias function:     \___/    _0   */
              if (fabs(dz)>slab.dz[k])
                En.pot+=slab.Kz2[k];
              else if (fabs(dz)<slab.dz[k]/2) {
                En.pot+=slab.Kz[k]*Sqr(dz);
                f[i][2]-=2*slab.Kz[k]*dz; }
              else {
                if (dz<0) dz=-slab.dz[k]-dz; else dz=slab.dz[k]-dz;
                En.pot+=slab.Kz2[k]-slab.Kz[k]*Sqr(dz);
                f[i][2]-=2*slab.Kz[k]*dz; } }
            else
//#    endif /*# SLAB & 1 */
            {
              /* harmonic or U-shaped forces */
              En.pot +=slab.Kz[k]*Sqr(dz);
              f[i][2]-=slab.Kz2[k]*dz;
            } } /* ns check */
          break; } }
#endif        /*# SLAB */
    } /* n */
#if defined(SLAB) && SLAB & 2
    En.Wcleave=locWcleave;
    StaSet(0,lag.err,2,lag.n);
    StaAdd("WPhi [Pa]",En.WPhi*Punit); /* PROBABLY NOT NEEDED... */
    StaAdd("Wcleave [Pa]",En.Wcleave*Punit);
#endif /*# defined(SLAB) && SLAB & 2 */
  }
}

#ifdef SLAB
void wallforces(ToIntPtr B, ToIntPtr A)  /*********************** wallforces */
/*
  1/ atom-wall integrated LJ potential
  2/ GOLD version: Coulomb potential of charge vs. its inverted image
*/
{
  vector *f;
  int n;
  molecule_t *mn;
  double U=0;

  if (!wall.is) return;

  wall.Punit_A=Punit/(box.L[0]*box.L[1]);

  loop (n,FROM,No.N) {
    mn=molec+n;
    f=rof(mn,B->rp);
    U += mol_wall(f,mn,A->rp); }

  En.pot += U;
}
#endif  /*# SLAB */

#ifdef POLAR
void elstforces(ToIntPtr B, ToIntPtr A)  /********************** elstforces */
/*
  homogeneous electrostatic field Eprogunits
*/
{
  vector *f,*fpol,*r,*rpol;
  int i,n;
  molecule_t *mn;
  siteinfo_t *si;
  double U=0;
  vector Et;

  if (Eext.isE) {

    loop (i,0,DIM) {
      Eext.arg[i]=2*PI*(Eext.f*t-Eext.phase[i]);
      Et[i]=Eext.E[i]*cos(Eext.arg[i]); }

    loop (n,FROM,No.N) {
      mn=molec+n;
      f=rof(mn,B->rp);
      fpol=polarrof(mn,B->rp);
      r=rof(mn,A->rp);
      rpol=polarrof(mn,A->rp);
      si=spec[mn->sp]->si;
      loop (i,0,mn->ns)
#  if POLAR&32
        if (si[i].qtype&32) {
          double qpol=rpol[i][0];
          VV(f[i],+=qpol*Et)
#error this line looks as nonsense:
          fpol[i][0]-=SCAL(Et,r[i]); /* doublecheck sign! */
          U-=qpol*SCAL(Et,r[i]);
        }
      else
#  endif /*# POLAR&32 */
        {
          double q=si[i].charge;
          double qpol=si[i].chargepol;
          vector rpi;

          VVV(rpi,=r[i],+rpol[i])

          VV(f[i],+=q*Et)
          VV(fpol[i],+=qpol*Et)
          /* DO NOT USE -- handled by the virial theorem
             En.vir+=SCAL(r[i],q*Et)
             En.vir+=SCAL(r[i],qpol*Et) */
#if PRESSURETENSOR&PT_VIR
          En.Pvir[0]+=r[i][0]*q*Et[0];
          En.Pvir[1]+=r[i][1]*q*Et[1];
          En.Pvir[2]+=r[i][2]*q*Et[2];
          En.Pvir[0]+=rpi[0]*qpol*Et[0];
          En.Pvir[1]+=rpi[1]*qpol*Et[1];
          En.Pvir[2]+=rpi[2]*qpol*Et[2];
#endif
#  if PRESSURETENSOR&PT_OFF
          /* doublecheck symmetry */
          En.Pvir[3]+=r[i][1]*q*Et[2];
          En.Pvir[4]+=r[i][2]*q*Et[0];
          En.Pvir[5]+=r[i][0]*q*Et[1];
          En.Pvir[3]+=rpi[1]*qpol*Et[2];
          En.Pvir[4]+=rpi[2]*qpol*Et[0];
          En.Pvir[5]+=rpi[0]*qpol*Et[1];
#endif
          /* doublecheck why not if (measure) but always */
          U-=q*SCAL(Et,r[i]);
          U-=qpol*SCAL(Et,rpi); } }
    En.el += U;
  }
}
#else /*# POLAR */
void elstforces(ToIntPtr B, ToIntPtr A)  /********************** elstforces */
/*
  external electric (homogeneous or oscillating) field
  external magnetic field
*/
{
  vector *f,*r;
  int i,n;
  molecule_t *mn;
  siteinfo_t *si;
  double U=0;
  vector Et;

  if (Eext.ism) {
    /* magnetic dipole interacting with the external field */
    /* NB: the magnetic dipoles do not interact mutually! */

    loop (n,FROM,No.N) {
      mn=molec+n;
      f=rof(mn,B->rp);
      r=rof(mn,A->rp);
      if (mn->sp==el.m.sp || el.m.sp<0 && el.m.sp+mn->sp>=0) {
        double rr=SQRD(r[el.m.plus],r[el.m.minus]);
        double qm=Eext.m/sqrt(rr);

        /* BUG: virial probably wrong*/
        VV(f[el.m.plus],+=qm*Eext.B)
#if PRESSURETENSOR&PT_VIR
        En.Pvir[0]+=r[el.m.plus][0]*qm*Eext.B[0];
        En.Pvir[1]+=r[el.m.plus][1]*qm*Eext.B[1];
        En.Pvir[2]+=r[el.m.plus][2]*qm*Eext.B[2];
#endif
#  if PRESSURETENSOR&PT_OFF
        /* doublecheck symmetry */
        En.Pvir[3]+=r[el.m.plus][1]*qm*Eext.B[2];
        En.Pvir[4]+=r[el.m.plus][2]*qm*Eext.B[0];
        En.Pvir[5]+=r[el.m.plus][0]*qm*Eext.B[1];
#endif
        U-=qm*SCAL(Eext.B,r[el.m.plus]);

        VV(f[el.m.minus],-=qm*Eext.B)
#if PRESSURETENSOR&PT_VIR
        En.Pvir[0]+=r[el.m.minus][0]*qm*Eext.B[0];
        En.Pvir[1]+=r[el.m.minus][1]*qm*Eext.B[1];
        En.Pvir[2]+=r[el.m.minus][2]*qm*Eext.B[2];
#endif
#  if PRESSURETENSOR&PT_OFF
        /* doublecheck symmetry */
        En.Pvir[3]+=r[el.m.minus][1]*qm*Eext.B[2];
        En.Pvir[4]+=r[el.m.minus][2]*qm*Eext.B[0];
        En.Pvir[5]+=r[el.m.minus][0]*qm*Eext.B[1];
#endif
        U+=qm*SCAL(Eext.B,r[el.m.minus]); } } }

  if (Eext.isE) {
    loop (i,0,DIM) {
      Eext.arg[i]=2*PI*(Eext.f*t-Eext.phase[i]);
      Et[i]=Eext.E[i]*cos(Eext.arg[i]); }

    loop (n,FROM,No.N) {
      mn=molec+n;
      f=rof(mn,B->rp);
      r=rof(mn,A->rp);
      si=spec[mn->sp]->si;
      loop (i,0,mn->ns) {
        double q=si[i].charge;

        VV(f[i],+=q*Et)
#if PRESSURETENSOR&PT_VIR
        En.Pvir[0]+=r[i][0]*q*Et[0];
        En.Pvir[1]+=r[i][1]*q*Et[1];
        En.Pvir[2]+=r[i][2]*q*Et[2];
#endif
#  if PRESSURETENSOR&PT_OFF
        /* doublecheck symmetry */
        En.Pvir[3]+=r[i][1]*q*Et[2];
        En.Pvir[4]+=r[i][2]*q*Et[0];
        En.Pvir[5]+=r[i][0]*q*Et[1];
#endif
        U-=q*SCAL(Et,r[i]); } } }
  En.el += U; /* also mag energy - will break virial check */
}
#endif /*#!POLAR */

#ifdef SLAB
void slabcutcor(ToIntPtr B, ToIntPtr A)  /********************** slabcutcor */
/*
  slab cutoff correction, called if slab.K
  wave numbers kz < slab.K used (kz=slab.K is NOT included)
  site type density profiles calculated HERE for all site types
  WARNING: includes inefficiently also sites with zero potential or not populated
*/
{
  vector *f,*r;
  int i,st,k,n;
  molecule_t *mn;
  siteinfo_t *si;
  double Resum=0; /* sum of atom energy corrections: debugging only */
  double U=0; /* as above, concise FT-based formula */
  double virsum=0; /* virial contribution, concise FT-based formula */
  complex **rho; /* [st][k] FT of z-profile (complex defined in ewald.h) */

  if (slab.K<2 || slab.K>1024) ERROR(("slabcutcor: slab.K=%d is invalid or out of reasonable range",slab.K))

  /* FT of the site z-profiles multiplied by V (ss->Skk are divided)*/
  alloc2Darrayzero(rho,nsites,slab.K);

  /* Fourier transform of the particle density multiplied by V (Skk is divided by V) */
  loop (n,FROM,No.N) {
    mn=molec+n;
    f=rof(mn,B->rp);
    r=rof(mn,A->rp);
    si=spec[mn->sp]->si;

    loop (i,0,mn->ns) if (si[i].isLJ) {
      double a=2*PI/box.L[2]*r[i][2];
      complex q={cos(a),sin(a)},x=q;

      st=si[i].st;

      rho[st][0].re+=1;
      rho[st][1].re+=x.re;
      rho[st][1].im-=x.im;

      loop (k,2,slab.K) {
        /* x=q^k */
        a=x.re*q.re-x.im*q.im;
        x.im=x.re*q.im+x.im*q.re;
        x.re=a;
        rho[st][k].re+=x.re;
        rho[st][k].im-=x.im; } } }

  /* corrected energy and forces */
  loop (n,FROM,No.N) {
    mn=molec+n;
    f=rof(mn,B->rp);
    r=rof(mn,A->rp);
    si=spec[mn->sp]->si;
    loop (i,0,mn->ns) {
      int st=si[i].st;
      int jst;
      double fsum=0;
      double a=2*PI/box.L[2]*r[i][2];
      complex q={cos(a),sin(a)},x=q;

      loop (jst,0,nsites) {
        sitesite_t *ss=&sstab[st][jst];

        if (ss->Skk) {
          if (measure) {
            Resum+=ss->Skk[0].E*rho[jst][0].re;
            Resum+=ss->Skk[1].E*(rho[jst][1].re*x.re-rho[jst][1].im*x.im); }
          fsum+=ss->Skk[1].E*(rho[jst][1].re*x.im+rho[jst][1].im*x.re);

          loop (k,2,slab.K) {
            a=x.re*q.re-x.im*q.im;
            x.im=x.re*q.im+x.im*q.re;
            x.re=a;
            if (measure) Resum+= ss->Skk[k].E*(rho[jst][k].re*x.re-rho[jst][k].im*x.im);
            fsum+=k*ss->Skk[k].E*(rho[jst][k].re*x.im+rho[jst][k].im*x.re); } } }

      //      if (option('v')&256) put2(r[i][2],Resum)
      f[i][2]+=4*PI/box.L[2]*fsum; } }

  if (measure) {
#  if 0
    /* verbose debug prints */
    if (option('v')&256) {
      loop (st,0,nsites) loop (k,0,slab.K)
        prt("rho[%d][%d] = %g %g",st,k,rho[st][k].re/box.V,rho[st][k].im/box.V); }
#  endif /*# 0 */

    /* concise energy and virial formulas */
    loop (k,0,slab.K) {
      int i,j;

      loop (i,0,nsites)
        loop (j,0,nsites) if (sstab[i][j].Skk) {
          U     +=(rho[i][k].re*rho[j][k].re+rho[i][k].im*rho[j][k].im)*sstab[i][j].Skk[k].E;
          virsum+=(rho[i][k].re*rho[j][k].re+rho[i][k].im*rho[j][k].im)*sstab[i][j].Skk[k].D; } }

    //prt("E Ecorr=%.9g Resum=%.8g\n",U,Resum);

    En.pot+=U;

#  if PRESSURETENSOR&PT_VIR
    /* does not make too much sense without PRESSURETENSOR&PT_VIR */
    VO(En.Pvir,+=U)
    En.Pvir[2]-=virsum;
#  endif /*# PRESSURETENSOR&PT_VIR */
    En.vir+=3*U-virsum;

    StaSet(0,2,2,lag.n);
    StaAdd("Eslabcorr [J/mol]",U*Eunit); /* no SI switch as in main.c */
    // debug old version: /box.V*Punit (total *1.1045188e-06)
    /* homogeneous contribution to the slab correction (does not affect surface tension) */
    StaAdd("Pslabcorr.xyz [Pa]",U/box.V*Punit);
    /* extra z-contribution to the slab correction */
    StaAdd("Pslabcorr.z [Pa]",-virsum/box.V*Punit); }

  //  if (option('v')&256) prt("D E=%g virsum=%g En.Pvir[2]=%g",U,virsum,U-virsum);

  free2Darray(rho);
}
#endif /*# SLAB */

#include "userforces.c"
