/* this module depends on #defines WORM and CUT */
/* PARALLEL=3 support removed - see old+misc/simpot-par3.c, par3.c, parinter.c */

#include "ground.h"
#include "simglob.h"
#include "simdef.h"
#include "siminit.h"
#include "simpot.h"
#include "interpot.h"

#ifndef LINKCELL
#  include "water.h"
#endif /*# LINKCELL */

#ifdef DUMP
/* WARNING: not correct for bonded */
#  define DEBUGU(X) { \
  prt("%9s = %10.5f  total=%10.5f kcal/mol",\
  #X,(U+En.el-U0)*0.0019872,(U+En.el)*0.0019872); U0=U+En.el; }
double U0=0;
#else /*# DUMP */
#  define DEBUGU(X)
#endif /*#!DUMP */

#ifdef LINKCELL

double intramol(vector f1[],molecule_t *m1,vector *rp) /*********** intramol */
/* only bonded forces (all pair forces solved by the link-cell method) */
{
  int i,j,k,l;
  double U=spec[m1->sp]->zero_energy;
  vector *r1=rof(m1,rp);
  bond_t *b;
  angle_t *an;
  torsion_t *t;
  specinfo_t *sp;

#  include "bonded.c"

  return U;
} /* intramol */

double intermol(vector f1[],vector f2[],molecule_t *m1,molecule_t *m2,vector *rp)
/* LINKCELL: only for configuration initializer ******************** intermol */
{
  int i,j;
  double U=0,q,qq;
  vector *r1=rof(m1,rp), *r2=rof(m2,rp);
  real *r1i,*f1i;
  int ns1=m1->ns,ns2=m2->ns;
  siteinfo_t *si1=spec[m1->sp]->si;
  siteinfo_t *si2=spec[m2->sp]->si;
  sitesite_t *ss;

#  define GLOB
#  include "cookmeas.c"
#  include "cookmm.c"
#  include "cookundf.c"
#  undef GLOB

  return U;
} /* intermol */

pot2_t *Potential(int i,int j) /********************************** Potential */
{
  i=spec[i]->pot; j=spec[j]->pot;
  if (i==5 && j==5) ERROR(("ST2 not supported for LINKCELL"))

  return intermol;
}

#else /*# LINKCELL */

double intermol(vector f1[],vector f2[],molecule_t *m1,molecule_t *m2,vector *rp)
/* not-LINKCELL general function for inter and intra molecular forces */
/* ? to be unrolled for m1==m2 and m1!=m2 */
{
  int i,j,k,l;

#  ifndef POLAR
  double q,qq;
#  endif /*# POLAR */

  vector *r1=rof(m1,rp);
  vector *r2=rof(m2,rp);

  real *r1i,*f1i;
  int ns1=m1->ns,ns2=m2->ns;
  bond_t *b;
  angle_t *an;
  torsion_t *t;
  specinfo_t *sp;
  sitesite_t *ss;
  siteinfo_t *si1=spec[m1->sp]->si;
  siteinfo_t *si2=spec[m2->sp]->si;
  double U=0;

  int j0,j1;
  exception_t *exc;

  if (m1!=m2) { /*** intermolecular interactions ***/

    if (measure) {

#  define GLOB /* used to be active with PARALLEL=3 - kept for compatibility */

#  include "cookmeas.c"
#  include "cookmm.c"
#  include "cookundf.c"
#ifdef LOG
      En.LJmn=U;
#endif
    }
#  undef GLOB
    else {
#  define GLOB
#  include "cooknmea.c"
#  include "cookmm.c"
#  include "cookundf.c"
#  undef GLOB
    }
  } /* intermolecular */

  else { /*** intramolecular interactions ***/
    /* m1 == m2 */
    /* all site-site interactions */

#  define GLOB
#  define I0 0
#  define I1 ns1

    if (measure) {
#  include "cookmeas.c"
#  include "cookm.c"
#  include "cookundf.c"
#ifdef LOG
      En.LJmn=U;
#endif
    }
    else {
#  include "cooknmea.c"
#  include "cookm.c"
#  include "cookundf.c"
    }

#  ifdef PERSUM
    /*** image-image interaction of the same molecule */
    if (measure) {
#    define GLOB /* used to be active with PARALLEL=3 - kept for compatibility */

#    include "cookmeasps.c"
#    include "cookmm.c"
#    include "cookundf.c"
#ifdef LOG
      En.LJmn=U;
#endif
    }
#    undef GLOB
    else {
#    define GLOB
#    include "cooknmeaps.c"
#    include "cookmm.c"
#    include "cookundf.c"
#    undef GLOB
    }
#  endif /*# PERSUM */

    DEBUGU(qq+LJ)
    {
      /* bonded interactions: dirty trick!!! */
      double Ubonded=0;
#  define U Ubonded
#  include "bonded.c"
#  undef U
      En.bonded+=Ubonded;
#ifdef LOG
      if (No.first) if (m1-molec<No.first) En.bonded0+=Ubonded;
#endif
      U+=Ubonded;
    }

    U += spec[m1->sp]->zero_energy;

  } /* intramolecular */

#  if 0
  prt("<0>: f1[0][0]=%g f1[%d][2]=%g U=%g En.el=%g En.vir=%g (simpot.c)",
      f1[0][0],ns1-1,f1[ns1-1][2],U,En.el,En.vir);
#  endif /*# 0 */

  return U;
} /* intermol */

double atomatom(vector f1[],vector f2[],molecule_t *m1,molecule_t *m2,vector *rp)
/* not LINKCELL optimized for monoatomic molecules */
{
#  ifndef POLAR
  double q,qq;
#  endif /*# POLAR */

  vector *r1=rof(m1,rp);
  vector *r2=rof(m2,rp);
  sitesite_t *ss;
  siteinfo_t *si1=spec[m1->sp]->si;
  siteinfo_t *si2=spec[m2->sp]->si;

  double U=0;

  if (m1!=m2) { /*** intermolecular interactions ***/

    if (measure) {
#  define GLOB
#  include "cookmeas.c"
#  include "cookaa.c"
#  include "cookundf.c"
#ifdef LOG
      En.LJmn=U;
#endif
    }
#  undef GLOB
    else {
#  define GLOB
#  include "cooknmea.c"
#  include "cookaa.c"
#  include "cookundf.c"
#  undef GLOB
    } }

  else {
#ifdef LOG
      En.LJmn=U;
#endif
    U += spec[m1->sp]->zero_energy;

#  ifdef PERSUM
    /*** image-image interaction of the same molecule (here atom) */
    if (measure) {
#    define GLOB /* used to be active with PARALLEL=3 - kept for compatibility */

#    include "cookmeasps.c"
#    include "cookaa.c"
#    include "cookundf.c"
#ifdef LOG
      En.LJmn=U;
#endif
  }
#    undef GLOB
    else {
#    define GLOB
#    include "cooknmeaps.c"
#    include "cookaa.c"
#    include "cookundf.c"
#    undef GLOB
    }
#  endif /*# PERSUM */
 }

  return U;
} /* atomatom */


pot2_t *Potential(int i,int j) /********************************** Potential */
{
  if (spec[i]->ns==1 && spec[j]->ns==1) return atomatom;

#  ifndef POLAR
  i=spec[i]->pot; j=spec[j]->pot;
  if (i==j && option('x')&1) switch (i) {
    case 3: return TIP3P;
    case 4: return TIP4P;
    case 5: return ST2;
#    ifdef WATERPLUS
#      if WATERPLUS==2
    case 6: return CO2;
    case 7: return MeOH;
    case 8: return acetone;
    case 9: return dljd;
#      endif /*# WATERPLUS==2 */
#      if WATERPLUS==3
    case 10: return ionion; /* see waterion.c and watercut.c */
#      endif /*# WATERPLUS==3 */
#    endif /*# WATERPLUS */
    }

#    if 0
  if (i==10 && j==3 || i==3 && j==10) return TIP3Pion; /* see waterion.c and watercut.c */
#    endif /*# 0 */

#  endif /*! POLAR */ /*# POLAR */

  return intermol;
}

#endif /*#!LINKCELL */


/******************************* two walls ******************************/

#ifdef SLAB
double mol_wall(vector f[],molecule_t *m,vector *rp) /************* mol_wall */
/*
  the molecule-wall potential, including gravity
  (could be optimized by preparing tables of LJ)
*/
{
  int ns=m->ns,i;
  vector *r=rof(m,rp);
  siteinfo_t *si=spec[m->sp]->si;
  double U=0,f1;

  loop (i,0,ns) {
    double z=r[i][2],zz=z*z;

#  ifdef GOLD
    double qq=Sqr(si[i].charge);

#    ifdef FREEBC
    f[i][2] -= qq/(4*zz);
    if (measure) En.el -= qq/(4*z);
#    else /*# FREEBC */
#      if COULOMB<0
#        error COULOMB<0 not supported
#      endif /*# COULOMB<0 */
    /* this is CUTELST version: cut-and-shift Coulomb needed */
    goldQQ(f[i],2*z,qq); /* updates En.el */
#    endif /*#!FREEBC */
#  endif /*# GOLD */

    if (wall.g) {
      /* gravity in the direction of z (the usual gravity is NEGATIVE) */
      f1=si[i].mass*wall.g;
      f[i][2]+=f1;
      if (measure) {
        U-=(si[i].mass*wall.g)*z;
        En.vir+=f1*(z-box.Lh[2]);
#if PRESSURETENSOR&PT_VIR
        En.Pvir[2]+=f1*(z-box.Lh[2]);
#endif
      } }

    if (wall.n) {
      struct wsstab_s *wss=wsstab+si[i].st;
      double z6;

      if (abs(wall.n)&1) {
        /* wall 0 */
        z-=box.L[2]*wall.z[0];
        zz=z*z; z6=zz*zz*zz;
        f1=(wss->AA/z6-wss->BB)/(zz*zz);
        f[i][2]+=f1;
        if (measure) {
          En.Pwall[0]+=f1*wall.Punit_A;
          U+=(wss->A/z6-wss->B)/(z*zz);
          En.vir+=f1*z;
#if PRESSURETENSOR&PT_VIR
          En.Pvir[2]+=f1*z;
#endif
        } }

      if (abs(wall.n)&2) {
        /* wall 1 */
        z=box.L[2]*wall.z[1]-r[i][2];
        zz=z*z; z6=zz*zz*zz;
        f1=(wss->AA/z6-wss->BB)/(zz*zz);
        f[i][2]-=f1;
        if (measure) {
          En.Pwall[1]+=f1*wall.Punit_A;
          U+=(wss->A/z6-wss->B)/(z*zz);
          En.vir+=f1*z;
#if PRESSURETENSOR&PT_VIR
          En.Pvir[2]+=f1*z;
#endif
        } }
    }
  } /* i */
  return U;
}
#endif /*# SLAB */
