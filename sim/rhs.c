#include "ground.h"
#include "statics.h"
#include "simglob.h"
#include "norm.h"
#include "forces.h" /* the same for LINKCELL */
#include "constrd.h"
#include "rhs.h"
#include "simgear.h"                                                            

void rhs(ToIntPtr B, ToIntPtr A, ToIntPtr V) /* ======================== rhs */
{
  int n,ns,sp,i;
  vector *v,p_mol;
  double m_mol,mi;
  molecule_t *mn;
  siteinfo_t *si;

  if (thermostat>=T_NPT) {
    /* virial pressure would be needed for the T_NPT algorithm at every step 
       (but it has not been implemented for Gear yet) */
    measure=1;
    if (gear.order>2) ERROR(("NPT/MTK not implemented for Gear")) }

#ifndef FREEBC
  VV(box.Lh,=0.5*box.L)
#endif /*# FREEBC */

  box.cq=Sqr(box.cutoff);
  En.kin=En.kin_tr=En.kin_in=0;

  /* Now all kinetic energy calculations via Gear and rhs moved here.
     In old versions they were integrated into constrd, which was
     slightly more efficient, but the inter and intramolecular energies
     were not always measured */

#  if PRESSURETENSOR&PT_KIN
  memset(En.Pkin,0,sizeof(En.Pkin));
#    if PRESSURETENSOR&PT_MOL
  memset(En.PKin,0,sizeof(En.PKin));
#      if PRESSURETENSOR&PT_MOM
  memset(En.PKin2,0,sizeof(En.PKin));
  memset(En.PKin3,0,sizeof(En.PKin));
#      endif /*# PRESSURETENSOR&PT_MOM */
#    endif /*# PRESSURETENSOR&PT_MOL */
#  endif /*# PRESSURETENSOR&PT_KIN */

  if (thermostat>=T_TR || measure) {
    /* also intra/inter Ekin needed */
    loop (n,FROM,No.N) {
      mn=molec+n;
      ns=mn->ns;
      sp=mn->sp;
      si=spec[sp]->si;
      v=rof(mn,V->rp);
      VO(p_mol,=0)
      m_mol=0;
      loop (i,0,ns) {
        mi=si[i].mass;
        En.kin += SQR(v[i])*mi;
#  if PRESSURETENSOR&PT_KIN
        /* site-based */
        En.Pkin[0]+=v[i][0]*v[i][0]*si[i].mass;
        En.Pkin[1]+=v[i][1]*v[i][1]*si[i].mass;
        En.Pkin[2]+=v[i][2]*v[i][2]*si[i].mass;
#    if PRESSURETENSOR&PT_OFF
        En.Pkin[3]+=v[i][1]*v[i][2]*si[i].mass;
        En.Pkin[4]+=v[i][2]*v[i][0]*si[i].mass;
        En.Pkin[5]+=v[i][0]*v[i][1]*si[i].mass;
#    endif /*# PRESSURETENSOR&PT_OFF */
#  endif /*# PRESSURETENSOR&PT_KIN */
        m_mol+=mi;
        VV(p_mol,+=mi*v[i]) }
#  if PRESSURETENSOR&PT_MOL
      /* CM-based */
      En.PKin[0]+=p_mol[0]*p_mol[0]/m_mol;
      En.PKin[1]+=p_mol[1]*p_mol[1]/m_mol;
      En.PKin[2]+=p_mol[2]*p_mol[2]/m_mol;

#    if PRESSURETENSOR&PT_MOM
      En.PKin2[0]+=Sqr(p_mol[0]*p_mol[0]/m_mol);
      En.PKin2[1]+=Sqr(p_mol[1]*p_mol[1]/m_mol);
      En.PKin2[2]+=Sqr(p_mol[2]*p_mol[2]/m_mol);
      En.PKin3[0]+=Cub(p_mol[0]*p_mol[0]/m_mol);
      En.PKin3[1]+=Cub(p_mol[1]*p_mol[1]/m_mol);
      En.PKin3[2]+=Cub(p_mol[2]*p_mol[2]/m_mol);
#    endif /*# PRESSURETENSOR&PT_MOM */

#    if PRESSURETENSOR&PT_OFF
      En.PKin[3]+=p_mol[1]*p_mol[2]/m_mol;
      En.PKin[4]+=p_mol[2]*p_mol[0]/m_mol;
      En.PKin[5]+=p_mol[0]*p_mol[1]/m_mol;
#    endif /*# PRESSURETENSOR&PT_OFF */
#  endif /*# PRESSURETENSOR&PT_MOL */
      En.kin_tr+=SQR(p_mol)/m_mol; /* twice intermolecular kin. energy */ }
    En.kin_in=En.kin-En.kin_tr; /* twice intramolecular energy */ }
  else {
    /* only total Ekin needed */
    loop (n,FROM,No.N) {
      mn=molec+n;
      ns=mn->ns;
      sp=mn->sp;
      si=spec[sp]->si;
      v=rof(mn,V->rp);
      loop (i,0,ns)
        En.kin += SQR(v[i])*si[i].mass; } }

  if (measure) {
    En.r1=constrainterror(A,V);
    En.v1=vconstrainterror;
  }

  depend_r(A,0);

#ifdef POLAR
  scf.nit=0;
  do { /* iterate self-field */
#endif /*# POLAR */

    zeroEn();
    forces(B,A);

#ifdef POLAR
  } while (!selffield(B,A,scf.eps,scf.omega,1));
  StaAdd("polar no of iter",scf.nit);
  StaAdd("polar Drude maxdr",scf.maxdr);
#endif /*# POLAR */

  depend_f(A,B);

/*.....if (measure) dipolemoments(A);*/

#ifdef SHEAR
  Shear(B,A,V,1.0);
  /*
    systematic flow kinetic energy removed
    NOTE: in some cases, the rest of the kinetic energy is calculated later
  */
  En.kin-=En.kinshear;
  En.kin_tr-=En.kinshear;
#endif /*# SHEAR */

  constraintdynamics(B,A,V);

  if (option('m')>2) {
    /*** Virtual volume change via <dU/dV> with Gear
         for SHAKE, this is in Shake(), module constrd.c 
         In V 3.6f moved here from maincps, fingers crossed... */
    measureP(1);
    measureP(2); }
} /* rhs */
