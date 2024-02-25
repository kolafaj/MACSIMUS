/*
  this module contains intramolecular (bonded) potential terms
  2020: Morse added (expansion to r^4 only, use #define MORSE)
  2016: bug in pressure tensor for Urey-Bradley fixed
  2008: pressure tensor support added
        sign convention now consistently uses f=-dU/dr
  1999: Urey-Bradley support added (use #define UREY_BRADLEY)
*/

#include "ground.h"
#include "vector.h" /* simglob.h not needed */
#include "intrapot.h"

// #define DUMP

#ifdef BLEND
#  ifdef COOK
#    error both COOK and BLEND #defined
#  endif /*# COOK */
#else /*# BLEND */
#  ifndef COOK
#    error none of COOK and BLEND #defined
#  endif /*# COOK */
#  include "simglob.h"
#  include "units.h"
#endif /*#!BLEND */

#if PRESSURETENSOR&PT_VIR
#  if PRESSURETENSOR&PT_OFF
#    define CALCULATEPRESSURETENSOR(DR,F) \
    En.Pvir[0]+=DR[0]*DR[0]*F; \
    En.Pvir[1]+=DR[1]*DR[1]*F; \
    En.Pvir[2]+=DR[2]*DR[2]*F; \
    En.Pvir[3]+=DR[1]*DR[2]*F; \
    En.Pvir[4]+=DR[2]*DR[0]*F; \
    En.Pvir[5]+=DR[0]*DR[1]*F;
#    define CALCULATEPRESSURETENSORG(G,F) \
    En.Pvir[0]+=G[0]*F[0]; \
    En.Pvir[1]+=G[1]*F[1]; \
    En.Pvir[2]+=G[2]*F[2]; \
    En.Pvir[3]+=G[1]*F[2]; \
    En.Pvir[4]+=G[2]*F[0]; \
    En.Pvir[5]+=G[0]*F[1];
#  else /*# PRESSURETENSOR&PT_OFF */
#    define CALCULATEPRESSURETENSOR(DR,F) \
    En.Pvir[0]+=DR[0]*DR[0]*F; \
    En.Pvir[1]+=DR[1]*DR[1]*F; \
    En.Pvir[2]+=DR[2]*DR[2]*F;
#    define CALCULATEPRESSURETENSORG(G,F) \
    En.Pvir[0]+=G[0]*F[0]; \
    En.Pvir[1]+=G[1]*F[1]; \
    En.Pvir[2]+=G[2]*F[2];
#  endif /*#!PRESSURETENSOR&PT_OFF */
#else /*# PRESSURETENSOR&PT_VIR */
#  define CALCULATEPRESSURETENSOR(DR,F) /* empty */
#  define CALCULATEPRESSURETENSORG(G,F) /* empty */
#endif /*#!PRESSURETENSOR&PT_VIR */

extern int measure; /* bad programming habit */

double phi; /* if measure==2 then phi is returned here
               improperpot: if measure==3 no U */

double bondpot( /*************************************************** bondpot */
  vector r0,vector r1,
  vector f0,vector f1,
  bondparm_t *parm)

/********************************************************************
  Bond potential of the form

                        2
  U = K (|r0-r1|-length)

  for bond    r0--r1   (length=0 allowed)

  #ifdef MORSE: expanded Morse potential (length>0 required)

  U =  D [ 1-exp(-a x)) ]^2   [original Morse]
    =  K x^2  * [1 - a x + 7/12*(a x)^2 ]
  where x=|r0-r1|-length, K = D a^2

  Generally, for potential U(x), x=|r0-r1|-length, it holds:
    rf = grad U(x) = U'(x)*(r0-r1)/|r0-r1|

  return values
    0 (if measure==0)
    U (if measure==1)

  if measure=2 then also phi:=bond length is returned
********************************************************************/

{
  vector dr;
  double r,rr,f,rf,U=0;

  VVV(dr,=r0,-r1)

  rr=SQR(dr);
  r=sqrt(rr);

  if (measure==2) phi=r;


#ifdef MORSE
  if (parm->a) {
    double x=r-parm->length;

    f=parm->Ki2*x/r*(1-parm->fa1*x+parm->fa2*Sqr(x));

    f0[0]+=rf=dr[0]*f; f1[0]-=rf;
    f0[1]+=rf=dr[1]*f; f1[1]-=rf;
    f0[2]+=rf=dr[2]*f; f1[2]-=rf;
    /*  VVO(f0,+=dr,*=f) VV(f1,-=dr) */

    if (measure) {
#  ifndef BLEND
      En.vir+=f*rr;
      CALCULATEPRESSURETENSOR(dr,f)
#  endif /*# BLEND */
      U=Sqr(x)*(parm->K-parm->Ka*x+parm->Ka2*Sqr(x));
      //      prt("%.9g %.9g %.9g MORSE r-r0 U f",x,U,f);
    }
  } else
#endif /*# MORSE */
  {
    /*  harmonic bond version
        f=(2*parm->K)*(parm->length/r-1);
        4/2012: optimized + parm->length=r=0 allowed */
    if (parm->length) f=parm->Ki2*(1-parm->length/r);
    else f=parm->Ki2;

    f0[0]+=rf=dr[0]*f; f1[0]-=rf;
    f0[1]+=rf=dr[1]*f; f1[1]-=rf;
    f0[2]+=rf=dr[2]*f; f1[2]-=rf;
    /*  VVO(f0,+=dr,*=f) VV(f1,-=dr) */

    if (measure) {
#ifndef BLEND
      En.vir+=f*rr;
      CALCULATEPRESSURETENSOR(dr,f)
#endif /*# BLEND */
      U=parm->K*Sqr(r-parm->length); }
  }

#ifdef DUMP
  prt("bond: %d-%d  r=%g u=%g  %.7f %.7f %.7f  %.7f %.7f %.7f",
      (vector*)r0-cfg[0]->rp,(vector*)r1-cfg[0]->rp,
      r,U,
      VARG(r0),VARG(r1));
#endif /*# DUMP */

  return U;
}


double anglepot( /************************************************ anglepot */
  vector r0,vector r1,vector r2,
  vector f0,vector f1,vector f2,
  angleparm_t *parm)

/********************************************************************

  angle>0: Bond angle potential of the form
                     2                       2
    U = K (phi-angle)  + Kub (|r0-r2|-length)

  angle<0: special form (close to phi=180 the same as above for
           angle=180, but not singular) (added 6/99)

    U = 2 K (1+cos phi)

  for angle       r1--r2
                  /
                 /  phi
                r0

  #ifdef UREY_BRADLEY (5/99): Kub and length are Urey-Bradley constants

  r0, r1, r2 = positions of atoms
  f0, f1, f2 = forces will be summed here
  parm = pointer to parameters
  return values
    0 (if measure==0)
    U (if measure==1)

  Notes:
    The virial (= trace of the pressure tensor) of the angle potential
      K (phi-angle)^2 is zero (the individual components are not).
    The virial of the Urey-Bradley term is not zero.
********************************************************************/

{
  vector a,b;
  vector grada,gradb;

  double ab,df,cf,norm,aa,bb,f;
  double U=0;

  VVV(a,=r0,-r1)
  VVV(b,=r2,-r1)
  norm=sqrt( (aa=SQR(a))*(bb=SQR(b)) );
  ab=SCAL(a,b);
  cf=ab/norm;

  /* to fix improbable numerical things like |cf|=1+1e-16 */
  if (cf<-1) cf=-1;
  if (cf>1) cf=1;

  if (measure==2) phi=acos(cf);
  if (parm->angle>0) {
    df=acos(cf)-parm->angle;
    norm=parm->K2/sqrt(1-Sqr(cf))*df/norm;
    if (measure) U=parm->K*Sqr(df); }
  else if (parm->angle<0) {
    /* angle==PI exception - marked by `negative angle' */
    norm=-parm->K2/norm;
    if (measure) U=parm->K2*(1+cf); }
  else {
    /* parm->angle=0 : needed in special cases only */
    norm=parm->K2/norm;
    if (measure) U=parm->K2*(1-cf); }

  f=norm*ab/aa; VVV(grada,=norm*b,-f*a)
  f=norm*ab/bb; VVV(gradb,=norm*a,-f*b)

  VV(f0,+=grada)
  VVV(f1,-=grada,+gradb)
  VV(f2,+=gradb)

#if PRESSURETENSOR&PT_VIR
  if (measure) {
    CALCULATEPRESSURETENSORG(grada,a)
    CALCULATEPRESSURETENSORG(gradb,b) }
#endif /*# PRESSURETENSOR&PT_VIR */

#ifdef UREY_BRADLEY
  if (parm->Kub) {
    VVV(a,=r0,-r2)

    aa=SQR(a);
    norm=sqrt(aa);

    /*.....if (measure==2) phi=r;*/

    f=parm->Kubi2*(1-parm->length/norm);
    f0[0]+=df=a[0]*f; f2[0]-=df;
    f0[1]+=df=a[1]*f; f2[1]-=df;
    f0[2]+=df=a[2]*f; f2[2]-=df;

    if (measure) {
      U+=parm->Kub*Sqr(norm-parm->length);
#  ifndef BLEND
      En.vir+=f*aa;
      CALCULATEPRESSURETENSOR(a,f)
#  endif /*# BLEND */
    }
  }
#endif /*# UREY_BRADLEY */

#ifdef DUMP
  prt("angle: %d-%d-%d  alpha=%g u=%g  %.7f %.7f %.7f  %.7f %.7f %.7f  %.7f %.7f %.7f",
      (vector*)r0-cfg[0]->rp,(vector*)r1-cfg[0]->rp,(vector*)r2-cfg[0]->rp,
      acos(cf)*(180/PI),U,
      VARG(r0),VARG(r1),VARG(r2));
#endif /*# DUMP */

  return U;
}

#ifdef DIHHIST
FILE *dihcp;
#endif /*# DIHHIST */

double dihedralpot( /***************************************** dihedralpot */
  vector r0,vector r1,vector r2,vector r3,
  vector f0,vector f1,vector f2,vector f3,
  torsionparm_t *parm)

/********************************************************************
  Torsion (dihedral) potential of the form

  n=0:           U=K[0]*(phi-K[1])^2       (K[1]=0 or PI only)
  n={1,2,3,4,6}: U=|K[0]|+K[0]*cos(n*phi)  [1 term]
  n=5:           K>0: U=K*cos^2(phi) for cos(phi)>0
                      U=0 for cos(phi)<0     (special for trans)
                 K<0: U=|K|*cos^2(phi) for cos(phi)<0
                      U=0 for cos(phi)>0     (special for cis)
                   |n|         i
  n<0:           U=SUM K[i] cos phi        [more terms]
                   i=0

  NOTE: U=K*(1+cos(n*phi-PI)) = |K'|+K'*cos(n*phi), where K'=-K

  for molecule      r1--r2
  (or part of       /    \
  molecule)       r0      r3

  phi is the dihedral angle

  return values
    0 (if measure==0)
    U (if measure==1)

  Note: virial of the dihedral potential is zero
********************************************************************/

{
  vector a,b,c,ap,cp;
  vector grada,gradb,gradc;

  double ab,bb,bc,cf,norm,apq,cpq;
  double abbb,bcbb;
  double U=0;
  double df,K;

  int i,n=parm->n;

  VVV(a,=r1,-r0)
  VVV(b,=r2,-r1)
  VVV(c,=r3,-r2)
  ab=SCAL(a,b);
  bc=SCAL(b,c);
  bb=SQR(b);
  abbb=ab/bb; VVV(ap,= -a,+abbb*b)
  bcbb=bc/bb; VVV(cp,=  c,-bcbb*b)
  /* vectors ap and cp are perpendicular to b */

  norm=1/sqrt( (apq=SQR(ap)) * (cpq=SQR(cp)) );

  /* cos phi */
  cf=SCAL(ap,cp)*norm;

  /* to fix improbable numerical things like |cf|=1+1e-16 */
  if (cf<-1) cf=-1;
  if (cf>1) cf=1;

#ifdef DIHHIST
  /* safe, but not fully under control: if (parm->dihhist) then grid!=0 ? */
  if (measure && parm->dihhist && parm->dihhist->grid) {
    vector perpab;
    double sinphi;

#  if 1
    /* NEW: the correct sign of phi. See also blend/dihedral.c */
    VECT(perpab,a,b) /* perpendicular to a,b */
    sinphi=SCAL(perpab,cp)/sqrt(SQR(perpab)*cpq);
    phi=atan2(sinphi,cf);
    phi=fmod(phi+2*PI,2*PI);
#  else /*# 1 */
    /* OLD: phi and -phi not distinguished */
    phi=acos(cf);
#  endif /*#!1 */


/*  prt("%s-%s-%s-%s %7.3f %d/%d",
        parm->dihhist->type[0],
	parm->dihhist->type[1],
	parm->dihhist->type[2],
	parm->dihhist->type[3],
	phi,
	(int)(phi*(parm->dihhist->grid/PI)),
	parm->dihhist->grid); */
    parm->dihhist->hist[(int)(phi*parm->dihhist->q)]++;
    if (phi>PI) phi-=2*PI;
    if (dihcp) fprintf(dihcp," %9.5f",phi*180/PI);
    if (fabs(phi)<parm->dihhist->gauche_trans) parm->dihhist->gauche++;
    else parm->dihhist->trans++; }
#else /*# DIHHIST */
  if (measure==2) phi=acos(cf);
#endif /*#!DIHHIST */

/* gradients: d cf/d r0 = -grada
              d cf/d r1 = grada-gradb
              d cf/d r2 = gradb-gradc
              d cf/d r3 = gradc       */
  apq=cf/apq; VVV(grada,= -norm*cp, +apq*ap)
  cpq=cf/cpq; VVV(gradc,=  norm*ap, -cpq*cp)
  VVV(gradb,= -abbb*grada, -bcbb*gradc)

  bb=Sqr(cf);
  K=parm->K[0];

  if (n==0) {
    /* n=0:           U=K[0]*(phi-K[1])^2       (K[1]=0 or PI only) */

    bb=sqrt(1-bb);
    if (bb==0) bb=1e-8; /* arcsin(1e-8)/1e-8 = 1-1e-16/6 */
    if (cf<0)
      if (parm->K[1]>0) df=-asin(bb);
      else df=acos(cf);
    else
      if (parm->K[1]>0) df=-acos(-cf);
      else df=asin(bb);

    if (measure) U=K*Sqr(df);
    bb=(-2*K)*df/bb; }

  else if (n>0) {

    if (measure) {

      /* U=cos(n*phi) first */
      switch (n) {
	case 1: U=cf; break;
	case 2: U=2*bb-1; break;
	case 3: U=(4*bb-3)*cf; break;
	case 4: U=8*(Sqr(bb)-bb)+1; break;
#ifdef BLEND
	case 5: if (K*cf<0) U=fabs(K)*bb; else U=0; goto pureU;
#endif /*# BLEND */
	case 6: U=((32*bb-48)*bb+18)*bb-1; break; }

      U=fabs(K)+K*U;
#ifdef BLEND
      pureU:;
#endif /*# BLEND */
      }

    switch (n) {
      case 1: bb=K; break;
      case 2: bb=(4*K)*cf; break;
      case 3: bb=(12*K)*(bb-0.25); break;
      case 4: bb=(32*K)*((bb-0.5)*cf); break;
#ifdef BLEND
      case 5: if (K*cf<0) bb=2*fabs(K)*cf; else bb=0; break;
#endif /*# BLEND */
      case 6: bb=(192*K)*( (Sqr(bb)-bb+0.1875)*cf ); break;
      default: ERROR(("dihedral: n=%d not supported",n)) }
    }

  else { /* n<0 */
    n=-n;
    if (measure) {
      U=parm->K[n];
      for (i=n-1; i>=0; i--) U=U*cf+parm->K[i]; }
    bb=n*parm->K[n];
    for (i=n-1; i>0; i--) bb=bb*cf+i*parm->K[i]; }

  VO(grada,*=bb) VO(gradb,*=bb) VO(gradc,*=bb)

  VV(f0,+=grada)
  VVV(f1,-=grada,-gradb)
  VVV(f2,-=gradb,-gradc)
  VV(f3,-=gradc)

#if PRESSURETENSOR&PT_VIR
  /* note that vector a=r1-r0 here, while a=r0-r1 in anglepot */
  if (measure) {
    CALCULATEPRESSURETENSORG(-grada,a)
    CALCULATEPRESSURETENSORG(-gradb,b)
    CALCULATEPRESSURETENSORG(-gradc,c) }
#endif /*# PRESSURETENSOR&PT_VIR */

#ifdef DUMP
  prt("dihedral: %d-%d-%d-%d  phi=%g u=%g  %.7f %.7f %.7f  %.7f %.7f %.7f  %.7f %.7f %.7f  %.7f %.7f %.7f",
      (vector*)r0-cfg[0]->rp,(vector*)r1-cfg[0]->rp,(vector*)r2-cfg[0]->rp,(vector*)r3-cfg[0]->rp,
      acos(cf)*(180/PI),U,
      VARG(r0),VARG(r1),VARG(r2),VARG(r3));
#endif /*# DUMP */

  return U;
} /* dihedralpot */


double improperpot( /****************************************** improperpot */
  vector r0,vector r1,vector r2,vector r3,
  vector f0,vector f1,vector f2,vector f3,
  torsionparm_t *parm)

/********************************************************************
  Improper torsion potential of the form

  n=0: U=K[0]*(phi-K[1])^2
  n!=0: not implemented

                            r2
  for molecule             /
  (or part of       r3---r0
  molecule)                \
                            r1

  Where phi is oriented angle of plane r1-r2-r3 from plane r0-r1-r2.
  In the above arrangement of atoms, phi is positive if r0 is below
  the plane r1-r2-r3; in other words, if vectors (r1-r0,r2-r0,r3-r0)
  are arranged as right-hand (x,y,z) Cartesian system.
  The chirality of r0, with respect to atoms r1,r2,r3 (in this order),
  is also positive.

  Good for bond angles |phi|<PI/2 (and not too close to PI/2).
  sin(phi) is used to calculate phi and hence phi and PI-phi are not
  distinguished.  This should not matter because the improper torsion potential
  is used for atoms 1, 2, and 3 bonded to 0 and they under normal conditions
  cannot bend so much.  The advantage of using sin(phi) is that -- unlike
  CHARMM guys -- we do not have problems for phi~0.  cos(phi) is calculated
  for testing purposes and error condition is set if cos(phi)<0.
  ... this has been fixed!

  If angle==0 then procedure torsion can be used instead (allowing any phi
  and a bit faster).

  return values
    0 (if measure==0)
    U (if measure==1)
  phi also returned if measure>=2
  no U calculated if measure==3

  Note: virial of the improper potential is zero
********************************************************************/

{
  vector v01,v12,v13, a,b,A,B, a0,a1,a2,b1,b2,b3;
  double s,t,q,sinphi,aa,bb,norm,ab,x,df;
  double U=0;

  VVV(v01,=r1,-r0)
  VVV(v12,=r2,-r1)
  VVV(v13,=r3,-r1)

  VECT(a,v01,v12) /* vector a is perpendicular to 012 plane */

  q=SQR(v12);
  s=SCAL(v12,v13)/q;
  VVV(b,=v13,-s*v12) /* vector b lies in 123 plane and is perpendicular to 12 */

  norm=sqrt( (aa=SQR(a))*(bb=SQR(b)) );
  ab=SCAL(a,b);
  sinphi=ab/norm;

  /* to fix improbable numerical things like |sinphi|=1+1e-16 */
  if (sinphi<-1) sinphi=-1;
  if (sinphi>1) sinphi=1;
  df=asin(sinphi);

  /* testing */
  VECT(B,v12,v13)
  if (SCAL(B,a)<0) {
    /* Error("|improper torsion angle|>PI/2"); */
    /* fix for |phi|>PI/2 (should not normally occur!) */
    if (df>0) df=PI-df;
    else df=-PI-df;
    norm=-norm; }

  if (measure>=2) {
    phi=df;
    if (measure==3) return 0; }

  df-=parm->K[1];

  norm=(-2*parm->K[0])/sqrt(1-Sqr(sinphi))*df/norm;
  x=norm*ab/aa; VVV(A,=norm*b,-x*a)
  x=norm*ab/bb; VVV(B,=norm*a,-x*b)

  VECT(a0,A,v12)
  VECT(a2,A,v01)
  VVV(a1,=a0,+a2) /* *-1 */

  t=SCAL(B,v12)/q;
  x=-2*s*t;
  VVVV(b2,=-s*B,+t*v13,+x*v12)

  VVV(b3,=B,-t*v12)

  VVV(b1,=b3,+b2) /* *-1 */

  if (measure) U=parm->K[0]*Sqr(df);

  VV(f0,+=a0)
  VVV(f1,-=a1,+b1)
  VVV(f2,+=a2,+b2)
  VV(f3,+=b3)

#if PRESSURETENSOR&PT_VIR
  if (measure) {
    CALCULATEPRESSURETENSORG(-a0,v01)
    CALCULATEPRESSURETENSORG(a2,v12)
    CALCULATEPRESSURETENSORG(b2,v12)
    CALCULATEPRESSURETENSORG(b3,v13) }
#endif /*# PRESSURETENSOR&PT_VIR   */

#ifdef DUMP
  prt("improper: %d-%d-%d-%d  phi=%g u=%g  %.7f %.7f %.7f  %.7f %.7f %.7f  %.7f %.7f %.7f  %.7f %.7f %.7f",
      (vector*)r0-cfg[0]->rp,(vector*)r1-cfg[0]->rp,(vector*)r2-cfg[0]->rp,(vector*)r3-cfg[0]->rp,
      acos(df)*(180/PI),U,
      VARG(r0),VARG(r1),VARG(r2),VARG(r3));
#endif /*# DUMP */

  return U;
}
