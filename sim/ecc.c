/* 
   since V3.4a, the ECC pressure is in En.P, but named as PECC in prints and CP
   for Pvir, see the paper
*/


#ifdef POLAR
#  error ECC + POLAR is nonsense
#endif /*# POLAR */

#ifdef GAUSSIANCHARGES
typedef double REAL;
#  include "gauss8.c"

static double extsigma;

static double toint(double r) /* ------------------------------------- toint */
{
  double efc=erfc(r/(sqrt(2)*extsigma));
  double ef=1-efc;
  double ee=exp(-Sqr(r)/(2*Sqr(extsigma)))/(sqrt(PI/2)*extsigma);

  return (ef+1)*efc/Sqr(r)-ee*(ee-2*ef/r);
}

/*
   Born energy of Gaussian(sigma) charge in an R-cavity
   see CGinsphere.mw
   without factor Q^2/(4 Pi eps0)*(er-1)/er
   inefficient but fool-proof
*/

static double UGCincavity(double sigma, double R) /* ----------- UGCincavity */
{
  double to=60*sigma;

  if (to<R) return 0;

  extsigma=sigma;

  /* - sign because of integration bounds swapped */
  return -(Gauss8(toint,to,R,160)+1/R)/2;
}
#endif /*# GAUSSIANCHARGES */

void ecc(int rescale) /************************************************* ecc */
/* EXPERIMENTAL,
   rescale=1: virtual volume change: calculation for rescaled box, see measureP
   rescale=-1: store original V
*/
{
  static double XA;

  if (rescale==-1) {
    /* virtual volume change: store factor 4*pi/3*SUM alpha_i for unscaled box */
    XA=(el.epsf-1)/(2+el.epsf)*box.V;
    return; }

  if (el.ecc) {
    double R=cbrt(3*box.V/(4*PI*No.N));
    double epsf; /* local epsf
                    =el.epsf if no virtual volume change
                    =newly calculated from XA for scaled box */
    double scale=1;

    /* in the notation of the Notes:
       En.ECC_U2 = En.el without self-energy = scaled pair elst energy
       En.ECC_U3 = self (Born/Lorentz) energy (does not depend on configuration
         except volume, but for simplicity treated here)
         independent noninteracting charges assumed (approximation) */

    if (rescale) epsf=(1+2*(XA/box.V))/(1-XA/box.V);
    else epsf=el.epsf;

    if (abs(el.ecc)==1) { /* ions */
      scale=el.epsf/epsf;

#ifdef GAUSSIANCHARGES
      /* self-energy U3, based on optimized code of ewald.c */
      {
        int i,n,ns,sp,nmol;
        molecule_t *mn;
        siteinfo_t *si;
        double e;

        En.ECC_U3=0;
        for (n=0;n<No.N;) {
          mn=molec+n;
          si=spec[mn->sp]->si;
          ns=mn->ns;

          sp=mn->sp;
          nmol=n;
          while (n<No.N && sp==molec[n].sp) n++; /* now n=next species */
          nmol=n-nmol; /* # of molecules of the same kind */
          loop (i,0,ns) {
            if ((e=si[i].charge)!=0)
              /* Sqr(e)*el.epsf = Q^2 */
              En.ECC_U3+=nmol*Sqr(e)*el.epsf * UGCincavity(sitedef[si[i].st].LJ[0].parm[SS_PARMS-1],R); }
        } /* n */
        En.ECC_U3*=(epsf-1)/epsf; }
#else /*# GAUSSIANCHARGES */
      /* point charges
         the formula with unscaled charges Q is:
           4 pi eps0 U3 = 1/2R * SUM Q^2 *(epsf-1)/epsf
         note that
           charges.molsumq = SUM q^2 = SUM Q^2/el.epsf
      */
      En.ECC_U3=-charges.molsumq*(epsf-1)/(2*R)*scale;
#endif /*#!GAUSSIANCHARGES */

      if (rescale) /* virtual volume change: rescale U2 */
        En.el *= scale;
      else {
        StaSet(0,lag.err,2,lag.n);
        StaAdd("ECC P2corr [Pa]",
               (2-epsf-Sqr(epsf))/epsf*En.el*((Punit/DIM)/box.V));
        StaAdd("ECC P2full [Pa]",
               (2-Sqr(epsf))/epsf*En.el*((Punit/DIM)/box.V));
        StaSet(0,2,2,0);
        StaAdd("ECC P3corr [Pa]",
               (1+2/epsf)*En.ECC_U3*((Punit/DIM)/box.V));
        En.ECC_Pcorr=
          ((2-epsf-Sqr(epsf))/epsf*En.el
           +(1+2/epsf)*En.ECC_U3)/(DIM*box.V);
        StaAdd("ECC P3full [Pa]",
               (2+2/epsf)*En.ECC_U3*((Punit/DIM)/box.V));
        En.ECC_P3ref=En.ECC_U3/(DIM*box.V);
        StaAdd("ECC P3ref [Pa]",En.ECC_P3ref*Punit);
        StaAdd("ECC U3 (incl. in Eel)",En.ECC_U3); } }

    else { /* dipoles, code for Gaussian charges not available */
      scale=el.epsf/Sqr((el.epsf+2)/3) / (epsf/(Sqr((epsf+2)/3)));
      En.ECC_U3=-charges.molsumdq/Cub(R)*(epsf-1)/(epsf+2)*scale;
      if (rescale) /* rescale U2 */
        En.el *= scale;
      else {
        StaSet(0,lag.err,2,lag.n);
        /* point-dipole-based formula (1/r^3) */
        StaAdd("ECC P2corr [Pa]",
               (Sqr(epsf)-3*epsf+2)/epsf*En.el*((Punit/DIM)/box.V));
        /* for finite-size dipoles, this is approximate */
        StaAdd("ECC P2full [Pa]",
               (Sqr(epsf)+2)/epsf*En.el*((Punit/DIM)/box.V));
        StaSet(0,2,2,0);
        StaAdd("ECC P3corr [Pa]",
               (Sqr(epsf)+2*epsf+2)/epsf*En.ECC_U3*((Punit/DIM)/box.V));
        En.ECC_Pcorr=
          ((Sqr(epsf)-3*epsf+2)/epsf*En.el+
           (Sqr(epsf)+2*epsf+2)/epsf*En.ECC_U3)/(DIM*box.V);
        StaAdd("ECC P3full [Pa]",
               (epsf+2)*(epsf+1)/epsf*En.ECC_U3*((Punit/DIM)/box.V));
        StaAdd("ECC P3ref [Pa]",
               En.ECC_P3ref=En.ECC_U3*((Punit/DIM)/box.V));
        StaAdd("ECC U3 (incl. in Eel)",En.ECC_U3); } }

    /* NOTE: the ECC self-energy (U3) is here included in the elst. energy, 
       it will be later used in  En.P as the virial */
    En.ECC_U2=En.el;
    En.el+=En.ECC_U3;
  }
}
