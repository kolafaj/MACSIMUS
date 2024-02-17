    /* external forces initialization and info */
    /* piece of code directly #included to main.c */


    /* harmonic central and similar (slab-keeping) forces */
    {
      int k;

      center.on=0;

      loop (k,0,DIM) {
        if (center.K[k]!=0) center.on |= 2;
        center.K2[k]=2*center.K[k];
        if (center.r0[k]!=0 && center.K[k]==0)
          WARNING(("center.r0[%d]=%g but center.K[k]=0 does not make sense\n\
*** (note that center.r0 was renamed to center.z0 from version 2.2g)",
                   k,center.r0[k])) }
      if (center.on&2) {
        if (SQR(center.r0))
          prt("atom-based force K=(%g %g %g) harmonic from r0=(%g %g %g)",
                                   VARG(center.K),             VARG(center.r0));
        else
          prt("atom-based harmonic force K=(%g %g %g) specified",
                                            VARG(center.K)); }

#  ifdef SLAB
#    ifdef FREEBC
#      error "combination SLAB+FREEBC is not recommended"
#    endif
#    if SLAB & 4
      if (slab.wall.epsrho) {
        center.on|=8;
        slab.wall.A=2*PI*slab.wall.epsrho*2./45*Cub(Cub(slab.wall.sig));
        slab.wall.B=2*PI*slab.wall.epsrho/3*Cub(slab.wall.sig);
        put2(slab.wall.A,slab.wall.B) }
#    endif /*# SLAB & 4 */
#    if SLAB & 2
      cleave.Ksigmah=cleave.K*cleave.sigma/2;
      if (!cleave.sigma) {
        cleave.n=0;
        prt("NOTE: cleave.n:=0 because cleave.sigma=0"); }

      if (cleave.n) {
        if (cleave.n<0 || cleave.n>2) {
          ERROR(("cleave.n=%d not in {0,1,2}",cleave.n))
          cleave.n=2; }

        center.on|=4;
        underline("cleaving walls");
        prt("cleave.sigma=%g AA  cleave.K=%g K/AA  cleave.n=%d  cleave.init=%d",
            cleave.sigma,cleave.K,cleave.n,cleave.init);

        if (cleave.i<0)
          prt("sites 0..%d (incl.) are susceptible to the cleaving potential",-cleave.i);
        else
          prt("site %d is susceptible to the cleaving potential",cleave.i);

        if (cleave.init&1) 
          prt("cleaving wall position(s) will be copied from data");

        if (cleave.init&2) {
          if (cleave.n==2) {
            /* barostat in both cleaved compartments */
            rescale |= RESCALE_CLEAVE;
            prt("In case of NPT in z-direction, the cleaving walls will be\n  %s", cleave.init&2?
                "moving according to pressures in compartments":
                "uniformly scaled with z-box"); } }
        else
          prt("only one cleaving wall, bit 2^1 of cleave.init ignored"); }
#    endif /*# SLAB & 2 */

      loop (k,0,NCENTER) {
        slab.nn[k]=slab.n[k];
        if (k) slab.nn[k]+=slab.nn[k-1];
        center.on |= slab.Kz[k]!=0;
        slab.Kz2[k]=2*slab.Kz[k];
        slab.dz[k]=slab.z1[k]-slab.z0[k]; }

      if (center.on&1)
        loop (k,0,NCENTER) if (slab.nn[k] && slab.Kz[k]) {
          if (slab.z1[k]) {
            prt_("slab step K=%g centered at z=L/2%+g, from z0=+-%g to z1=+-%g\n  for n in [%d,%d)",
                slab.Kz[k],slab.z[k],slab.z0[k],slab.z1[k],
                k==0?0:slab.nn[k-1],slab.nn[k]);
            slab.Kz2[k]=slab.Kz[k]/2*Sqr(slab.dz[k]);
            prt_("  step height=%g K  width=%g A",slab.Kz2[k],slab.dz[k]);
            if (thermostat) prt_("  Boltzmann factor=%g",exp(-slab.Kz2[k]/T));
            _n }
          else
            prt("z-force K=%g slabed at z=L/2%+g, from z0=+-%g for n in [%d,%d)",
                slab.Kz[k],slab.z[k],slab.z0[k],
                k==0?0:slab.nn[k-1],slab.nn[k]); }
#  endif /*# SLAB */
    }

    if (center.cmn<0) center.cmn=0;

    if (center.cmn==0)
      if (center.cmK[0]!=0 || center.cmK[1]!=0 || center.cmK[2]!=0)
        prt("Nonzero center.cmK specified with center.cmn=0 ==>\n\
  center.cmn=%d (=No.N=all molecules) assigned",center.cmn=No.N);

    center.cmmass=0;
    if (center.cmn) {
      int n,i;

      loop (n,FROM,center.cmn) {
        molecule_t *mn=molec+n;
        siteinfo_t *si=spec[mn->sp]->si;

        loop (i,0,mn->ns)
        center.cmmass+=si[i].mass; }

      prt("Whole-configuration harmonic force K*(r_CM-r0)^2 to center-of-mass specified\n\
  for %d molecules of total mass %g (prog.units):",center.cmn,center.cmmass);
      loop (i,0,DIM) if (center.cmK[i]!=0) {
        double tau=sqrt(center.cmmass/(2*center.cmK[i]));

        prt("K%c = %g  ==> tau = %g ps (oscillation period=%g ps)",i+'x',
            center.cmK[i],tau,2*PI*tau); } }

#ifndef FREEBC
    if (center.cmn || center.on)
      prt("WARNING: the above forces violate periodic b.c.\n");
#endif /*# FREEBC */

    /* electric field */

    /* conversion of V/m to internal prog. units */
    /*  NB: (chargeunit/forceunit)=1./352259391.) */
    VVO(Eext.E,=el.E,*(chargeunit/forceunit))
    Eext.f=el.f*timeunit; /* el.f is in Hz */
    VV(Eext.phase,=el.phase)
#ifdef POLAR
    Eext.Eepspol=scf.E*(chargeunit/forceunit);
#endif /*# POLAR */
    Eext.isE=el.E[0]!=0 || el.E[1]!=0 || el.E[2]!=0;
    if (Eext.E[0]!=0 && Eext.E[1]==0 && Eext.E[2]==0) Eext.isE=-3;
    if (Eext.E[1]!=0 && Eext.E[0]==0 && Eext.E[2]==0) Eext.isE=-2;
    if (Eext.E[2]!=0 && Eext.E[1]==0 && Eext.E[0]==0) Eext.isE=-1;
    _n
    if (Eext.isE) prt("Electrostatic field el.E = [ %g %g %g ] V/m,\n\
* WARNING: energy in ionic system may not be conserved",VARG(el.E));
    if (Eext.isE && lag.J==0) prt("\
* lag.J=0 so that the current density is not measured and conductivity\n\
  cannot be calculated\n\
* program unit of elst field (used in .cp) = %g V/m",1/(chargeunit/forceunit));
#ifdef COULOMB
    if (Eext.isE<0)
      prt(No.ion?"J%c will be monitored instead of |J|":"M%c will be monitored instead of |M|",Eext.isE+3+'x');
#endif /*# COULOMB */

    Eext.isB=el.B[0]!=0 || el.B[1]!=0 || el.B[2]!=0;
    if (Eext.isB) {
      if (gear.order>2) ERROR(("magnetic field not implemented for Gear, check -m"))
      prt("\
Magnetic field = [ %g %g %g ] T,\n\
               applies to ions and magnetic dipoles (if any)",VARG(el.B));
      /* static magnetic field transformed to p.u. and divided by h */ 
      VV(Eext.B,=(chargeunit*timeunit/massunit)*el.B)
      VVO(Eext.Bh,=Eext.B,/h) }

    Eext.ism=Eext.isB && (el.m.m!=0);
    if (Eext.ism) {
      if (el.m.plus==el.m.minus) ERROR(("el.m.plus==el.m.minus"))
      Eext.m=el.m.m/(chargeunit*sqr(lengthunit)/timeunit);
      prt("Magnetic dipoles sp=%d bond %d(+)--%d(-) m=%g A m2",
          el.m.sp,el.m.plus,el.m.minus,el.m.m); }

  
