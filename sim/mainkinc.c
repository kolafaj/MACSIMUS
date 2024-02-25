/* info about the CoM sampling or the kinetic pressure correction */

    if (tau.CM) {
      underline("CoM thermostat");
      prt("\
The %s (tau=%g ps) thermostat for the center of mass is applied.\n\
Double check variables conserved, norm, corr&16, fixed sites, walls etc.\n",tau.CM<0?"Langevin":"Maxwell-Boltzmann",fabs(tau.CM));
      if (thermostat==T_ANDERSEN || thermostat==T_MAXWELL || thermostat==T_LANGEVIN || thermostat==T_LANGEVIN_CM) ERROR(("Center-of-mass thermostat is incompatible with stochastic thermostats"))
      if (thermostat==0) WARNING(("Center-of-mass thermostat and no particle thermostat does not make sense."))
      if (option('m')>2) WARNING(("Center-of-mass thermostat not implemented for Gear."))
      if (corr&4) ERROR(("Center-of-mass thermostat and kinetic pressure correction."))
      if ( (rescale&7) != 7 || (rescale&32) ) WARNING(("Center-of-mass thermostat may not be compatible with rescale=%d",rescale))
#ifdef SLAB
    if (SLAB&4) ERROR(("Center-of-mass thermostat for x,y-periodic system not implemented\n\
*** Replace vector operations in sim/thermocm.c by x,y-only versions."))
#endif
      if (No.free_s!=No.s) WARNING(("Center-of-mass thermostat and fixed sites."))
      if (conserved) WARNING(("Center-of-mass thermostat and conserved=%d.",conserved))
      if (drift&63) ERROR(("Center-of-mass thermostat and drift=%d affecting CoM are incompatible.",drift))
      if (drift) ERROR(("Center-of-mass thermostat and nonzero drift=%d, pls double check.",drift))

      /* no "kinetic pressure correction */
      isPkincorr=0;
      No.Pkinq=No.PkinqCM=1; }

    else {
      underline("kinetic pressure correction");
      isPkincorr=thermostat<=T_NOSE || (thermostat>=T_BUSSI && thermostat<T_DUMMY);
      prt("\
In 3D periodic boundary conditions and NVE or some thermostats, momenta are\n\
conserved. Hence, correction kB*T/V should be added to pressure, especially\n\
for gaseous systems. On the other hand, this is not suitable for the slab\n\
geometry (P_zz = saturated vapor pressure above a slab is more accurate without\n\
this correction.) With T derived from the kinetic temperature, this is\n\
equivalent to multiplying the kinetic part of pressure by factor 1+3/Nf, where\n\
Nf is the number of not-conserved degrees of freedom. See also variable corr.\n\
WARNING: applied to Ptens components in 3D, but does not allow separate\n\
         factors in 1D/2D symmetries.");
      if (corr&16) {
        /* automatic setup of the kinetic pressure correction */
        corr&=0x7fffffe3;
        if (isPkincorr) corr|=12;
        prt("\
>>> Automatic setup of the kinetic pressure correction because corr&16:\n\
>>> corr=%d was set to be consistent with the thermostat and barostat.\n\
>>> WARNING: may fail for fixed atoms, walls, etc. Check the factors below!\n\
",corr); }

      if (corr&4) {
        No.Pkinq=(double)No.f0/No.f;
        if (No.Pkinq!=1)
          prt("corr&4 set: correction factor in site-frame = %.8f",No.Pkinq); }
      else {
        No.Pkinq=1;
        if (No.f0!=No.f && thermostat<T_NPT)
          prt("corr&4 unset: no kinetic pressure correction in the site-frame"); }

      if (corr&8) {
        No.PkinqCM=(double)No.f0_tr/No.f_tr;
        if (No.PkinqCM!=1)
          prt("corr&8 set: correction factor in molecule-frame = %.8f",No.PkinqCM); }
      else {
        No.PkinqCM=1;
        if (No.f0!=No.f && rescale&8 && dV!=0)
          prt("corr&8 unset: no kinetic (id. gas) pressure correction in the molecule frame\n\
         (used in the virtual volume change method)"); } }
      
    underline("Correction factors");
    prt("\
No.Pkinq = %.9g site-based kinetic pressure correction factor,\n\
           for NVT and Berensen barostat; must be 1 for MTTK NPT",No.Pkinq);
    prt("\
No.PkinqCM = %.9g as above but molecule (center-of-mass) based,\n\
             for NVT and Berensen barostat with pressure calculated by the\n\
             virtual volume change method, not used in MTTK NPT",No.PkinqCM);
    prt("\
No.NPT = %.9g for the Hoover/MTK NPT method, factor appearing\n\
         in the equation for scaling factor Vfactor",No.NPT);
    prt("\
No.invf = %.9g for the Hoover/MTK NPT method, factor appearing\n\
          in the equation for d2lambda/dt2",No.invf);

    if (isPkincorr ^ (No.Pkinq!=1)) if ((corr&128)==0) WARNING (("\
The kinetic pressure correction seems to be inconsistent with the chosen\n\
*** thermostat and barostat. See the manual, Sec. Pressure tensor,\n\
*** variables corr, drift, thermostat.\n\
*** Note that the automatic setup (corr&16):\n\
*** - may fail in special cases as fixed atoms, modified drift etc.\n\
*** - is applied at start and may fail on thermostat/barostat change\n\
*** - the correction is likely not suitable in the slab geometry\n\
*** - the correction is important for a homogeneous gas or fluid\n\
*** This warning can be suppressed by corr|128."))
