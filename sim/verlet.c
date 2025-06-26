/*
  Verlet step (for 1 molecule) with posible thermostats:
  - atom-based Andersen/Maxwell thermostat
  - CM-based Andersen/Maxwell thermostat
  - Nose with predicted velocity
  #included from constrd.c (before each version of SHAKE)
  NB: friction (Berendsen) thermostats are out of the VERLET/SHAKE loop
  input:  p[i]     = force on atom i at time t
          r[i]     = position of atom i at time t
          r1[i]    = r[i](t)-r[i](t-h) = h*v(t-h/2)
          v[i]     = predicted velocity_i at time t = rof(mn,vpred)
          Vpred[0] = h*d logs(t)/dt, predicted
          vpred[i] = h*velocity of atom i, predicted
          Vfactor  = v_f (diagonal tensor or scalar)
  output: p[i] = r[i](t+1) (before SHAKE)
          r1 is not changed

  Andersen/Maxwell thermostat for models with constraints is not recommended!
  (previous version of verlet.c without Nose = shakemax.c)
*/

  /* Nose-Hoover thermostat with predicted velocities, MTK barostat */
  if (thermostat==T_NOSE || thermostat>=T_NPT || tau.rho<0) {
    if (thermostat>=T_NPT && rescale&RESCALE_PT) {
      /* separate x,y,z */
      loop (i,0,ns) { // @6
        VV(p[i],=hh*si[i].imass*p[i])
        p[i][0]-=(Vfactor+Vpred[1])*v[i][0];
        p[i][1]-=(Vfactor+Vpred[2])*v[i][1];
        p[i][2]-=(Vfactor+Vpred[3])*v[i][2];
        VVV(p[i],+=r[i],+r1[i]) } }
    else {
      /* plain Nose or isotropic MTK thermostat+barostat */
      loop (i,0,ns) { // @6
        VVV(p[i],=hh*si[i].imass*p[i],-Vfactor*v[i])
        VVV(p[i],+=r[i],+r1[i]) } } }

  /* Andersen/Maxwell thermostat */
  else if (prob) {
    if (CMbasedthermostat) {
      if (rnd()<prob) { /* center-of-mass-based thermostat */
        vector cmv; /* velocity of the center-of-mass */
        double im=1./spec[mn->sp]->mass;
        double sigma=sqrt(Thh*im);

        VO(cmv,=0)
#ifdef SLAB
        if (slab.qT) {
          /* z-modulated thermostat */
          double cmz=0; /* center-of-mass z */

          loop (i,0,ns) {
            VV(cmv,-=si[i].mass*r1[i])
            cmz+=si[i].mass*r[i][2]; }
          VO(cmv,*=im)
          cmz*=im/box.L[2];
#  ifndef SLIT
          cmz=fmod(cmz,1);
#  endif  /*# SLIT */
          if (cmz>=slab.Tz0 && cmz<slab.Tz1) sigma*=slab.qT;
          VO(cmv,+=rndgauss()*sigma)
          loop (i,0,ns) {
            VO(p[i],*=hhh*si[i].imass)
            VVVV(p[i],+=r[i],+r1[i],+cmv) } }
        else
#endif /*# SLAB */
         {
          loop (i,0,ns) VV(cmv,-=si[i].mass*r1[i])
          VVO(cmv,=im*cmv,+rndgauss()*sigma)
          loop (i,0,ns) {
            VO(p[i],*=hhh*si[i].imass)
            VVVV(p[i],+=r[i],+r1[i],+cmv) } } }
      else loop (i,0,ns) { /* no thermostat now */
        VO(p[i],*=hh*si[i].imass)
        VVV(p[i],+=r[i],+r1[i]) } }

    else loop (i,0,ns) /* atom-based thermostat */
      if (rnd()<prob) {
        double sigma=sqrt(Thh*si[i].imass);

#ifdef SLAB
        if (slab.qT) {
          /* z-modulated thermostat */
          double z=r[i][2]/box.L[2];
#  ifndef SLIT
          z=fmod(z,1);
#  endif  /*# SLIT */
          if (z>=slab.Tz0 && z<slab.Tz1) sigma*=slab.qT;
          VO(p[i],*=hhh*si[i].imass)
            VVO(p[i],+=r[i],+rndgauss()*sigma)
        } else
#endif /*# SLAB */
      {
        VO(p[i],*=hhh*si[i].imass)
        VVO(p[i],+=r[i],+rndgauss()*sigma) } }

      else { /* no thermostat now */
        VO(p[i],*=hh*si[i].imass)
        VVV(p[i],+=r[i],+r1[i]) } }

  else /* no thermostat */
    loop (i,0,ns) {
      VO(p[i],*=hh*si[i].imass)
      VVV(p[i],+=r[i],+r1[i]) }
