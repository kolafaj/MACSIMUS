/*
  #included from constrd.c, Shake():
  thermostats PASS 2
  final rescaling in MTTK barostat/thermostat   P
*/

  if (tau.rho<0)
    /* @Eb scale[] set in trvpscale */
    box.V=rescalecfg(cfg[0],(rescale&RESCALE_XYZ)|RESCALE_L,0,scale);

  if (thermostat>=T_NPT) {
    vector oldlambda,ddlambda={0,0,0};
    real ldlogs=cfg[1]->logs;

    if (rescale&RESCALE_PT) {
      /* NPT: 3 general (dependent or independent) scalings */
      En.ext = No.fx*T*cfg[0]->logs;
      En.logs = cfg[0]->logs; /* at time t (before leap-frog) */

      /* Vpred[0] could be improved using newly calc. vel. */
      En.kinT=No.M_Th/hh*Sqr(Vpred[0]);
      En.kinP=No.M_Ph/hh*(Sqr(Vpred[1])+Sqr(Vpred[2])+Sqr(Vpred[3]));

      /* Verlet (=leap-frog here) for logs; NB: En.kin is doubled */
      cfg[1]->logs += hh*(En.kin + 2*En.kinP - No.fx*T)/No.M_T; // @A3
      cfg[0]->logs += cfg[1]->logs; // @B3, only good for the total Hamiltonian

      /* to be exported via SIMNAME.cpi with En.logs above */
      En.dlogs = (ldlogs+cfg[1]->logs)/(2*h); /* velocity Verlet estimate */
      VV(En.lambda,=cfg[0]->lambda)
      VVO(En.dlambda,=(Vpred+1),/h) /* predictor estimate */

      En.ext += En.kinT + En.kinP + box.V*No.P;

      /* NB: No.Pkinq here for consistency, however, should be 1 for NPT/MTK */
      /* En.Ptens is in Pa and without cutoff corrections */

      VV(oldlambda,=cfg[0]->lambda)
#if PRESSURETENSOR&PT_VIR
      loop (j,0,3) if (rescale&(1<<j))
        /* @93 & @C3; similarly below for other RESCALE options */
        ddlambda[j]=hh*(box.V*(En.Ptens[j]-No.P)+En.kin*No.invf+En.corr/box.V)/No.M_P - Vpred[0]*Vpred[j+1];
#endif /*# PRESSURETENSOR&PT_VIR */
      if (rescale&RESCALE_XisYisZ) {
        ddlambda[0]=ddlambda[1]=ddlambda[2]=(ddlambda[0]+ddlambda[1]+ddlambda[2])/3;
        /* not needed, but increases precision */
        cfg[1]->lambda[0]=cfg[1]->lambda[1]=cfg[1]->lambda[2]=(cfg[1]->lambda[0]+cfg[1]->lambda[1]+cfg[1]->lambda[2])/3; }
      else if (rescale&RESCALE_XisY) {
        ddlambda[0]=ddlambda[1]=(ddlambda[0]+ddlambda[1])/2;
        /* not needed, but increases precision */
        cfg[1]->lambda[0]=cfg[1]->lambda[1]=(cfg[1]->lambda[0]+cfg[1]->lambda[1])/2; }
      loop (j,0,3) if (rescale&(1<<j)) {
        cfg[1]->lambda[j] += ddlambda[j]; // @D3
        cfg[0]->lambda[j] += cfg[1]->lambda[j];
        scale[j]=exp(cfg[0]->lambda[j]-oldlambda[j]); }
      //      else cfg[1]->lambda[j]=0,scale[j]=1;

      // DO NOT USE (unfinished attempt of CM-based version)
      // box.V=rescalecfg(cfg[0],(rescale&(RESCALE_XYZ|RESCALE_CM))|RESCALE_L,0,scale); }

      // @E3 (see scale[j]= above); NB: also box.L rescaled
      box.V=rescalecfg(cfg[0],(rescale&RESCALE_XYZ)|RESCALE_L,0,scale); }

    else { /* NPT: isotropic scaling */
      En.ext = No.fx*T*cfg[0]->logs;
      En.logs = cfg[0]->logs; /* at time t (before leap-frog) */

      /* Vpred[0] could be improved using newly calc. vel. */
      En.kinT=No.M_Th/hh*Sqr(Vpred[0]);
      En.kinP=No.M_Ph/hh*Sqr(Vpred[1]);

      /* Verlet (=leap-frog here) for logs; NB: En.kin is doubled */
      cfg[1]->logs += hh*(En.kin + 6*En.kinP - No.fx*T)/No.M_T; // @A1
      cfg[0]->logs += cfg[1]->logs; /* @B1, only good for the total Hamiltonian */

      /* to be exported via SIMNAME.cpi with En.logs above */
      En.dlogs = (ldlogs+cfg[1]->logs)/(2*h); /* velocity Verlet estimate */
      En.lambda[0]=cfg[0]->lambda[0];
      En.dlambda[0]=Vpred[1]/h; /* predictor estimate */

      En.ext += En.kinT + En.kinP*3 + box.V*No.P;

      z=cfg[0]->lambda[0]; // remember lambda[0]

      // @C1
      cfg[1]->lambda[1]=
      cfg[1]->lambda[2]=
      cfg[1]->lambda[0] += hh*(box.V*(En.Pcfg-No.P)+En.kin*No.invf)/No.M_P - Vpred[0]*Vpred[1];

      // @D1
      cfg[0]->lambda[1]=
      cfg[0]->lambda[2]=
      cfg[0]->lambda[0] += cfg[1]->lambda[0];

      scale[2]=
      scale[1]=
      scale[0]=exp(cfg[0]->lambda[0]-z); /* z = old lambda[0] */

      // DO NOT USE (unfinished attempt of CM-based version)
      // box.V=rescalecfg(cfg[0],RESCALE_XYZ|rescale&RESCALE_CM|RESCALE_L,0,scale);

      /* @E1 (see scale[]= above); NB: also box.L rescaled */
      box.V=rescalecfg(cfg[0],RESCALE_XYZ|RESCALE_L,0,scale); } }

  else if (thermostat>=T_TR) {
    /* decoupled friction thermostats */
    z=ze=0;

    switch (thermostat) {
      case T_TR: /* friction thermostat - intermolecular only */
        ze=blogT(En.kin_tr/No.f_tr,T)*(h/tau.T);
        break;
      case T_IN:
        z=blogT(En.kin_in/No.f_in,T)*(h/tau.T);
        break;
      case T_FRICTIONS: {
        double T_tr=(No.f_in+No.f_tr)*T/(T_tr_in*No.f_in+No.f_tr);
        double T_in=(No.f_in+No.f_tr)*T/(No.f_in+No.f_tr/T_tr_in);

        z=blogT(En.kin_in/No.f_in,T_in)*(h/tau.T);
        ze=blogT(En.kin_tr/No.f_tr,T_tr)*(h/tau.T); } }

    ze-=z;
    loop (n,FROM,No.N) {
      si=spec[molec[n].sp]->si;
      r1=rof(&molec[n],cfg[1]->rp);
      ns=molec[n].ns;
      VO(v_mol,=0)
      m_mol=0;
      loop (i,0,ns) {
        mi=si[i].mass;
        m_mol+=mi;
        VV(v_mol,+=mi*r1[i]) }

      VO(v_mol,/=m_mol) /* h*velocity of CM of a molecule */

      loop (i,0,ns) VVV(r1[i], -=z*r1[i], +ze*v_mol) } }

  else if (thermostat==T_BERENDSEN) {
#ifdef SLAB
    if (slab.qT) {
      /* separate Berendsen (friction) thermostat for a slab */
      z=exp(blogT(En.T,T)*-h/tau.T);
      ze=exp(blogT(En.T,slab.T)*-h/tau.T);
      loop (n,FROM,No.N) {
        r1=rof(&molec[n],cfg[1]->rp);
        loop (i,0,molec[n].ns) {
          double Z=rof(&molec[n],cfg[0]->rp)[i][2]/box.L[2];
#  ifndef SLIT
          Z=fmod(Z,1);
#  endif /*# SLIT */
          if (Z>=slab.Tz0 && Z<slab.Tz1) VO(r1[i],*=ze)
          else VO(r1[i],*=z) } } }
    else
#endif /*# SLAB */
    {
      /* standard Berendsen (friction) thermostat */
      z=exp(blogT(En.T,T)*-h/tau.T); /* z=pow(En.T/T,-h/2/tau.T); */
      loop (n,FROM,No.N) {
        r1=rof(&molec[n],cfg[1]->rp);
        loop (i,0,molec[n].ns) VO(r1[i],*=z) } } }

  else if (thermostat==T_BUSSI) {
    double e=exp(-h/tau.T);
    double u1=rndgauss();
    double uf=0;

    loop (i,1,No.f) uf+=sqr(rndgauss()); /* inefficient */

    z=T/(No.f*En.T);

    z=sqrt(e+z*(1-e)*(u1*u1+uf)+2*e*sqrt(z*(1-e))*u1);

    loop (n,FROM,No.N) {
      r1=rof(&molec[n],cfg[1]->rp);
      loop (i,0,molec[n].ns) VO(r1[i],*=z) } }

  else if (thermostat==T_NOSE) {
    real ldlogs=cfg[1]->logs;

    /* Verlet (=leap-frog here) for logs */
    En.ext = cfg[0]->logs;
    En.logs = cfg[0]->logs; /* at time t (before leap-frog) */

    cfg[1]->logs += hh*(En.T/T-1)/Sqr(tau.T);
    cfg[0]->logs += cfg[1]->logs; /* only good for the total Hamiltonian */

    /* to be exported via SIMNAME.cpi with En.logs above */
    En.dlogs = (ldlogs+cfg[1]->logs)/(2*h); /* velocity Verlet estimate */

#if VERLET==0
    En.ext += Sqr(tau.T*ldlogs/h)/2;
#elif VERLET==1
    En.ext += Sqr(tau.T*(ldlogs+cfg[1]->logs)/h)/8;
#elif VERLET==2
    En.ext += Sqr(tau.T/h)*ldlogs*cfg[1]->logs/2;
#elif VERLET==3
    En.ext += Sqr(tau.T/h)*(Sqr(ldlogs)+Sqr(cfg[1]->logs))/4;
#elif VERLET==4
    En.ext += Sqr(tau.T*Vpred[0]/h)/2;
#elif VERLET==5
    En.ext += Sqr(tau.T/h)*(5./24*(Sqr(ldlogs)+Sqr(cfg[1]->logs))+1./12*ldlogs*cfg[1]->logs);
#elif VERLET==9
    En.ext += Sqr(tau.T*cfg[1]->logs/h)/2;
#endif /*#!VERLET==0!VERLET==1!VERLET==2!VERLET==3!VERLET==4!VERLET==5!VERLET==9 */
    En.ext *= No.f*T; }
