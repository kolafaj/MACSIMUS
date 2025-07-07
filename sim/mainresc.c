/* rescaling stuff for main.c, periodic b.c., also saturation */

#ifdef FREEBC
      if (tau.rho || tau.P || dV) ERROR(("tau.rho, tau.P, dV not supported with FREEBC"))
#else /*# FREEBC */

      /* tau.rho<0 solved elsewhere */

      if (tau.rho>0) {
        /* rescale to reach given density or L */
        if (rhorescale(rescale,tau.rho,L,halftcyc)) {
          if (rhoreached) prt("^^^^^^^^^^^^^^^^ the reference box (density) has been reached ^^^^^^^^^^^^^^^^");
          rhoreached=0; } }

      /* Berendsen (friction) barostat or adjusting molecule size */

      if (rescale & RESCALE_PT) {
        /* rescaling based on pressure tensor components */
        if (thermostat>=0 && thermostat<T_NPT && tau.P!=0) {
#  if (PRESSURETENSOR&PT_ANY) == PT_ANY
          int ii;
          double beta_T; /* compressibility in 1 direction */

          if (No.bulkmodulus) beta_T=1/(3*No.bulkmodulus);
          else beta_T=box.V/(No.f*En.T); /* ideal gas */

          if (dV)
            ERROR(("implementation limitation:\n\
*** pressure tensor components calculated by virtual volume (area, shape)\n\
*** change (with dV set) cannot be used for shape change (rescale & %d)",RESCALE_PT))

          loop (ii,0,DIM) {
            double Ploc=En.Ptens[ii]+((corr&2)==2)*En.corr/Sqr(box.V);
#ifdef SLAB
            if (abs(wall.n)==3 && (constrd.mode&7)==4) Ploc=(En.Pwall[0]+En.Pwall[1])/2;
#endif
            Pscale[ii]=scaling(tau.P,noint,h*beta_T*(No.P-Ploc),maxscale); }
/*  Pscale=exp((tau.P+tau.R<0?noint:-1)*h*box.V/No.f/En.T/(tau.P+tau.R)*(No.P-locP)); */

          if (rescale & RESCALE_XisY) Pscale[0]=Pscale[1]=sqrt(Pscale[0]*Pscale[1]);
          if (rescale & RESCALE_XisYisZ) Pscale[0]=Pscale[1]=Pscale[2]=cbrt(PROD(Pscale));


#    ifdef SLAB
#      if SLAB & 1
          if (tau.L) {
            double Lx=slab.Lx[0]+En.T*(slab.Lx[1]+En.T*(slab.Lx[2]+En.T*slab.Lx[3]));
            double Ly=slab.Ly[0]+En.T*(slab.Ly[1]+En.T*(slab.Ly[2]+En.T*slab.Ly[3]));

            if (option('v')&4)
              prt("T=%.3f target Lx=%.5f Ly=%.5f, now Lx=%.5f Ly=%.5f",
                  En.T,Lx,Ly,box.L[0],box.L[1]);

            Pscale[0]=scaling(tau.L,noint, h*(box.L[0]-Lx), maxscale);
            Pscale[1]=scaling(tau.L,noint, h*(box.L[1]-Ly), maxscale); }
#      endif
#      if SLAB & 2
          if (rescale & RESCALE_CLEAVE)
            cleavescaling(Pscale,tau.P,noint,
                          h*beta_T*(En.fcleave[1]-En.fcleave[0])/(box.L[0]*box.L[1]),
                          maxscale);
#      endif /*# SLAB & 2 */
#    endif

          if (tau.P<0) /* rescale now (=every cycle), otherwise every step */
            rescalecfg(cfg[0],rescale|RESCALE_L,0,Pscale);
#  else /*# (PRESSURETENSOR&PT_ANY) == PT_ANY */
          ERROR(("\
rescale & RESCALE_PT (%d) was specified but this version\n\
*** of cook was not compiled with needed #define PRESSURETENSOR 3 (or more)",RESCALE_PT))
#  endif /*#!(PRESSURETENSOR&PT_ANY) == PT_ANY */
        } }
      else { /* !(rescale & RESCALE_PT) */
        /* isotropic (virial-based) rescaling */
        if (thermostat>=0 && thermostat<T_NPT && (tau.P!=0 || tau.sig!=0)) {
          double tauloc=tau.P;
          /* tau.P<0: rescaled every cycle
             tau.P>0: rescaled every step (keeps scaling factor through whole cycle)
             note: P and Pvir are in Pa */
          double Ploc;
          double beta_T; /* compressibility */

          if (No.bulkmodulus) beta_T=1/(3*No.bulkmodulus);
          else beta_T=box.V/(No.f*En.T);

          if (tau.sig) tauloc=tau.sig;

          if (dV) Ploc=corr&2?En.PdV.c:En.PdV.n;
          else Ploc=En.Pref; /* corr&2 already taken into account */

#ifdef SLAB
          if (abs(wall.n)==3 && (constrd.mode&7)==4)
            Ploc=(En.Pwall[0]+En.Pwall[1])/2;
#endif
          Pscale[0]=Pscale[1]=Pscale[2]=scaling(tauloc,noint,h*beta_T*(No.P-Ploc),maxscale);
          if (rescale & RESCALE_CLEAVE) ERROR(("cleaving rescaling incompatible with virial-based pressure"))

          if (tau.P<0) /* rescale now (otherwise every step) */
            rescalecfg(cfg[0],rescale|RESCALE_L,0,Pscale);

          if (tau.sig<0) { /* adjusting sigvdW (LJsigma) now */
            *sigvdWptr/=Pscale[0];
            setss(&sstab[tau.i][tau.j],tau.i,tau.j,LJcutoff,0);
            if (tau.i!=tau.j) sstab[tau.j][tau.i]=sstab[tau.i][tau.j]; } } }

#  ifdef COULOMB
/*** saturation control ***/
      if (tau.sat!=0) {
        double sat,q;

        if (!Q) ERROR(("internal: Q (Ewald charge structure) not allocated"))

        if (Eext.isE<0) {
          /* adjust external electrostatic field to saturation */
          sat=Q->M[Eext.isE+3]/el.sumM;
          q=exp((1-sat/el.sat)*No.dt/fabs(tau.sat));
          Min(q,1.01) Max(q,0.99)
          Eext.E[Eext.isE+3]*=q;
          if (option('v')&4) prt("%c: sat=%7.4f E=%.5e q=%7.4f SAT",'x'+Eext.isE+3,sat,Eext.E[Eext.isE+3]/(chargeunit/forceunit),q);
          StaAdd(string("E%c [V/m]",'x'+Eext.isE+3),
                 Eext.E[Eext.isE+3]/(chargeunit/forceunit)); }
        else if (Eext.isE>0)
          ERROR(("elst field perpendicular to one of axes required with tau.sat,el.sat"))
        else {
          /* adjust el.xinf=1/(2*el.epsinf+1) to saturation */
          sat=sqrt(SQR(Q->M))/el.sumM;
          q=log(sat/el.sat)*No.dt/fabs(tau.sat);
          Min(q,0.01) Max(q,-0.01)
          el.xinf=1/(2*el.epsinf+1)+q; /* el.epsinf is the primary stored quantity */
          el.epsinf=0.5/el.xinf-0.5;
          if (option('v')&4) prt("sat=%7.4f el.epsinf=%8.4f q=%9.6f SAT",sat,el.epsinf,q);
          StaAdd("1/(2*epsinf+1)",el.xinf); }
      }

#  endif /*# COULOMB */

#endif /*#!FREEBC */
