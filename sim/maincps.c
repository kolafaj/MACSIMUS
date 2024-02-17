/*
   directly #included from main.c:
   - called at the end of a cycle (by noint steps)
   - more measurements added
   - convergence profile recorded here
   - statistics of (some) variables added here
*/
#ifdef ANCHOR
loop (i,0,anchor.col) {
  anchor.rec[i].var=anchor.rec[i].sum/noint;
  anchor.rec[i].sum=0; }
#endif /*# ANCHOR */
      {
        struct ssd_s *ssd;
#ifdef FREEBC
        double Epot=En.pot*Eunit;
#else /*# FREEBC */
        double Epot=(En.pot+En.corr/box.V)*Eunit;
#endif /*#!FREEBC */
#ifdef RGYR
        double endtoend,Rgyr;
#endif /*# RGYR */
#ifdef WIDOM
        double Wid;
#endif /*# WIDOM */
#ifdef XSECTION
        double Xsec=0;
#endif /*# XSECTION */

#ifdef CLUSTERS
        analyzeclusters(icyc); /* WARNING - stride? */
#endif /*# CLUSTERS */
#ifdef XSECTION
        if (xs.freq<0 || (xs.freq>0 && icyc%xs.freq==0)) {
          if (!xs.grid) {
            xs.grid=1;
            WARNING(("xs.freq specified but xs.grid=0: xs.grid:=1.0 set")) }
          StaSet(DT*(xs.freq<0?1:xs.freq),lag.err,2,lag.n);
#  ifdef CLUSTERS
          if (xs.mode&16)
            Xsec=clXsection(xs.A);
          else
#  endif /*# CLUSTERS */
            if (xs.mode&8)
              Xsec=cXsection(cfg[0]);
            else
              Xsec=mXsection(cfg[0]);
          StaSet(DT,lag.err,2,lag.n); }
#endif  /*# XSECTION  */

        calculateSF();
        if (diff.mode&3) calculatediff(reread.frame,reread.to);

#ifdef FREEBC
        VO(box.center,=0)
#else /*# FREEBC */
        VV(box.center,=box.Lh)
#endif /*#!FREEBC */
#ifdef SLAB
        /* NB: center for CoM and angular momentum calculated here */
        measuredpr(slab.mode);
#endif /*# SLAB */

        /* remove CoM/momentum/angM if requested */
        if ((drift&DRIFT_WHEN)==DRIFT_CYCLE) removedriftssta(option('v')&64);

#ifdef RGYR
        Rgyr=measureplus(&endtoend,rg.end,rg.cp);
#else /*# RGYR */
        if (lag.M || lag.J)  measureplus();
#endif /*#!RGYR */

        if (lag.CM || lag.LM || lag.AM) measuredrifts();

        /* prior V3.6f, the following code:
           if (gear.order>2) { measureP(1); measureP(2);}
           was here, now moved to rhs */
        if (gear.order<2 && option('f')) {
          measureP(1);
          /*** Virtual volume change via <dU/dV> in the playback mode
               for NVT only: The kinetic part of pressure based on T is used
               NB: measureP(2) accepts twice the kinetic energy and adds the
               kinetic pressure correction. */
          if (constrd.mode & RESCALE_CM) En.kin_tr=No.f_tr*T;
          else En.kin=No.f*T;
          measureP(2); }

#ifdef POLAR
        if (scf.E) measureepspol();
#endif /*# POLAR */

        ssdistance();

#ifdef WIDOM
        if (widom.spreal>=0) Wid=XWidom(widom.spreal,widom.sp);
#  ifdef SLAB
        else Wid=Widom(widom.sp,widom.n,widom.z0,widom.z1,widom.dz,widom.mode);
#  else /*# SLAB */
        else Wid=Widom(widom.sp,widom.n,widom.corr);
#  endif /*#!SLAB */
#endif /*# WIDOM */

#ifdef DIHHIST
        StaSet(DT,lag.err,2,lag.n);
        gauchetrans(-2);
#endif /*# DIHHIST */

#ifdef HARMONICS
        measureharmonics();
#endif /*# HARMONICS */

#ifdef BJERRUM
#  include "mainbjerrum.c"
#endif /*# BJERRUM */


        /*** record convergence profiles ***/
        if (CPbuf) {
          float *CPrec=CPbuf+icyc*NCP;
          int count=0;

          memset(CPrec,0,NCP*sizeof(CPrec[0]));
#define recordCP(X) { if (count<NCP) CPrec[count++]=(float)(X); }
          /* The order of recordCP must match CPinfo in simmeas.c */
          /* The column (index in CPrec[] is by 1 less) */
          /*1*/ recordCP(En.tot)
          /*2*/ recordCP(En.T)
          /*3*/ recordCP(Epot)

          /* columns 4,5 of the convergence profile */
          if (tau.P==0 && tau.sig==0 && tau.rho==0) {
            if (dV==0) {
              if (thermostat==T_NOSE || thermostat>=T_NPT)
                recordCP(En.U*Eunit) /* col 4 */
              else
#ifdef LOG
                {
                  if (No.first) recordCP(En.pot0*Eunit) /* col 4 */
                  else recordCP(En.intra*Eunit) /* col 4 */
                }
#else /*# LOG */
                recordCP(rhounit*No.free_mass/box.V) /* col 4 */
#endif /*#!LOG */
              }
            else /* dV!=0 */
              recordCP(En.PdV*Punit) /* col 4 */
            recordCP(En.Pcfg*Punit) /* col 5, En.P prior 3.6l, see variable virial */ }
          else {
            if (tau.sig) recordCP(*sigvdWptr) /* col 4 */
            else recordCP(rhounit*No.free_mass/box.V); /* col 4 */
            if (dV==0) recordCP(En.Pcfg*Punit) /* col 5, En.P prior 3.6l, see variable virial */
            else recordCP(En.PdV*Punit) /* col 5 */ }

          /*6*/ recordCP(En.T_in)
          /*7*/ recordCP(En.T_tr)

          /* additional columns over "standard" 7 columns */
          /*8*/ /* used to be for COULOMB only, now also for FREE */
            if (Eext.isE<0) {
              if (No.ion) recordCP(En.J[Eext.isE+3])
#ifndef FREEBC
              else if (Q) recordCP(Q->M[Eext.isE+3])
#endif /*# FREEBC */
              else recordCP(En.M[Eext.isE+3]) }
            else {
              if (No.ion) recordCP(sqrt(SQR(En.J)))
#ifndef FREEBC
              else if (Q) recordCP(sqrt(SQR(Q->M)))
#endif /*# FREEBC */
              else recordCP(sqrt(SQR(En.M))) }
#ifdef SHEAR
          /*9*/ recordCP(En.Cv)
#endif /*# SHEAR */
#ifdef POLAR
          /*10*/ recordCP(En.self)
#endif /*# POLAR */
#ifdef DIHHIST
          /*11*/ recordCP(gauchetrans(dih.cp))
          if (dihcp) fprintf(dihcp,"\n");
#endif /*# DIHHIST */
#ifdef RGYR
          /*12*/ recordCP(Rgyr)
          /*13*/ recordCP(endtoend)
#endif /*# RGYR */
#ifdef XSECTION
          /*14*/ recordCP(Xsec) /* warning: need not be at the same frequency! */
#endif /*# XSECTION */
#ifdef WIDOM
          /*15*/ recordCP(Wid)
#endif /*# WIDOM */

          /* bug fixes:
             V3.6r: bad units forgotten from previous change (P real, En.P p.u.)
             V3.4l: cutoff correction was added twice; use P for NPT 
             V2.7f: was En.pot instead of En.U */
          En.Hnc=En.Unc+(tau.P?P/Punit:En.Pnc)*box.V;
          En.H=En.U+(tau.P?P/Punit:En.P)*box.V;

#if (PRESSURETENSOR&PT_ANY) == PT_ANY
          En.Hz=En.U+(tau.P?P/Punit:En.Ptens[2])*box.V;
#endif /*# (PRESSURETENSOR&PT_ANY) == PT_ANY */

          /* variables recorded in SIMNAME.cpi (incl. distances) */
          looplist (ssd,ssd0) {
            if (ssd->indx[0]==-1) recordCP(*ssd->u.q*Eunit)
            else if (ssd->indx[0]==-2) recordCP(*ssd->u.q*massunit*lengthunit/Sqr(timeunit)/chargeunit)
            else if (ssd->indx[0]==-3) recordCP(*ssd->u.q*Punit)
            else if (ssd->indx[0]==-9) recordCP(*ssd->u.q)
            else recordCP(ssd->u.dist); } }

        StaSet(DT,lag.err,2,lag.n);
        StaAdd("Tkin",En.T); issta=1;
        StaAdd("Tin",En.T_in);
        StaAdd("Ttr",En.T_tr);
#ifdef LOG
        StaAdd("Ein [J/mol]",En.intra*Eunit);
        if (No.first) StaAdd(string("Epot0 (<%d) [J/mol]",No.first),En.pot0*Eunit);
#endif /*# LOG */
        StaSet(DT,2,2,0);
        StaAdd("virial",En.vir);
#ifdef ECC
        StaAdd("ECC Pvir [Pa]",En.ECC_Pvir*Punit); /* "virial pressure" via Uel, former P (new in V3.4a) */
        StaAdd("ECC Pnel [Pa]",((En.kin*2*No.Pkinq+En.ECC_virnel)/DIM+En.corr/box.V)/box.V*Punit);
        StaSet(DT,lag.err,2,lag.n);
        StaAdd("ECC P [Pa]",En.P*Punit); /* ECC pressure (new in V3.4a) */
        StaAdd("ECC Pscaled [Pa]",En.ECC_Pscaled*Punit); /* conventional pressure (ef=const, new in V3.4a) */
        StaAdd("ECC Uscaled [J/mol]",Epot-En.ECC_U3*Eunit);
        StaAdd("ECC Epot [J/mol]",Epot);
#else /*# ECC */
        StaSet(DT,lag.err,2,lag.n);
        StaAdd("Epot [J/mol]",Epot); /* Epot is in J/mol */
        StaSet(DT,lag.err,2,lag.n);
        StaAdd("Pcfg [Pa]",En.Pcfg*Punit);
        StaAdd("P(vir) [Pa]",En.P*Punit);
#endif /*#!ECC */
        if (dV)
          StaAdd("Pcfg-P(dV) [Pa]",(En.Pcfg-En.PdV)*Punit);
        if (tau.P) {
          StaAdd("V",box.V);
          StaAdd("(Pcfg-P)*V",(En.Pcfg-P/Punit)*box.V); /* in K */
          if (rescale&RESCALE_PT) {
            if (rescale&RESCALE_X) StaAdd("Lx",box.L[0]);
            if (rescale&RESCALE_Y) StaAdd("Ly",box.L[1]);
            if (rescale&RESCALE_Z) StaAdd("Lz",box.L[2]); }
          StaSet(DT,2,2,0);
          StaAdd("Ecorr [J/mol]",En.corr/box.V*Eunit);
          StaAdd("Pcorr [Pa]",En.corr/Sqr(box.V)*Punit); }
        StaSet(DT,lag.err,2,lag.n);
        if (tau.sig)
          StaAdd("sigvdW",*sigvdWptr);

#ifdef SLAB
        if (wall.n) {
          if (abs(wall.n)&1) StaAdd("Pwall[0] [Pa]",En.Pwall[0]*Punit);
          if (abs(wall.n)&2) StaAdd("Pwall[1] [Pa]",En.Pwall[1]*Punit); }
#endif /*# SLAB */

#if PRESSURETENSOR&PT_OFF
        if (lag.Pt) {
          StaSet(DT,lag.Pt,2,0);
          StaAdd("Ptyz [Pa]",En.Ptens[3]*Punit);
          StaAdd("Ptzx [Pa]",En.Ptens[4]*Punit);
          StaAdd("Ptxy [Pa]",En.Ptens[5]*Punit);
          /* Daivis, Evans: JCP 100 541 (1994) */
          StaAdd("Ptxx-tr [Pa]",(En.Ptens[0]-En.trPt)*Punit);
          StaAdd("Ptyy-tr [Pa]",(En.Ptens[1]-En.trPt)*Punit);
          StaAdd("Ptzz-tr [Pa]",(En.Ptens[2]-En.trPt)*Punit); }
#endif /*# PRESSURETENSOR&PT_OFF */

        StaSet(DT,lag.err,2,lag.n);
        StaAdd("Etot",En.tot);
        StaAdd("internal energy [J/mol]",En.U*Eunit);
        StaAdd("enthalpy [J/mol]",En.H*Eunit);
        if ((thermostat==T_BERENDSEN || thermostat>=T_ANDERSEN) && tau.E==0)
          LRAdd("Etot",1,t,En.tot);
        if (lastEtot<2.999e33) {
          lastEtot -= En.tot;
          StaSet(DT,2,2,0);
          StaAdd("|dEtot|",fabs(lastEtot));
          StaAdd("dEtot^2",lastEtot*lastEtot); }
#if (PRESSURETENSOR&PT_ANY) == PT_ANY
        StaSet(DT,lag.err,2,lag.n);
        En.trPt=SUM(En.Ptens)/3+En.corr/(box.V*box.V);
        StaAdd("P(tens) [Pa]",En.trPt*Punit);
        StaSet(DT,2,2,0);
        StaAdd("P(tens)-P(vir) [Pa]",(En.trPt-En.P)*Punit);
        if (dV) StaAdd("P(tens)-P(dV) [Pa]",(En.trPt-En.PdV)*Punit);
#endif /*# (PRESSURETENSOR&PT_ANY) == PT_ANY */

        StaSet(DT,lag.err,2,lag.n);
        looplist (ssd,ssd0) if (ssd->stat) {
          if (ssd->indx[0]==-1) StaAdd(string("%s [J/mol]",ssd->name),*ssd->u.q*Eunit);
          else if (ssd->indx[0]==-2) StaAdd(string("%s [V/m]",ssd->name),*ssd->u.q*massunit*lengthunit/Sqr(timeunit)/chargeunit);
          else if (ssd->indx[0]==-3) StaAdd(string("%s [Pa]",ssd->name),*ssd->u.q*Punit);
          else if (ssd->indx[0]==-9) StaAdd(ssd->name,*ssd->u.q);
          else StaAdd(ssd->name,ssd->u.dist); }
      } /* convergence profile and statistics */
