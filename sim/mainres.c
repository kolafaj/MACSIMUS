#ifdef LINKCELL
    if (box.max14)
      prt("fixed max 1-4 distance used: box.max14=%g",box.max14);
    else {
      prt("automatic max 1-4 distance setup:");
      put3(box.over14,box.Max14,box.max14v)
      box.Max14=0; }
    //    box.max14v=0; NO, compromises the 14 predictor, (problem with -w0)
#endif /*# LINKCELL */

#ifdef XSECTION
    put(xs.maxsize)
#endif /*# XSECTION */

    underline("statistics");
    if (option('v')&2) StaPrintAll(NULL);
    if (option('v')&4) StaPrintAll("");
    else prt("HINT: use `staprt %s' for detailed time-series analysis",simils.simname);
if (option('m')<2) {
      if (option('f')) prt("WARNING: T, Ptr, enthalpy.. may be wrong in the reread mode");
      else prt("WARNING: forces and energy not calculated because -f0"); }
    prt("V = %g m3 (last value)",box.V*1e-30);
    _n

#ifdef WIDOM
#  include "widomres.c"
#endif /*# WIDOM */

#ifdef DIHHIST
    printdihhist();
#endif /*# DIHHIST */

  if (tau.P) {
    double dV=StaStdErr("V");
    double V=StaMean("V");

    prt("rho = mass/<V> = %.3f %.3f [kg/m3]",rhounit*No.free_mass/V,rhounit*No.free_mass/Sqr(V)*dV);
    prt("<L> = %.5f %.5f [AA] (cube%s)",cbrt(V),dV/V*cbrt(V),iscube()?"":"-equivalent"); }

#if 0 /* energy per molecule */
    prt("<Epot>/molecule  = %12.2f [J/mol] =%11.5f [kcal/mol] = %g [K]",
        StaMean("Epot [*")/No.N,
        StaMean("Epot [*")/No.N/kcal,
        StaMean("Epot [*")/No.N/Eunit);
    prt("<Eintra>/molecule= %12.2f [J/mol] =%11.5f [kcal/mol] = %g [K]",
        StaMean("Ein [*")/No.N,
        StaMean("Ein [*")/No.N/kcal,
        StaMean("Ein [*")/No.N/Eunit);
#endif /*# 0 */

#if 0 /* mean quadratic dipole moments */
    {
      int n;
      double D;
      char s[16];

      loop (n,0,nspec) if (spec[n]->N) {
        sprintf(s,"dipmom[%d]^2",n);
        if ( (D=StaMean(s)) ) {
          D=sqrt(D);
          prt("sqrt(%s) = %f = %.4f D",s,D,D*Debye); } } }
#endif /*# 0 */

#if (PRESSURETENSOR&PT_ANY) == PT_ANY
    { /*
         this actually duplicates Pvirerr, but is more fool-proof
         because is based on the final pressure tensor
      */
      double x=StaMean("Ptr [Pa]")-StaMean("Pevir [Pa]");

      prt("\n\
%sressure accuracy: tr(<Ptens>)/3+corr-<P>=%g Pa\n\
EXPLANATION: The standard virial pressure P is calculated assuming that\n\
  the electrostatic energy equals electrostatic virial whereas the pressure\n\
  tensor components Ptens are calculated directly from forces; both pressures\n\
  include standard homogeneous (nonelectrostatic) cutoff correction is applied.\n"
#  if COULOMB<0
"  The difference thus assesses the Ewald summation accuracy.\n"
#  else /*# COULOMB<0 */
"  The Ptens-based value is probably more accurate.\n"
"  Consider using rescale|=64 (rescale=\"XXX\" or \"XXXCM\") with the barostat.\n"
#  endif /*#!COULOMB<0 */
"  Both values may differ also while fixing some atoms (e.g., via option -j).\n"
          ,fabs(x)>el.Perr?"WARNING: Low p":"P",x);
      if (fabs(x)>el.Perr*10) WARNING(("Too low el. virial pressure vs. tensor accuracy (%g Pa)\n\
*** This is normal for Guassian charges or cutoff electrostatics.\n\
*** For Ewald with point charges, check el.alpha, el.kappa, el.grid!",x))
    }
#endif /*# (PRESSURETENSOR&PT_ANY) == PT_ANY */

#if PRESSURETENSOR&PT_MOM
    {
      char x;
      underline("kurtosis and normalized 6-th central moment of v_CM");
      prt_("defined so that Gauss distribution gives 0 (kurtosis-3)");
      header("DIR  kurtosis 6-th moment ");
      loop (x,'x','x'+DIM) {
        double E1=StaMean(string("PKin.%c%c [Pa]",x,x))/No.N;
        double E2=StaMean(string("PKin2.%c%c [Pa^2]",x,x))/No.N;
        double E3=StaMean(string("PKin3.%c%c [Pa^3]",x,x))/No.N;

        prt(" %c  %9.6f  %9.6f",x,E2/Sqr(E1)-3,E3/Cub(E1)-15*E2/Sqr(E1)+30); }
      header("");
    }
#endif /*# PRESSURETENSOR&PT_MOM */

#ifdef COULOMB
    {
      int k;
      double V;

      if (tau.P) V=StaMean("V");
      else V=box.V;

      if (tau.sat!=0) {
        double x;
        char *E;

        /* saturation: write back (possibly) changed el.E: convert to V/m */
        VVO(el.E,=Eext.E,/(chargeunit/forceunit))
        underline("saturation autoset");
        prt("final values after a sweep:\n\
  el.epsinf = %g\n\
  el.E = [ %g %g %g ] V/m",el.epsinf,VARG(el.E));
        prt("NB: to restart simulation from this POINT, use in the data:");
        if (Eext.isE)
          prt("  el.E[%d]=%g",3+Eext.isE,el.E[3+Eext.isE]);
        else
          prt("  el.epsinf=%g",el.epsinf);
        prt("averaged values with standard errors:");
        if (Eext.isE) {
          E=string("E%c [V/m]",'x'+3+Eext.isE);
          x=StaMean(E);
          prt("  <%s> = %g %g",E,x,StaStdErr(E));
          if (tau.sat<0) {
            if (3+Eext.isE<0 || 3+Eext.isE>=DIM) ERROR(("Eext.isE=%g out of range",Eext.isE))
            el.E[3+Eext.isE]=x; } }
        else {
          double dx=StaStdErr("1/(2*epsinf+1)");

          x=StaMean("1/(2*epsinf+1)");
          E="epsinf";
          prt("  <1/(2*epsinf+1)> =  %g %g",x,dx);
          prt("  el.epsinf = %g  67\%: [%g %g]",
              (1/x-1)/2,(1/(x+dx)-1)/2,(1/(x-dx)-1)/2);
          if (tau.sat<0) el.epsinf=(1/x-1)/2; }
        if (tau.sat<0) {
          init=1;
          tau.sat=0;
          prt("\n\
====== You specified AUTOMATIC START by tau.sat<0. Therefore: ======\n\
* %s was assigned to the above AVERAGED value\n\
* init=1 (init=\"append\") was assigned to be used for the next sweep\n\
         (you may override init at next sweep start)\n\
? BUG: init=2 cannot be used with tau.sat<0 directly\n\
* tau.sat=0 was assigned\n\
NB: to restart simulation from this AVERAGE, use the data:",E);
          if (Eext.isE)
            prt("  el.E[%d]=%g",3+Eext.isE,x);
          else
            prt("  el.epsinf=%g",el.epsinf); }
      } /* tau.sat!=0 */

      /*** dielectric constant ***/
      if (No.ion && tau.sat==0) {
#  ifdef POLAR
        double a=No.A*(4*PI/3/V);
        prt("\
box polarizability No.A = %.8g\n\
polarizability density No.A/V = %.8f\n\
a = No.A*(4*PI/3/V) = (epsf-1)/(epsf+2) = %.8f\n\
Clausius-Mossotti high-frequency epsf = (1+2*a)/(1-a) = %.8f",
            No.A,No.A/V,a,(1+2*a)/(1-a));
#  endif /*# POLAR */
        prt("%d ions detected, dielectric constant calculation suppressed",No.ion);
      }
      else
        prtdielconst(V,Q!=NULL); /* new: also if elst. field */

      prt("\nEXTERNAL ELECTRIC FIELD = [ %g %g %g ] V/m", VARG(el.E));
      if (Eext.isE)
        prt("%d ions detected, so I will calculate %s",No.ion,
            No.ion?"conductivity from current":"dielectric constant from dipole moment");

      loop (k,0,DIM) if (el.E[k]) {
        char xyz=k+'x';

        if (No.ion) {
          char *name=string("J%c [A/m2]",xyz);
          double J=StaMean(name),dJ;
          int n;

          if (StaError) {
            prt("missing lag.J?");
            continue; }
          dJ=StaStdErr(name);
          prt("conductivity_%c = %g %g  S/m (mean stderr)",
              xyz,J/el.E[k],dJ/el.E[k]);
          if (lag.J>0)
            loop (n,0,nspec)
              prt("transference number [%d] = %.6f",n,StaMean(string("J%c[%d] [A/m2]",k+'x',n))/J); }
        else {
          char *name=string("M%c",xyz);
          double M=StaMean(name),dM;
          double eps[5],deps;
          int j;

          if (StaError) {
            name=string("M[%d]",k);
            M=StaMean(name); }
          if (StaError) {
            name=string("xM%c",xyz);
            M=StaMean(name); }
          
          if (StaError) continue;
          deps=M/Eext.sumM;
          prt("saturation <M%c>/M_tot = %g (M_tot=%g p.u.)",xyz,deps,Eext.sumM);
          if (fabs(deps)>0.1) prt("WARNING: The dipole moment of the cell is %.3f of the saturated one.\n\
*** The dielectric constant may be affected by nonlinear response.\n\
*** Consider decreasing the field and/or el.epsinf.",deps);

          dM=StaStdErr(name);
          M*=4.4266218e+09/(V*el.E[k]);
          dM*=4.4266218e+09/(V*el.E[k]);
          prt("chi_%c = %9.7f %9.7f",xyz,M,dM);
          loopto (j,-2,2) eps[j+2]=1+1/(1/(M+j*dM)-1/(el.epsinf*2+1));
          deps=Sqr((eps[2]-1)/M)*dM;
          prt("epsilon_%c = %9.7f %9.7f mean stderr",xyz,eps[2],deps);
          prt("68%% range: [%7.4f %7.4f]  95.45%% range: [%7.4f %7.4f]",
              eps[1],eps[3],eps[0],eps[2]); } }

      if (lag.J) {
        underline("EMD conductivity hints");
        prt("(applies only in directions where el.E_x,y,z=0)\n\
* calculate the time autocovariances of el. current density by (bash):\n\
    for f in x y z ; do staprt -t2 -nJ$f'[A/m2]' %s ; done\n\
* check that the data have sufficiently fine grid (otherwise decrease noint)\n\
* check that the data have sufficiently long lag (otherwise increase lag.J)\n\
* calculate their running integrals by commands:\n\
    for f in x y z ; do runint -e -m2 < %s.J${f}Am2.cov > %s.J$f.run ; done\n\
* calculate the average by:\n\
  mergetab %s.JxAm2.run:1:2 %s.JyAm2.run:1=1:2 %s.JzAm2.run:1=1:2 | tabproc A '(B+C+D)/3' > %s.Jsum.run\n\
* plot these curves and determine the limiting averaged value:\n\
    plot %s.J?.run\n\
* the conductivity in S/m is this value multiplied by V/kT=%g",
            simils.simname,
            simils.simname,simils.simname,
            simils.simname,simils.simname,simils.simname, simils.simname,
            simils.simname,
            V*Cub(lengthunit)/energyunit/StaMean("Tkin"));
      }
      prt("\nEXTERNAL MAGNETIC FIELD = [ %g %g %g ] T", VARG(el.B));
    }
#  if (COULOMB<=-1) && (PRESSURETENSOR&PT_VIR) && defined(SLAB)
    {
      int k;
      double V;

      if (tau.P) V=StaMean("V");
      else V=box.V;

      header("Yeh-Berkowitz correction of pressure tensor components");
      loop (k,0,3) {
        char xyz=k+'x';
        char *name=string("EwM%c",xyz);
        double MM=StaVar(name);

        if (StaError) {
          name=string("M%c",xyz);
          MM=StaVar(name); }
        if (StaError) { MM=0; continue; }

        En.YB[k]=2*PI*MM/Sqr(V)*Punit;
        prt("Delta Ptens.%c%c = %g Pa (%s)",
            xyz,xyz,
            En.YB[k],
            el.corr&1<<k
             ? el.corr&8 ? "added, forces were not corrected":"fully included incl. force correction"
             : "not included"); }
      header("");
    }
#  endif /*# (COULOMB<=-1) && (PRESSURETENSOR&PT_VIR) && defined(SLAB) */
#endif /*# COULOMB */

#if PRESSURETENSOR&PT_OFF
    if (lag.Pt) {
      double V; /* unnecessarily repeated! */

      if (tau.P) V=StaMean("V");
      else V=box.V;

      /* NB: ?\?- used not to activate trigraph ??- = ~ */

      underline("EMD viscosity hints");
#  ifdef SHEAR
      if (shear) prt("WARNING: shear is set, EMD invalid or inaccurate except Ptxy")
#  endif /*# SHEAR */
      prt("\
* calculate the time autocovariances of pressure tensor components (in Pa^2):\n\
    staprt -t2 -nPtxy* %s\n\
    staprt -t2 -nPtyz* %s\n\
    staprt -t2 -nPtzx* %s\n\
    staprt -t2 -nPtxx-tr* %s\n\
    staprt -t2 -nPtyy-tr* %s\n\
    staprt -t2 -nPtzz-tr* %s\n\
* check that the data have sufficiently fine grid exp. near t=0 (otherwise\n\
  decrease noint and long enough lag (otherwise increase lag.Pt):\n\
    plot :A:B:2P- %s.Pt?\?.cov \"@:A:B*.75:3P-\" %s.Pt?\?-tr.cov\n\
* calculate their running integrals:\n\
    runint -o -m2 < %s.Ptxy.cov > %s.Ptxy.run\n\
    runint -o -m2 < %s.Ptyz.cov > %s.Ptyz.run\n\
    runint -o -m2 < %s.Ptzx.cov > %s.Ptzx.run\n\
    runint -o -m2 -q.75 < %s.Ptxx-tr.cov > %s.Ptxx-tr.run\n\
    runint -o -m2 -q.75 < %s.Ptyy-tr.cov > %s.Ptyy-tr.run\n\
    runint -o -m2 -q.75 < %s.Ptzz-tr.cov > %s.Ptzz-tr.run\n\
* plot these curves and determine the limiting averaged value:\n\
    plot :A:B:2P- %s.Pt?\?.run \"@:A:B:3P-\" %s.Pt?\?-tr.run\n\
* NOTE THAT:\n\
  - the curves are correlated\n\
  - the cyan traceless diagonal components are less accurate (rel.weight=2/3)\n\
* in case the limit has not been reached, the long-t part of the curves may\n\
  be approximated by the hydrodynamic tail, a+b/sqrt(t)\n\
* the total value is\n\
     eta=(eta_xy+eta_yz+eta_zx)/5+(eta_xx+eta_yy+eta_zz)/7.5\n\
  where eta_xy is based on Ptxy, eta_xx on Ptxx-tr, etc.\n\
* the viscosity (in Pa.s, for dimensionless eta) is eta*V/kT=eta*%.7g\n\
* the diffusivity correction (in m2/s) in this simulation is +%.7g/eta",
          simils.simname,simils.simname,simils.simname,
          simils.simname,simils.simname,simils.simname,
          simils.simname,simils.simname,
          simils.simname,simils.simname,simils.simname,
          simils.simname,simils.simname,simils.simname,
          simils.simname,simils.simname,simils.simname,
          simils.simname,simils.simname,simils.simname,
          simils.simname,simils.simname,
          V*Cub(lengthunit)/energyunit/StaMean("Tkin"),
          2.837297*energyunit/6/PI*sqr(StaMean("Tkin"))*energyunit/(pow(V,4./3)*Pow4(lengthunit))); }
#endif /*# PRESSURETENSOR&PT_OFF */

#ifdef SHEAR
    /*** shear viscosity ***/
    prtviscosity();
#endif /*# SHEAR */

    if (LRRes("Etot","N")>1) { LRPrint(' '); _n }

#ifdef ANCHOR
    if (anchor.f) {
      fclose(anchor.f);
      anchor.f=NULL; }
#endif /*# ANCHOR */

    sortmolecules(sort);
    sort=0;

    if (option('w')>1) /* obsolete, deprecated */
      waitfordiskspace(No.s
#ifdef POLAR
                           *2
#endif /*# POLAR */
                             /7 + no*NCP/200 + 256*(1+!!option('l'))
#ifdef SLAB
                                                                      +128
#endif /*# SLAB */
                                                                          );

     /* closing sta,plb here so that in case of a removed molecule, 
        the plb-file still contains the molecule being removed */
     if (icyc) {
      backup("sta");
      StaKeySave(Fn("sta"),stoptime);
      if (option('y')) {
        if (dt.plb==0) writeplayback();
        closeplayback();
        plbopened=0; }
      saveCP(icyc,stoptime); }

    /*** save cfg, rdf ***/
    if (option('w')) {
      if (tau.sig) sigvdW=*sigvdWptr; /* otherwise sigvdW=0 already set */

      /* remove a molecule */
      if (removemol.n>=0) {
        /* removedrifts(1), depend_r(cfg[0],0) moved to remove1mol */
        remove1mol(removemol.n);
        savecfg(-2,(int4)stoptime,&sigvdW); }
      else
        /* removedrifts(1), depend_r(cfg[0],0) in savecfg */
        savecfg(-1,(int4)stoptime,&sigvdW); }

    writeasc(); /* according to option -l */

#ifdef DIHHIST
    if (dihcp) fclose(dihcp);
#endif /*# DIHHIST */

    if (icyc || option('w')==1) {
#ifdef DIHHIST
      savedih(dih.grid);
#endif /*# DIHHIST */
#ifdef HARMONICS
      saveharmonics();
#endif /*# HARMONICS */
#ifdef SLAB
      /* density profiles and surface tension */
      savedpr(slab.grid);
      printdpr(slab.mode,slab.prt,dV,constrd.PdVname);
      if (slab.ext.span) {
        if (slab.sp>=nspec) WARNING(("stacking correction (slab.ext.span set) requires normally\n\
*** autocenter (slab.sp). The results may be wrong."))
        extenddpr(slab.mode,slab.prt,dV,constrd.PdVname,
                  slab.ext.zero,slab.ext.center,slab.ext.span); }
#  if SLAB & 2
      if (cleave.n>0)
        prt("cleave.z0=%.12f (specified)  cleave.z[0]=%.12f (final)",cleave.z0,cleave.z[0]);
      if (cleave.n>1)
        prt("cleave.z1=%.12f (specified)  cleave.z[1]=%.12f (final)",cleave.z1,cleave.z[1]);
#  endif /*# SLAB & 2 */
#endif /*# SLAB */
#if PARALLEL==1
      parallelizerdf(0);
#endif /*# PARALLEL==1 */
      saverdf(rdf.grid); }

    /* original place of closing sta,plb,cp */

#ifdef POLAR
#  include "polarres.c"
#endif /*# POLAR */
    if (nm.dr) {
      if (nm.method&3) normalmodesc(); /* with constraints */
      else normalmodes(); /* flexible */ }

    /* release the heap */
    StaFree();
    release(mark);
#ifdef WIDOM
    widomrdf(-1);
#endif /*# WIDOM */

    if (No.c) {
      int sp;

      loop (sp,0,nspec)
        if (constrit[sp].nit[COR_LAST])
          prt("spec=%d: omegac=%f (%s) last niter/molec=%f\n(cf. \"constrit%d spec%d\")",
            sp,constrit[sp].omegac,
            omegac<0?"optimized":"fixed",
            (double)constrit[sp].nit[COR_LAST]/spec[sp]->N,sp,0); }

    prt("sweep finished: L=[%.9f %.9f %.9f] V=%.11g\n\
En.pot=%.9g En.kin=%.9g En.tot=%.12g",
        VARG(box.L),
        box.L[0]*box.L[1]*box.L[2],
        En.pot,En.kin,En.tot);

    fflush(out);

    if (En.kin && option('v')&4) put2(En.bonded,En.bonded/En.kin)

    if (T) if (En.bonded>En.kin*10) {
      static int pass=0;
      if (!pass)
        WARNING(("Bonded energy is %g times greater than kinetic energy\n\
(more warnings will be suppressed)",En.bonded/En.kin))
      pass=1; }

    printfSF();
    if (MSD.mode&3) printdiff();
#ifdef CLUSTERS
    printclusters();
#endif /*# CLUSTERS */
