/* 
  printing results of Widom insertion particle calculations
  directly #included from main.c
*/

#ifdef SLAB
    if (widom.n) {
      double z,zi,w;
      int m;
      int nz=(widom.z1-widom.z0)/widom.dz+0.9999,iz;
      FILE *f=fopen(Fn("wid"),"wt");

      prt(">>> Widom (density profile) results written to file %s",lastFn);

      for (m=2; m<=8; m<<=1) if (widom.mode&m) {

        fprintf(f,"\n# WIDOM mode=%d (%s)\n\
# molecule=%d %d insertions/slab:\n\
#  z/AA      prob            mu/K     [prob_sym     mu_sym]\n",
                m,m==2?"total mu":m==4?"w/o wall":"only wall",widom.sp,widom.n);

        loopto (iz,0,nz) {
          z=widom.z0+iz*widom.dz;
          zi=widom.z0+(nz-iz)*widom.dz;

          w=StaMean(string("Widom%d z=%g",m,z));

          fprintf(f,"%7.3f %12.7g %12.7g",z,w,-T*log(w));
          if (widom.mode&1) {
            if (fabs(zi-z)>widom.dz/2)
              w=StaMean(string("Widom%d z=%g",m+1,fmin(z,zi)));
            fprintf(f,"  %12.7g %12.7g",w,-T*log(w)); }
          fprintf(f,"\n"); } }
      fclose(f); }
    if (widom.spreal>=0) {
      int n;
      char *s;
      double w,dw;

      prt_("XWidom, raw chem.pot. for species (identity) change (scaling):");
      header(" mol   <exp[-DU/T]>    stderr        Dmu/K    stderr/K ");
      loop (n,0,No.N) if (cfg[n].sp==widom.spreal) {
        s=string("XWidom%d",n);
        w=StaMean(s);
        dw=StaStdErr(s);

        prt("%4d %12.6g %11.6g  %12.5f %10.5f",n,w,dw,-T*log(w),T*dw/w); }
      header(""); }
#else /*# SLAB */
    if (widom.n) {
      double w=StaMean("Widom"),dw=StaStdErr("Widom");
      double V,rho,Vm;

      if (tau.P) V=StaMean("V");
      else V=L[0]*L[1]*L[2];
      rho=rhounit*No.free_mass/V;
      Vm=V/No.N;

      prt("Widom residual chemical potential:\n\
  mu = %g +- %g kT (corrected)\n\
  mu = %g +- %g K (corrected)\n\
  mu = %g +- %g J/mol (corrected)\n\
  mu = %g +- %g K (raw value without cutoff correction)\n\
  corr = %g/V (standard cutoff correction)",
          -log(w),dw/w,
          -log(w)*T,dw/w*T,
          -log(w)*T*Eunit,dw/w*T*Eunit,
          -log(StaMean("Widom raw"))*T,StaStdErr("Widom raw")/StaMean("Widom raw")*T,
          widom.corr);
      prt("Henry law constants:\n\
  concentration -> pressure, p=c*Kc ([c]=mol/m^3, [p]=Pa):\n\
    Kc = %g +- %g Pa.m^3.mol^-1",
          Eunit*T/w,Eunit*T*dw/Sqr(w));
      prt("  pressure -> molality, mm=p*Km ([p]=bar, [mm]=mol/kg):\n\
    Km = %g +- %g mol.kg^-1.bar^-1",
          1e5/(Eunit*T/w*rho),dw*1e5/(Eunit*T*rho));
      prt("  molar fraction -> pressure, p=x*Kx ([p]=MPa, [x]=1):\n\
    Kx = %g +- %g MPa\n",
          1e-6*(Punit*T/w/Vm),dw*1e-6*(Punit*T/Vm/Sqr(w))); }

    if (widom.spreal>=0) {

      double w=StaMean("XWidom"),dw=StaStdErr("XWidom");
      prt("Widom chemical potential difference:\n\
  mu = %g +- %g kT (raw value without cutoff correction)\n\
  mu = %g +- %g K (raw value without cutoff correction)\n\
  mu = %g +- %g J/mol (raw value without cutoff correction)\n",
          -log(w),dw/w,
          -log(w)*T,dw/w*T,
          -log(w)*T*Eunit,dw/w*T*Eunit); }
#endif /*# SLAB */
