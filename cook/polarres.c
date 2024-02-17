if (scf.eps>0) {
  int N=StaN("selffield iter");
  int iter=(int)(N*StaMean("selffield iter"));

  if (option('p')/10%10>=2) {
    if (iter>N+2)
      WARNING(("ASPC and by %d more iterations than MD steps (%d)!\n\
*** Small surplus at start, esp. if far from equilibrium, is acceptable.\n\
*** Check variables =polar one-step maxerr=, =polar one-step stderr=\n\
*** =selffield iter=, and =selffield maxdr=.\n\
*** To monitor, use keywords pmax and/or pstd in %s.\n\
*** Consider decreasing scf.omega=%g or (if stable) increase scf.eps=%g.",
               iter-N,N,scf.omega,scf.eps,Fn("cpi")))
    else if (iter>N)
      prt("\
WARNING: ASPC and by %d more iterations than MD steps (%d)\n\
         scf.omega=%g scf.eps=%g (OK at start)",iter-N,N,scf.omega,scf.eps); }
}

/* special #include modul for POLAR: selected results to file */
/* Probably quite obsolete... */
if (option('v')&16) {
  long unsigned Nsta;
  FILE *f=NULL;
  double
    Ekin=No.f*StaMean("T")/2*1e-4, /* tab value (rel.err) multiplied by 1e4 */
    muperm=1e-4/Debye,             /* reduced by units of 1e4 Debye */
/*  muperm=240.26974*1e-4,          this is acetone : tab value * 1e4 */
    B=0,dB=0,dy=0;

  f=fopen(Fn("polres"),"at");
  if (!f) ERROR(("cannot open %s for writing/appending",lastFn))
  Nsta=StaN("selffield maxerr");

  if (LRRes("Etot","?")) {
    B=LRRes("Etot","B")/Ekin;
    dB=LRRes("Etot","dB")/Ekin;
    dy=LRRes("Etot","dy")/Ekin; }

  fprintf(f,"# %s -p%d -m%d -a%d -_%d scf.eps=%g h=%g\n",
	  simils.sysname,
          option('p'),option('m'),option('a'),option('_'),
          scf.eps,h);

  fprintf(f,
	  "%4.2f %d %3d %5.3f  %6.2f  %10.7f %9.7f  %7.5f %7.5f  %9.4f %6.4f %6.4f %.0f/%lu\n",
	  option('a')/optionscaling, option('p')/100, option('^')/option('_'), scf.omega,
	  StaMean("T"),
	  StaMean("selffield maxerr")/muperm, StaMean("selffield stderr")/muperm,
	  sqrt(StaMean("dEtot^2"))/Ekin, StaMean("|dEtot|")/Ekin,
	  B,dB,dy,
	  StaMean("selffield iter")*Nsta, Nsta );
  fclose(f);
}
