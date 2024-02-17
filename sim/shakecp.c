/* 
  #include module to constrd.c, function Shake()
  extended Lagrangian (Car-Parrinello-like) method 
  WARNING: no separate thermostating of these degrees of freedom
*/
  En.kinpol=0;
  if (option('p')/10%10==0) {
    /* Lagrange (Car-Parrinello) method: Verlet integration of the aux site */
    FILE *masspol=NULL;
    static int first=1;

    if (first && option('v')&32) {
      masspol=fopen(Fn("mass.pol"),"wt");
      first=0;
      fprintf(masspol,"#mol si mass[prog.u] mass[g/mol]\n"); }

    En.self=0;

    loop (n,FROM,No.N) {
      mn=molec+n;
      sp=mn->sp;
      ns=mn->ns;
      r=polarrof(mn,rp); /* at time t */
      r1=polarrof(mn,cfg[1]->rp); /* at time t-h */
      p=polarrof(mn,rp2); /* forces */
      si=spec[sp]->si;
      loop (i,0,ns) if (si[i].chargepol) {
        double K=1./si[i].alpha_qq; /* = Q*Q/alpha */
        double M=Sqr(tau.dip)*K,IM=1./M;

        if (masspol) fprintf(masspol,"%3d %3d  %g  %g\n",n,i,M,M*Munit);

        VV(p[i],-=K*r[i]) /* harm. oscillator force */
        VO(p[i],*=hh*IM) /* acceleration */

        loop (al,0,3) {
          z=r[i][al]-r1[i][al]; /* old vel. */
          IM=z+p[i][al]; /* new vel. */
#if VERLET==2
          En.kinpol+=M*(z*IM); /* 2*hh*kin.energy */
#elif VERLET==3 || VERLET==30
          En.kinpol+=M*(z*z+IM*IM); /* 4*hh*kin.energy */
#else /*#!VERLET==2!VERLET==3 || VERLET==30 */
          ERROR(("Ekin with VERLET = 0,1 not supported for Car-Parrinello"))
#endif /*#!VERLET==2!VERLET==3 || VERLET==30 */
          En.self+=K*Sqr(r[i][al]); /* 2*pot.en. */

          /* Verlet: */
          z=r[i][al]+IM;
          r1[i][al]=r[i][al];
          r[i][al]=z; } } }

    if (masspol) fclose(masspol);

    En.self/=2;
    En.kin+=En.kinpol; /* scaled by 4*hh or 2*hh */
    En.Tpol=En.kinpol/(No.free_pol*hh*(DIM*2));
  } /* Car-Parrinello */
