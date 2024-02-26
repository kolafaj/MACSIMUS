/*
  #included from constrd.c
  Maxwell/Boltzmann (tau.T>0) and simplified Langevin (tau.T<0)
  thermostat for the center of mass
  WARNING: does not recalculate the kinetic energy using the new v(t+h/2)
*/
  vector vCM,dvCM;
  double sigma;

  if (tau.CM<0)
    /* Langevin; NB: different factor than for LANGEVIN thermostat */
    sigma=h*sqrt(2*T/No.mass*(h/fabs(tau.CM)));
  else {
    /* Maxwell-Boltzmann */
    sigma=0;
    static int it=0;
    double when=(t+h/2)/tau.CM-it;

    if (when>1) {
      it++;
      sigma=h*sqrt(T/No.mass); } }

  /* calculate the velocity of the CoM */
  if (sigma) {
    VO(vCM,=0)
    loop (n,FROM,No.N) {
      mn=molec+n;
      v=rof(mn,cfg[1]->rp);
      si=spec[mn->sp]->si;
      ns=mn->ns;

      /* cf. norm.c */
      loop (i,0,ns) VV(vCM,+=si[i].mass*v[i]) }
    VO(vCM,/=No.mass)

    if (tau.CM<0) /* Langevin */
      VVO(dvCM,=(h/tau.CM)*vCM,+sigma*rndgauss())
    else /* Maxwell-Boltzmann */
      VVO(dvCM,=-vCM,+sigma*rndgauss())

    loop (n,FROM,No.N) {
      mn=molec+n;
      v=rof(mn,cfg[1]->rp);
      si=spec[mn->sp]->si;
      ns=mn->ns;

      /* cf. norm.c */
      loop (i,0,ns) VV(v[i],+=dvCM) }

    sigma=SQR(vCM);
    VV(vCM,+=dvCM)
    sigma+=SQR(vCM);

    StaAdd("TkCM",No.mass*sigma/6/Sqr(h)); }
