/* #included from  ewald.c */

if (charges.sumq!=0 && el.kappa>0) {
  double V=PROD(box.L);
  double acq=Sqr(el.alpha*box.cutoff);
  double dEr=charges.sumq*sqrt(box.cutoff/2/V)/acq/exp(acq);
  double dfr=2*sqrt(charges.sumq/box.cutoff/V)/exp(acq);
  double eaL=exp(-Sqr(PI/el.alpha*el.kappa));
  double dEk=charges.sumq*el.alpha/(PI*PI*el.kappa*sqrt(V*el.kappa))*eaL;
  double dfk=sqrt(8*charges.sumq/el.kappa/V)*el.alpha/PI*eaL;
  double rat;

  charges.bgE=-PI/2*Sqr(charges.sum/el.alpha);
  
  /* PROBLEMS - alleviated in V3.3b, see ewald.c */
  if (!el.sf && el.alpha>=0 && Sqr(charges.sum)/charges.sumq>Sqr(el.epsq)) {
    if (fabs(el.epsinf)<1e16) 
      ERROR(("The system is charged (%g e) and not tin-foil (el.epsinf=%g),\n\
*** I don\'t know how to handle the background",charges.sum/electron, el.epsinf))
    else 
      prt("The system is charged (%g e),\n\
the \"background\" term %g/V will be added to the energy",
          charges.sum/electron, charges.bgE); }

  if (el.epsinf<0) WARNING(("\
The dielectric constant of the surrounding continuum is negative:\n\
***   el.epsinf=%g",el.epsinf))

  if (fabs(el.epsinf+0.5)<1e-5)
    ERROR(("el.epsinf=%g out of range (too large dipolar response)",el.epsinf))

  prt("\n++++++ Ewald initialized ++++++  alpha = %.6f/AA   kappa = %.5f/AA",
      el.alpha,el.kappa);

  if ((option('v')&1) && el.alpha>0) { /* <<<<<<<<<<<<<<<<<<<<<<<<<<<< */

    prt(
        "Q->N = #{k,|k|<=K,k=-k} = %i   "
        "Nkk = #{k,|k|<=K,kx>=0,ky>=0,kz>=0} = %i",
        iQ,ikk);
    prt("%d site charges  SUM_i(q_i^2) = %g K*A = %g e^2\n\
SUM_i(q_i) = %g  SUM_i((q_i/m_i)^2) = %g",
        Nq,charges.sumq,charges.sumq/Sqr(electron),
        charges.sum,charges.summq);

#ifndef POLAR
    prt("sums over molecules: SUM_n(q_n^2) = %g  SUM_n,uncharged(mu_n^2) = %g",
      charges.molsumq,charges.molsumdq);
#endif

    /*** error estimates ***/

    /* standard deviations of expected errors caused by cutoffs */
    prt("\n\
The following estimates are based on the assumption of random distribution of\n\
charges. They work quite well for ionic systems. For dipolar systems without\n\
free charges (e.g., water), the k-space estimates are too pessimistic and\n\
actual errors are roughly 10x smaller; thus, epsk=10*epsr is recommended.");
    /* diagonal (always positive) terms omitted by |k|<K cutoff */
    if (el.diag>=0)
      prt("* diagonal correction to k-space cutoff error of total energy (%sincluded):\n\
  %g K",el.diag?"":"not ", Ediag);
    prt("* expected non-diagonal cutoff error of total energy:\n\
  r: %.2e K        k: %.2e K", dEr,dEk);
    prt("* expected cutoff error of force on maximum charge %g = %g e:\n\
  r: %.2e K/AA     k: %.2e K/AA",
        charges.max, charges.max/electron,
        charges.max*dfr, charges.max*dfk);

    rat=charges.max*dfr/el.epsr;
    if (rat>1.03)
      prt("WARNING: estimated r-space error is %.3fx larger than limit epsr=%g K/AA",rat,el.epsr);
    rat=charges.max*dfk/el.epsk;
    if (rat>1.03)
      prt("WARNING: estimated k-space error is %.3fx larger than limit epsk=%g K/AA",rat,el.epsk);

    /* square roots of sums of squared errors over all charges */
    e=sqrt(charges.sumq);
    prt("* expected summarized standard cutoff error of forces:\n\
  r: %.2e K/AA     k: %.2e K/AA", estimr=e*dfr,estimk=e*dfk); // Bug fixed by T.Trnka

    e=sqrt(charges.summq);
    prt("* expected summarized standard cutoff error of accelerations:\n\
  r: %.2e AA/ps^2  k: %.2e AA/ps^2",
        //        estimr=e*dfr,estimk=e*dfk); BUG: removed by T.Trnka
        e*dfr,e*dfk);

    prt("NB: summarized standard errors are defined as sqrt[SUM_i(f_i-f_iref)^2]\n\
    - there is no division by N because of problematic definition in mixtures");

  } /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

  loop (i,0,3)
    if (box.cutoff>box.L[i]*0.50000001)
      prt("WARNING cutoff/box.L[%d] = %.7f > 1/2",i,box.cutoff/box.L[i]); }
