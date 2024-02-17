    VV(dr0,=dr)

    /* center-center */
    if ( (qq=si1->charge*si2->charge) ) {

#if POLAR&1
      if (anion1) f+=si1->qkappa*frep;
      if (anion2) f+=si2->qkappa*frep;
#endif /*# POLAR&1 */

#ifdef COULOMB
      f+=erd(rr)*qq;
#else /*# COULOMB */
      f+=qq/(sqrt(rr)*rr);
#endif /*#!COULOMB */
    }

    ADDFORCES(f1,f2,f,dr) /* charge-charge + LJ */

    /* pol-center */
    if ( (qq=si1->chargepol*si2->charge) ) {
      VVV(dr,=dr0,+r1pol)
      rr=SQR(dr);

#ifdef COULOMB
      f=erd(rr)*qq;
#else /*# COULOMB */
      f=qq/(sqrt(rr)*rr);
#endif /*#!COULOMB */

#if POLAR&1
      if (anion1 && rr < ss->C2q) {
        double fff;
#  define f fff

        SS_NOMEASURE_rep
        if (rr < ss->C1q)
          SS_NOMEASURE(=)
#  if 0 /* this applied when f==frep, now U,f not needed */
        else {
          x=rr-ss->C2q; f=ss->A4*x; }
#  endif /*# 0 */
#  undef f
        f-=frep*si1->qkappa; }
#endif /*# POLAR&1 */

      ADDFORCES(f1pol,f2,f,dr) }

    if (si2->chargepol) {

      VVV(dr,=dr0,-=r2pol)

      /* center-pol */
      if (si1->charge) {
        qq=si1->charge*si2->chargepol;
        rr=SQR(dr);

#ifdef COULOMB
        f=erd(rr)*qq;
#else /*# COULOMB */
        f=qq/(rr*sqrt(rr));
#endif /*#!COULOMB */

#if POLAR&1
        if (anion2 && rr < ss->C2q) {
          double fff;
#  define f fff

          SS_NOMEASURE_rep
          if (rr < ss->C1q)
            SS_NOMEASURE(=)
#  if 0 /* this applied when f==frep, now U,f not needed */
          else {
            x=rr-ss->C2q; f=ss->A4*x; }
#  endif /*# 0 */
#  undef f
          f-=frep*si2->qkappa; }
#endif /*# POLAR&1 */

        ADDFORCES(f1,f2pol,f,dr) }

      /* pol-pol */
      if (si1->chargepol) {
        qq=si1->chargepol*si2->chargepol;
        VV(dr0,+=r1pol)
        rr=SQR(dr0);
#ifdef COULOMB
        f=erd(rr)*qq;
#else /*# COULOMB */
        f=qq/(rr*sqrt(rr));
#endif /*#!COULOMB */
        ADDFORCES(f1pol,f2pol,f,dr0) } }
