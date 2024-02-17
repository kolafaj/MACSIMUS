    VV(dr0,=dr)

    /* center-center */
    if ( (qq=si1->charge*si2->charge) ) {

#if POLAR&1
      if (anion1)
        U+=si1->qkappa*Urep,ENVIR+=si1->qkappa*frep*rr,f+=si1->qkappa*frep;
      if (anion2)
        U+=si2->qkappa*Urep,ENVIR+=si2->qkappa*frep*rr,f+=si2->qkappa*frep;
#endif /*# POLAR&1 */

#ifdef COULOMB
      ENEL+=erud(rr)*qq; f+=byerd*qq;
#else /*# COULOMB */
      ENEL+=y=qq/sqrt(rr); f+=y/rr;
#endif /*#!COULOMB */
    }

    MADDFORCES(f1,f2,f,dr) /* charge-charge + LJ */

    /* pol-center */
    if ( (qq=si1->chargepol*si2->charge) ) {
      VVV(dr,=dr0,+r1pol)
      rr=SQR(dr);

#ifdef COULOMB
      ENEL+=erud(rr)*qq; f=byerd*qq;
#else /*# COULOMB */
      ENEL+=y=qq/sqrt(rr); f=y/rr;
#endif  /*#!COULOMB */

#if POLAR&1
      if (anion1 && rr < ss->C2q) {
        double z,frep=0,Urep=0,fff,UUU;
#  define f fff
#  define U UUU
        SS_MEASURE_rep
        if (rr < ss->C1q)
          SS_MEASURE
#  if 0 /* this applied when f==frep, now U,f not needed */
        else {
          x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
#  endif /*# 0 */
        frep*=si1->qkappa;
        ENVIR-=rr*frep;
#  undef U
#  undef f
        U-=si1->qkappa*Urep;
        f-=frep; }
#endif /*# POLAR&1 */

      MADDFORCES(f1pol,f2,f,dr) }

    if (si2->chargepol) {

      VVV(dr,=dr0,-=r2pol)

      /* center-pol */
      if (si1->charge) {
        qq=si1->charge*si2->chargepol;
        rr=SQR(dr);

#ifdef COULOMB
        ENEL+=erud(rr)*qq; f=byerd*qq;
#else /*# COULOMB */
        ENEL+=y=qq/sqrt(rr); f=y/rr;
#endif /*#!COULOMB */

#if POLAR&1
        if (anion2 && rr < ss->C2q) {
          double fff,UUU;
#  define f fff
#  define U UUU

          SS_MEASURE_rep
          if (rr < ss->C1q)
            SS_MEASURE
#  if 0 /* this applied when f==frep, now U,f not needed */
          else {
            x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
#  endif /*# 0 */
          frep*=si2->qkappa;
          ENVIR-=rr*frep;
#  undef U
#  undef f
          U-=si2->qkappa*Urep;
          f-=frep; }
#endif /*# POLAR&1 */

          MADDFORCES(f1,f2pol,f,dr) }

      /* pol-pol */
      if (si1->chargepol) {
        qq=si1->chargepol*si2->chargepol;
        VV(dr0,+=r1pol)
        rr=SQR(dr0);

#ifdef COULOMB
        ENEL+=erud(rr)*qq; f=byerd*qq;
#else /*# COULOMB */
        ENEL+=y=qq/sqrt(rr); f=y/rr;
#endif /*#!COULOMB */
        MADDFORCES(f1pol,f2pol,f,dr0) } }
