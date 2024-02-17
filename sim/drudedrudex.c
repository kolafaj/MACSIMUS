#ifdef QQTAB
#  error "QQTAB not supported (yet?)"
#endif

#if defined(FREEBC) || defined(NIBC) || COULOMB>=0
#  if POLAR&64
     /* ADIM without Ewald correction */
#    define ADIMONLY
#  else /*# POLAR&64 */
     /* no contribution: should be caught before #include "drudedrudex.c" */
#    error internal: drudedrudex.c should not be #included in this case
#  endif /*#!POLAR&64 */
#endif /*# defined(FREEBC) || defined(NIBC) || COULOMB>=0 */

  /* Drude-Drude for:
     - either POLAR&64 (intramolecular induced dipole-induced dipole)
     - or Ewald corrections
     - or both
  */
    VV(dr0,=dr)

#if POLAR&64
    qqpp=si1->chargepol*si2->chargepol;
#endif /*# POLAR&64 */
    /* center-center */
    if ( (qq=si1->charge*si2->charge) ) {

#ifdef ADIMONLY
      ENEL+=y=qqpp/sqrt(rr); f=y/rr;
#elif COULOMB==-2
      ENEL+=exacterud_sqrt(rr,&y)*qq; f=y*qq;
#  if POLAR&64
      ENEL+=y=qqpp/sqrt(rr); f+=y/rr;
#  endif /*# POLAR&64 */
#else /*#!ADIMONLY!COULOMB==-2 */
      ENEL+=(erud(rr)-(y=1/sqrt(rr)))*qq; f=(byerd-y/rr)*qq;
#  if POLAR&64
      ENEL+=qqpp*y; f+=qqpp*y/rr;
#  endif /*# POLAR&64 */
#endif /*#!ADIMONLY!COULOMB==-2 */

      MADDFORCES(f1,f2,f,dr) }

    /* pol-center */
    if ( (qq=si1->chargepol*si2->charge) ) {
      VVV(dr,=dr0,+r1pol)
      rr=SQR(dr);
#ifdef ADIMONLY
      ENEL-=y=qqpp/sqrt(rr); f=-y/rr;
#elif COULOMB==-2
      ENEL+=exacterud_sqrt(rr,&y)*qq; f=y*qq;
#  if POLAR&64
      ENEL-=y=qqpp/sqrt(rr); f-=y/rr;
#  endif /*# POLAR&64 */
#else /*#!ADIMONLY!COULOMB==-2 */
      ENEL+=(erud(rr)-(y=1/sqrt(rr)))*qq; f=(byerd-y/rr)*qq;
#  if POLAR&64
      ENEL-=qqpp*y; f-=qqpp*y/rr;
#  endif /*# POLAR&64 */
#endif /*#!ADIMONLY!COULOMB==-2 */

      MADDFORCES(f1pol,f2,f,dr) }

    if (si2->chargepol) {

      VVV(dr,=dr0,-=r2pol)

      /* center-pol */
      if (si1->charge) {
        qq=si1->charge*si2->chargepol;
        rr=SQR(dr);
#ifdef ADIMONLY
        ENEL-=y=qqpp/sqrt(rr); f=-y/rr;
#elif COULOMB==-2
        ENEL+=exacterud_sqrt(rr,&y)*qq; f=y*qq;
#  if POLAR&64
        ENEL-=y=qqpp/sqrt(rr); f-=y/rr;
#  endif /*# POLAR&64 */
#else /*#!ADIMONLY!COULOMB==-2 */
        ENEL+=(erud(rr)-(y=1/sqrt(rr)))*qq; f=(byerd-y/rr)*qq;
#  if POLAR&64
        ENEL-=qqpp*y; f-=qqpp*y/rr;
#  endif /*# POLAR&64 */
#endif  /*#!ADIMONLY!COULOMB==-2 */

        MADDFORCES(f1,f2pol,f,dr) }

      /* pol-pol */
      if (si1->chargepol) {
        qq=si1->chargepol*si2->chargepol;
        VV(dr0,+=r1pol)
        rr=SQR(dr0);
#ifdef ADIMONLY
        ENEL+=y=qqpp/sqrt(rr); f=y/rr;
#elif COULOMB==-2
        /* note: qq=qqpp */
        ENEL+=exacterud_sqrt(rr,&y)*qq; f=y*qq;
#  if POLAR&64
        ENEL+=y=qqpp/sqrt(rr); f+=y/rr;
#  endif /*# POLAR&64 */
#else /*#!ADIMONLY!COULOMB==-2 */
#  if POLAR&64
        ENEL+=erud(rr)*qq; f=byerd*qq;
#  else /*# POLAR&64 */
        ENEL+=(erud(rr)-(y=1/sqrt(rr)))*qq; f=(byerd-y/rr)*qq;
#  endif /*#!POLAR&64 */
#endif  /*#!ADIMONLY!COULOMB==-2 */

        MADDFORCES(f1pol,f2pol,f,dr0) } }
