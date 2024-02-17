/* #included only once (from interpol.c)
   14 forces are not supporded by:
     fluctuating charge (interfq.c)
     qq-based splines
     gaussian charges (they need qq-based splines)
*/
    VV(dr0,=dr)

#if POLAR&64
    qqpp=si1->chargepol*si2->chargepol;
#endif /*# POLAR&64 */
    /* center-center */
    if ( (qq=si1->charge*si2->charge) ) {

#if defined(FREEBC) || defined(NIBC)
#  if POLAR&64
      ENEL+=y=(qq*qfactor-qqpp*qfactor_1)/sqrt(rr); f+=y/rr;
#  else /*# POLAR&64 */
      ENEL+=y=qq*qfactor/sqrt(rr); f+=y/rr;
#  endif /*#!POLAR&64 */
#elif COULOMB>=0
#  if POLAR&64
      ENEL+=(y=qq*qfactor-qqpp*qfactor_1)*erud(rr); f+=y*byerd;
#  else /*# POLAR&64 */
      ENEL+=(y=qfactor*qq)*erud(rr); f+=y*byerd;
#  endif /*#!POLAR&64 */
#else /*#!defined(FREEBC) || defined(NIBC)!COULOMB>=0 */
#  if POLAR&64
      y=(qq-qqpp)*qfactor_1/sqrt(rr);
      ENEL+=qq*erud(rr)+y; f+=qq*byerd+y/rr;
#  else /*# POLAR&64 */
      ENEL+=(erud(rr)+(y=qfactor_1/sqrt(rr)))*qq; f+=(byerd+y/rr)*qq;
#  endif /*#!POLAR&64 */
#endif /*#!defined(FREEBC) || defined(NIBC)!COULOMB>=0 */
    }

    MADDFORCES(f1,f2,f,dr) /* charge-charge + LJ */

    /* pol-center */
    if ( (qq=si1->chargepol*si2->charge) ) {
      VVV(dr,=dr0,+r1pol)
      rr=SQR(dr);
#if defined(FREEBC) || defined(NIBC)
#  if POLAR&64
      ENEL+=y=(qq*qfactor+qqpp*qfactor_1)/sqrt(rr); f=y/rr;
#  else /*# POLAR&64 */
      ENEL+=y=qq*qfactor/sqrt(rr); f=y/rr;
#  endif /*#!POLAR&64 */
#elif COULOMB>=0
#  if POLAR&64
      ENEL+=(y=qq*qfactor+qqpp*qfactor_1)*erud(rr); f=y*byerd;
#  else /*# POLAR&64 */
      ENEL+=qfactor*qq*erud(rr); f=qfactor*qq*byerd;
#  endif /*#!POLAR&64 */
#else /*#!defined(FREEBC) || defined(NIBC)!COULOMB>=0 */
#  if POLAR&64
      y=(qq+qqpp)*qfactor_1/sqrt(rr);
      ENEL+=qq*erud(rr)+y; f=qq*byerd+y/rr;
#  else /*# POLAR&64 */
      ENEL+=(erud(rr)+(y=qfactor_1/sqrt(rr)))*qq; f=(byerd+y/rr)*qq;
#  endif /*#!POLAR&64 */
#endif /*#!defined(FREEBC) || defined(NIBC)!COULOMB>=0 */

      MADDFORCES(f1pol,f2,f,dr) }

    if (si2->chargepol) {

      VVV(dr,=dr0,-=r2pol)

      /* center-pol */
      if (si1->charge) {
        qq=si1->charge*si2->chargepol;
        rr=SQR(dr);
#if defined(FREEBC) || defined(NIBC)
#  if POLAR&64
        ENEL+=y=(qq*qfactor+qqpp*qfactor_1)/sqrt(rr); f=y/rr;
#  else /*# POLAR&64 */
        ENEL+=y=qq*qfactor/sqrt(rr); f=y/rr;
#  endif /*#!POLAR&64 */
#elif COULOMB>=0
#  if POLAR&64
        ENEL+=(y=qq*qfactor+qqpp*qfactor_1)*erud(rr); f=y*byerd;
#  else /*# POLAR&64 */
        ENEL+=qfactor*qq*erud(rr); f=qfactor*qq*byerd;
#  endif /*#!POLAR&64 */
#else /*#!defined(FREEBC) || defined(NIBC)!COULOMB>=0 */
#  if POLAR&64
        y=(qq+qqpp)*qfactor_1/sqrt(rr);
        ENEL+=qq*erud(rr)+y; f=qq*byerd+y/rr;
#  else /*# POLAR&64 */
        ENEL+=(erud(rr)+(y=qfactor_1/sqrt(rr)))*qq; f=(byerd+y/rr)*qq;
#  endif /*#!POLAR&64 */
#endif  /*#!defined(FREEBC) || defined(NIBC)!COULOMB>=0 */

        MADDFORCES(f1,f2pol,f,dr) }

      /* pol-pol */
      if (si1->chargepol) {
#if !(POLAR&64)
        qq=si1->chargepol*si2->chargepol;
#endif /*# !(POLAR&64) */
        VV(dr0,+=r1pol)
        rr=SQR(dr0);
#if defined(FREEBC) || defined(NIBC)
#  if POLAR&64
        ENEL+=y=qqpp/sqrt(rr); f=y/rr;
#  else /*# POLAR&64 */
        ENEL+=y=qq*qfactor/sqrt(rr); f=y/rr;
#  endif /*#!POLAR&64 */
#elif COULOMB>=0
#  if POLAR&64
        ENEL+=erud(rr)*qqpp; f=byerd*qqpp;
#  else /*# POLAR&64 */
        ENEL+=qfactor*qq*erud(rr); f=qfactor*qq*byerd;
#  endif /*#!POLAR&64 */
#else  /*#!defined(FREEBC) || defined(NIBC)!COULOMB>=0 */
#  if POLAR&64
        ENEL+=erud(rr)*qqpp; f=byerd*qqpp;
#  else /*# POLAR&64 */
        ENEL+=(erud(rr)+(y=qfactor_1/sqrt(rr)))*qq; f=(byerd+y/rr)*qq;
#  endif /*#!POLAR&64 */
#endif  /*#!defined(FREEBC) || defined(NIBC)!COULOMB>=0 */
        MADDFORCES(f1pol,f2pol,f,dr0) } }
