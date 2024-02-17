/*
  Ewald corrections with QQTAB support
*/
#ifndef QQTAB
#  error "QQTAB required"
#endif /*# QQTAB */
#if POLAR&(1+32+64)
#  error "unsupported POLAR version"
#endif /*# POLAR&(1+32+64) */
#if !defined(COULOMB) || COULOMB>=0
#  error "unsupported COULOMB version"
#endif /*# !defined(COULOMB) || COULOMB>=0 */

#if QQTAB==0

#undef ERFC
#  define ERFC ss->qqtabx[1][1]
    /* pol-pol */
    VVV(dr,+=r1pol,-r2pol)
    rr=SQR(dr);
    ENEL+=erud(rr); f=byerd;
    MADDFORCES(f1pol,f2pol,f,dr)

#else /*# QQTAB==1 */
    /* not optimized as for QQTAB==2 */
    VV(dr0,=dr)

#undef ERFC
#  define ERFC ss->qqtabx[0][0]
    /* center-center */
    ENEL+=erud(rr); f=byerd;
    MADDFORCES(f1,f2,f,dr)

#undef ERFC
#  define ERFC ss->qqtabx[1][0]
    /* pol-center */
    VVV(dr,=dr0,+r1pol)
    rr=SQR(dr);
    ENEL+=erud(rr); f=byerd;
    MADDFORCES(f1pol,f2,f,dr)

#undef ERFC
#  define ERFC ss->qqtabx[0][1]
    /* center-pol */
    VVV(dr,=dr0,-=r2pol)
    rr=SQR(dr);
    ENEL+=erud(rr); f=byerd;
    MADDFORCES(f1,f2pol,f,dr) 

#undef ERFC
#  define ERFC ss->qqtabx[1][1]
    /* pol-pol */
    VV(dr0,+=r1pol)
    rr=SQR(dr0);
    ENEL+=erud(rr); f=byerd;
    MADDFORCES(f1pol,f2pol,f,dr0)

#endif /*#!QQTAB==1 */
