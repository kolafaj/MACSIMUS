/*
  with QQTAB support
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

    ADDFORCES(f1,f2,f,dr) /* LJ */

#  undef ERFC
#  define ERFC ss->qqtab[1][1]
    /* pol-pol only (central atom not charged = COS) */
    VVV(dr,+=r1pol,-r2pol)
    rr=SQR(dr);
    f=erd(rr);
    ADDFORCES(f1pol,f2pol,f,dr)

#elif QQTAB==1

    /* both charges, all 4 terms */
    VV(dr0,=dr)
#  undef ERFC
#  define ERFC ss->qqtab[0][0]
    /* center-center */
    f+=erd(rr);
    ADDFORCES(f1,f2,f,dr) /* charge-charge + LJ */

#  undef ERFC
#  define ERFC ss->qqtab[1][0]
    /* pol-center */
    VVV(dr,=dr0,+r1pol)
    rr=SQR(dr);
    f=erd(rr);
    ADDFORCES(f1pol,f2,f,dr)

#  undef ERFC
#  define ERFC ss->qqtab[0][1]
    /* center-pol */
    VVV(dr,=dr0,-=r2pol)
    rr=SQR(dr);
    f=erd(rr);
    ADDFORCES(f1,f2pol,f,dr)

#  undef ERFC
#  define ERFC ss->qqtab[1][1]
    /* pol-pol */
    VV(dr0,+=r1pol)
    rr=SQR(dr0);
    f=erd(rr);
    ADDFORCES(f1pol,f2pol,f,dr0)

#elif QQTAB==2

    /* all possibilities optimized */
    if (si1->qtype&QTYPE_CENTER) {
      if (si2->qtype&QTYPE_CENTER) {
        /* both centers charged */
        VV(dr0,=dr)

#  undef ERFC
#  define ERFC ss->qqtab[0][0]
        /* center-center: charge-charge+LJ */
        f+=erd(rr);
        ADDFORCES(f1,f2,f,dr)

#  undef ERFC
#  define ERFC ss->qqtab[1][0]
        /* pol-center */
        VVV(dr,=dr0,+r1pol)
        rr=SQR(dr);
        f=erd(rr);
        ADDFORCES(f1pol,f2,f,dr)

#  undef ERFC
#  define ERFC ss->qqtab[0][1]
        /* center-pol */
        VVV(dr,=dr0,-=r2pol)
        rr=SQR(dr);
        f=erd(rr);
        ADDFORCES(f1,f2pol,f,dr)

#  undef ERFC
#  define ERFC ss->qqtab[1][1]
        /* pol-pol */
        VV(dr0,+=r1pol)
        rr=SQR(dr0);
        f=erd(rr);
        ADDFORCES(f1pol,f2pol,f,dr0) }
      else {
        /* center1 charged, center2 uncharged (COS) */

        /* center-center: LJ only */
        ADDFORCES(f1,f2,f,dr)

#  undef ERFC
#  define ERFC ss->qqtab[0][1]
        /* center-pol */
        VV(dr,-=r2pol)
        rr=SQR(dr);
        f=erd(rr);
        ADDFORCES(f1,f2pol,f,dr)

#  undef ERFC
#  define ERFC ss->qqtab[1][1]
        /* pol-pol */
        VV(dr,+=r1pol)
        rr=SQR(dr);
        f=erd(rr);
        ADDFORCES(f1pol,f2pol,f,dr) } }
    else {
      if (si2->qtype&QTYPE_CENTER) {
        /* center1 uncharged (COS), center2 charged */

        /* center-center: LJ */
        ADDFORCES(f1,f2,f,dr)

#  undef ERFC
#  define ERFC ss->qqtab[1][0]
        /* pol-center */
        VV(dr,+=r1pol)
        rr=SQR(dr);
        f=erd(rr);
        ADDFORCES(f1pol,f2,f,dr)

#  undef ERFC
#  define ERFC ss->qqtab[1][1]
        /* pol-pol */
        VV(dr,-=r2pol)
        rr=SQR(dr);
        f=erd(rr);
        ADDFORCES(f1pol,f2pol,f,dr) }
      else {
        /* pol-pol only (none central atom charged = COS) */
        ADDFORCES(f1,f2,f,dr) /* LJ */

#  undef ERFC
#  define ERFC ss->qqtab[1][1]
        /* pol-pol only (central atom not charged = COS) */
        VVV(dr,+=r1pol,-r2pol)
        rr=SQR(dr);
        f=erd(rr);
        ADDFORCES(f1pol,f2pol,f,dr)
      } }
#else /*#!QQTAB==0!QQTAB==1!QQTAB==2 */
#  error "QQTAB must be one of 0,1,2"
#endif  /*#!QQTAB==0!QQTAB==1!QQTAB==2 */

#if 0
    prt("%d %g %g %g  %g %g %g %g qqtab",QQTAB,
        rr,ENEL,f,
        ss->qqtab[0][0].qq,
        ss->qqtab[1][0].qq,
        ss->qqtab[0][1].qq,
        ss->qqtab[1][1].qq);
#endif /*# 0 */
