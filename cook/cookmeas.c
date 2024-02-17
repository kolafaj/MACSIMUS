/* defines macros for cookm.c, cookmm.c, water.c - measurements included */

#ifndef POLAR

#  ifdef QQTAB
/* omit qq from the parameter list (qq in cookm.c, cookmm.c should be optimized out!!!) */
#    if PARALLEL==3
#      define LJQQX(ss,r1,r2,f1,f2,qq,G) U+=LJQQM(ss,r1,r2,f1,f2,G)
#      define LJQQ14X(ss,r1,r2,f1,f2,qq,G) U+=LJQQ14M(ss,r1,r2,f1,f2,G)
#      define QQX(ss,r1,r2,f1,f2,qq,G) QQM(ss,r1,r2,f1,f2,G)
#    else /*# PARALLEL==3 */
#      define LJQQX(ss,r1,r2,f1,f2,qq) U+=LJQQM(ss,r1,r2,f1,f2)
#      define LJQQ14X(ss,r1,r2,f1,f2,qq) U+=LJQQ14M(ss,r1,r2,f1,f2)
#      define QQX(ss,r1,r2,f1,f2,qq) QQM(ss,r1,r2,f1,f2)
#    endif /*#!PARALLEL==3 */
#    define LJX U+=LJM
#  else /*# QQTAB */
#    define LJQQX U+=LJQQM
#    define LJQQ14X U+=LJQQ14M
#    define QQX QQM
#    define LJX U+=LJM
#  endif /*#!QQTAB */

#  if !defined(GOLD) && (defined(FREEBC) || defined(NIBC))
#    define XQQX(R1,R2,F1,F2,QQ) /* no fix for 1-2,1-3 */
#  else /*# !defined(GOLD) && (defined(FREEBC) || defined(NIBC)) */
#    define XQQX U+=XQQM
#  endif /*#!!defined(GOLD) && (defined(FREEBC) || defined(NIBC)) */

#  ifdef QQTAB
#    define QQPARM /* no qq needed */
#  else /*# QQTAB */
#    define QQPARM ,qq
#  endif /*#!QQTAB */

#else /* POLAR */ /*# POLAR */

#  define polarLJQQX U+=polarLJQQM
#  define polarLJQQ14X U+=polarLJQQ14M

#  if defined(GOLD) || defined(NIBC)
#    error GOLD || NIBC not supported
#  endif /*# defined(GOLD) || defined(NIBC) */

#  if (POLAR&64) || (COULOMB<0)
/* intramolecular polarizability (ADIM) and/or Ewald correction */
#    define polarXQQX U+=polarXQQM
#  else /*# (POLAR&64) || (COULOMB<0)  */
/* no fix for 1-2,1-3 */
#    define polarXQQX(SS,R1,R2,F1,F2,S1,S2)
#  endif /*#!(POLAR&64) || (COULOMB<0)  */

#endif /*#!POLAR */
