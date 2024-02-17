/* defines macros for cookm.c and cookmm.c - no measurement */

#ifndef POLAR

#  ifdef QQTAB
/* omit qq from the parameter list (qq in cookm.c, cookmm,c should be optimized out!!!) */
#    if PARALLEL==3
#      define LJQQX(ss,r1,r2,f1,f2,qq,G) LJQQ(ss,r1,r2,f1,f2,G)
#      define LJQQ14X(ss,r1,r2,f1,f2,qq,G) LJQQ14(ss,r1,r2,f1,f2,G)
       /* QQX not used if PARALLEL==3 */
#    else /*# PARALLEL==3 */
#      define LJQQX(ss,r1,r2,f1,f2,qq) LJQQ(ss,r1,r2,f1,f2)
#      define LJQQ14X(ss,r1,r2,f1,f2,qq) LJQQ14(ss,r1,r2,f1,f2)
#      define QQX(ss,r1,r2,f1,f2,qq) QQ(ss,r1,r2,f1,f2)
#    endif /*#!PARALLEL==3 */
#  else /*# QQTAB */
#    define LJQQX LJQQ
#    define LJQQ14X LJQQ14
#    define QQX(ss,r1,r2,f1,f2,qq) QQ(r1,r2,f1,f2,qq)
#  endif /*#!QQTAB */

#  define LJX LJ

#  if !defined(GOLD) && (defined(FREEBC) || defined(NIBC))
#    define XQQX(R1,R2,F1,F2,QQ) /* no fix for 1-2,1-3 */
#  else /*# !defined(GOLD) && (defined(FREEBC) || defined(NIBC)) */
#    define XQQX XQQ
#  endif /*#!!defined(GOLD) && (defined(FREEBC) || defined(NIBC)) */


#else /* POLAR */ /*# POLAR */

#  define polarLJQQX polarLJQQ
#  define polarLJQQ14X polarLJQQ14M /* not available in no-measurement version */

#  if defined(GOLD) || defined(NIBC)
#    error GOLD || NIBC not supported
#  endif /*# defined(GOLD) || defined(NIBC) */

#  if (POLAR&64) || (COULOMB<0) 
/* intramolecular polarizability (ADIM) and/or Ewald correction */
#    define polarXQQX polarXQQM
#  else /*# (POLAR&64) || (COULOMB<0)  */
/* no fix for 1-2,1-3 */
#    define polarXQQX(SS,R1,R2,F1,F2,S1,S2)
#  endif /*#!(POLAR&64) || (COULOMB<0)  */

#endif /*#!POLAR */
