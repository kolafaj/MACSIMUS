/* defines macros for cookm.c, cookmm.c, water.c - measurements included */

#ifndef POLAR

#  ifdef QQTAB
/* omit qq from the parameter list (qq in cookm.c, cookmm.c should be optimized out!!!) */
#    if PARALLEL==3
#      define LJQQX(ss,r1,r2,f1,f2,qq,G) U+=LJQQM_MM(ss,r1,r2,f1,f2,G)
#      define LJQQ14X(ss,r1,r2,f1,f2,qq,G) U+=LJQQM_MM(ss,r1,r2,f1,f2,G)
#      define QQX(ss,r1,r2,f1,f2,qq,G) QQM_MM(ss,r1,r2,f1,f2,G)
#    else /*# PARALLEL==3 */
#      define LJQQX(ss,r1,r2,f1,f2,qq) U+=LJQQM_MM(ss,r1,r2,f1,f2)
#      define LJQQ14X(ss,r1,r2,f1,f2,qq) U+=LJQQM_MM(ss,r1,r2,f1,f2)
#      define QQX(ss,r1,r2,f1,f2,qq) QQM_MM(ss,r1,r2,f1,f2)
#    endif /*#!PARALLEL==3 */
#    define LJX U+=LJM_MM
#  else /*# QQTAB */
#    define LJQQX U+=LJQQM_MM
#    define LJQQ14X U+=LJQQM_MM
#    define QQX QQM_MM
#    define LJX U+=LJM_MM
#  endif /*#!QQTAB */

#  ifdef QQTAB
#    define QQPARM /* no qq needed */
#  else /*# QQTAB */
#    define QQPARM ,qq
#  endif /*#!QQTAB */

#else

#  define polarLJQQX U+=polarLJQQM_MM
#  define polarLJQQ14X U+=polarLJQQ14M_MM

#  if (POLAR&64) || (COULOMB<0)
# error very likely not compatible combination
/* intramolecular polarizability (ADIM) and/or Ewald correction */
#    define polarXQQX U+=polarXQQM
#  else /*# (POLAR&64) || (COULOMB<0)  */
/* no fix for 1-2,1-3 */
#    define polarXQQX(SS,R1,R2,F1,F2,S1,S2)
#  endif /*#!(POLAR&64) || (COULOMB<0)  */

#endif /*#!POLAR */
