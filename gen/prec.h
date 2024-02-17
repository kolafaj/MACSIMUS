#ifndef PREC_INCLUDED
#  define PREC_INCLUDED

/*
  Interface to PRECISION for directly #included modules; for separately
  compiled modules, use the ground.h interface.
  Must be #included after all system header files and after func.c
  (see, e.g., 2dfit.c).

  PRECISION==0 : float version (not tested)
  PRECISION==1 : double version (default if PRECISION is not #defined)
  PRECISION==2 : long double version
                 DOS/TC warning: check longprec.h and #define MATHLL !
  PRECISION>2  : emulated high precision numbers (HIGHPREC=PRECISION)
                 in C++ (see modules highprec.h, highprec.cc)
  (HIGHPREC==2 not available, but it almost does not make sense)
*/

#  ifndef PRECISION
#    define PRECISION 1
#  endif /*# PRECISION */

#  if PRECISION<=2
#    define Val(X) (double)(X)
#    define Int(X) (int)(X)
#  endif /*# PRECISION<=2 */

#  if PRECISION>2
#    define HIGHPREC PRECISION
#    include "highprec.h"
#  elif PRECISION==2
#    include "longprec.h"
#  elif PRECISION==1
#    define REAL double
#  elif PRECISION==0
#    include "fltprec.h"
#  else /*#!PRECISION>2!PRECISION==2!PRECISION==1!PRECISION==0 */
#    error bad PRECISION
#  endif /*#!PRECISION>2!PRECISION==2!PRECISION==1!PRECISION==0 */

#endif /*# PREC_INCLUDED */
