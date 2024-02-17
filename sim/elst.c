#include "ground.h"
#include "elst.h"

#ifdef COULOMB

#  if COULOMB<0
#    ifdef QQTAB
#      ifdef GAUSSIANCHARGES
#        include "gcelst.c"
#      else /*# GAUSSIANCHARGES */
#        include "erfcqq.c"
#      endif /*#!GAUSSIANCHARGES */
#    else /*# QQTAB */
#      include "erfc.c"
#    endif /*#!QQTAB */
#  elif COULOMB==0
#    include "cutelstd.c"
#  elif COULOMB==2
#    include "cutelst.c"
#  elif COULOMB==3
#    include "fgelst.c"
#  else /*#!COULOMB<0!COULOMB==0!COULOMB==2!COULOMB==3 */
#    error this COULOMB value not supported
#  endif /*#!COULOMB<0!COULOMB==0!COULOMB==2!COULOMB==3 */

#endif /*# COULOMB */
