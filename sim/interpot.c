/*
  control module for definitions of site-site potentials
  (#included from sim/DIR/sitesite.c and sim/DIR/sitesite.h)

  this module depends on the following #defines:
     WORM, CUT, TWODIM, POLAR, PARALLEL, STARS, PERSUM

  NOTE: initcombrule() and combrule() are in sitesite.c
*/

#include "ground.h"
#include "simglob.h"

#ifdef COULOMB
#  include "elst.h"
#else /*# COULOMB */
#  define ermacrodcl /* empty */
#  undef ermacro
#endif /*#!COULOMB */

#include "interpot.h"

#ifdef TWODIM
#  error TWODIM out of date
#  include "inter2d.h"
#else /*# TWODIM */
#  include "intermac.h"
#endif /*#!TWODIM */
#include "simdef.h"
#include "units.h"

/* to move ??? -- and also LJQQ14* ? */
sitesite_t **sstab,**sstab14;
real factor14=0.5,factor14_1=-0.5;

/* shared memory PARALLEL support */
#if PARALLEL==3
#  error NOT TESTED RECENTLY
#  define GLOB ,pll_global_t *glob
#  define ENVIR glob->Envir
#  define ENPVIR glob->Pvir
#  define ENEL glob->Enel
#else  /*# PARALLEL==3 */
#  define GLOB
#  define ENVIR En.vir
#  define ENPVIR En.Pvir
#  define ENEL En.el
#endif /*#!PARALLEL==3 */

#ifdef MARK
#  include "mark.h"
#else /*# MARK */
#  define MARK_PATCH /*empty*/
#endif /*#!MARK */

/*
  this ugly patch is for simulating gravitational forces instead of Coulomb ones
*/
#ifdef STARS
#  define SG qq=-qq; /* gravity, FREEBC recommended */
#else /*# STARS */
#  define SG         /* (empty) the default charge-charge */
#endif /*#!STARS */

#ifdef POLAR

#  if POLAR&32
#    include "interfq.c"
#  else /*# POLAR&32 */
#    include "interpol.c"
#    ifdef PERSUM
#      define PERSUM_MM
#      error NOT IMPLEMENTED
#      include "interpol.c"
#    endif /*# PERSUM */
#  endif /*#!POLAR&32 */

#else /*# POLAR */

#  ifdef METAL
#    include "metalnp.c"
#    ifdef PERSUM
#      error NOT IMPLEMENTED
#    endif /*# PERSUM */
#  else /*# METAL */
#    include "internp.c"
#    ifdef PERSUM
#      define PERSUM_MM
#      include "internp.c"
#    endif /*# PERSUM */
#  endif /*#!METAL */

#endif /*#!POLAR */
