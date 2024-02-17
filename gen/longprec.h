#ifndef LONGPREC_INCLUDED
#  define LONGPREC_INCLUDED

#  define REAL long double

#  define copysign copysignl
#  define nextafter nextafterl
#  define acos acosl
#  define asin asinl
#  define atan atanl
#  define atan2 atan2l
#  define ceil ceill
#  define cos cosl
#  define cosh coshl
#  define exp expl
#  define fabs fabsl
#  define floor floorl
#  define fmod fmodl
#  define frexp frexpl
#  define hypot hypotl
#  define ldexp ldexpl
#  define log logl
#  define log10 log10l
#  define modf modfl
#  define pow powl
#  define rint rintl
#  define sin sinl
#  define sinh sinhl
#  define sqrt sqrtl
#  define tan tanl
#  define tanh tanhl
#  define trunc truncl

#  define sqr  sqrl
#  define cub  cubl
#  define pow4 pow4l
#  define pow5 pow5l
#  define pow6 pow6l
#  define powi powil
#  define fmax fmaxl
#  define fmin fminl
#  define frac fracl
#  define cbrt cbrtl

#  ifdef Sign
#    undef Sign
#    define Sign sign
#  endif /*# Sign */

#  define PI_LONG 3.141592653589793238462643L

#  ifdef PI
#    undef PI
#    define PI PI_LONG
#  endif /*# PI */

#endif /*# LONGPREC_INCLUDED */
