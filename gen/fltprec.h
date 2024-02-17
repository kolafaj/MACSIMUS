/* WARNING: never tested */
#ifndef FLTPREC_INCLUDED

#  define FLTPREC_INCLUDED

#  define REAL float

#  define copysign copysignf
#  define nextafter nextafterf
#  define acos acosf
#  define asin asinf
#  define atan atanf
#  define atan2 atan2f
#  define ceil ceilf
#  define cos cosf
#  define cosh coshf
#  define exp expf
#  define fabs fabsf
#  define floor floorf
#  define fmod fmodf
#  define frexp frexpf
#  define hypot hypotf
#  define ldexp ldexpf
#  define log logf
#  define log10 log10f
#  define modf modff
#  define pow powf
#  define rint rintf
#  define sin sinf
#  define sinh sinhf
#  define sqrt sqrtf
#  define tan tanf
#  define tanh tanhf
#  define trunc truncf

#  define sqr  sqrf
#  define cub  cubf
#  define pow4 pow4f
#  define pow5 pow5f
#  define pow6 pow6f
#  define powi powif
#  define fmax fmaxf
#  define fmin fminf
#  define frac fracf
#  define cbrt cbrtf

#  ifdef Sign
#    undef Sign
#    define Sign sign
#  endif /*# Sign */

#  define PI_FLOAT 3.14159265358F

#  ifdef PI
#    undef PI
#    define PI PI_FLOAT
#  endif /*# PI */

#endif /*# FLTPREC_INCLUDED */
