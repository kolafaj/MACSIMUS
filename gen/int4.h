/* legacy code.. */
#ifndef _INT4
#  define _INT4

typedef int int4;
typedef unsigned unsigned4;

#if 0
/* looks dirty, but prevents unnecessary int -> long conversion */
#  ifndef sizeof
#    define sizeof (int)sizeof /* to suppress some warnings */
#  endif /*# sizeof */
#  ifndef strlen
#    define strlen (int)strlen /* to suppress some warnings */
#  endif /*# strlen */
#endif

#endif /*# _INT4 */
