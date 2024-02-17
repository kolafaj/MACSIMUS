/*** options parsing and using support
 *** usage: ***

defining:
^^^^^^^^^
  #define x badoption
  int optionlist[32] =
  (* ` a b c d e f g h i j k l m n o p q r s t u v w x y z { | } ~   *)
    {x,x,x,x,x,x,x,x,x,0,x,x,x,0,x,x,0,0,x,x,x,x,0,x,x,x,x,x,x,x,x,x};
  (* @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z [ \ ] ^ _ *)
  #undef x

parsing:
^^^^^^^^
  loop (i,0,narg)
    if (arg[i][0]=='-') getoption(arg[i])

value of option -x :
^^^^^^^^^^^^^^^^^^^^
  option('x')

list option in the form " -a1"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  prtoption('a')

warning:
^^^^^^^^
macro option is #defined also in ground.h
*/

/* ` a b c d e f g h i j k l m n o p q r s t u v w x y z { | } ~   */
/* @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z [ \ ] ^ _ */

/*** OOPS!
the following smart trick had to be abandoned because the damned Turbo C 2.0
preprocessor was not able to compile it correctly !
 #define option(X) optionlist[* #X & 31] ***/


#define badoption -32333
extern int optionlist[32];

#define option(X) optionlist[X & 31]  /* #defined also in ground.h */

#define getoption(A) { \
  int val; char *o=A+2; \
  if ( (*o=='+' && o[1]==0) || *o==0) val=1; \
  else if (*o=='-' && o[1]==0) val=0; \
  else val=strtol(o,NULL,0); /* accepts oct,hex */ \
  if (optionlist[o[-1]&31]==badoption) { \
    prts(A); Error("unknown option"); } \
  optionlist[o[-1]&31]=val; }

#define prtoption(X,OL) { \
  prt_(" -%c","@abcdefghijklmnopqrstuvwxyz[|]^_"[X&31]); \
  if (option(X)!=badoption) prt_("%d",OL[X&31]); \
  else prtc('X'); }
