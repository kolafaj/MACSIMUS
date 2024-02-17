/* 
   this is long version of options.h
   to be used e.g. for DOS when ints are too short
*/


#define badoption -2147483333L
extern long int optionlist[32];

#define option(X) optionlist[X & 31]  /* #defined also in ground.h */

#define getoption(A) { \
  long int val; char *o=A+2; \
  if ( (*o=='+' && o[1]==0) || *o==0) val=1; \
  else if (*o=='-' && o[1]==0) val=0; \
  else val=atol(o); \
  if (optionlist[o[-1]&31]==badoption) { \
    prts(A); Error("unknown option"); } \
  optionlist[o[-1]&31]=val; }

#define prtoption(X,OL) { \
  prt_(" -%c","@abcdefghijklmnopqrstuvwxyz[|]^_"[X&31]); \
  if (option(X)!=badoption) prt_("%ld",OL[X&31]); \
  else prtc('X'); }
