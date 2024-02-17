#include "ground.h"
#include "simglob.h"
#include "simgear.h"

/* MAXGEARORDER=9 assumed, see simglob.h */
struct gear_s gear = {
  {-9e9,-9e9,-9e9,-9e9,-9e9,-9e9,-9e9,-9e9,-9e9}, /* C[] */
  3, /* default = Janek */
  0 /* lastorder */
};

#define D(X) ((real*)&((X)->logs))

void rhs(ToIntPtr B, ToIntPtr A, ToIntPtr V);

/* 
   gear2.c is in ../gen/, since V3.6g for both POLAR and nonpolar version 
   #define VELOCITY (is in simgear.h) is required
   edit gear2.C and use C2c !
*/

#include "gear2.c"
