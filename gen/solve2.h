#include "prec.h"

typedef REAL func_2(REAL,REAL);

int Solve2(func_2 f,func_2 g,
  REAL *x,REAL minx,REAL maxx,REAL epsx,
  REAL *y,REAL miny,REAL maxy,REAL epsy, int maxit,int method);

#if 0
const char *Solve2Name(int method);
#endif

extern struct keys_s Solve2Key[];
