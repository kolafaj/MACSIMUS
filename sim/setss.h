/* WARNING: watch the directory of sitesite.c, sitesite.h (see metamake) */
#include "sitesite.h" 

/*
module with definitions of site-site potentials
This is the LJ version:

  uLJ(r) = 4*eps*((sigma/r)^12-(sigma/r)^6)

This version is the `default' version in `sim' directory:
interpot.c can be overriden by a system-specific version - interpot.h should
remain the same in `sim'
*/

void setss(sitesite_t *ss, int i,int j, double C2,int onefour);

extern double poteps;

#ifdef WIDOM
void widomrdf(int mode);
#endif
