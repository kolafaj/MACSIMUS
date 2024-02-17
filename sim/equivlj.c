/*
  To be #included by sim/XXX/sitesite.c to initssaux() if SLAB.
  See equivlj.h for more info.
  NOTE: may be overwritten by wall.LJ[].* data, see initwall in simdef.c
*/

#include "ms.h"

{
  double x,y,z,f,rr=30,U=0;
  struct locaux_s { pairaux_t a; } locss,*ss=&locss;

  a->LJsig=1;
  a->LJeps=0;
  if (a->Aij) {
    ss->a=*a;
    do {
      if ( (rr*=0.99) < 0.1 ) break;
      SS_NOMEASURE(=)
    } while (f<0);

    MS_BEGIN(rr,1e-10)
      SS_NOMEASURE(=)
      MS_f=f;
    MS_END(rr,1)
    SS_MEASURE

    a->LJsig=sqrt(rr/cbrt(2));
    a->LJeps=-U; }

  a->Sq=Sqr(a->LJsig);
  a->E4=4*a->LJeps;

  if (rr<0.1) 
    WARNING(("initssaux (%g %g): cannot determine potential minimum (f=0)",
             pp->eps,pp->sig))
  else
    prt("equivalent LJeps=%g LJsig=%g will be used for integrated atom-wall potential",
        a->LJeps,a->LJsig);
}
