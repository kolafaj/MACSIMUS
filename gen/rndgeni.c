/* 
  random number generator
  ^^^^^^^^^^^^^^^^^^^^^^^
  to be directly #included
  for a module, see rndgen.c and rndgen.h

  four-tap register-shift generator:
  - grid of rnd() is at least 2^-53 (as fine as possible)
  - good quality/speed ratio

  initialization: rndinit(tablen,seed)
  where tablen<0: query whether already initialized
        tablen=0: use default
        tablen>0: shuffle table lenght for the initializer
                  (the generator itself does not use shuffle table)
        seed = any integer (0 uses seed from time)

  RNDMETHOD is the sum of (applies for IEEE numbers):
    rnd():
    0: rnd in [0,1] 
       if cast to IEEE double, can be exactly: 
         0 with probability 2^-64=5.42e-20
         1 with probability 2^-54=5.55e-17
       if kept in register (10 bytes long double), can be exactly:
         0 with probability 2^-64=5.42e-20
         1 never
    1: rnd in [0,1) = exact (double)1 is never reached
      
    rndcos():
    0 = rndcos in [-1,1] (can be exactly +1 or -1 with prob. < 1e-16)
    1 = rndcos in [-1,1), IEEE numbers assumed
    1+2 = rndcos in (-1,1), IEEE numbers assumed
  
    irnd(N):
    0: based on 32 bits (simpler and faster, good enough for most cases)
    4: based on 64 bits (rarely important)
    NB: irnd(N) uses integer aritmetic and NEVER gives N

  Note that generally code like "double x,y; x=y; if (x==y)" is NOT guaranteed
  to pass because the first x may be kept in a register (10 bytes long) and
  the second not.
*/

#ifndef RNDMETHOD 
#  define RNDMETHOD 0
#endif /*# RNDMETHOD  */

#include "rndgen.h"
#include "rndseed.c"

#define A 471  
#define B 1586 
#define C 6988 
#define D 9689 
#define M 16383
#define trnd (++nd,ra[nd&M]=ra[(nd-A)&M]^ra[(nd-B)&M]^ra[(nd-C)&M]^ra[(nd-D)&M])

static unsigned4 ra[M+1],nd;

/* now a higher-quality generator (rndgen1) used as initializer */
#define rndX2 (int4)113L      /* any not too big prime */
#define rndM2 (int4)19004269L /* max prime < 2^31/rndX2 */
#define rndX1 (int4)78125L    /* so that the cycle mod 2^32 is max (2^29) */

struct irndgen_s { int4 seed1,seed2,fac,tab[1/*var len*/]; };

unsigned4 initIrnd(struct irndgen_s *rndgen) /********************* initIrnd */
{
  int4 *p,r;

  p = rndgen->tab + (rndgen->seed2 = (rndgen->seed2*rndX2)%rndM2) / rndgen->fac;
  r = *p;
  *p = rndgen->seed1 *= rndX1;

  return (unsigned4)r;
}

unsigned4 rndinit(int4 tablen,int4 seed) /************************** rndinit */
{
  struct irndgen_s *rndgen;
  static int initialized=0;
  int i;

  if (tablen<0) return initialized;

  if (seed==0) seed=seed_from_time();

  if (tablen==0) tablen=47;
  if (tablen>1023) tablen=1023;
#ifdef alloc  
  alloc(rndgen,sizeof(struct irndgen_s)+tablen*sizeof(int4));
#else   /*# alloc   */
  rndgen=(struct irndgen_s*)malloc(sizeof(struct irndgen_s)+tablen*sizeof(int4));
#endif /*#!alloc   */

  initialized=1;

  rndgen->seed1 = seed | 1; /* to be odd */
  rndgen->seed2 = abs(seed)%(rndX2-1)+1; /* to be in [1,rndX2-1] */
  rndgen->fac = rndM2/tablen+1;
  
  for (i=0; i<7; i++) rndgen->seed1 *= rndX1;
  
  for (i=0; i<tablen; i++) rndgen->tab[i] = rndgen->seed1 *= rndX1;

  for (nd=0; nd<=M; nd++) ra[nd]=initIrnd(rndgen);
  /* to have all bits sufficiently independent on seed */
  for (nd=0; nd<=M; nd++) ra[nd]^=initIrnd(rndgen)>>7;
  for (nd=0; nd<=M; nd++) ra[nd]^=initIrnd(rndgen)>>13;
  for (nd=0; nd<=M; nd++) ra[nd]^=initIrnd(rndgen)>>19;

  free(rndgen); /* must be the compatible free(0 to either alloc() or malloc() */

  return 0xffffffffU;
}

#undef rndX2
#undef rndX1
#undef rndM2

unsigned4 Irnd(void) /************************************************* Irnd */
{
  return trnd;
}

unsigned4 irnd(unsigned N)
{
#if RNDMETHOD&4
  unsigned4 t0=((long long unsigned)(trnd)*N)>>32;
  return ((long long unsigned)(trnd)*N+t0)>>32;
#else   /*# RNDMETHOD&4 */
  return ((long long unsigned)(trnd)*N)>>32;
#endif   /*#!RNDMETHOD&4 */
}
  
double rnd(void) /****************************************************** rnd */
{
#if RNDMETHOD&1
  /* random number in [0,1) */
  double r=(trnd&0xfffff800U)*(1./4294967296e0);
  return (r+(double)trnd)*(1./4294967296e0);
#else   /*# RNDMETHOD&1 */
  /* random number in [0,1] */
  double r=trnd*(1./4294967296e0);
  return (r+(double)trnd)*(1./4294967296e0);
#endif /*#!RNDMETHOD&1 */
}

double rndcos(void) /************************************************ rndcos */
{
#if RNDMETHOD&2
  /* random number in (-1,1) */
  double r=(trnd&0xfffff800U|0x3ff)*(1./4294967296e0);
  return (r+(double)trnd)*(1./2147483648e0)-1;
#elif RNDMETHOD&1
  /* random number in [-1,1) */
  double r=(trnd&0xfffffc00U)*(1./4294967296e0); 
  return (r+(double)trnd)*(1./2147483648e0)-1;
#else   /*#!RNDMETHOD&2!RNDMETHOD&1 */
  /* random number in [-1,1] */
  double r=trnd*(1./4294967296e0);
  return (r+trnd)*(1./2147483648e0)-1;
#endif /*#!RNDMETHOD&2!RNDMETHOD&1 */
}

#undef A
#undef B
#undef C
#undef D
#undef M

#include "rndetc.c"
