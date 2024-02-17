/* Spline preparation of two functions (eru and erd) at once.
   To be used for electrostatic/Ewald.
   Robust code replacing the old versions, still backward compatible
   Needs function
     erfcs(i,&eru,&erd) -> eru(i*Erfc.h), erd(i*Erfc.h)
   which is called with non-increasing i's (therefore compatible with
   the old code using numerical integration from long i).
   SUBGRID accepted (as for old numerical integration code)

SPLINE==2: A+x*(B+x*C) - fastest, least accurate
SPLINE==-2: A+B/(x+C) - good for r-space Ewald, cannot use for nonmonotonous functions
SPLINE==3: A+x*(B+x*(C+x*D)) - good for Gaussian charges + Ewald
*/

static void makesplines(int ito,double alpha) /********* makesplines */
/*
   ito in units of step
   alpha = scaling factor (only old module erfc.c)
*/
{
  int i,n;
  double Du,Dd,x,h,xh;
  double Duh,Ddh;
  ertab_p p;
  struct {
    double u,d;
  } *f;

  h=1/Erfc.sgrid;
  Erfc.h=h/SUBGRID;
  n=ito*SUBGRID;

  /* precalculate table of function values at grid points */
  allocarray(f,n+3); /* doublecheck possible overflow with SUBGRID>1 */
  for (i=n+2; i>=0; i--) erfcs(i,&f[i].u,&f[i].d);

  // loopto(i,-2,4) prt("%g %g %g",i*Erfc.h,f[i].u,f[i].d);

#  if SUBGRID==1
  /* 4th order formula - needed with SUBGRID=1 */
  Du=(-f[n+2].u+8*f[n+1].u-8*f[n-1].u+f[n-2].u)/(12*Erfc.h);
  Dd=(-f[n+2].d+8*f[n+1].d-8*f[n-1].d+f[n-2].d)/(12*Erfc.h);
  /* NB: quadratic splines used to have Du=Dd=0 */
#  else /*# SUBGRID==1 */
  /* derivatives by 2nd order formula, enough with SUBGRID=4 (original version) */
  Du=(f[n+1].u-f[n-1].u)/(2*Erfc.h);
  Dd=(f[n+1].d-f[n-1].d)/(2*Erfc.h);
#  endif /*#!SUBGRID==1 */

  x=n*Erfc.h;

  while (ito--) {
    Duh=Du;
    Ddh=Dd;
    xh=x; /* x+h (h=SUBGRID*Erfc.h) */
    i=ito*SUBGRID;
    x=i*Erfc.h;

    p=Erfc.tab+ito;

#if SPLINE==2
    p->Cu=(h*Du+f[i].u-f[i+SUBGRID].u)/Sqr(h);
    p->Bu=Du-2*xh*p->Cu;
    p->Au=f[i].u-x*p->Bu-x*x*p->Cu;
    Du= p->Bu+2*p->Cu*x;

    p->Cd=(h*Dd+f[i].d-f[i+SUBGRID].d)/Sqr(h);
    p->Bd=Dd-2*xh*p->Cd;
    p->Ad=f[i].d-x*p->Bd-x*x*p->Cd;
    Dd= p->Bd+2*p->Cd*x;

    p->Au *= alpha;       p->Ad *= Cub(alpha);
    p->Bu *= Cub(alpha);  p->Bd *= Pow5(alpha);
    p->Cu *= Pow5(alpha); p->Cd *= Pow4(alpha)*Cub(alpha);

#elif SPLINE==-2
    p->Cu=(x*(f[i+SUBGRID].u-f[i].u)-xh*Du*h) / (Du*h-(f[i+SUBGRID].u-f[i].u));
    p->Bu= -Du*Sqr(p->Cu+xh);
    p->Au=f[i+SUBGRID].u+Du*(p->Cu+xh);
    Du= -p->Bu/Sqr(p->Cu+x);

    p->Cd=(x*(f[i+SUBGRID].d-f[i].d)-xh*Dd*h) / (Dd*h-(f[i+SUBGRID].d-f[i].d));
    p->Bd= -Dd*Sqr(p->Cd+xh);
    p->Ad=f[i+SUBGRID].d+Dd*(p->Cd+xh);
    Dd= -p->Bd/Sqr(p->Cd+x);

    p->Au *= alpha;     p->Ad *= Cub(alpha);
    p->Bu /= alpha;     p->Bd *= alpha;
    p->Cu /= Sqr(alpha);p->Cd /= Sqr(alpha);

#elif SPLINE==3
#  if SUBGRID>1
    /* derivatives by 2nd order formula, good with SUBGRID==4 */
    if (i>0) {
      Du=(f[i+1].u-f[i-1].u)/(2*Erfc.h);
      Dd=(f[i+1].d-f[i-1].d)/(2*Erfc.h); }
    else {
      Du=(-1.5*f[i].u+2*f[i+1].u-0.5*f[i+2].u)/Erfc.h;
      Dd=(-1.5*f[i].d+2*f[i+1].d-0.5*f[i+2].d)/Erfc.h; }
#  else /*# SUBGRID>1 */
    /* 4th order formula - needed with SUBGRID=1 (NB: 6th order only slightly better) */
    if (i>1) {
      Du=(-f[i+2].u+8*f[i+1].u-8*f[i-1].u+f[i-2].u)/(12*Erfc.h);
      Dd=(-f[i+2].d+8*f[i+1].d-8*f[i-1].d+f[i-2].d)/(12*Erfc.h); }
    else if (i>0) {
      // 4th order:  a = -1/4, b = -5/6, c = 3/2, d = -1/2, e = 1/12
      Du=(-1./4*f[i-1].u-5./6*f[i].u+3./2*f[i+1].u-1./2*f[i+2].u+1./12*f[i+3].u)/Erfc.h;
      Dd=(-1./4*f[i-1].d-5./6*f[i].d+3./2*f[i+1].d-1./2*f[i+2].d+1./12*f[i+3].d)/Erfc.h; }
    else {
      Du=(-25./12*f[i].u+4*f[i+1].u-3*f[i+2].u+4./3*f[i+3].u-1./4*f[i+4].u)/Erfc.h;
      Dd=(-25./12*f[i].d+4*f[i+1].d-3*f[i+2].d+4./3*f[i+3].d-1./4*f[i+4].d)/Erfc.h; }
#  endif /*#!SUBGRID>1 */
    p->Au=(h-2*x)*Sqr(xh)*f[i].u + Sqr(x)*(h+2*xh)*f[i+SUBGRID].u - x*h*Sqr(xh)*Du -xh*h*Sqr(x)*Duh;
    p->Bu=6*x*xh*(f[i].u-f[i+SUBGRID].u) +  h*(2*x+xh)*xh*Du +h*x*(2*xh+x)*Duh;
    p->Cu=-3*(x+xh)*(f[i].u-f[i+SUBGRID].u) -  h*(x+2*xh)*Du -h*(2*x+xh)*Duh;
    p->Du=2*(f[i].u-f[i+SUBGRID].u)+h*(Du+Duh);

    p->Ad=(h-2*x)*Sqr(xh)*f[i].d + Sqr(x)*(h+2*xh)*f[i+SUBGRID].d - x*h*Sqr(xh)*Dd -xh*h*Sqr(x)*Ddh;
    p->Bd=6*x*xh*(f[i].d-f[i+SUBGRID].d) +  h*(2*x+xh)*xh*Dd +h*x*(2*xh+x)*Ddh;
    p->Cd=-3*(x+xh)*(f[i].d-f[i+SUBGRID].d) -  h*(x+2*xh)*Dd -h*(2*x+xh)*Ddh;
    p->Dd=2*(f[i].d-f[i+SUBGRID].d)+h*(Dd+Ddh);

    p->Au *= alpha/Cub(h);       p->Ad *= Cub(alpha)/Cub(h);
    p->Bu *= Cub(alpha)/Cub(h);  p->Bd *= Pow5(alpha)/Cub(h);
    p->Cu *= Pow5(alpha)/Cub(h); p->Cd *= Pow4(alpha)*Cub(alpha)/Cub(h);
    p->Du *= Pow4(alpha)*Cub(alpha)/Cub(h);
    p->Dd *= Pow4(alpha)*Cub(alpha)*Sqr(alpha)/Cub(h);

#else /*#!SPLINE==2!SPLINE==-2!SPLINE==3 */
#  error "unsupported SPLINE type, expected: -2=hyperbolic, 2=quadratic, 3=cubic"
#endif /*#!SPLINE==2!SPLINE==-2!SPLINE==3 */
  }

  Erfc.sgrid *= Sqr(alpha);

  free(f);
} /* eerfcs */
