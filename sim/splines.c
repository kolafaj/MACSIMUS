/* Spline preparation of two functions (eru and erd) at once, to be used in
   the (r-space part of) the electrostatic forces.
   The electrostatic energy/force between charges qi,qj is then
     uij = qi*qj*eru(rr)
     fij = qi*qj*eru(rr)*rij (vector)
   where rij=rj-ri (vector) and rr=rij^2.

   Robust code replacing the old versions, still backward compatible.

   Needed:
     erfcs(i,&eru,&erd) should return eru(i*Erfc.h), erd(i*Erfc.h)
   where Erfc.h=1/grid/SUBGRID
   This function is be called with non-increasing i's 
   (therefore it remains compatible with the old code using 
   numerical integration from long i).

   SPLINE values:
   2  A+x*(B+x*C)       - fastest, least accurate
   -2 A+B/(x+C)         - good monotonously decreasing functions as erfc
   3  A+x*(B+x*(C+x*D)) - good for Gaussian charges + Ewald
*/

static void makesplines(int ito,double alpha) /***************** makesplines */
/*
   ito = max. spline interval index
   alpha = scaling factor (with erfc.c only: USE 1 OTHERWISE)
*/
{
  int i,n;
#if SPLINE==3  
  /* long double improves spline precision for SPLINE=3 */
  long
#endif
    double Du,Dd,x,h,xh, Duh,Ddh;
  ertab_p p;
  struct f_e {
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
    /* See splines.mw and maple2spline.sh
       In V3.6w and older, numerically unstable formulas were used which made
       problems for el.grid>64 and/or demanding applications like virtual
       volume change method.  For standard MD, el.grid=32 worked fine.      
     */
    p->Au=-Du*x*h*(xh+x)*xh+Sqr(x)*h*(Du-Duh)*xh + (f[i].u*Cub(h)+(f[i+SUBGRID].u-f[i].u)*Sqr(x)*(h + 2*xh));
    p->Bu=Sqr(h)*(x+xh)*Du + h*x*(x+2*xh)*(Duh+Du) + 6*x*xh*(f[i].u-f[i+SUBGRID].u);
    p->Cu=-Sqr(h)*Du + (-Sqr(h) - 3*h*x)*(Duh+Du) + 3*(xh+x)*(f[i+SUBGRID].u-f[i].u);
    p->Du=2*(f[i].u-f[i+SUBGRID].u)+h*(Du+Duh);

    p->Ad=-Dd*x*h*(xh+x)*xh+Sqr(x)*h*(Dd-Ddh)*xh + (f[i].d*Cub(h)+(f[i+SUBGRID].d-f[i].d)*Sqr(x)*(h + 2*xh));
    p->Bd=Sqr(h)*(x+xh)*Dd + h*x*(x+2*xh)*(Ddh+Dd) + 6*x*xh*(f[i].d-f[i+SUBGRID].d);
    p->Cd=-Sqr(h)*Dd + (-Sqr(h) - 3*h*x)*(Ddh+Dd) + 3*(xh+x)*(f[i+SUBGRID].d-f[i].d);
    p->Dd=2*(f[i].d-f[i+SUBGRID].d)+h*(Dd+Ddh);

    p->Au *= alpha/Cub(h);
    p->Ad *= Cub(alpha)/Cub(h);
    p->Bu *= Cub(alpha)/Cub(h);
    p->Bd *= Pow5(alpha)/Cub(h);
    p->Cu *= Pow5(alpha)/Cub(h);
    p->Cd *= Pow4(alpha)*Cub(alpha)/Cub(h);
    p->Du *= Pow4(alpha)*Cub(alpha)/Cub(h);
    p->Dd *= Pow4(alpha)*Cub(alpha)*Sqr(alpha)/Cub(h);

#else /*#!SPLINE==2!SPLINE==-2!SPLINE==3 */
#  error "unsupported SPLINE type, expected: -2=hyperbolic, 2=quadratic, 3=cubic"
#endif /*#!SPLINE==2!SPLINE==-2!SPLINE==3 */
  }

  Erfc.sgrid *= Sqr(alpha);

  prt("Splines calculated with compile-time switches: SPLINE=%d, SUBGRID=%d",SPLINE,SUBGRID);

  free(f);
} /* makesplines */
