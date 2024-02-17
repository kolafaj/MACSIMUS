/***
    Exact but slow variant of erud without direct elst term, used for
    "Ewald intramolecular correction" for excluded 1-2 and 1-3 charges,
    intramolecular charges close together (like excluded 1-2 and 1-3
    charges); exacterud_sqrt_1 for POLAR: Drude charges.  
    - value of eru minus direct elst term returned, 
    - erd minus direct elst term passed as *erd.

    exacterud_sqrt(rr) = -erf(alpha*sqrt(rr))/sqrt(rr)
                       = alpha*[e(sqrt(rr)*alpha)-1/(sqrt(rr)*alpha)]
    exacterud_sqrt_1(rr) = exacterud_sqrt(rr) + alpha*2/SQRTPI
    erd = -erf(alpha*sqrt(rr))/sqrt(rr)/rr+alpha*2/SQRTPI/rr*exp(-alpha^2*rr)
        = -alpha^3*[e'(sqrt(rr)*alpha)+1/(sqrt(rr)*alpha)^2]/(sqrt(rr)*alpha)

    see above for e, e'

    SUBGRID=1: uses library erf() for larger args, Taylor for smaller;
               more accurate than below
    SUBGRID>1: uses the Taylor series, efficient for small rr only
               (which is the case of the Ewald correction indeed).
***/

#  if SUBGRID==1

double exacterud_sqrt(double rr,erreal *erd)
{
  double err=Erfc.alphaq*rr;

  if (err<0.01) {
    int n=1;
    double su=0,sd=0,q,r;

    q=1;

    for (;;) {
      r=q/(2*n+1);
      sd+=r*n;
      su+=r;
      if (fabs(q)<1e-16) break;
      q*=-err/++n; }

    *erd=Erfc.B*sd;

    return Erfc.A*(1-su*err); }
  else {

    double r=sqrt(rr);
    double y=-erf(Erfc.alpha*r)/r; 

    *erd=(y-Erfc.A*exp(-err))/rr;

    return y; }
}

double exacterud_sqrt_1(double rr,erreal *erd)
{
  double err=Erfc.alphaq*rr;

  if (err<0.01) {
    int n=1;
    double su=0,sd=0,q,r;

    q=1;

    for (;;) {
      r=q/(2*n+1);
      sd+=r*n;
      su+=r;
      if (fabs(q)<1e-17) break;
    q*=-err/++n; }

    *erd=Erfc.B*sd; 

    return Erfc.A*(-su*err); }
  else {
    double r=sqrt(rr);
    double y=-erf(Erfc.alpha*r)/r; 
    
    *erd=(y-Erfc.A*exp(-err))/rr;

  return y-Erfc.A; }
}

#  else /*# SUBGRID==1 */

double exacterud_sqrt(double rr,erreal *erd)
{
  int n=1;
  double su=0,sd=0,q,r;

  rr*=Erfc.alphaq;

  q=1;

  for (;;) {
    r=q/(2*n+1);
    sd+=r*n;
    su+=r;
    if (fabs(q)<erexact_eps) break;
    q*=-rr/++n; }

  *erd=Erfc.B*sd;

  return Erfc.A*(1-su*rr);
}

double exacterud_sqrt_1(double rr,erreal *erd)
{
  int n=1;
  double su=0,sd=0,q,r;

  rr*=Erfc.alphaq;

  q=1;

  for (;;) {
    r=q/(2*n+1);
    sd+=r*n;
    su+=r;
    if (fabs(q)<erexact_eps) break;
    q*=-rr/++n; }

  *erd=Erfc.B*sd;

  return Erfc.A*(-su*rr);
}
#  endif /*#!SUBGRID==1 */
