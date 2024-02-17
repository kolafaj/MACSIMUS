/* 
  oriented angle of plane b x r from plane a x r
  = oriented angle of vector bp from vector ap, where ap and bp are
    projections of a and b, respectively, to a plane perpendicular 
    to r (in the direction of r)
      _
      /|
   a /
    /       
   +------------------>
            r          \ b
                        \
                         \|
                         -

*/

#define ORIENTEDANGLE(angle,a,b,r,rr,absr) { \
  double ra=SCAL(r,a)/rr, rb=SCAL(r,b)/rr; \
  vector ap,bp,axr; \
  VVV(ap,=a,-ra*r) VVV(bp,=b,-rb*r) \
  VECT(axr,ap,r) \
  angle=atan2(SCAL(axr,bp),SCAL(ap,bp)*sqrt(rr)); \
  }
