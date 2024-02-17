/* powers and min,max for double and long double */
/* this file may be #included twice, once with REAL=double 
   and once more with REAL=according to PRECISION
*/   

REAL FN(sqr)(REAL x) { return x*x; }

REAL FN(cub)(REAL x) { return x*x*x; }

REAL FN(pow4)(REAL x) { REAL y=x*x; return y*y; }

REAL FN(pow5)(REAL x) { REAL y=x*x; return y*y*x; }

REAL FN(pow6)(REAL x) { REAL y=x*x*x; return y*y; }

REAL FN(powi)(REAL x,int i)
/* the optimum square-and-multiply algorithm */
{
  int j;
  REAL p;

  if (!i) return 1;
  if (i<0) { i= -i; x=1/x; };
  p=x;
  for (j=1; j<=i; j<<=1);
  j>>=1;
  while ((j>>=1)) { p*=p; if (j & i) p*=x; };

  return p;
}

REAL FN(frac)(REAL x)
{
  return x-(int)x;
} 

#ifndef C99
/* fmin,fmax are part of C99 standard and in some implementations may
   interfere with these definitions */
REAL FN(fmax)(REAL x,REAL y)
{
  if (x<y) return y; else return x;
}

REAL FN(fmin)(REAL x,REAL y)
{
  if (x<y) return x; else return y;
}
#endif
