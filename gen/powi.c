#ifndef REAL
#define REAL double
#endif

REAL powi(REAL x,int i)
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
