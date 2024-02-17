/*
  Gauss 2 point formula for numerical integration, error ~ n^-4
  4n function calls
*/

REAL Gauss4(REAL (*f)(REAL),REAL a,REAL b,int n)
{
REAL h=(b-a)/(n*2),sum=0,q=sqrt((REAL)1/3);
int i;

if (n<=0) return 0;

for (i=n*2-1; i>0; i-=2) sum += (*f)(a+(i+q)*h) + (*f)(a+(i-q)*h);

return h*sum;
}
