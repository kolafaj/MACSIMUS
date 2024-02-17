/* log(x+1), without loss of precision around x=0 */

static double logseries(double x)
{
int n=2;
double s=0,xx=-x*x,q=xx;

do {
  q*=-x;
  n++;
  s+=q/n;
  } while (fabs(q)>1e-17);
 
return (s+xx/2)+x;
}

double log1(double x)
{
if (fabs(x)<0.25) return logseries(x);
else return log(1+x);
}
