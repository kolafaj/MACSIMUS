double rndball(double *v) /***************************************** rndball */
/* random vector in 1-ball */
{
  double rr;
  
  do {
    v[0]=rndcos(); v[1]=rndcos(); v[2]=rndcos();
  } while ((rr=v[0]*v[0]+v[1]*v[1]+v[2]*v[2])>=1);
  
 return rr;
}

#if 0
double rnddisk(double *v)
/* random vector in 1-disk */
{
  double rr;
 
  do {
    v[0]=rndcos(); v[1]=rndcos();
  } while ((rr=v[0]*v[0]+v[1]*v[1])>=1);
  
  return rr;
}
#endif /*# 0 */

void rndsphere(double *v) /*************************************** rndsphere */
/* random unit vector in 3D */
{
  double rr;

  do {
    v[0]=rndcos(); v[1]=rndcos(); v[2]=rndcos();
    rr=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
  } while (rr>1 || rr<0.0625);
 
  rr=sqrt(rr);
  v[0]/=rr; v[1]/=rr; v[2]/=rr;
}

#ifndef PI_DOUBLE
#  define PI_DOUBLE 3.14159265358979323846
#endif /*# PI_DOUBLE */

double rndgauss(void) /******************************************** rndgauss */
/* normal Gaussian distribution, exact version */
{
  static unsigned phase;
  static double gauss2;

  if (phase++&1)
    return gauss2;
  else {
    double x=rndcos()*PI_DOUBLE,y=sqrt(-2*log(rnd()));
    gauss2=sin(x)*y;
    return cos(x)*y; }
}
