/* cc -O2 -o lattice lattice.c -lm 
 */
#include "../gen/include.h"

int check110(double r,double x, double y,double z)
{
  if (r==0) return 1;
  r-=1e-5;

  if (fabs(x)+fabs(y)>r) return 0;
  if (fabs(z)+fabs(y)>r) return 0;
  if (fabs(x)+fabs(z)>r) return 0;

  return 1;
}

int check111(double r,double x, double y,double z)
{
  if (r==0) return 1;
  r-=1e-5;

  if (fabs(x)+fabs(y)+fabs(z)>r) return 0;

  return 1;
}

int main(int narg,char **arg)
{
  double x110=0,x111=0;
  int N,n,k,kmax;
  double a;
  int x,y,z,i;
  enum { NDEF, SC, BCC, FCC } lattice=NDEF;

  if (narg<3) {
    fprintf(stderr,"Make cubic lattice or crystal (x,y,z). Call by:\n\
Auto select the lattice matching exactly number of atoms N:\n\
  lattice N a\n\
Make given lattice of k times the basic unit in each direction:\n\
  lattice k a lattice [r110 r111]\n\
Options:\n\
  a       lattice constant\n\
  lattice 1=sc,2=bcc,3=fcc\n\
  r110    distance (from center, in the units of a) to limit the 110 planes\n\
  r111    distance (from center, in the units of a) to limit the 111 planes\n\
");
  exit(0); }
    
  N=atoi(arg[1]);
  a=atof(arg[2]);
  if (narg==3) {
    kmax=(int)(pow((double)N,1.0/3))+2;
    loop (i,1,kmax) {
      if      (N==  i*i*i) { lattice=SC;  k=i; break; }
      else if (N==2*i*i*i) { lattice=BCC; k=i*2; a/=2; break; }
      else if (N==4*i*i*i) { lattice=FCC; k=i*2; a/=2; break; } } }
  else {
    k=N;    
    if (narg>3) lattice=atoi(arg[3]);
    if (narg>4) x110=atof(arg[4]);
    if (narg>5) x111=atof(arg[5]); }

  if (lattice<=NDEF || lattice>FCC) Error("bad lattice or no lattice matching N");

  i=0; n=0;
  loop (x,0,k) loop (y,0,k) loop (z,0,k) {
    if (lattice==BCC) if ( (x&1)!=(y&1) || (x&1)!=(z&1) ) continue;
    if (lattice==FCC) if ( (x+y+z)&1 ) continue;
    i++;
    if (check110(x110,x-(k-1)/2.,y-(k-1)/2.,z-(k-1)/2.) 
     && check111(x111,x-(k-1)/2.,y-(k-1)/2.,z-(k-1)/2.)) {
      n++;
      printf("%8.5f %8.5f %8.5f\n",a*x,a*y,a*z); } }
  fprintf(stderr,"%d atoms (in %d)\n",n,i);
    
  return 0;
}
