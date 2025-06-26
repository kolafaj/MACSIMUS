/* cc -Wall -O2 -o gro2atm gro2atm.c
 */
#include "../gen/include.h"
#include "../gen/forscanf.c"

int main(int narg,char **arg)
{
  int N;
  char line[1024];
  struct cfg_s {
    int n;
    int i;
    char molid[8];
    char atid[8];
    double x,y,z; 
  } *cfg;
  double Lx,Ly,Lz;
  int i,ns,molN;

  if (narg<2) {
    fprintf(stderr,"\
Convert Gromacs .gro configuration to .atm. Call by:\n\
  gro2atm MOLNUMBER < FILE.gro > FILE.atm\n\
MOLNUMBER >=1 : extract only given molecule (as given in columns 0--4)\n\
MOLNUMBER ==0 : extract whole configuration\n");
    exit(0); }

  molN=atoi(arg[1]);

  gets(line);
  fprintf(stderr,"%s\n",line);

  gets(line);
  N=atoi(line);
  if (N<1 || N>1000000) Error("N out of range");

  allocarray(cfg,N);

  loop (i,0,N) {
/*
Groningen Machine for Chemical Simulation
31910
    1MOL     O1    1   0.162   5.065  10.775  0.2239  0.1997  0.1681
    1MOL     C1    2   0.057   5.017  10.815  0.6393 -0.7165  0.1863
    1pXY     C1    1   0.234   4.000   0.561
nnnnnMMMMMaaaaaIIIIIxxxxxxxxyyyyyyyy
012345678901234567890123456789012345678901234567890123456789
*/
    if (!gets(line)) Error("unexpected EOF");
    if (7!=sforscanf(line,"%5d%5s%5s%5s%8lf%8lf%8lf",
                    &cfg[i].n,
                    cfg[i].molid,
                    cfg[i].atid,
                    &cfg[i].i,
                    &cfg[i].x,&cfg[i].y,&cfg[i].z)) Error("line format"); }

  gets(line);
  sscanf(line,"%lf%lf%lf",&Lx,&Ly,&Lz);
  
  ns=0;
  if (molN==0)
    ns=N;
  else {
    loop (i,0,N) if (cfg[i].n==molN) ns++; }
  if (ns==0) Error("no such molecule");
  
  printf("%d\n %g %g %g box\n",ns,Lx*10,Ly*10,Lz*10);
  
  loop (i,0,N)
    if (!molN || cfg[i].n==molN)
      printf("%6s %g %g %g\n",cfg[i].atid,
             cfg[i].x*10,cfg[i].y*10,cfg[i].z*10);

  return 0;
}
