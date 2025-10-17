/* cc -O2 -o makedep makedep.c -lm
 */
#include "../gen/include.h"
typedef double real;
#include "../gen/vector3d.h"
#include "../gen/rndgeni.c"

int N;
vector r0[4];
vector r;
double w[4],ww[4];
double f,ff,sw;
double rad;

double dist(double *w) /*********************************************** dist */
{
  int i;
  vector rx;

  VO(rx,=0)
  loop (i,0,N) VV(rx,+=w[i]*r0[i])

  return SQRD(rx,r);
}

int main(int narg,char **arg) /**************************************** main */
{
  char line[256];
  int i;

  if (narg<2) {
    fprintf(stderr,"\
Make dependants for ble-file. Call by:\n\
  makedep N\n\
DATA (on stdin):\n\
  x y z of base atom 1\n\
  ...\n\
  x y z of base atom N\n\
  x y z of atom, weights will be printed\n\
  ...\n\
See also:\n\
  waterdep\n\
");
    exit(0); }

  rndinit(0,0);

  N=atoi(arg[1]);

  getsbufsize=256;
  if (N<2 || N>4) Error("N must be 2,3,4");
  loop (i,0,N) {
    do (gets(line)); while (line[0]=='#');
    sscanf(line,"%lf%lf%lf",r0[i],r0[i]+1,r0[i]+2); }

  while (gets(line)) if (line[0]!='#') {
    sscanf(line,"%lf%lf%lf",r,r+1,r+2);

    loop (i,0,N) w[i]=1./N;
    rad=1;
    f=dist(w);

    while (rad>1e-14) {

      sw=0;
      loop (i,0,N) sw+=ww[i]=w[i]+(rnd()-rnd())*rad;
      loop (i,0,N) ww[i]/=sw;
      ff=dist(ww);

      if (ff<f) {
	rad*=1.7;
	f=ff;
	loop (i,0,N) w[i]=ww[i]; }
      else
	rad*=.99; }

    loop (i,0,N) printf(" %.8f",w[i]);
    _n }

  return 0;
}
