/* cc -O2 -o coordn coordn.c -lm
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
char line[128];
int col;
int i;
double g,r0=0,r,r1,n=0,dr;
char *tok;
double rho=1,g0=0;

if (narg<3) {
  fprintf(stderr,"Coordination number from r g(r) file.  Call by:\n\
  coord COLUMN_OF_G DR [RHO G0 CONST] < INFILE > OUTFILE\n\
calculates CONST + RHO integral_0^R [g(r)-G0] dV\n\
default RHO=1, CONST=0\n\
where the integral over dV = 4 PI r^2 dr is replaced by the sum over\n\
dV = 4 PI/3*[(r+DR/2)^3-(r-DR/2)^3]\n");
  exit(0); }

col=atoi(arg[1]);
dr=atof(arg[2]);
if (narg>3) rho=atof(arg[3]);
if (narg>4) g0=atof(arg[4]);
if (narg>5) n=atof(arg[5]);

rho*=4*PI/3;

while (gets(line)) if (!strchr("!#",line[0])) {
  r=atof(strtok(line," \t\n"));
  loop (i,1,col) {
    tok=strtok(NULL," \t\n");
    if (!tok) Error("no data in specified column"); }
  g=atof(tok)-g0;
  r1=r+dr/2;
  n+=rho*(Cub(r1)-Cub(r0))*g;
  printf("%6.3f %8.3f %8.4f\n",r1,n,g); 
  r0=r1; }

return 0;
}
