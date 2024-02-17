/* cc -o esstest esstest.c -lm ; esstest 300

TEST of essential dynamics -- three independent perpendicular motions
to be used with test.3db, where test.che is:

test

HA HA HA
  \| /
   CT
   |
   CT
  /| \
HA HA HA

essential dynamics is called by:
blend -EA -P3 -F11 -r4 test

*/

#include "../gen/include.h"

#define N 8
float hdr[2]={N,0};
float r0[N][3],r[N][3];

int main(int narg,char **arg)
{
FILE *c=fopen("test.3db","rb");
FILE *f=fopen("test.plb","wb");
int i,k,n;

if (narg<2) {
  fprintf(stderr,"Call by:\n\
  esstest #_OF_FRAMES\n");
  exit(0); }

fread(r0,N,12,c);
fwrite(hdr,2,4,f);

loop (n,0,atoi(arg[1])) {
  loop (i,0,N) loop (k,0,3)
    r[i][k]=r0[i][k]
/*.....      +(1-2*(i<4))*(k==2)*0.5*((n%20)-9.5)/10*/
      +(1-2*(i<4))*(k==2)*0.5*cos((double)n)
      +(i==7)*(k==1)*0.1*cos(0.618*n)
/*.....      +(i<2)*(1-2*i)*(k==0)*0.2*cos(1.23457*n)*/
      ;
  fwrite(r,N,12,f); }

fclose(f);
fclose(c);

return 0;
}
