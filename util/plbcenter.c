/* cc -O2 -Wall -o plbcenter plbcenter.c -lm
*/
#include "../gen/include.h"
typedef double real;
#include "../gen/vector3d.h"

int main(int narg,char **arg)
{
  char *fn1=NULL,*fn2=NULL;
  int iarg,i,k,n,f=0,ns,varL=0,box=0,cluster=0;
  float h[2],L[3];
  float (*r)[3];
  double cm[3];
  FILE *plb1,*plb2;

  if (narg<2) {
    fprintf(stderr,"\
Center a configuration (in a plb-file). Call by:\n\
  plbcenter [OPTIONS] INPUT.plb [OUTPUT.plb]\n\
OPTIONS:\n\
  -b       center to the box center (and keep the box)\n\
           default: omit box info (put in vacuum) and center to (0,0,0)\n\
  -l       try to move the configuration (cluster) to minimize split over\n\
           periodic boundary conditions\n\
Missing OUTPUT.plb -> stdout\n\
See also:\n\
  plbcluster plbinfo plbtran plbbox\n");
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-') switch (arg[iarg][1]) {
      case 'b': box++;  break;
      case 'l': cluster++;  break;
      default: Error("plbcenter: bad option"); }
    else {
      if (fn1) fn2=arg[iarg];
      else fn1=arg[iarg]; }

  if (fn1 && fn2 && !strcmp(fn1,fn2)) Error("plbcenter: infile=outfile");
  if (fn1) plb1=fopen(fn1,"rb"); else plb1=stdin;
  if (fn2) plb2=fopen(fn2,"wb"); else plb2=stdout;

  if (fread(h,4,2,plb1)!=2) Error("plbcenter: too short plb file");
  varL=h[1]<0;
  VO(L,=h[1])
  ns=h[0];
  if (ns<1 || ns>16777216) Error("plbcenter: wrong ns - endian?");
  allocarray(r,ns);

  h[1]=-3;
  fwrite(h,4,2,plb2);

  for (f=1;;f++) {

    if (varL) if (3!=fread(L,4,3,plb1)) goto end;

    if (ns!=fread(r,12,ns,plb1)) goto end;

    if (cluster)
      loop (n,1,ns) {
        VO(cm,=0)
        loop (i,0,n) VV(cm,+=r[i])
        VO(cm,/=n)
        loop (k,0,3) {
          while (r[n][k]>cm[k]+L[k]/2) r[n][k]-=L[k];
          while (r[n][k]<cm[k]-L[k]/2) r[n][k]+=L[k]; } }

    VO(cm,=0)
    loop (i,0,ns) VV(cm,+=r[i])
    VO(cm,/=ns)

    if (box) VV(cm,+=0.5*L)
    else VO(L,=0) /* written as not periodic */
    fwrite(L,4,3,plb2);

    loop (i,0,ns) VV(r[i],-=cm)

    if (ns!=fwrite(r,12,ns,plb2)) goto end;
  }

 end:
  fclose(plb1);
  fclose(plb2);

  return 0;
}
