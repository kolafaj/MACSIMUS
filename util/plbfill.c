/* cc -O2 -Wall -o plbfill plbfill.c -lm
*/

#include "../gen/include.h"
#define real float
#include "../gen/vector3d.h"

struct mol_s {
  vector r;
  int oldnm,newnm,ns;
  int oldpos,newpos;
} mol[128];

int main(int narg,char **arg) /**************************************** main */
{
  float header[2];
  vector *oldr,*newr;
  vector L;
  int iarg,f,i,j,n,ns,oldns,newns,varL,append=0;

  if (narg<4) {
    fprintf(stderr,"\
Fill or remove molecules from a plb file. Call by:\n\
  plbfill [OPTIONS] NM NEWNM NS [OPTIONS] [NM NEWNM NS ..]\n\
    < INPUT-PLB > OUTPUT-PLB\n\
OPTIONS (before triplets):\n\
  -a  = omit header = append (use >> OUTPUT-PLB)\n\
  -x# = x-position of added atoms\n\
  -y# = y-position of added atoms\n\
  -z# = z-position of added atoms\n\
TRIPLETS:\n\
  NM = number of molecules on input\n\
  NEWNM = new number of molecules\n\
  NS = number of sites in given molecule\n\
Example:\n\
  plbfill -x40 -y80 -z80 123 200 4 498 600 3 < evap601.plb > filled.plb\n\
See also:\n\
  m2m wat2wat plb2plb plbmerge mol2mol ice\n");
    exit(0); }

  n=0;
  oldns=newns=0;

  loop (iarg,1,narg) {
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'a': append++; break;
        case 'x': mol[n/3].r[0]=atof(arg[iarg]+2); break;
        case 'y': mol[n/3].r[1]=atof(arg[iarg]+2); break;
        case 'z': mol[n/3].r[2]=atof(arg[iarg]+2); break;
        default: Error("plbfill: unknown option"); }
    else {
      switch (n%3) {
        case 0:
          mol[n/3].oldnm=atoi(arg[iarg]); break;
        case 1:
          mol[n/3].newnm=atoi(arg[iarg]); break;
        case 2:
          mol[n/3].ns=atoi(arg[iarg]);
          mol[n/3].oldpos=oldns;
          oldns+=mol[n/3].ns*mol[n/3].oldnm;
          mol[n/3].newpos=newns;
          newns+=mol[n/3].ns*mol[n/3].newnm; }
      n++; } }

  if (n%3) Error("plbfill: the number of arguments except -x,-y,-z must be a multiple of 3");
  n/=3;

  fprintf(stderr,"plbfill: nspec=%d oldns=%d newns=%d\n",n,oldns,newns);

  if (!fread(header,4,2,stdin)) Error("plbfill: no data on input");
  ns=header[0];
  if (ns!=oldns) Error("plbfill: the number of sites on input does not match arguments");
  varL=header[1]<0;
  if (!varL) L[0]=L[1]=L[2]=header[1];

  if (!append) {
    header[0]=newns;
    header[1]=-3;
    fwrite(header,4,2,stdout); }
  allocarrayzero(oldr,oldns);
  allocarrayzero(newr,newns);

  for (f=0;;f++) {
    if (varL) {
      if (fread(L,4,3,stdin)!=3) break; }
    i=fread(oldr,12,oldns,stdin);
    if (i==0) break;
    if (i!=ns) Error("plbfill: plb is truncated");
    loop (i,0,n) {
      loop (j,0,mol[i].newnm*mol[i].ns) {
        if (j<mol[i].oldnm*mol[i].ns) VV(newr[j+mol[i].newpos],=oldr[j+mol[i].oldpos])
        else VV(newr[j+mol[i].newpos],=mol[i].r) } }

    fwrite(L,4,3,stdout);
    fwrite(newr,12,newns,stdout);
  }

  fprintf(stderr,"plbfill: %d frames\n",f);

  return 0;
}
