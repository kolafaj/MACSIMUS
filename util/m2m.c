/* cc -O2 -Wall -o m2m m2m.c -lm
*/

#include "../gen/include.h"
#define real float
#include "../gen/vector3d.h"

int main(int narg,char **arg) /**************************************** main */
{
  static float header[2];
  vector *r,*ro;
  vector L;
  int l1,l2;
  int f,i,io,m,n,no,ns,nso;
  int NM=0,HEAD=0,TAIL,COPYHEAD=1,COPYTAIL=1,iarg,varL;
  double Q=1;
  int *ren;

  if (narg<3) {
    fprintf(stderr,"\
Convert configurations by various models. Call by:\n\
  m2m INPUT-SITES OUTPUT-SITES [OPTIONS] < INPUT-PLB > OUTPUT-PLB\n\
SITES = string of chars denoting sites\n\
BUG: in case of multiple chars on input, the first one applies\n\
OPTIONS:\n\
  -qQ    scale coordinates+box by Q [default=1]\n\
  -cN    copy the first N sites unchanged\n\
  -sN    skip the first N sites\n\
  -mM    process M molecules (as denoted by SITES), copy the rest\n\
  -oM    process M molecules (as denoted by SITES), omit the rest\n\
Example:\n\
  m2m -c64 HOMX HOX < TIP4Pbrine.plb > SPCEbrine.plb\n\
See also:\n\
  wat2wat plb2plb plbmerge mol2mol ice\n");
    exit(0); }

  loop (iarg,3,narg) switch (arg[iarg][1]) {
    case 'q': Q=atof(arg[iarg]+2); break;
    case 's': COPYHEAD=0;
    case 'c': HEAD=atoi(arg[iarg]+2); break;
    case 'o': COPYTAIL=0;
    case 'm': NM=atoi(arg[iarg]+2); break;
    default: Error("m2m: bad option"); }

  if (!fread(header,4,2,stdin)) Error("m2m: no data on input");
  ns=header[0];
  varL=header[1]<0;
  if (!varL) L[0]=L[1]=L[2]=header[1]*Q;

  l1=strlen(arg[1]);
  l2=strlen(arg[2]);

  if (NM==0) {
    NM=(ns-HEAD)/l1;
    if (HEAD+NM*l1!=ns)
      Error("m2m: inconsistent input number of sites (not integer multiple)"); }

  TAIL=ns-HEAD-NM*l1;
  if (TAIL<0) Error("m2m: too many sites on input is requested");
  TAIL*=COPYTAIL;
  nso=HEAD*COPYHEAD+NM*l2+TAIL;

  fprintf(stderr,"m2m: ns_in=%d ns_out=%d nm=%d head=%d tail=%d\n",ns,nso,NM,HEAD,TAIL);
  
  header[0]=nso;
  header[1]=-3;
  fwrite(header,4,2,stdout);

  allocarray(r,ns);
  allocarray(ro,nso);

  allocarrayzero(ren,l2);

  loop (io,0,l2) {
    int c=arg[2][io];
    loop (i,0,l1)
      if (c==arg[1][i]) {
        ren[io]=i;
        break; } }

  for (f=0;;f++) {
    if (varL) {
      if (fread(L,4,3,stdin)!=3) break;
      VO(L,*=Q) }
    i=fread(r,12,ns,stdin);
    loop (i,0,ns) VO(r[i],*=Q)
    if (i==0) break;
    if (i!=ns) Error("plb is truncated");

    no=0;
    if (HEAD) {
      if (COPYHEAD) {
        loop (i,0,HEAD) VV(ro[i],=r[i])
        no=HEAD; } }
    n=HEAD;

    loop (m,0,NM) {

      loop (i,0,l2) {
        if (no+i>=nso) Error("m2m: internal: overflow on output");
        if (n+ren[i]>ns) Error("m2m: internal: overflow on input");
        VV(ro[no+i],=r[n+ren[i]]) }

      no+=l2;
      n+=l1; }

    if (TAIL) {
      if (no+TAIL!=nso) Error("m2m: internal: tail mismatch on output");
      if (n+TAIL!=ns) Error("m2m: internal: tail mismatch on intput");
      loop (i,0,TAIL) VV(ro[no+i],=r[n+i]) }  
    fwrite(L,4,3,stdout);
    fwrite(ro,12,nso,stdout);
  }

  fprintf(stderr,"m2m: %d frames\n",f);
  
  return 0;
}
