/* cc -o mergeplb -O2 mergeplb.c
(c) J.Kolafa 11/95
updated 8/97: header[1] printed as t0 (10/00: changed into L)
              bug fixed (FROM did not work correctly)
              -BY added
updated 10/00: arbitrary selection added
updated 11/04: more frames allowed
*/

#include "../gen/include.h"

int main(int narg,char **arg)
{
  FILE **f;
  float r[3],hdr[2];
  int iarg,*ns,i,frame=0;
  float L=0;
 
  if (narg<2) {
    fprintf(stderr,"\
Merge several plb files into another plb file, frames synchronously.\n\
Stops when any plb-file reaches EOF.\n\
merge plb files into another one-frame plb file\n\
  mergeplb PLB-FILE [PLB-FILE ...] > MERGED-PLB-FILE\n\
WARNING: L set to maximum L\n\
use molcfg to create mol and gol-files\n\
see also: plbinfo plbmsd cutplb filtplb rotplb plb2diff plbasc\n");
    exit(0); }
  
  hdr[0]=hdr[1]=0;
  
  allocarray(f,narg);
  allocarray(ns,narg);
  
  loop (iarg,1,narg) {
    f[iarg]=fopen(arg[iarg],"rb");
    if (!f[iarg]) Error(arg[iarg]);
    if (2!=fread(r,4,2,f[iarg])) Error(arg[iarg]);
    fprintf(stderr,"%20s: ns=%-5.0f L=%f\n",arg[iarg],r[0],r[1]);
    hdr[0]+=r[0];
    ns[iarg]=r[0];
    Max(hdr[1],r[1]) }

  fprintf(stderr,"              OUTPUT: ns=%-5.0f L=%f\n",hdr[0],hdr[1]);

  fwrite(hdr,4,2,stdout);

  for (;;) {
    loop (iarg,1,narg)
      loop (i,0,ns[iarg]) {
	if (3!=fread(r,4,3,f[iarg])) {
	  fprintf(stderr,"%d frames written\n",frame);
	  return 0; }
	fwrite(r,4,3,stdout); }
	frame++; }

}
