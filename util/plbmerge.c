/* cc -Wall -o plbmerge -O2 plbmerge.c
   part of MACSIMUS
History:
  10/21: bug fixed (-c of the old format)
  03/21: bugs fixed
  02/21: concatenation changed (files not opened at once)
  05/08: append option added
  06/05: new plb format (variable L[3])
  11/04: more frames allowed
  10/00: arbitrary selection added
  08/97: header[1] printed as t0 (10/00: changed into L)
         bug fixed (FROM did not work correctly)
         -BY added
  11/95: first version
*/

#include "../gen/include.h"

typedef float fvector[3];

int main(int narg,char **arg) /**************************************** main */
{
  float hdr[2];
  fvector L3,r;
  int i,iarg,iargfile=0,NS=0,CAT=0;
  unsigned NF=0xffffffffu,frame=0;

  if (narg<2) {
    fprintf(stderr,"\
Merge binary playback files.  Call by:\n\
  plbmerge [OPTIONS] PLB-FILE [PLB-FILE ...] > MERGED-PLB-FILE\n\
OPTIONs:\n\
  -m\tmerge frames synchronously [default]:\n\
\t- stops when any plb-file reaches EOF\n\
\t- output box size L is set to MAXIMUM of all plb-files\n\
\t- hint: use molcfg to create merged mol- and gol-files\n\
  -cNS\ttruncate playback files to NS sites and concatenate\n\
  -c\tas concatenate, the first file is the shortest (or all the same length)\n\
  -fNF\tmerge/concatenate at most NF frames from each file [default=all]\n\
The output plb file contains box info (cf. plbconv or mergeplb for old format)\n\
Examples:\n\
  plbmerge -c run1.plb run2.plb > runs1+2.plb\n\
  plbmerge -f1 sim.p00 sim.p01 > sim.plb\n\
See also:\n\
  plbinfo plbmsd plbcut plbfilt plbtran plbrot plb2diff plbasc\n");
    exit(0); }

  hdr[1]=-3; /* always variable-box format on output */

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-') {
      switch (arg[iarg][1]) {
        case 'c':
          CAT=1;
          NS=atoi(arg[iarg]+2);
          break;
        case 'm': 
          CAT=0;
          break;
        case 'f':
          NF=atoi(arg[iarg]+2);
          break;
        default: 
          Error("unknown option"); } }
    else {
      iargfile=iarg;
      break; }

  if (iargfile==0) Error("no file");
  
  if (CAT) {
    /* concatenate */
    fvector *R;
    int varL;

    fprintf(stderr,"concatenating plb-files (serially):\n");
  
    loop (iarg,iargfile,narg) {
      FILE *f=fopen(arg[iarg],"rb");

      fprintf(stderr,"%20s:",arg[iarg]);
      
      if (!f) Error(arg[iarg]);

      if (2!=fread(r,4,2,f)) Error(arg[iarg]);
      fprintf(stderr," ns=%-5.0f L/key=%f\n",r[0],r[1]);
      varL=r[1]<0;
      L3[0]=L3[1]=L3[2]=r[1];

      if (iarg==iargfile) {
        if (NS==0) NS=r[0];
        allocarray(R,NS);
        if (r[0]<NS) Error("1st file shorter than -cNS"); }

      if (r[0]<NS) Error("file shorter than -cNS");
      if (r[0]>NS) fprintf(stderr,"*** file will be truncated\n");

      if (iarg==iargfile) {
        /* write header */
        hdr[0]=NS;
        fwrite(hdr,4,2,stdout); }

      loop (frame,0,NF) {
        if (varL) if (3!=fread(L3,4,3,f)) break;
        if (NS!=fread(R,sizeof(fvector),NS,f)) break;
        fwrite(L3,4,3,stdout);
        fwrite(R,sizeof(fvector),NS,stdout);
        loop (i,NS,(int)r[0]) fread(R[0],4,3,f); }
    
      fclose (f); } }

  else {
    /* merge files synchronously */
    FILE **f;
    int *ns,*varL,nsmin=0x7fffffff;
    float L=0,L3[3],*LL;

    hdr[0]=0;

    allocarray(f,narg);
    allocarray(ns,narg);
    allocarray(varL,narg);
    allocarray(LL,narg);

    fprintf(stderr,"merging plb-files synchronously (in parallel):\n");

    loop (iarg,iargfile,narg) {
      f[iarg]=fopen(arg[iarg],"rb");
      if (!f[iarg]) Error(arg[iarg]);
      if (2!=fread(r,4,2,f[iarg])) Error(arg[iarg]);
      fprintf(stderr,"%20s: ns=%-5.0f L/key=%f\n",arg[iarg],r[0],r[1]);
      varL[iarg]=r[1]<0;
      Max(L,r[1])
        LL[iarg]=r[1];
      Min(nsmin,(int)r[0])
        hdr[0]+=r[0];
      ns[iarg]=r[0]; }

    fprintf(stderr,"              OUTPUT: ns=%-5.0f key=%g (variable box)\n",hdr[0],hdr[1]);

    fwrite(hdr,4,2,stdout);

    loop (frame,0,NF) {
      L3[0]=L3[1]=L3[2]=L;
      loop (iarg,iargfile,narg)
        if (varL[iarg]) {
          if (3!=fread(r,4,3,f[iarg])) {
            fprintf(stderr,"%d frame(s) written\n",frame);
            return 0; }
          Max(L3[0],r[0])
          Max(L3[1],r[1])
          Max(L3[2],r[2]) }

      fwrite(L3,4,3,stdout);

      loop (iarg,iargfile,narg) {
        loop (i,0,ns[iarg]) {
          if (3!=fread(r,4,3,f[iarg])) {
            fprintf(stderr,"EOF on %s: %d frame(s) written\n",arg[iarg],frame);
            return 0; }
          fwrite(r,4,3,stdout); } }
    } }

  return 0;
}
