/* cc -O2 -o plbframe plbframe.c
 */

#include "../gen/include.h"

typedef float fvector[3];

int main(int narg,char **arg)
{
  FILE *in,*out;
  float hdr[2];
  fvector *r0,*r1,*sw,L;
  int f0=0,f1=0,last=0;
  double frame=1;
  int ns,nr,i,varL;

  if (narg<3) {
    fprintf(stderr,"\
Extract/interpolate one frame from a plb file. Call by:\n\
  plbframe {INPUT.plb|-} {OUTPUT.plb|-} [FRAME]\n\
where\n\
  - = stdin or stdout\n\
  FRAME = frame (1st=1), fractional FRAME is interpolated\n\
  missing FRAME = last frame\n\
See also:\n\
  plbcut plbinfo plb2asc asc2plb plbconv plb2plb\n");
    exit(0); }

  if (!strcmp(arg[1],arg[2])) Error("plbframe: INPUT.plb=OUTPUT.plb");

  if (strcmp(arg[1],"-")) in=fopen(arg[1],"rb");
  else in=stdin;
  if (!in) Error(arg[1]);

  if (strcmp(arg[2],"-")) out=fopen(arg[2],"wb");
  else out=stdout;
  if (!out) Error(arg[2]);

  if (narg>3) frame=atof(arg[3]);
  else last++;

  if (2!=fread(hdr,4,2,in)) Error("plbframe: no header");
  varL=hdr[1]<0;
  ns=hdr[0]+varL;
  if (!varL) L[0]=L[1]=L[2]=hdr[1];

  fwrite(hdr,4,2,out);

  allocarray(r0,ns);

  if (last) {
    for (;;) {
      nr=fread(r0,sizeof(fvector),ns,in);
      if (nr==0) {
        fwrite(r0,sizeof(fvector),ns,out);
        break; }
      if (nr<ns) Error("plbframe: (last) file truncated"); } }
  else {
    allocarray(r1,ns);

    for (;;) {
      nr=fread(r0,sizeof(fvector),ns,in);
      //      fprintf(stderr,"nr=%d ns=%d\n",nr,ns);
      if (nr<ns) Error("plbframe: file truncated");
      else f0++;

      if (f0==frame) {
        fwrite(r0,sizeof(fvector),ns,out);
        break; }        
      
      if (f1 && frame>=f1 && frame<=f0 || nr==0) {
        if (!f1) {
          /* this may happen for f0=1 only */
          if (frame!=f0) Error("plbframe: at least 2 frames needed to interpolate/extrapolate");
          else fwrite(r0,sizeof(fvector),ns,out); }
        else {
          /* interpolate */
          loop (i,0,ns) {
            fvector av;
            int k;

            if (varL && i==0) loop (k,0,3) L[k]=r0[0][k];

            loop (k,0,3) {
              if (L[k]) {
                while (r0[i][k]-r1[i][k]>L[k]/2) r0[i][k]-=L[k];
                while (r0[i][k]-r1[i][k]<-L[k]/2) r0[i][k]+=L[k]; }
              av[k]=(frame-f1)*r0[i][k]+(f0-frame)*r1[i][k]; }
            fwrite(av,sizeof(fvector),1,out); } }
        break; }  
      f1=f0;
      sw=r1; r1=r0; r0=sw; } }

  fclose(out);

  return 0;
}
