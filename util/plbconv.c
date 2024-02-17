/* cc -o plbconv -O2 plbconv.c
*/

#include "../gen/include.h"

typedef float vector[3];

int main(int narg,char **arg)
{
  FILE *f1,*f2;
  int ns,i,n=0;
  static float header[2];
  vector *r;
  int varL=0,outvarL,isoutL=0;
  vector L,outL;
  float oldL;

  if (narg<4) {
    fprintf(stderr,"Converts old (fixed L) and new (variable L[3]) plb formats\n\
  plbconv INPUT_FILE OUTPUT_FILE new [Lx Ly Lz]\n\
  plbconv INPUT_FILE OUTPUT_FILE old [L]\n\
{old|new} selects the output file format\n\
optional L or Lx,Ly,Lz will replace box sizes\n\
Example:\n\
  plbconv sim.plb 4blend.plb old 0\n\
See also: plbmerge plbcut plbinfo plbfilt plbbox plb2asc\n");
    exit(0); }

  if (!strcmp(arg[1],arg[2])) Error("INPUT_FILE = OUTPUT_FILE");

  f1=fopen(arg[1],"rb");
  if (f1==NULL) Error("no INPUT_FILE");

  f2=fopen(arg[2],"wb");
  if (f2==NULL) Error("cannot open OUTPUT_FILE");

  fread(header,sizeof(header),1,f1);
  ns=header[0];

  if (header[1]<0)
    if (header[1]<0) varL=1;
    else Error("wrong header (L<0, L!=-3)");

  printf("# of sites=%d  size of 1 frame=%d  L=%f (%s)\n",
	 ns,(int)(ns*sizeof(vector)),header[1],varL?"variable box":"old format");
  if ((header[0]-ns)!=0 || ns<1 || header[0]>1e6) Error("bad header");

  if (varL) {
    fread(L,sizeof(L),1,f1);
    oldL=fmax(fmax(L[0],L[1]),L[2]);
    rewind(f1);
    fread(header,sizeof(header),1,f1); }

  outL[0]=outL[1]=outL[2]=L[0]=L[1]=L[2]=header[1];
 
  allocarray(r,ns);

  if (arg[3][0]=='n') outvarL=1;
  else if (arg[3][0]=='o') outvarL=0;
  else Error("3dr parm must be new or old");

  if (narg>4) {
    isoutL=1;
    outL[0]=outL[1]=outL[2]=atof(arg[4]); }    

  if (narg>5) outL[1]=outL[2]=atof(arg[5]);
  if (narg>6) outL[2]=atof(arg[6]);

  if (outvarL) header[1]=-3;
  else {
    if (varL) header[1]=oldL;
    if (isoutL) header[1]=outL[0]; }

  fwrite(header,sizeof(header),1,f2);

  for (;;) {
    if (varL) {
      if (!fread(L,sizeof(vector),1,f1)) break;
      if (!outvarL && !isoutL)
	if (L[0]!=oldL || L[1]!=oldL || L[2]!=oldL)
	  printf("frame %d: L=%g %g %g does not match L=%g in the old plb header\n",n+1,L[0],L[1],L[2],oldL); }

    if (ns!=fread(r,sizeof(vector),ns,f1)) break;
    if (outvarL) {
      if (isoutL) L[0]=outL[0],L[1]=outL[1],L[2]=outL[2];
      fwrite(L,sizeof(vector),1,f2); }
    if (ns!=fwrite(r,sizeof(vector),ns,f2))
      Error("cannot write output plb file");
    n++; }

  fclose(f2);
  fclose(f1);

  printf("%d frames converted\n",n);
  return 0;
}
