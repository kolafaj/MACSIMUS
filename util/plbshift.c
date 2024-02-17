/* cc -o plbshift -O2 plbshift.c -lm
*/

#include "../gen/include.h"

typedef float vector[3];

int main(int narg,char **arg)
{
  FILE *f1,*f2;
  int NS=0,ns,i,j,k,n=0;
  static float header[2];
  vector *r;
  vector DR;
  vector L;
  vector unit={1,1,1};

  if (narg<5) {
    fprintf(stderr,"\
Shift a plb-file periodically. Call by:\n\
  plbshift INPUT_FILE OUTPUT_FILE DX[L] DY[L] DZ[L] [NS]\n\
where\n\
  DX,DY,DZ = shifts (in AA)\n\
  DXL,DYL,DZL = shifts (in box size)\n\
  NS = number of sites in a molecule [default = 0 = not periodic]\n\
Bug:\n\
  not suitable for mixtures of polyatomic molecules\n\
See also:\n\
  plbstack plbtran plbrot\n");
    exit(0); }

  loop (i,3,6)
    DR[i-3]=atof(arg[i]);
  if (narg>6) NS=atoi(arg[6]);

  if (!strcmp(arg[1],arg[2])) Error("plbshift: INPUT_FILE = OUTPUT_FILE");

  f1=fopen(arg[1],"rb");
  if (f1==NULL) Error("plbshift: no INPUT_FILE");

  if (1!=fread(header,sizeof(header),1,f1)) Error("plbshift: empty input file");
  ns=header[0];

  if (header[1]!=-3)
    Error("plbshift: wrong format: only variable-L plb-file supported (use plbconv)");

  printf("plbshift: %s: of sites=%d  size of 1 frame=%d\n",
	 arg[1],ns,(int)(ns*sizeof(vector)));
  if ((header[0]-ns)!=0 || ns<1 || header[0]>16777216.) Error("plbshift: bad header");

  f2=fopen(arg[2],"wb");
  if (f2==NULL) Error("plbshift: cannot open OUTPUT_FILE");

  fwrite(header,sizeof(header),1,f2);
  allocarray(r,ns);

  for (;;) {
    if (1!=fread(L,sizeof(L),1,f1)) break;
    loop (k,0,3) if (strchr(arg[k+3],'L')) unit[k]=L[k];
    n++;
    fwrite(L,sizeof(L),1,f2);
    if (ns!=fread(r,sizeof(vector),ns,f1)) Error("plbshift: frame truncated");
    loop (i,0,ns) loop (k,0,3) r[i][k]+=DR[k]*unit[k];
    
    if (NS) for (i=0; i<ns; i+=NS) {
      vector minr;
      
      loop (k,0,3) minr[k]=3e33;
      loop (j,i,i+NS) 
        loop (k,0,3) Min(minr[k],r[j][k])
      loop (k,0,3) if (L[k]) minr[k]=L[k]*floor(minr[k]/L[k]);
      loop (j,i,i+NS)
        loop (k,0,3) if (L[k]) r[j][k]-=minr[k]; }
    
    fwrite(r,sizeof(vector),ns,f2); }

  fclose(f2);
  fclose(f1);

  return 0;
}
