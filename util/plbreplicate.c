/* cc -O2 -Wall -o plbreplicate plbreplicate.c
*/

#include "../gen/include.h"

typedef float vector[3];

int perm=0;

void permute(vector v) /******************************************** permute */
{
  float x;

  switch (perm) {
    case 0: break;
    case 2: x=v[0]; v[0]=v[1]; v[1]=v[2]; v[2]=x; break;
    case 4: x=v[0]; v[0]=v[2]; v[2]=v[1]; v[1]=x; break;
    default: Error("plbreplicate: bad permutation"); }
}

int main(int narg,char **arg) /**************************************** main */
{
  FILE *f1,*f2;
  int ns,i,NX,NY,NZ,ix,iy,iz;
  static float header[2];
  vector *r;
  vector L,newL;
  double XSQUEEZE_Y=0,XSQUEEZE_Z=0,YSQUEEZE_Z=0;

  if (narg<5) {
    fprintf(stderr,"\
Replicate a cell. Call by:\n\
  plbreplicate INPUT.plb OUTPUT.plb NX NY NZ XSQUEZE_Y [XSQUEEZE_Z [YSQUEEZE_Z]]\n\
ARGUMENTS:\n\
  INPUT.plb = plb-file with box info (new format)\n\
  OUTPUT.plb = the same format\n\
  NX NY NZ = replication factors; negative values determine input shuffling:\n\
    NY<0 = cyclic permutation: y -> x -> z -> y (input y becomes x, etc.)\n\
    NZ<0 = cyclic permutation: z -> x -> y -> z (input z becomes x, etc.)\n\
  XSQUEEZE_Y = add XSQUEEZE_y*y/newLx*newLy to x\n\
  XSQUEEZE_Z = add XSQUEEZE*z/newLz*newLx to x\n\
  YSQUEEZE_Z = add Ysqueeze*z/newLz*newLy to y\n\
WARNING:\n\
  not suitable as cook* input if the cell contains more species\n\
See also:\n\
  cook* (load.n[] in data)   molcfg plbbox crystal ice naclcryst\n");
    exit(0); }

  if (!strcmp(arg[1],arg[2])) Error("plbreplicate: INPUT_FILE = OUTPUT_FILE");

  f1=fopen(arg[1],"rb");
  if (f1==NULL) Error("plbreplicate: no INPUT_FILE");

  f2=fopen(arg[2],"wb");
  if (f2==NULL) Error("plbreplicate: cannot open OUTPUT_FILE");

  if (1!=fread(header,sizeof(header),1,f1)) Error("plbreplicate: no header");
  ns=header[0];

  if (header[1]!=-3) Error("plbreplicate: wrong header (variable box required), try plbconv");

  allocarray(r,ns);

  perm=0;
  NX=atoi(arg[3]); if (NX<0) perm+=1,NX=-NX;
  NY=atoi(arg[4]); if (NY<0) perm+=2,NY=-NY;
  NZ=atoi(arg[5]); if (NZ<0) perm+=4,NZ=-NZ;

  if (narg>6) XSQUEEZE_Y=atof(arg[6]);
  if (narg>7) XSQUEEZE_Z=atof(arg[7]);
  if (narg>8) YSQUEEZE_Z=atof(arg[8]);

  header[0]=ns*NX*NY*NZ;
  if (header[0]>16777216.||header[0]<1) Error("plbreplicate: bad number of sites (max. 16777216)");

  fwrite(header,sizeof(header),1,f2);
  if (1!=fread(L,sizeof(L),1,f1)) Error("plbreplicate: plb too short");
  permute(L);

  newL[0]=L[0]*NX;
  newL[1]=L[1]*NY;
  newL[2]=L[2]*NZ;
  if (ns!=fread(r,sizeof(vector),ns,f1)) Error("plbreplicate: plb file truncated");
  loop (ix,0,ns) permute(r[ix]);

  fwrite(newL,sizeof(newL),1,f2);

  loop (ix,0,NX)
    loop (iy,0,NY)
      loop (iz,0,NZ)
        loop (i,0,ns) {
          vector newr;

          newr[2]=L[2]*iz+r[i][2];
          newr[0]=L[0]*ix+r[i][0]+XSQUEEZE_Z*newL[0]/newL[2]*newr[2];
          newr[1]=L[1]*iy+r[i][1]+YSQUEEZE_Z*newL[1]/newL[2]*newr[2]
                                 +XSQUEEZE_Y*newL[0]/newL[1]*newr[0];
          fwrite(newr,sizeof(newr),1,f2); }

  fclose(f2);
  fclose(f1);

  printf("! plbreplicate %s %s %d %d %d\n\
N[0]=%.0f ! the total number of sites\n\
L[0]=%.7f\n\
L[1]=%.7f\n\
L[2]=%.7f\n",
         arg[1],arg[2],NX,NY,NZ,
         header[0],newL[0],newL[1],newL[2]);

  return 0;
}
