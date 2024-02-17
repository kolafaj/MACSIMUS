/* make plbmatch
*/

#include "ground.h"

/* common with show.c: */
typedef float fvector[3];

double drot[3][3]={
  {1,0,0},
  {0,1,0},
  {0,0,1} };

fvector rot[3]={
  {1,0,0},
  {0,1,0},
  {0,0,1} };

static fvector center;

int bodyframe; /* rotations WRT fixed frame/body frame */

void rotate(int axis,double angle) /******************************** rotate */
{
  double sa=sin(angle), ca=cos(angle);
  double o[3][3], rot2[3][3];
  int i=(axis+1)%3,j=(i+1)%3,k=(j+1)%3,a,b;

  if (bodyframe) sa=-sa;

  o[i][i]=ca; o[i][j]=-sa; o[i][k]=0;
  o[j][i]=sa; o[j][j]= ca; o[j][k]=0;
  o[k][i]= 0; o[k][j]=  0; o[k][k]=1;

  loop (a,0,3) loop (b,0,3)
    rot2[a][b]=o[a][0]*drot[0][b]+o[a][1]*drot[1][b]+o[a][2]*drot[2][b];

  loop (a,0,3) loop (b,0,3) {
    drot[a][b]=rot2[a][b];
    if (bodyframe) rot[b][a]=drot[a][b];
    else rot[a][b]=rot2[a][b]; }
}

int OPTION_C=0;

#include "match.c"

int main(int narg,char **arg)
{
  FILE *plbin,*plbout;
  fvector hdr;
  int varL,ns,i;
  fvector *r;

  if (narg<3) {
    fprintf(stderr,"\
Match (rotate and move) configurations according to selected sites. Call by:\n\
  plbmatch INPUT OUTPUT [MATCH[:FRAME[:FIXED]]]\n\
INPUT  plb-file\n\
OUTPUT plb-file; - = stdout (the same format is INPUT)\n\
MATCH  plb-file to be used as the reference, the same number of sites required\n\
       missing MATCH = INPUT (frame 1, all sites)\n\
FRAME  the frame to be used as the reference [default=1]\n\
FIXED  file of sites (one per line) to match, [default=all]\n\
Example:\n\
  plbmatch cluster-gauss.plb cluster-norot.plb cluster-gauss.plb:1:cluster.fix\n\
See also:\n\
  show (option -m and hot keys m M)\n");
    exit(0); }

  initscroll(0);
  plbin=fopen(arg[1],"rb");
  if (!plbin) Error(arg[1]);
  fread(hdr,4,2,plbin);
  varL=hdr[1]<0;
  ns=hdr[0];
  
  plbout=fopen(arg[2],"wb");
  if (plbout==NULL) Error("cannot open OUTPUT");
  fwrite(hdr,4,2,plbout);

  if (narg>3) matchfile=strdup(arg[3]);
  else matchfile=strdup(arg[1]);

  allocarrayzero(r,ns);
  
  for (;;) {
    if (varL) fread(hdr,4,3,plbin);
    if (ns!=fread(r,sizeof(fvector),ns,plbin)) break;
    match(r,ns,1e-6);
    if (varL) fwrite(hdr,4,3,plbout);
    loop (i,0,ns) {
      fvector cr,R;
      
      /* see showsphere() */
      cr[0]=r[i][0]-center[0];
      cr[1]=r[i][1]-center[1];
      cr[2]=r[i][2]-center[2];
      R[0]=rot[0][0]*cr[0]+rot[0][1]*cr[1]+rot[0][2]*cr[2];
      R[1]=rot[1][0]*cr[0]+rot[1][1]*cr[1]+rot[1][2]*cr[2];
      R[2]=rot[2][0]*cr[0]+rot[2][1]*cr[1]+rot[2][2]*cr[2];
            
      fwrite(R,sizeof(fvector),1,plbout); } }

  fclose(plbout);        
  fclose(plbin);        
  
  return 0;
}
