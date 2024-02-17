/* cc -o plbstack -O2 plbstack.c
*/

#include "../gen/include.h"

typedef float vector[3];

#define NMOL 20 /* max no. of species */

struct file_s {
  char *fn;
  FILE *f;
  int ns;
  vector L;
  vector *r; /* [ns] */
  int chunk[NMOL];
  int key;
  float z; } *plb;

vector L, LL;

void static swapxyz(float *x,int swapload,int inv) /**************** swapxyz */
{
  float sw;

  switch (swapload%10) {
    case 0: break;
    case 1: sw=x[0]; x[0]=x[1]; x[1]=sw; break; /* x <-> y */
    case 2: sw=x[2]; x[2]=x[1]; x[1]=sw; break; /* y <-> z */
    case 3: sw=x[2]; x[2]=x[0]; x[0]=sw; break; /* x <-> z */
    case 4: sw=x[2]; x[2]=x[1]; x[1]=x[0]; x[0]=sw; break; /* z -> y -> x -> z */
    case 5: sw=x[2]; x[2]=x[0]; x[0]=x[1]; x[1]=sw; break; /* z -> x -> y -> z */
    default: Error("bad key"); }

  //  fprintf(stderr,"%g %g %g\n",L[0],L[1],L[2]); sleep(1);

  if (inv) {
    if ((swapload/10)&1) x[0]=L[0]-x[0];
    if ((swapload/10)&2) x[1]=L[1]-x[1];
    if ((swapload/10)&4) x[2]=L[2]-x[2]; }
}

int main(int narg,char **arg)
{
  int iplb,nplb,ns=0,i,varL,frame;
  int nmol,imol,maxnmol=0,nssum;
  float hdr[2],dz;
  char *sep;
  
  if (narg<2) {
    fprintf(stderr,"\
Stack simulation boxes (binary playback files) in z-direction.  Call by:\n\
  plbstack KEY PLB-FILE[:FRAME] KEY PLB-FILE[:FRAME] ... > MERGED-PLB-FILE\n\
KEY = DZ[:TR][,NS0,NS1,...]]\n\
  DZ = additional shift in z-direction (separation of boxes, to avoid overlaps)\n\
       (note that the first DZ is the last in the periodic b.c.)\n\
  TR = transformation of coordinates (on read, if possible sum of):\n\
    0 = none\n\
    1 = x <-> y             10 = x <-> Lx-x\n\
    2 = y <-> z             20 = y <-> Ly-y\n\
    3 = x <-> z             40 = z <-> Lz-z\n\
    4 = z -> y -> x -> z\n\
    5 = z -> x -> y -> z\n\
  NSi = number of sites in all molecules of species i (product N[i]*ns)\n\
        (default = total number of sites, or the remaining sites)\n\
PLB-FILE = playback file\n\
FRAME = frame number (1=first etc., default=0=last)\n\
Output MERGED-PLB-FILE is in variable-L format, number of sites is the sum,\n\
  it is ordered by species as required by cook (if NSi given),\n\
  output box = max in x,y-directions, sum of all Lz and DZ in the z-direction\n\
Example (stack ice and water by 1A, rotate):\n\
  plbstack 1:4 Ih533.plb 1:4 water540.plb > Ih533-z.plb\n\
Example (stack nacl.plb with 108*nacl and sol.plb with 27*nacl + 206*spce):\n\
  plbstack 1,108,108 nacl.plb 1,27,27 sol.plb > naclsol.plb\n\
  molcfg -135 na cl -206 spce naclsol\n\
  show naclsol\n\
Example (select last frame of sim.plb)\n\
  plbstack 0 sim.plb > lastframe.plb\n\
See also:\n\
  cook -[\n\
  plbmerge plbinfo plbmsd plbcut plbfilt plbtran plbrot plb2diff\n\
  plb2asc asc2plb cfg2plb plb2cfg\n");
    exit(0); }
  
  nplb=narg/2;
  allocarrayzero(plb,nplb+1);
  if (narg!=nplb*2+1) Error("arguments not in pairs KEY PLB-FILE");
  
  loopto (iplb,1,nplb) {

    plb[iplb].key=atoi(arg[iplb*2-1]);

    dz=atof(arg[iplb*2-1]);
    sep=strchr(arg[iplb*2-1],':');
    plb[iplb].key = sep ? atoi(sep+1) : 0;
    
    plb[iplb].fn=strdup(arg[iplb*2]);
    sep=strchr(plb[iplb].fn,':');
    if (sep) {
      *sep=0;
      frame=atoi(sep+1); }
    else
      frame=0;
    
    fprintf(stderr,"opening %s...\n",plb[iplb].fn);
    if (!(plb[iplb].f=fopen(plb[iplb].fn,"rb"))) Error("cannot open");
    if (2!=fread(hdr,4,2,plb[iplb].f)) Error("no header");
    varL=hdr[1]<0;
    ns += plb[iplb].ns = hdr[0];
    
    fseek(plb[iplb].f,0,SEEK_END);
    if (frame<=0) frame=ftell(plb[iplb].f)/((plb[iplb].ns+varL)*12);
    
    if (fseek(plb[iplb].f,8+(frame-1)*(plb[iplb].ns+varL)*12,SEEK_SET))
      Error("no such frame or file truncated");
    if (varL) {
      if (3!=fread(L,4,3,plb[iplb].f)) Error("too short file");
      swapxyz(L,plb[iplb].key,0); }
    else 
      L[0]=L[1]=L[2]=hdr[1];
    loop (i,0,3) plb[iplb].L[i]=L[i];
    Max(LL[0],L[0])
    Max(LL[1],L[1])
    LL[2]+=dz;
    plb[iplb].z=LL[2];
    LL[2]+=L[2]; 
    fprintf(stderr,"  frame=%d, L=%g %g %g, ns=%d\n",frame,L[0],L[1],L[2],plb[iplb].ns);
    
    sep=strchr(arg[iplb*2-1],',');
    nssum=0;
    nmol=0;
    while (sep) {
      sep++;
      if (nmol>NMOL) Error("too many species, increase NMOL and recompile");
      nssum+=atoi(sep);
      if (nssum>plb[iplb].ns) Error("too many sites specified");
      plb[iplb].chunk[nmol]=nssum;
      sep=strchr(sep,','); 
      nmol++; }
    if (nssum<plb[iplb].ns) {
      if (nmol>NMOL) Error("too many species, increase NMOL and recompile");
      plb[iplb].chunk[nmol]=plb[iplb].ns;
      nmol++; }
            
    allocarray(plb[iplb].r,plb[iplb].ns);
    if (plb[iplb].ns!=fread(plb[iplb].r,sizeof(vector),plb[iplb].ns,plb[iplb].f)) Error("file truncated"); 
    loop (i,0,plb[iplb].ns) swapxyz(plb[iplb].r[i],plb[iplb].key,1);

    fprintf(stderr,"  chunks: 0");
    loop (imol,0,nmol) fprintf(stderr," %d",plb[iplb].chunk[imol]);
    fprintf(stderr," sites\n");
    Max(maxnmol,nmol) }

  fprintf(stderr,"output box=%g %g %g, ns=%d\n",LL[0],LL[1],LL[2],ns);

  hdr[0]=ns; hdr[1]=-3;
  fwrite(hdr,4,2,stdout);
  fwrite(LL,4,3,stdout);

  loopto (imol,0,maxnmol) {
    loopto (iplb,1,nplb) {
      int from=0,to=plb[iplb].chunk[imol];
      if (imol) from=plb[iplb].chunk[imol-1];
      
      loop (i,from,to) {
        plb[iplb].r[i][2]+=plb[iplb].z;
        fwrite(plb[iplb].r[i],4,3,stdout); } } }
      
  fclose(stdout);
            
  return 0;
}            
