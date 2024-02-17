/* cc -O2 -o concprof concprof.c -lm 
 */
#include "../gen/include.h"
typedef double real;
#include "../gen/vector3d.h"

typedef float fvector[3];

int grid=100;
double *aux;
void dosmooth(double *c) /***************************************** dosmooth */
{
  int i;

  loop (i,0,grid) aux[(i+1)%grid]=(c[i]+c[(i+1)%grid]+c[(i+2)%grid])/3;
  loop (i,0,grid) c[i]=aux[i];
}

int main(int narg,char **arg) /**************************************** main */
{
  FILE *plb,*out=NULL;
  char name[256];
  int iarg,i,ns,ns1=0,smooth=2,frame,ii,from,to,verbose=0,layer=7;
  float hdr[2];
  fvector L,r;
  unsigned *hist1,*hist2;
  double *c1,*c2,c2max,sum;
  double drop=0.5;

  if (narg<3) {
    fprintf(stderr,"\
Calculates the concentration of (atoms of) the solute in the region of\n\
a z-slab where bulk solution is detected.\n\
The solution region is the region on the smoothed z-density profile of the\n\
solvent where its concentration/max.concentration > DROP.\n\
A layer LAYER thick at each side is not included in the bulk solution region.\n\
\n\
  |          *                                     ***\n\
  | **   **** **                                  *   ***\n\
  |   ***       *     oo         oo  oooo <--MAX *\n\
  |              **  o  o o  oooo  o     o   ****     * = solute\n\
  |                *o                     o**         o = solvent\n\
  |   DROP*MAX --> o*            *        *o\n\
  |               o| **** * ***** ** ***** |o\n\
  |              o |   | * *        *  |   | o\n\
0 +-ooooooooooooo--+---+---------------+---+--ooooooooo\n\
                   LAYER bulk solution LAYER\n\
\n\
Call by:\n\
  concprof [OPTIONS] NAME [OPTIONS]\n\
FILES:\n\
  NAME.plb input plb-file (new format with box size required)\n\
OPTIONS:\n\
  -gGRID   z-grid, in histogram bins per z-box (Lz) [default=%d]\n\
  -nNM     number of atoms of the solute (e.g., both ions)\n\
  -lLAYER  width of the surface layer (from slab end), in histogram bins [%d]\n\
  -sSMOOTH filter the z-profiles SMOOTH times through the 3-point formula [%d]\n\
  -tDROP   rel. concentration of atoms of the solvent (e.g., water) with\n\
           respect to the maximum (of the smoothed profile) which defines the\n\
           region of the solution [%g]\n\
  -v       verbose (print all profiles)\n\
Algorithm:\n\
  * calculates the z-profiles in [0,Lz] divided into GRID bins, periodic b.c.\n\
  * smoothes the profiles [if verbose, prints them all]\n\
  * selects the region where the solvent concentration > DROP*MAX\n\
  * removes the surface layers (L histogram bins thick) from this region\n\
  * calculates averaged concentrations of the solute in this region\n\
Example:\n\
  concprof cbo4-3-5 -g100 -n250 -l8 -t0.5 -s5\n\
See also:\n\
  densprof plb2dens\n",grid,layer,smooth,drop);
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'l': layer=atof(arg[iarg]+2); break;
        case 'g': grid=atoi(arg[iarg]+2); break;
        case 'n': ns1=atoi(arg[iarg]+2); break;
        case 't': drop=atof(arg[iarg]+2); break;
        case 's': smooth=atoi(arg[iarg]+2); break;
        case 'v': verbose++; break;
        default: Error("unknown option"); }
    else 
      strcpy(name,arg[iarg]);

  if (!(plb=fopen(string("%s.plb",name),"rb"))) Error(name);
 
  allocarray(hist1,grid);
  allocarray(hist2,grid);
  allocarray(c1,grid);
  allocarray(c2,grid);
  allocarray(aux,grid);
 
  fread(hdr,4,2,plb);
  ns=hdr[0];
  if (hdr[1]!=-3) Error("only new (L3) format of plb-files is accepted");

  printf("#frame conc[1/AA^3] from to\n");

  for (frame=1;;frame++) {
    if (3!=fread(L,4,3,plb)) break;

    loop (i,0,grid) hist1[i]=hist2[i]=0;
    
    loop (i,0,ns) {
      if (3!=fread(r,4,3,plb)) Error("plb truncated");
      r[2]+=L[2];
      ii=(int)(r[2]/L[2]*grid)%grid;
      if (i<ns1) hist1[ii]++;
      else hist2[ii]++; }
      
    loop (i,0,grid) {
      c1[i]=hist1[i];
      c2[i]=hist2[i]; }

    loop (i,0,smooth) { dosmooth(c1); dosmooth(c2); }
   
    c2max=0;
    loop (i,0,grid) Max(c2max,c2[i])
    
    if (verbose) {
      out=fopen(string("%s.f%05d.z",name,frame),"wt");
      loop (i,0,grid) fprintf(out,"%g %g %g\n",(double)i/grid,c1[i],c2[i]); }
    
    from=-1,to=-1;
    c2max*=drop;
    loop (i,0,grid) {
      if (c2[i]<c2max && c2[(i+1)%grid]>=c2max) {
        if (from>=0) Error("not smooth enough: decrease -g or increase -s"); 
        from=i; }
      if (c2[i]>c2max && c2[(i+1)%grid]<=c2max) {
        if (to>=0) Error("not smooth enough: decrease -g or increase -s");
        to=(i+1)%grid; } }
    from+=layer;
    to-=layer;
    if (from>=to) to+=grid;
    if (from<0 || to<0) Error("cannot determine range");
    
    sum=0;
    if (out) fprintf(out,"\n");
    loop (i,from,to) {
      if (out) fprintf(out,"%g 0\n",(double)i/grid);
      sum+=c1[i%grid]; }
    if (out) fclose(out);
    printf("%d %g %d %d\n",frame,sum/(L[0]*L[1]*L[2]/grid*(to-from)),from,to); }

  fclose(plb);

  return 0;
}
