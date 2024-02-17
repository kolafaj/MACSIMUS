/* \make plbmsd
 */

#include "ground.h"
#include "alloc.h"
#include "linregr.h"

/*
  For given series of sites (given on the command line), plbmsd calculates:
  
  (1) square displacement from the 1st frame, as the function of frame
      number t: (r(t)-r(1))^2
  (2) msd over "time" period 1..t: (r(t)-r(1))^2/(t-1)
  (3) linear regression of (1) in interval FROM..maxframe.
      This msd should be more reliable in Brownian region: set FROM
      to the beginning of the linear region of (1)
  (4) individual one-frame square displacements, (r(t)-r(t-1))^2
  
  Note that in (1)--(3), the base configuration is always the 1st frame
  -- do not confuse with FROM which is barely the first datum for linear
  regression!  If reliable autodiffusion coefficients are to be
  obtained, unequilibrated configurations have to be cut out first.
  
  Motion of waters around protein is probably far from Brownian, then,
  (1) and (4) are of interest.
  
  The autodiffusion coefficient D = msd/(2*DIM)
  
  05/2008:
   - DIM added
   - D printed
   - average of all sites added (if >1 site)
   - weight changed to 1/n
   - output of .msd1 (consecutive frames) removed
  06/2005:
   - updated for new plb format
   - reverse endian removed
   - variable L allowed
   - because of variable L, uses reduced r/L variables internally
   - hence, the last L is used to convert back to real units
*/

typedef double vector[3];
typedef float fvector[3];

FILE *plb;
int ns,reverse=0,varL;
double argL;
float hdr[2];

int readframe(vector L,vector *r) /******************************* readframe */
{
  int i,k;
  fvector fL,*fr;

  allocarray(fr,ns);

  if (varL) {
    if (fread(fL,sizeof(fL),1,plb)!=1) return 0; }
  else 
    fL[0]=fL[1]=fL[2]=hdr[1];

  if (argL>0) L[0]=L[1]=L[2]=argL;
  else L[0]=fL[0],L[1]=fL[1],L[2]=fL[2];

  if (ns!=fread(fr,sizeof(fvector),ns,plb)) return 0;

  loop (k,0,3) loop (i,0,ns) r[i][k]=fr[i][k];

  if (argL<0)
    L[0]=L[1]=L[2]=0; /* marks free b.c. */
  else
    loop (k,0,3) if (L[k])
      loop (i,0,ns) r[i][k]=fr[i][k]/L[k];

  return 1;
}

int main(int narg,char **arg) /*************************************** main */
{
  vector *r0,*r1,*r,L,L0,L1;
  int i,n,k,is,from=2,DIM=3;
  char *fn,*end,*c;
  double *msd;
  int *site;
  double sum,x;

  initscroll(0);

  if (narg<3) {
  fprintf(stderr,"\
(Mean) square displacement of selected sites.  Call by:\n\
  plbmsd NAME[:FROM] L[:DIM] SITE [SITE ...]\n\
where\n\
  input playback file = NAME.p00 or NAME.plb (check in this order)\n\
  output file=NAME.msd\n\
  FROM: start of linear regression, default=2=first datum\n\
  L=0: read L from plb-file\n\
  L>0: use cubic box L^3 (ignore L from file)\n\
  L<0: force FREE boundary conditions\n\
  DIM: 1=1D diffusion (x only), 2=2D (x+y), 3=3D (x+y+z, default)\n\
  SITE=site #\n\
See also:\n\
  plbmsd.sh - extension to make averages over sites and blocks\n\
  plb2diff  - automatic calculation of diffusivity and more via cook\n\
  plb2plb plbsites plb2hist testplbmsd\n\
");
  exit(0); }

  alloc(fn,strlen(arg[1])+5);
  strcpy(fn,arg[1]);
  if ((end=strchr(fn,':'))) from=atoi(end+1);
  else end=strend(fn);
  if (from<2) ERROR(("wrong FROM=%d: must be at least 2\n\
first square displacement is available between t=1 and 2",from))
  argL=atof(arg[2]);
  if ((c=strchr(arg[2],':'))) DIM=atoi(c+1);  
  if (DIM<1 || DIM>3) Error("wrong DIM");

  strcpy(end,".p00");
  plb=fopen(fn,"rb");
  fprintf(stderr,"opening %s\n",fn);
  if (!plb) {
    strcpy(end,".plb");
    fprintf(stderr,"no such file, trying %s\n",fn);
    plb=fopen(fn,"rb");
    if (!plb) ERROR(("no such file")) }

  if (fread(hdr,sizeof(hdr),1,plb)!=1) ERROR(("file too short"))

  ns=hdr[0];
  varL=hdr[1]<0;
  if (varL) fprintf(stderr,"%d sites, variable L (new format)\n",ns);
  else fprintf(stderr,"%d sites, L=%f\n",ns,hdr[1]);
  if (ns<=0 || ns>16777216) ERROR(("wrong number of sites %d: endian?",ns))

  strcpy(end,".msd");
  out=fopen(fn,"wt");

  alloc(r0,ns*sizeof(vector));
  alloc(r1,ns*sizeof(vector));
  alloc(r,ns*sizeof(vector));
  alloc(site,narg*sizeof(site[0]));
  alloc(msd,narg*sizeof(msd[0]));

  prt("# square displacements from configuration 1: (r(t)-r(1))^2, DIM=%d",DIM);

  prt_("# t");
  loop (i,3,narg) {
    site[i]=atoi(arg[i]);
    prt_(" %7d  ",site[i]);
    if (site[i]<0 || site[i]>=ns) ERROR(("site %d out of range",site[i])) }
  if (narg>3) prt_("    aver");
  _n

  readframe(L0,r0); n=1;
  memcpy(L1,L0,sizeof(vector));
  memcpy(r1,r0,ns*sizeof(vector));

  while (readframe(L,r)) {
    n++;  

    loop (k,0,DIM) {
      if (L[k]) loop (i,0,ns) {
        double D=r[i][k]-r1[i][k];
  
        while (D>0.5) { r[i][k]-=1; D-=1; }
        while (D<-0.5) { r[i][k]+=1; D+=1; } } }

    prt_("%3d",n);

    sum=0;
    loop (i,3,narg) {
      is=site[i];
      msd[i]=0;
      loop (k,0,DIM) 
        if (L[k]) msd[i]+=Sqr(L[k]*(r[is][k]-r0[is][k]));
        else msd[i]+=Sqr(r[is][k]-r0[is][k]);
      if (n>=from) LRAdd(arg[i],1./n,n,msd[i]);
      prt_(" %9.5f",msd[i]); 
      sum+=msd[i]; }

    if (narg>3) {
      prt_(" %9.5f",sum/(narg-3));
      if (n>=from) LRAdd("aver",1./n,n,sum/(narg-3)); }

    _n
    memcpy(r1,r,ns*sizeof(vector)); }

  prt_("\n# D [AA^2/dt.plb] between t=1 and t=%d:\n   ",n);
  loop (i,3,narg) {
    x=msd[i]/((n-1)*2*DIM);
    sum+=x;
    prt_(" %9.5f",x); }
  if (narg>3) prt_("  %9.5f",sum/(narg-3));
  
  prt_("\n\n# D [AA^2/dt.plb] from linear regression between %d and %d:\n   ",from,n);
  loop (i,3,narg) prt_(" %9.5f",LRRes(arg[i],"B")/2/DIM);
  if (narg>3) prt_("  %9.5f",LRRes("aver","B")/2/DIM);
  
  prt("  D FINAL");
  prt("\n# linearly regressed data:");
  for (is=from; is<n; is+=(n-from)-1) {
    prt_("%3d",is);
    loop (i,3,narg) prt_(" %9.5f",LRRes(arg[i],"A")+is*LRRes(arg[i],"B"));
    if (narg>3) prt(" %9.5f",LRRes("aver","A")+is*LRRes("aver","B")); }
    prt("\n\n# linear regression protocol:"); 
    LRPrint('l');
    fclose(out);

  fclose(plb);

  return 0;
}
