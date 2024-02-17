/* make sfourier
   2016 - update (fft.c+fft.h instead of #include fourier.c)
   2009 - update (L added)
   1999 - first version
 */

#include "ground.h"
#include "fft.h"
/* float if -DFLOAT, otherwise double
   CAVEAT: not conforming to -DPRECISION */
typedef fftcomplex complex;
typedef fftreal real;

struct {
  char name[12];
  char *check;
  double b; 
} *site;

int ns;

int findsite(char *s)
{
  int i;
 
  if (!ns) return 0;
  if (!s || !s[0]) ERROR(("# site1:site2 expected"))
  loop (i,0,ns) if (!strcmp(site[i].name,s)) return i;
  ERROR(("site %s found in g-file but not specified",s));

  return 0;
}

int main(int narg,char **arg)
{
  double GRID,hlast,rlast,CUTOFF;
  int NG,N,i,nzero,iarg,site1,site2,npairs,npairstot=0;
  char line[256];
  double *g,*sf,w,wtot=0,V,rho,x;
  complex *s;
  
  if (narg<4) {
    fprintf(stderr,"Calculate structure factor from radial distribution functions.\n\
\n\
Partial SF from site-site RDF:\n\
  sfourier GRID CUTOFF NS < NAME.X.X.g > NAME.sfg\n\
\n\
Total SF for a mixture from a set of RDFs:\n\
  sfourier GRID CUTOFF NS SITE1:b1 SITE2:b2 ... < NAME.g > NAME.sfg\n\
\n\
ARGUMENTS:\n\
  NAME.X.X.g = site-site RDF obtained by `rdfg -g NAME' (see rdfg.c for format)\n\
  NAME.g     = concatenated RDFs (use `cat NAME.*.*.g')\n\
  GRID       = grid points/A, must match the grid of the RDF files\n\
  CUTOFF     = rdf.cutoff (see cook input data)\n\
               if CUTOFF<rdf.cutoff then the data are truncated\n\
               if CUTOFF>rdf.cutoff then the data are padded by g(r)=1\n\
               NOTE: a check is made for incomplete last histogram bin\n\
  NS         = total # of sites (atoms; not types of sites)\n\
  SITE#      = name of site #. The sites may be listed in any order\n\
  b#         = scattering length of site #\n\
FFT is used: faster if trunc(GRID*CUTOFF) does not contain large primes\n\
NOTE: the structure factor calculated directly by cook is more accurate\n\
See also:\n\
  rdfg smoothg spectrum\n\
");
   exit(0); }

  initscroll(0);

  GRID=atof(arg[1]);
  CUTOFF=atof(arg[2]);
  NG=(int)(GRID*CUTOFF+0.5);
  prt("# %d grid points",NG);
  if (fabs(NG/(double)GRID-CUTOFF)>1e-5) 
    WARNING(("CUTOFF is not an integer multiple of 1/GRID: rounded"))
  N=atoi(arg[3]);
  ns=narg-4;

  if (ns) {
    allocarrayzero(site,ns);
    loop (iarg,4,narg) {
      char *c=strchr(arg[iarg],':');

      if (!c) ERROR(("`:' expected in argument SITE:b"))
      site[iarg-4].b=atof(c+1);
      copy(site[iarg-4].name,arg[iarg],c-arg[iarg]);
      alloczero(site[iarg-4].check,ns); }
    alloczero(sf,sizeof(sf[0])*NG); }
  allocarrayzero(g,NG);
  allocarrayzero(s,4*NG);
    
  for (;;) {
    char *c;
    double r,gg;
 
    for (;;) {
      if (!gets(line)) goto end;
      if (line[0]=='#' && line[1]!='#') break; }
    memset(g,0,sizeof(g[0])*NG);
 
    printf("%-8s ",line);
    site1=findsite(strtok(line+1,"\n :"));
    site2=findsite(strtok(NULL,"\n :"));
    
    if (!gets(line)) ERROR(("unexpected EOF"))
    if ( !(c=strchr(line,']')) ) ERROR(("line %s: two `]' expected",line))
    if ( !(c=strchr(c+1,']')) ) ERROR(("line %s: two `]' expected",line))
    npairs=atoi(strtok(c+1," \n"));
    npairstot+=npairs;
    if (ns) w=site[site1].b*site[site2].b*npairs;
    else w=1;
    wtot+=w;
 
    if (!gets(line)) ERROR(("unexpected EOF"))
    if ( !(c=strchr(line,'=')) ) ERROR(("line %s: two `=' expected",line))
    V=atof(c+1);
    if ( !(c=strchr(c+1,'=')) ) ERROR(("line %s: two `=' expected",line))
    x=atof(strtok(c+1," \n"));
    if (fabs(GRID-x)>1e-6)
      ERROR(("grid in file=%g does not match argument=%g\n\
*** rerun sfourier with the correct grid!",x,GRID))
    rho=N/V;
    prt(" %5d pairs  weight=%-12.6g  number density=%g",npairs,w,rho);
 
    do
      if (!gets(line)) ERROR(("unexpected EOF"))
      while (strlen(line)<2 || line[0]=='#');
 
    for (;;) {
      if (!gets(line) || strlen(line)<2 || line[0]=='#') break;
      sscanf(line,"%lf%lf",&r,&gg);
      i=r*GRID;
      if (i>=0 && i<NG) g[i]=gg; }
 
    nzero=0;
    hlast=rlast=0; 
    loop (i,0,NG) {
      /* to fix some bad ends... */
      if (nzero>10 && i<NG-1 && g[i+1]==0 && g[i]<0.7) g[i]=1;
      if (g[i]) nzero++;
      else if (nzero>10) g[i]=1;
      r=(i+0.5)/GRID;
      s[2*i].im=((g[i]-1)*r+hlast*rlast)/2;
      //printf("%f %f\n",(r+rlast)/2,s[2*i].im/(r+rlast)*2);
      hlast=g[i]-1;
      s[2*i+1].im=hlast*r;
      //printf("%f %f\n",r,s[2*i+1].im/r);
      rlast=r; }
 
    loop (i,1,2*NG) s[4*NG-i].im=-s[i].im;
    
    Fourier(s,4*NG);
 
    prt("#k/AA^-1    S(k)");
 
    loop (i,1,NG) {
      x=-rho/GRID*s[i].re/(GRID*i/NG)+1;
      if (!ns)
        printf("%7.4f %7.4f\n",0.5*GRID*i/NG,x);
      else
        sf[i]+=w*x; } }
 
 end:
  if (ns) {
    if (npairstot!=N*(N-1)/2)
      ERROR(("%d pairs found, expected N*(N-1)/2=%d (wrong N or missing data)",npairstot,N*N))
    loop (i,1,NG) 
      printf("%7.4f %7.4f\n",0.5*GRID*i/NG,sf[i]/wtot); }
 
  return 0;
}
