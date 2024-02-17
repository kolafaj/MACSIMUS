/* cc -O2 -o histogr histogr.c -lm
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
  double from,to,grid,x,x0=9e99,x9=-9e99;
  int *hist;
  int nhist,i,sumhist=0,before=0,after=0;
  char line[1024];

  if (narg<4) {
    fprintf(stderr,"\
Calculate histogram of 1 column data. Call by:\n\
  histogr FROM TO { DX | 1/GRID } < INPUT > OUTPUT\n\
where\n\
  DX = histogram bin width\n\
  GRID = 1/DX = # of bins/unit (not necessarily integer)\n\
  FROM = first bin starts at FROM\n\
  TO = range end rounded up to FROM+INTEGER*DX\n\
  INPUT = ASCII stream of real numbers, comments or invalid lines ignored\n\
  OUTPUT = x,y-file to plot\n\
The output histogram (incl. data outside) is normalized to integral=1\n\
Example:\n\
  tab 0 99999 | tabproc '\\rnd(0)' | histogr 0 1 1/20 | hist2hist 0.05 | plot -\n\
See also:\n\
  histogram # dynamic allocation - slower, does not need FROM, TO\n\
  hist2hist # transform data to rectangle-like plot\n");
   exit(0); }

  from=atof(arg[1]);
  to=atof(arg[2]);
  if (arg[3][0]=='1' && arg[3][1]=='/') grid=atof(arg[3]+2);
  else grid=1/atof(arg[3]);
  nhist=(to-from)*grid+(1-1e-14);
  allocarrayzero(hist,nhist);
 
  while (gets(line))
    if (sscanf(line,"%lf",&x)==1) {
      Min(x0,x) Max(x9,x)
      sumhist++;
      if (x>=from) {
        int i=(x-from)*grid;
 
        if (i<0) before++;
        else if (i>=nhist) after++;
        else hist[i]++; }
      else before++; }
 
  sprintf(line,"# %s %s %s nhist=%d sumhist=%d before=%g after=%g\n",
 	 arg[1],arg[2],arg[3],
 	 nhist,sumhist,
 	 before/(double)sumhist,after/(double)sumhist);
  fputs(line,stderr);

  fputs(line,stdout);
  loop (i,0,nhist) 
    printf("%15.7g %15.7g %d\n",
 	  from+(i+0.5)/grid,
 	  hist[i]/(double)sumhist*grid,
 	  hist[i]);
 	  
  { 
    double dx=(x9-x0)/200,d=1;
    
    while (d<dx) d*=10;
    while (d/10>dx) d/=10;

    if (d*5>dx) d/=5;
    if (d*2>dx) d/=2;
    
//    fprintf(stderr,"# d=%g dx=%g %g bins\n",d,dx,(x9-x0)/d);
    
    x0=floor(x0/d)*d;
    x9=ceil(x9/d)*d;
        
    fprintf(stderr,"histogr %.10g %.10g %.2g\n",x0,x9,d);
//    fprintf(stderr,"# %g bins\n",(x9-x0)/d);
  }
 
  return 0;
}
