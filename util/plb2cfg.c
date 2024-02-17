/* \make plb2cfg
*/

#include "ground.h"
#include "sds.h"
#include "varfile.h"
#include "rndgen.h"
#include "options.h"
#include "optlist.c"

typedef double vector[3];

int main(int narg,char **arg)
{
  FILE *f;
  float header[2];
  vector L,rmin;
  vector scale={1,1,1};
  float Lf[3];
  char *fn;
  struct a_s {
    int  size;    /* whole struct in bytes */
    int  dep;     /* set if dependants calculated, unset if cfg changes */
    double logs;  /* log of the Nose variable s */
    vector lambda;/* log(L): NEW since V2.7a 
                     active for derivatives, cfg[0].lambda is always derived from box.L*/
    vector offdiag;/* reserved, NEW since V2.7a */
    vector rp[1]; /* contiguous array of No.s vectors (POLAR: 2*No.s) */
  } *cfg;  
  struct rec_s {
    int key;
    int intval;
    vector vecval;
  } rec;
  int i,k,cfgkey,N=0,ns,NS,varL,frame,norm=0,nspec=0;
  double h=0.001;
  char *c;
  double T=0;

  if (narg<3) {
    fprintf(stderr,"Convert plb-file into cfg-file. Call by:\n\
  plb2cfg PLB-FILE[:FRAME] CFG-FILE [[-]N [T/M [NS [XS YS ZS]]]\n\
where\n\
  NAME.plb will be converted into NAME.cfg\n\
  FRAME = frame (1st=1, default=0=last)\n\
  old NAME.cfg is renamed to NAME.cfg~\n\
  N = number of molecules (not needed for loading to cook, cf. load.N)\n\
 -N = as above and try to shift periodically to avoid unnormalized molecules\n\
  T = temperature in K [default=0]\n\
  M = \"typical\" molar mass of atom in g/mol [default=1]\n\
  NS = number of sites (to truncate the configuration) [default=0=from plb]\n\
  XS,YS,ZS = if box missing in plb: the box (default=1)\n\
             if box present: scale box and coordinates by given factors\n\
Notes:\n\
  nspec=0 used, table of sites and molecules not written (cook can read this)\n\
  The CFG-FILE should be used with init=2\n\
  T and M serve for setting initial velocities. All atoms are assumed to have\n\
    the same mass M which is imprecise for molecules composed of various atoms;\n\
    in this case, specifying  initvel=99999  in cook data gives better results\n\
  If T=0 or not given, the cfg-file contains only positions; cook then assumes\n\
    that velocities are zero (unless initvel=99999 is given)\n\
Example:\n\
  plb2cfg ice.plb ice.cfg 360 273/6\n\
See also:\n\
  plbinfo plb2asc asc2plb plbconv plb2plb cfg2plb cfg2atm\n");
    exit(0); }

  initscroll(0);
  
  fn=strdup(arg[1]);
  if ( (c=strchr(fn,':'))) frame=atoi(c+1),*c=0;
  else frame=0;

  f=fopen(fn,"rb");
  if (!f) ERROR(("cannot open %s",fn))
  else prt("reading %s",fn);

  if (2!=fread(header,4,2,f)) ERROR(("too short %s",fn))
  ns=header[0];
  varL=(header[1]<0);
  fseek(f,0,SEEK_END);
  if (frame<=0) {
    frame=ftell(f)/((ns+varL)*12);
    prt("finding the last frame"); }
  put2(ns,frame);
  
  if (narg>3) N=atoi(arg[3]);
  else WARNING(("N not specified (zero used, OK for loading to cook)"))
  if (N<0) norm++,N=abs(N);
  
  if (narg>4) {
    char *c=strchr(arg[4],'/');
    double M=1;
    
    T=atof(arg[4]);
    if (c) M=atof(c+1);
    if (M==0) M=1;
    /* v=sqrt(R*T/M); note units. Cook accepts a[1] multiplied by h */
    T=h*sqrt(0.8314472*T/M); }    

  if (narg>5) NS=atoi(arg[5]);
  else NS=ns;
  if (NS==0) NS=ns;
  if (NS>ns) prt("WARNING: NS>ns specified, NS=ns=%d assumed",NS=ns);
  if (NS<ns) prt("WARNING: smaller NS=%d specified, configuration will be truncated",NS);

  if (narg>6) scale[0]=atof(arg[6]);
  if (narg>7) scale[1]=atof(arg[7]);
  if (narg>8) scale[2]=atof(arg[8]);
  
  if (fseek(f,8+(frame-1)*(ns+varL)*12,SEEK_SET))
    ERROR(("no such frame or %s truncated",fn))
  if (varL) {
    if (3!=fread(Lf,4,3,f)) ERROR(("no such frame or %s truncated",fn)) }
  else 
    Lf[0]=Lf[1]=Lf[2]=header[1];

  loop (i,0,3) {
    L[i]=Lf[i]*scale[i];
    if (L[i]<0) {
      L[i]=Lf[i]=scale[i]; 
      prt("WARNING: zero L%c replaced by scale=%g",i+'x',scale[i]); } }

  rename(arg[2],string("%s~",arg[2]));
  VarOpen(arg[2],"w");
  prt("writing %s",arg[2]);
  cfgkey=1; /* version 2.7 */
  VarPut(&cfgkey,sizeof(cfgkey));
  optionlist['m'&31]=1+(T!=0); /* 1=just positions,2=also velocities */
  VarPut(optionlist,sizeof(optionlist));
  put2(cfgkey,N)

  VarPut(&nspec,sizeof(nspec));

  memset(&rec,0,sizeof(rec));

  rec.key=1;
  rec.intval=N;
  rec.vecval[1]=h;
  VarPut(&rec,sizeof(rec));

  rec.key=2;
  rec.intval=ns; 
  loop (k,0,3) rec.vecval[k]=L[k];
  VarPut(&rec,sizeof(rec));
  fprintf(stderr,"%g %g %g\n",rec.vecval[0],rec.vecval[1],rec.vecval[2]);

  sdsalloc(cfg,sizeof(*cfg)+(NS-1)*sizeof(vector));

  loop (k,0,3) rmin[k]=0;
  
  loop (i,0,NS) {
    if (3!=fread(Lf,4,3,f)) ERROR(("no such frame or %s truncated",fn))
    loop (k,0,3) {
      Min(rmin[k],Lf[k])
        cfg->rp[i][k]=(double)Lf[k]*scale[k]; } }
  fclose(f);

  loop (k,0,3) if (rmin[k]<0) 
    prt("WARNING: negative %c-coordinate=%g: %s",
                           k+'x',        rmin[k],
        norm?"cfg will be shifted":"some versions of cook may fail");
    
  if (norm) loop (i,0,NS) loop (k,0,3) cfg->rp[i][k]-=rmin[k];

  VarPut(cfg,cfg->size);
  if (T) {
    rndinit(0,0);
    loop (i,0,NS) loop (k,0,3) cfg->rp[i][k]=rndgauss()*T;
    VarPut(cfg,cfg->size); }
    
  KeyClose(0);

  return 0;
}
