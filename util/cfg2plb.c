/* \make cfg2plb
   convert MACSIMUS cfg-file to plb-file
*/

#include "ground.h"
#include "sds.h"
#include "varfile.h"

typedef double vector[3];

int main(int narg,char **arg) /**************************************** main */
{
  FILE *f;
  float header[2];
  vector L={0,0,0};
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
  } *a;  
  int nspec;
  struct spec_s {
    int N;           /* number of molecules */
    int ns;          /* number of sites */
  } *spec; /* [nspec] */
  struct rec_s {
    int key;
    int intval;
    vector vecval;
  } rec;
  int i,k,cfgkey,N=0,ns=0,nsx,NS,oldsize;
  int optionlist[32];
  double t,h,Entot;

  if (narg<2) {
    fprintf(stderr,"\
Convert cfg-file into plb-file. Call by:\n\
  cfg2plb NAME [NS]\n\
where\n\
  NAME.cfg will be converted into NAME.plb (one frame)\n\
  (asks for confirmation if NAME.plb would be overwritten)\n\
  NS is optional number of sites (to truncate the configuration)\n\
Limitations:\n\
  some information is lost on conversion\n\
  cannot handle cfg-files written in float or with reversed endianness\n\
See also:\n\
  cfg2drude cfg2flt6 cfg2atm plbinfo plb2asc asc2plb plbconv plb2plb plb2cfg\n");
    exit(0); }

  initscroll(0);

  VarOpen(fn=string("%s.cfg",arg[1]),"r");
  prt("reading %s",fn);
  VarRead(&cfgkey,sizeof(cfgkey));
  VarRead(optionlist,sizeof(optionlist));
  VarRead(&nspec,sizeof(nspec));

  if (nspec) {
    allocarray(spec,nspec);
    VarRead(spec,nspec*sizeof(spec[0])); }

  while (VarFile.size==sizeof(struct rec_s)) {
    VarRead(&rec,sizeof(rec));
    switch (rec.key) {
      case 1: 
        N=rec.intval;
        put2(cfgkey,N)
        t=rec.vecval[0];
        h=rec.vecval[1];
        Entot=rec.vecval[2];
        put3(t,h,Entot)
        break;

      case 2:
        ns=rec.intval;
        loop (k,0,3) L[k]=rec.vecval[k];
        prt("Ns=%d L=%g %g %g loaded",ns,L[0],L[1],L[2]);
        break;
    } }


  nsx=(VarFile.size-sizeof(*a))/sizeof(vector)+1;
  if (cfgkey&2) {
    prt("POLAR configuration detected");
    nsx/=2; }
  if (ns && nsx!=ns) prt("WARNING: ns=%d from record length, %d from key=2: %d used",nsx,ns,nsx);
  ns=nsx;
  put(ns)
  if (ns<N) ERROR(("ns<N read from file"))
  if (narg>2) NS=atoi(arg[2]);
  else NS=ns;
  if (NS>ns) prt("WARNING: NS>ns specified, NS=ns=%d assumed",NS=ns);
  if (NS<ns) prt("WARNING: smaller NS=%d specified, configuration will be truncated",NS);

  sdsalloc(a,VarFile.size);
  oldsize=VarFile.size;
  VarRead(a,VarFile.size);
  VarClose();
  if (a->size!=oldsize)
    ERROR(("size of SDS in file (%d) does not match VarFile.size (%d)",a->size,oldsize))

  fn=string("%s.plb",arg[1]);
  f=fopen(fn,"r");
  if (f) {
    WARNING(("%s exists and is about to be overwritten",fn))
    fclose(f); }

  f=fopen(fn,"wb");
  if (!f) ERROR(("cannot open %s for writing",fn))

  header[0]=NS;
  header[1]=-3;

  if (fwrite(header,4,2,f)!=2) ERROR(("cannot write to %s",fn))
  loop (k,0,3) Lf[k]=L[k];
  if (fwrite(Lf,4,3,f)!=3) ERROR(("cannot write to %s",fn))
  loop (i,0,NS) {
    loop (k,0,3) Lf[k]=a->rp[i][k];
    if (fwrite(Lf,4,3,f)!=3) ERROR(("cannot write to %s",fn)) }
  fclose(f);

  prt("%s written",fn);

  return 0;
}
