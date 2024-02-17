/* \make cfg2atm
*/

#include "ground.h"
#include "sds.h"
#include "varfile.h"

typedef double vector[3];

int main(int narg,char **arg)
{
  FILE *f;
  vector L={0,0,0};
  char *fn;
  struct a_s {
    int  size;    /* whole struct in bytes */
    int  padd;    /* padded to 8 bytes */
    double logs;  /* log of the Nose variable s */
    vector lambda;/* log(box) */
    vector shape; /* (reserved for box shape) */
    vector rp[1]; /* contiguous array of all r[ns] (and p[ns]) */
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
  double BOX=1;

  if (getenv("BOX")) BOX=atof(getenv("BOX"));

  if (narg<3) {
    fprintf(stderr,"Convert cfg-file into atm-file. Call by:\n\
  cfg2atm NAME ATOM1 [ATOM2 ...]\n\
where\n\
  NAME.cfg will be converted into NAME.atm (one frame)\n\
  (asks for confirmation if NAME.atm would be overwritten)\n\
  If the configuration is longer, ATOMs are repeated\n\
  Empty ATOM or -: atom is omitted\n\
Environment:\n\
  NS=number of sites read from input (to truncate the configuration)\n\
  BOX=box multiplication factor (default=1); 0=omit box line\n\
Limitations:\n\
  some information is lost on conversion\n\
  cannot handle cfg-files written in float or with reversed endianness\n\
Example (tip4p256.cfg contains 256 TIP4P molecules in order HOMH):\n\
  cfg2atm tip4p256 H O - H\n\
See also:\n\
  cfg2asc cfg2plb plb2asc\n");
    exit(0); }

  initscroll(0);
  
  VarOpen(fn=string("%s.cfg",arg[1]),"r");
  prt("reading %s",fn);
  VarRead(&cfgkey,sizeof(cfgkey));
  VarRead(optionlist,sizeof(optionlist));
  VarRead(&nspec,sizeof(nspec));
  allocarray(spec,nspec);
  VarRead(spec,nspec*sizeof(spec[0]));

  while (VarFile.size==sizeof(struct rec_s)) {
    VarRead(&rec,sizeof(rec));
    switch (rec.key) {
      case 1: N=rec.intval;
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
  NS=ns;
  if (getenv("NS")) NS=atoi(getenv("NS"));
  if (NS>ns) prt("WARNING: NS>ns specified, NS=ns=%d assumed",NS=ns);
  if (NS<ns) prt("WARNING: smaller NS=%d specified, configuration will be truncated",NS);

  sdsalloc(a,VarFile.size);
  oldsize=VarFile.size;
  VarRead(a,VarFile.size);
  VarClose();
  if (a->size!=oldsize)
    ERROR(("size of SDS in file (%d) does not match VarFile.size (%d)",a->size,oldsize))


  fn=string("%s.atm",arg[1]);
  f=fopen(fn,"r");
  if (f) {
    WARNING(("%s exists and is about to be overwritten",fn))
    fclose(f); }
            
  f=fopen(fn,"wb");
  if (!f) ERROR(("cannot open %s for writing",fn))

  ns=0;
  loop (i,0,NS) if (arg[2+i%(narg-2)][0]>'-') ns++;

  fprintf(f,"%d\n",ns);
  if (BOX) fprintf(f," box %9.6f %9.6f %9.6f\n",L[0]*BOX,L[1]*BOX,L[2]*BOX);
  else fprintf(f,"\n");
  loop (i,0,NS) if (arg[2+i%(narg-2)][0]>'-')
    fprintf(f,"%2s %9.6f %9.6f %9.6f\n",arg[2+i%(narg-2)],a->rp[i][0],a->rp[i][1],a->rp[i][2]);
  fclose(f);

  prt("%s written",fn);

  return 0;
}
