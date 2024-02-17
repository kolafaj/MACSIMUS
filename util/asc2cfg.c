/* \make asc2cfg

*/

#include "ground.h"
#include "sds.h"
#include "varfile.h"

typedef double vector[3];

int main(int narg,char **arg)
{
  FILE *asc;
  char *ascfn,*cfgfn;
  struct a_s {
    int  size;    /* whole struct in bytes */
    int  padd;    /* padded to 8 bytes */
    double logs;  /* log of the Nose variable s */
    vector lambda;/* log(box) */
    vector shape; /* (reserved for box shape) */
    vector rp[1]; /* contiguous array of all r[ns] (and p[ns]) */
  } *a[3];
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
  int i,cfgkey,ns=0,size;
  int np; /* 1:nonpolar, 2:polar */
  int optionlist[32];
  char line[1024];

  if (narg<2) {
    fprintf(stderr,"Convert ascii-image back into cfg-file. Call by:\n\
  cfg2asc ASCNAME CFGNAME\n\
  cfg2asc NAME\n\
where\n\
  CFGNAME is MACSIMUS cfg-file and ASCNAME is its ascii-image\n\
  if only one parameter is given, ASCNAME=NAME.cfa and CFGNAME=NAME.cfg\n\
  Note that missing coordinates in ASCNAME are replaced by zero\n\
See also:\n\
  cfg2asc cfg2atm plb2cfg cfgconv\n");
    exit(0); }

  initscroll(0);

  if (narg==3) {
    cfgfn=arg[2];
    ascfn=arg[1]; }
  else if (narg==2) {
    cfgfn=string("%s.cfg",arg[1]);
    ascfn=string("%s.cfa",arg[1]); }
  else
    ERROR(("too many arguments"))

  asc=fopen(ascfn,"rt");
  if (!asc) ERROR(("cannot open \"%s\"",ascfn))
  VarOpen(cfgfn,"w");

  prt("reading %s",ascfn);
  if (!fgets(line,1024,asc)) ERROR(("\"%s\" truncated (no header)",ascfn))
  cfgkey=atoi(line);
  VarPut(&cfgkey,sizeof(cfgkey));
  np = cfgkey&2 ? 2 : 1;

  loop (i,0,32) {
    if (!fgets(line,1024,asc)) ERROR(("\"%s\" truncated (optionlist)",ascfn))
    optionlist[i]=atoi(line+2);
    if (i!=(line[1]&31)) ERROR(("in file \"%s\", line:\n%s\n*** option -%c expected",ascfn,line,tolower(i+'@'))) }
  VarPut(optionlist,sizeof(optionlist));

  if (!fgets(line,1024,asc)) ERROR(("\"%s\" truncated (nspec)",ascfn))
  nspec=atoi(line);
  VarPut(&nspec,sizeof(nspec));

  allocarray(spec,nspec);
  loop (i,0,nspec) {
    if (!fgets(line,1024,asc)) ERROR(("\"%s\" truncated (species %i)",ascfn,i))
    sscanf(line,"%d%d",&spec[i].N,&spec[i].ns); }
  VarPut(spec,nspec*sizeof(spec[0]));

  for (;;) {
    if (!fgets(line,1024,asc)) ERROR(("\"%s\" truncated (keys)",ascfn))
    if (line[0]=='.') break;

    sscanf(line,"%d%d%lf%lf%lf",&rec.key,&rec.intval,&rec.vecval[0],&rec.vecval[1],&rec.vecval[2]);
    if (rec.key==2) { ns=rec.intval; put(ns) }
    VarPut(&rec,sizeof(rec)); }

  if (!ns) ERROR(("record key=2 missing or ns=0"))

  size=sizeof(*a[0])+sizeof(vector)*(ns*np-1);
  put(size)
  loop (i,0,3) sdsalloczero(a[i],size);
  if (!fgets(line,1024,asc)) ERROR(("unexpected EOF while reading logs"))
  if (line[0]=='!') fgets(line,1024,asc);
  if (!strstr(line,"logs")) ERROR(("%s: logs expected",line))
  sscanf(line,"%lf%lf%lf",&a[0]->logs, &a[1]->logs, &a[2]->logs);

  if (!fgets(line,1024,asc)) ERROR(("unexpected EOF while reading ln(L)=lambda"))
  sscanf(line,"%lf%lf%lf%lf%lf%lf%lf%lf%lf",
          &a[0]->lambda[0], &a[0]->lambda[1], &a[0]->lambda[2],
          &a[1]->lambda[0], &a[1]->lambda[1], &a[1]->lambda[2],
          &a[2]->lambda[0], &a[2]->lambda[1], &a[2]->lambda[2]);

  loop (i,0,ns*np) {
    if (!fgets(line,1024,asc)) {
      WARNING(("\"%s\" truncated (reading record %d)",ascfn,i))
      break; }
    while (line[0]=='!') fgets(line,1024,asc);
    sscanf(line,"%lf%lf%lf%lf%lf%lf%lf%lf%lf",
           &a[0]->rp[i][0], &a[0]->rp[i][1], &a[0]->rp[i][2],
           &a[1]->rp[i][0], &a[1]->rp[i][1], &a[1]->rp[i][2],
           &a[2]->rp[i][0], &a[2]->rp[i][1], &a[2]->rp[i][2]);
  }
  VarPut(a[0],size);
  VarPut(a[1],size);
  VarPut(a[2],size);
  VarClose();

  fclose(asc);

  prt("%s written",cfgfn);

  return 0;
}
