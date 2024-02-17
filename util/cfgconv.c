/* \make cfgconv
*/

#include "ground.h"
#include "sds.h"
#include "varfile.h"

typedef double vector[3];

int main(int narg,char **arg)
{
  vector L;
  char *fn,*name,*end;
  struct oldcfg_s {
    int  size;    /* whole struct in bytes */
    int  padd;    /* padded to 8 bytes */
    double logs;  /* log of the Nose variable s */
    vector rp[1]; /* contiguous array of all r[ns] (and p[ns]) */
  } *oldcfg[9];  
  struct newcfg_s {
    int  size;    /* whole struct in bytes */
    int  padd;    /* padded to 8 bytes */
    double logs;  /* log of the Nose variable s */
    vector lambda;/* log(box) */
    vector shape; /* (reserved for box shape) */
    vector rp[1]; /* contiguous array of all r[ns] (and p[ns]) */
  } *newcfg;
  /* new key-based record; must be less than 7 doubles */
  struct rec_s {
    int key;
    int intval;
    vector vecval;
  } rec;
  int i,k,cfgkey,N,ns,NS,oldsize,newsize,m=0,nspec=0,polar=0;
  int optionlist[32];
  double t,h,Entot,RvdW,zero,taudip;

  if (narg<2) {
    fprintf(stderr,"Convert old cfg-file (<V2.6m) into new one (>=V2.7). Call by:\n\
  cfgconv NAME[.cfg] [NS]\n\
NAME.cfg will be renamed to V2.6~NAME.cfg, NAME.cfg in new format written\n\
See also:\n\
  cfg2plb cfg2atm plbinfo plb2asc asc2plb plbconv plb2plb plb2cfg\n");
    exit(0); }

  initscroll(0);

  name=strdup(arg[1]);
  end=strend(name);
  if (end>name && !strcmp(end-4,".cfg")) end[-4]=0;

  fn=string("%s.cfg",name);
  VarOpen(fn,"r");
  prt("reading %s",fn);
  VarRead(&cfgkey,sizeof(cfgkey));
  VarRead(optionlist,sizeof(optionlist));
  VarRead(&N,sizeof(N));
  put2(cfgkey,N)
    if (cfgkey&1) ERROR(("%s.cfg this is already the new version",name))
  if (!(cfgkey&8)) ERROR(("sorry, such very old version of cfg-file is not supported"))
  if (VarFile.size!=8) ERROR(("cfg-file is not in double or has bad format"))
  VarRead(L,sizeof(L[0]));
  VarRead(&t,sizeof(t));
  VarRead(&h,sizeof(h));
  VarRead(&Entot,sizeof(Entot));
  put3(t,h,Entot)

  i=12-2;
  /* added in 2.4b */
  VarRead(&taudip,sizeof(taudip)); i--;
  VarRead(&RvdW,sizeof(RvdW)); i--;

  while (i--) VarRead(&zero,sizeof(zero));
  VarRead(L+1,sizeof(L[1]));
  VarRead(L+2,sizeof(L[2]));
  putv(L)
  
  ns=(VarFile.size-sizeof(struct oldcfg_s))/sizeof(vector)+1;
  if (cfgkey&2) {
    prt("POLAR configuration detected");
    polar++;
    ns/=2; }
  else
    RvdW=taudip;
  put(ns)
  if (ns<N) ERROR(("ns<N read from file"))
  if (narg>2) NS=atoi(arg[2]);
  else NS=ns;
  if (NS>ns) prt("WARNING: NS>ns specified, NS=ns=%d assumed",NS=ns);
  if (NS<ns) prt("WARNING: smaller NS=%d specified, configuration will be truncated",NS);
  
  oldsize=VarFile.size;
  while (oldsize==VarFile.size) {
    if (m>=9) ERROR(("too many cfg records (likely a bad format)"))
    sdsalloc(oldcfg[m],VarFile.size);
    VarRead(oldcfg[m],VarFile.size);
    m++; }

  VarClose();

  rename(fn,string("V2.6~%s.cfg",name));
  prt("%s.cfg renamed to V2.6~%s.cfg",name,name);

  VarOpen(fn,"w");

  cfgkey&=0xffffffe7;
  cfgkey|=1;
  put(cfgkey)
  VarPut(&cfgkey,sizeof(cfgkey));
  VarPut(optionlist,sizeof(optionlist));
  VarPut(&nspec,sizeof(nspec));
  if (nspec) {
    ERROR(("nspec: not finished...")) }
  //  loop (i,0,nspec) {
  //    simils.spec[i].N=spec[i]->N;
  //    simils.spec[i].ns=spec[i]->ns; }
  //  VarPut(simils.spec,sizeof(struct spec_s)*nspec);

  memset(&rec,0,sizeof(rec));

  rec.key=1;
  rec.intval=N;
  rec.vecval[0]=t; rec.vecval[1]=h; rec.vecval[2]=Entot;
  VarPut(&rec,sizeof(rec));

  rec.key=2;
  rec.intval=ns;
  loop (k,0,3) rec.vecval[k]=L[k];
  VarPut(&rec,sizeof(rec));

  if (polar) {
    rec.key=17;
    rec.vecval[0]=taudip;
    VarPut(&rec,sizeof(rec)); }

  newsize=oldsize+2*sizeof(vector);
  sdsalloczero(newcfg,newsize);
  put3(oldsize,newsize,m)

  loop (i,0,m) {
    newcfg->logs=oldcfg[i]->logs;
    memcpy(newcfg->rp,oldcfg[i]->rp,ns*sizeof(vector));
    loop (k,0,3) {
      newcfg->lambda[k]=0;
      if (i==0 && L[k]!=0) newcfg->lambda[k]=log(L[k]); }
    VarPut(newcfg,newcfg->size); }

  VarClose();

  prt("%s in format of cook* V2.7 written",fn);

  return 0;
}
