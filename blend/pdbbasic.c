#include "ground.h"
#include "options.h"

#include "pdbbasic.h"

vector undefvector={UNDEFATOM,UNDEFATOM,UNDEFATOM};

residue_t *reshead;

char pdbheader[82],pdbcompnd[82],optioninfo[82];
char *molname; /* mol name without extensions */
char *pdbname; /* pdb name without extension */

FILE *openfile(char *ext,char *mode,int chk) /********************* openfile */
/*
  open file molname.ext with given mode, (disaster if chk and no file)
*/
{
  char fn[FNLEN];
  FILE *f;

  if (!molname || !molname[0]) return 0;

  strcpy(fn,molname);
  strcat(fn,".");
  strcat(fn,ext);
  if (option('v')) prt("%s %s",mode[0]=='r'?"opening":"writing",fn);
  f=fopen(fn,mode);
  if (!f && chk) DISASTER(("cannot open %s",fn))

  return f;
}     

/* SSBOND command from the PDB file */
ssbond_t *sshead;

void includeSS(ssbond_t *ss) /************************************ includeSS */
/* includes ss into list sshead->... if not already there */
{
  ssbond_t *s;
  
  looplist (s,sshead)
    if (
      ( s->chain[0]==ss->chain[0] && s->chain[1]==ss->chain[1]
     && s->resno[0]==ss->resno[0] && s->resno[1]==ss->resno[1] )
     ||
      ( s->chain[0]==ss->chain[1] && s->chain[1]==ss->chain[0]
     && s->resno[0]==ss->resno[1] && s->resno[1]==ss->resno[0] ) ) {
      if (option('v')) prt("double defined SS bond: CYS %c %d - CYS %c %d",
                            s->chain[0],s->resno[0],s->chain[1],s->resno[1]);
      free(ss);
      return; }
      
  ss->site[0]=ss->site[1]=NULL;
  ss->next=sshead;
  sshead=ss;
}

/* CONECT command from the PDB file */
connect_t *connecthead;

/* reconnect patches */
reconnect_t *reconnecthead;

char *removedigits(char *id) /********************************* removedigits */
/*
  returns string id with leading and trailing digits removed
  1HA2 -> HA
  1C -> C
  C2 -> C
  TWO calls to removedigits supported without static string clash
  ASCII code assumed !!!
  NEW (8/87) does not do this:
  C11 -> C
*/
{
  static char str1[STRLEN],str2[STRLEN],*w;
  char *c,*ret;

  if (w==str1) w=str2; else w=str1;

  if (strlen(id)>=STRLEN) ERROR(("too long string"))

  strcpy(w,id); c=w;

  while (*c>='0' && *c<='9') c++;
  ret=c;
  while (*c>='A') c++;
  /* NEW: */
  if (strlen(c)<2) *c=0;
  
  return ret;
}

/* omitted sites (are in *.pdb but not supported by *.rsd), see also -i */
omit_t *omithead;


void prtres(residue_t *res) /**************************************** prtres */
{
  if (res)
    prt_("residue %4s %c %d %c",res->resnm,res->chain,res->resno,res->resins);
  else
    prts_("NULL residue");
}


site_t *findsite(int n) /****************************************** findsite */
/*
  finds site with given number n, returns NULL if not found
*/
{
  residue_t *res;
  site_t *s;
  
  looplist (res,reshead)
    looplist (s,res->site)
      if (s->n==n) return s;
  
  return NULL;
}

void checkneighbors(void) /******************************** checkneighbors */
/*
  check of consistency of neighbor tables
*/
{
  residue_t *res;
  site_t *s,*si;
  int i,j,n;

  looplist (res,reshead)
    looplist (s,res->site)
      loop (i,0,s->nnbr) {
        n=s->nbr[i];
        si=findsite(n);
        if (!si) {
          prtres(res);
          ERROR(("site %d %s : neighbor %d not found",s->n,s->ids->id,n))
          continue; }
        if (si->n!=n) ERROR(("internal"))
        loop (j,0,si->nnbr) if (s->n==si->nbr[j]) goto OK;
        prtres(res);
        ERROR(("site %d %s : neighbor %d %s does not have back reference",
               s->n,s->ids->id,n,si->ids->id))
        OK:; }
  if (option('v')) prt("neighbors checked");
}

#if 1
/* DEBUGGING TOOLS */
void prtsitelist(char *key,site_t *s0) /************************ prtsitelist */
{
  site_t *s;
  struct id_s *id;
  int i;
 
  prt("--- site list %s ---",key);
  looplist (s,s0) {
    prt_("%3d",s->n);
    looplist (id,s->ids) prt_("%c%4s",", "[id==s->ids],id->id);
    if (s->patch[0]) prt_(":%-4s",s->patch);
    else prt_("     ");
    prt_(" %4s %7.4f %d  %d",s->type,s->charge,s->chir,s->nnbr);
    loop (i,0,s->nnbr) prt_(" %d",s->nbr[i]);
    prt(s->unique?" U":" NU"); }
}
#endif

char *blendpath;
FILE *file;
char line[LINELEN];
char fn[FNLEN];
const char *rsdext=".rsd"; /* now constant - see subdir */
char *rsddir;
char *parameter_set;

char *getrsddir(void) /****************************************** getrsddir */
{
  FILE *par;
  char *dir,*tok;

  if (!(blendpath=getenv("BLENDPATH"))) blendpath="";
  if (*blendpath) if (*strlast(blendpath) != *SLASH) strcat(blendpath,SLASH);

  strcpy(fn,blendpath); strcat(fn,parameter_set); strcat(fn,".par");
  par=fopen(fn,"rt");

  if (par)
    while (fgets(line,LINELEN,par)) {
      tok=strtok(line," \t\n=");
     
      if (tok && !strcmp(tok,"rsddir")) {
        dir=strtok(NULL," \t\n="); 
        if (!dir) ERROR(("wrong rsddir statement in file %s",fn))
                    return dir; } }
  else
    WARNING(("no parameter file %s (check BLENDPATH)",fn))

  prt("no rsddir statement in parameter file %s",fn);
  
  return parameter_set;
}

void freeid(site_t *s) /********************************************* freeid */
/*
  free the s->ids list
*/
{
  struct id_s *id,*nextid;

  for (id=s->ids; id; id=nextid) {
    nextid=id->next;
    //    fprintf(stderr,"free(%p) \"%s\" %p\n",id,id->id,nextid);
    free(id);
  }
  s->ids=NULL;
}

void appendid(site_t *s, char *id) /******************************* appendid */
/*
  appends site id to the s->ids list (to the end)
  these id's are in the rsd-files separated by , 
  and represent equivalent names (aliases)
*/
{
  struct id_s **idp,*ide;

  for (idp=&(s->ids); ; idp=&(ide->next)) {
    ide=*idp;
    if (ide==NULL) {
      alloconezero(ide);
      *idp=ide;
      strcpy(ide->id,id);
      break; } }
}
