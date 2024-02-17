/*
  %%%%%%%%%%%%%%%%%%%%%%% reading/writing rsd/mol/cfg-files %%%%%%%%%%%%%%%%%%%
*/

#include "ground.h"
#include "options.h"
#include "pdbbasic.h"
#include "pdbmol.h"
#include "pdbrep.h"

static char *mygetline(void) /************************************** getline */
/***
  gets one `line' form `file', skipping comments (lines beginning with !)
  returns line on success, otherwise NULL
  line length is tested (max LINELEN) and no LF is appended
  comment lines may be longer than LINELEN
  characters after the last LF until EOF are ignored
***/
{
  char *c=line, *cmax=line+(LINELEN-1);
  int i;

  line[0]=0;

 nextline:

  /* 1st character in the line */
  if ( (i=fgetc(file))<0 ) return NULL;

  if (i=='!') {
    /* skip comment */
    do {
      if ( (i=fgetc(file))<0 ) return NULL;
      } while (i!='\n');
    goto nextline; }

  while (i!='\n') {
    if (c>=cmax)
      DISASTER(("line too long while reading %s (check, or increase LINELEN)",fn))
    *c++=i;
    if ( (i=fgetc(file))<0 ) return NULL; }

  *c=0;

  return line;
}

static site_t *allocsite(int ns,site_t *template) /*************** allocsite */
/*
  allocates array[ns] of sites incl. nbr[option('l')] and list structure
    and initializes them by 0
  if (template) then copies all from the template (of the same structure)
*/
{
  site_t *site;
  int i;

  if (!ns) return NULL;
  allocarrayzero(site,ns);
  loop (i,0,ns) site[i].line=-1; /* undef. line of pdb */
  if (template) memcpy(site,template,ns*sizeof(site_t));
  alloc(site[0].nbr,ns*sizeof(int)*option('l'));
  if (template) memcpy(site[0].nbr,template[0].nbr,ns*sizeof(int)*option('l'));
  loop (i,0,ns) {
    site[i].next=site+i+1;
    site[i].nbr=site[0].nbr+i*option('l'); }
  site[ns-1].next=NULL;

  return site;
}


/* template for storing read *.rsd files */

typedef struct templs_s {
  struct templs_s *next;
  site_t *site;
  int ns;
  enum moltype_e moltype;
  char name[8]; } templs_t;

static templs_t *tshead;

enum moltype_e test_cterp,test_nterp,test_patch;

site_t *readRSD(char *name,enum moltype_e *moltype_p,int resno,char chain)
/*                                                       *********** readRSD *
  reads one residue in `strict mol-format', namely:

    <type> <NAME(s)>

    number_of_atoms = <number>

    sites
    <table of sites (see blend manual)>
  
  exception: XX:YY in id-field means replacement (patch):
    XX is new atom id, YY is id from *.rsd file to be replaced
  <type> is one of {aminoacid,water,molecule,patch,nter,cter}
  if moltype_p=is one of &test_{patch,nter,cter}, type is tested
  line `parameter_set = <name>' before `number_of_atoms = <number>'
    is accepted but ignored (this format is produced by `blend')

  resno and chain are used to make name replacement if given by select0
    (see sel-file)
  resno=-1 means no replacement
*/
{
  char lcname[8];
  int i,j,ns,indx;
  char *tok,*c,*lc;
  site_t *site,*s,*s0;
  enum moltype_e moltype;
  templs_t *ts;
  select_t *select;
  struct id_s *id,*idj;

  /* name to lowercase */
  for (c=name,lc=lcname; *c; c++,lc++) *lc=tolower(*c);
  *lc=0;

  if (option('v')>1) prt("reading %s",name);

  if (option('n')) {
    /* choose neutral residues (option('n')<0 => with counterions) */
    c=NULL;
    if (!strcmp(name,"ASP")) c = option('n')>0 ? "asph" : "aspna";
    if (!strcmp(name,"GLU")) c = option('n')>0 ? "gluh" : "gluna";
    if (!strcmp(name,"HIS")) c = option('n')>0 ? "hisn" : "hiscl";
    if (!strcmp(name,"ARG")) c = option('n')>0 ? "argn" : "argcl";
    if (!strcmp(name,"LYS")) c = option('n')>0 ? "lysn" : "lyscl";
    if (c) strcpy(lcname,c); }

  /* this was default. now try use `select' info */
  if (resno>=0) 
    looplist (select,select0)
      if (!strcmp(name,select->resnm)
          && (select->resno==-1 || select->resno==resno)
          && (select->chain=='*' || select->chain==chain) ) {
        strcpy(lcname,select->resfn);
        prt("%s %d %c replaced by %s",name,resno,chain,lcname);
        break; }

  /* now lcname is the name in lowercase, correct charged state aminoacid */

  if (option('\\')) {

    /* looking in the list of residue templates */
    looplist (ts,tshead)
      if (!strcmp(ts->name,lcname)) {
        /* template found: copied from it */
        ns=ts->ns;
        site=allocsite(ns,ts->site);
        *moltype_p=ts->moltype;
        return site; }

    /* template not found */
    alloc(ts,sizeof(templs_t));
    ts->next=tshead; tshead=ts;
    strcpy(ts->name,lcname);
    /* ... and continue to read it from file */ }

  strcpy(fn,blendpath); strcat(fn,rsddir); strcat(fn,SLASH);
  strcat(fn,lcname); strcat(fn,rsdext);
  file=fopen(fn,"rt");
  if (!file) {
    ERROR(("residue file %s not found (check BLENDPATH)",fn))
    return NULL; }

  mygetline();
  tok=strtok(line," \t\n");
  moltype = !strcmp(tok,"aminoacid") ? AMINOACID :
            !strcmp(tok,"patch") ? PATCH :
            !strcmp(tok,"water") ? WATER :
            !strcmp(tok,"molecule") ? MOLECULE :
            !strcmp(tok,"nter") ? NTER :
            !strcmp(tok,"cter") ? CTER :
            UNKNOWN;
  if (moltype==UNKNOWN) ERROR(("%s: unknown type of rsd-file",fn))
  *moltype_p=moltype;

  /* patch now not used, only cter and nter */
  if (moltype!=PATCH && moltype_p==&test_patch) ERROR(("%s: patch expected",fn))

  if (moltype_p==&test_cterp) {
    if (moltype==PATCH) WARNING(("%s: rsd file - type `cter' expected\n\
*** (type `patch' deprecated)",fn))
    else if (moltype!=CTER) ERROR(("%s: bad rsd file - type `cter' expected",fn)) }

  if (moltype_p==&test_nterp) {
    if (moltype==PATCH) WARNING(("%s: rsd file - type `nter' expected\n\
*** (type `patch' deprecated)",fn))
    else if (moltype!=NTER) ERROR(("%s: bad rsd file - type `nter' expected",fn)) }

  do {
    if (!mygetline()) { ERROR(("%s format: unexpected EOF",fn)) return NULL; }
    } while (line[0]==0 || line[0]=='\n' || line[0]=='!');

  /* ignore line parameter_set=<name> */
  if (!strcmp(c=strtok(line," \t="),"parameter_set")) {
    mygetline(); c=strtok(line," \t="); }

  if (strcmp(c,"number_of_atoms")) {
    ERROR(("%s format: `number_of_atoms\' expected",fn)) return NULL; }
  ns=atoi(strtok(NULL," \n\t="));

  do {
    if (!mygetline()) { ERROR(("%s format: unexpected EOF",fn)) return NULL; }
    } while (line[0]==0 || line[0]=='\n');

  if (strcmp(strtok(line," \t\n"),"atoms")) {
    ERROR(("%s format: `atoms\' keyword expected",fn)) return NULL; }

  if (ns==0) ERROR(("%s: number_of_atoms=0",fn))

  site=allocsite(ns,NULL);

  loop (i,0,ns) {
    s=site+i;
    if (!mygetline()) {
      ERROR(("%s format: `atoms\' table too short",fn))
      break; }
    tok=line;
    indx=atoi(strtok(tok," \t"));
    if (i!=indx) {
      ERROR(("%s `atoms\' table format: numbering",fn))
      break; }
    s->n=i;
    tok=strtok(NULL," \t");
    if ( (c=strchr(tok,':')) ) {
      *c=0; strcpy(s->patch,c+1); }
    else
      s->patch[0]=0;

    /* the id field: now a comma-separated list allowed */
    s->ids=NULL;
    do {
      for (c=tok; *c && *c!=','; c++);
      if (*c==',') *c++=0;
      
      if (strlen(tok)>=STRLEN-1) {
        ERROR(("%s `atoms\' table format: id \"%s\" too long",fn,tok))
        break; }
      appendid(s,tok);
      tok=c;
      } while (*c);

    strcpy(s->type,strtok(NULL," \t"));
    s->charge=atof(strtok(NULL," \t"))*(option('q')/100.0);
    s->chir=atoi(strtok(NULL," \t"));
    if ( (s->nnbr=atoi(strtok(NULL," \t"))) > option('l') ) {
      ERROR(("%s: atom %s %s: too many bonds - check option -l!",
           fn,s->ids->id,s->type))
      s->nnbr=option('l'); }
    loop (j,0,s->nnbr) {
      c=strtok(NULL," \t");
      if (!c) ERROR(("%s: atom %s %s: wrong neighbor table", fn,s->ids->id,s->type))
      s->nbr[j]=atoi(c); }
    s->next=site+i+1; /* next in the list */
    s->r=undefvector; }

  site[ns-1].next=NULL; /* end of list */
  fclose(file);

  /* element uniqueness (i.e., after stripping digits, there is such) */
  loop (i,0,ns) {
    s0=site+i;
    s0->unique=1;
    loop (j,0,ns) if (i!=j) 
      looplist (id,s0->ids)
        looplist (idj,site[j].ids)
          if (!strcmp(removedigits(id->id),removedigits(idj->id)))
            s0->unique=0;
    }

  if (option('\\')) {
    /* to store the read residue as a template */
    ts->ns=ns;
    ts->site=allocsite(ns,site);
    ts->moltype=moltype; }

  return site;
}

static double totalcharge; /* returned by countsites */
static int maxval;         /* returned by countsites */
static float hdr[2]={0,-3};/* since V1.4a */
static float L[3];

static int countsites(void) /************************************ countsites */
{
  residue_t *res;
  site_t *s;
  int ns=0;

  totalcharge=0;
  maxval=0;

  looplist (res,reshead)
    looplist (s,res->site) {
      Max(maxval,s->nnbr)
      if (ns!=s->n) {
        ERROR(("numbering (pdb internal error): site->n=%d ns=%d",s->n,ns))
      ns=s->n; }
      ns++;
      totalcharge+=s->charge; }

  return ns;
}

void writeMOL(void) /********************************************** writeMOL */
/*
  Writes the molecule file *.mol
*/
{
  FILE *f;
  static char id[32];
  residue_t *res;
  site_t *s;
  int ns,j;

  f=openfile("mol","wt",1);

  if (f==NULL) {
    ERROR(("cannot open %s.mol for output",molname)) return; }

  /* compound name without leading spaces */
  if (pdbcompnd[0]) fprintf(f,"%s\n",pdbcompnd);
  else fprintf(f,"%s\n",molname);
  if (pdbheader[0]) fprintf(f,"! %s\n",pdbheader);
  fputs("\n! generated by `pdb\' V" VERSION "\n\n",f);
  fputs(optioninfo,f);

  ns=countsites();

  fprintf(f,"! total charge = %.3f\n\n",totalcharge);
  if (option('v')) {
    prt("total charge = %.3f%s",
        totalcharge,
        fabs(totalcharge)>0.0005 && option('n')==0 ? " - check option -n!" : "");
    prt("max valence = %d%s",maxval,maxval>4?" (?)":""); }

  if (fabs(totalcharge-(int)(totalcharge+16000.5)+16000)>0.0005)
    WARNING(("fractional total charge %f",totalcharge))

  fprintf(f,"parameter_set = %s\n",parameter_set);
  fprintf(f,"number_of_atoms = %d\n\n",ns);

  fputs("atoms\n! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n",f);

  looplist (res,reshead) 
    looplist (s,res->site) {
      if (res->chain==' ')
        sprintf(id,"%s%d%s",res->resnm,res->resno,s->ids->id);
      else
        sprintf(id,"%s%d%c%s",res->resnm,res->resno,tolower(res->chain),s->ids->id);
      fprintf(f,"%3i %-10s %-4s %7.4f %2i %4i ",
              s->n,id,s->type,s->charge,s->chir,s->nnbr);
      loop (j,0,s->nnbr) fprintf(f," %d",s->nbr[j]);
      if (fputc('\n',f)<0) ERROR(("writing %s.mol",molname)) }

  fputc('\n',f);
  if (fclose(f)) ERROR(("closing %s.mol",molname))
}

static void endian(void *c) /**************************************** endian */
/* changes endianess of 4 bytes pointed to by c */
{
  char x;

  x=((char*)c)[0]; ((char*)c)[0]=((char*)c)[3]; ((char*)c)[3]=x;
  x=((char*)c)[1]; ((char*)c)[1]=((char*)c)[2]; ((char*)c)[2]=x;
}

void write3D(char *mode, int opts) /******************************** write3D */
/*
  mode=="wt": write configuration in the 3 column (x y z) text format
  mode=="wb": write configuration in the float[ns][3] (binary) format
  (for reads, all structures must be prepared in advance)
*/
{
  char ext[4];
  FILE *f;
  residue_t *res;
  site_t *s;
  int i,mask,opt;

  strcpy(ext,"3d?");

  loopto (mask,1,2) if ( (opt=(abs(opts)&mask)*opts/abs(opts)) ) {

    if (abs(opt)==2) memcpy(ext,"pl",2);

    ext[2]=mode[1]; /* 'b' or 't' */
    f=openfile(ext,mode,1);

    if (abs(opt)==2) {
      hdr[0]=(float)(i=countsites());
      if (mode[1]=='b') {
        if (opt<0) endian(hdr);
        fwrite(hdr,sizeof(float),2,f); 
        fwrite(L,sizeof(float),3,f); }
      else {
        fprintf(f,"%d -3\n",i); 
        fprintf(f,"0 0 0\n"); } }

    looplist (res,reshead)
      looplist (s,res->site)
        if (mode[1]=='b') {
          if (opt<0) loop (i,0,3) endian(s->r+i);
          fwrite(s->r,sizeof(vector),1,f); }
        else
          fprintf(f,
            s->r==undefvector ? "%8.0f%8.0f%8.0f\n" : "%8.3f%8.3f%8.3f\n",
            s->r[0],s->r[1],s->r[2]);
    fclose(f); }
}

int retonerror=0;

int read3D(char *mode, int opt, int no) /**************************** read3D */
/*
  mode=="rt": read configuration in the 3 column (x y z) text format
  mode=="rb": read configuration in the float[ns][3] (binary) format
  no-th cfg is read (no=1 means 1st, no=-1 means the last)
  all structures must be prepared in advance!
  if retonerror then returns 1 on reading  error "no such frame"
*/
{
  char ext[4]="3d?";
  FILE *f;
  residue_t *res;
  site_t *s;
  int i=0,ns=countsites(),z,varL=0;

  if (abs(opt)>2) ERROR(("-b%d or -t%d ambiguous for reading cfg",opt,opt))
  if (mode[1]=='t' && no!=1) ERROR(("bad -p for text input - only -p1 supported"))

  if (abs(opt)==2) memcpy(ext,"pl",2);

  ext[2]=mode[1]; /* 'b' or 't' */
  f=openfile(ext,mode,1);

  if (abs(opt)==2) {
    if (mode[1]=='b') {
      if (fread(hdr,sizeof(float),2,f)!=2)
        ERROR((".%s file too short",ext))
      if (opt<0) endian(hdr);
      if (hdr[0]!=(float)ns)
        ERROR(("bad ns=%f in .%s file (expected %d)",hdr[0],ext,ns))
      varL=hdr[1]==-3; }
    else {
      char line[64];

      if (!fgets(line,64,f))
        ERROR(("reading .%s, line %d (header)",ext,i))
      i++;
      if (sscanf(line,"%d%d",&ns,&z)!=2)
        ERROR((".%s line %d (header):format",ext,i));
      if (ns!=countsites())
        ERROR(("bad ns=%d in .%s file (expected %d)",ns,ext,countsites()))
      varL=z==-3; } }

  if (mode[1]=='b') {
    if (option('v')) prt("reading %d%s%s configuration",
                         abs(no),
                         abs(no)==1?"st":abs(no)==2?"nd":abs(no)==3?"rd":"th",
                         no>0?"":" last");
    if (fseek(f,(long)(ns+varL)*sizeof(vector)*(no-(no>0)),no<0?SEEK_END:SEEK_CUR)) {
      if (retonerror) return 1;
      ERROR(("cannot find configuration - too short file?")) } }

  looplist (res,reshead)
    looplist (s,res->site) {
      if (mode[1]=='b') {
        if (varL)
          if (!fread(L,sizeof(vector),1,f))
            ERROR(("reading .%s, missing L",ext))
        if (fread(s->r,sizeof(vector),1,f)!=1) {
	  if (retonerror) return 1;
          ERROR(("reading .%s, rec %d",ext,i)) }
        if (opt<0) loop (i,0,3) endian(s->r+i); }
      else {
        char line[64];
        
        if (varL)
          if (!fgets(line,64,f))
            ERROR(("reading .%s, missing L",ext))  
        if (!fgets(line,64,f)) {
	  if (retonerror) return 1;
	  ERROR(("reading .%s, line %d",ext,i)) }
        if (sscanf(line,"%f%f%f",&s->r[0],&s->r[1],&s->r[2])!=3) {
	  if (retonerror) return 1;
          ERROR((".%s line %d:format",ext,i)); } }
      i++; }
  fclose(f);

  return 0;
}
