/* cc -O2 -o makemake makemake.c
*/

#define METAMAKE "metamake.mmk" /* source file */

/*** makemake.c *** (c) J.Kolafa 07/1993--2024 ***

Generate unix-style `makefile' from `metamake.mmk'.
Simplified version with all DOS/EMX legacy features removed.

History:
  02/2024: !private added because of easier GitHub compatibility
           `metamake' renamed to `metamake.mmk'
  06/2022: DOS support removed, spaces in names allowed
  02/2017: multiply defined files reported, the LAST one is active
           (in previous versions, the FIRST one was active)
           thus, also the order of -I changed
  10/2010: linux mode made the default
  11/2009: print PWD added
  10/2009: getline renamed to getlongline (because of getline() in stdio)
  05/2009: target-specific options OPT,LIBOPT added
  02/2004: linking via CPLUS if at least 1 C++ module
           (because new gcc needs this)
           some string sizes increased
  07/2001: !include statements added (because of MACSIMUS needs), can be nested
  09/2000: bug fixed (conditionals with | &)
  11/1999: small bug fixed
  11/1997: conditionals (!if !else !endif ...) added,
           turboc.cfg automatically generated (because of NSK needs)
  09/1996: $(OPT) now after -I (can use -I- in gcc)
  01/1996: C++ support
  09/1995: Turbo make format, all formats integrated
  06/1995: aesthetical bug (too many empty lines printed) fixed
  01/1995: bug in getline() fixed
  11/1993: update
  07/1993: 1st version

USING `makemake'
^^^^^^^^^^^^^^^^
Synopsis:
  makemake DEF [DEF ...]
where DEF is name (identifier) to define.
These DEFines may be used in the `metamake' file in conditional
processing (in the similar way as #ifdef NAME ... in the C preprocessor,
see below for details).

Example:
  makemake high

In addition:
  filelist = will generate file "filelist" with a list of all files
             (relatively to !home, alphabetically)
  makefile = all modules will depend also on the makefile

Arrangement of work:
^^^^^^^^^^^^^^^^^^^^

* There are several directories containing source files (*.c), header files
  (*.h), and other files #included while compilation.  These directories are
  subdirectories (sub-subdirectories...) of one `home' directory (specified by
  statement `!home').  These directories (without the `home' part) are listed
  in statement `!dir'.  The order of directories corresponds to the order in
  which #include files are looked for.  Spaces (and some other special
  characters) are now allowed in `home'.

* The `#include "NAME"' statements in the source files should contain only file
  names, not / nor .  `makemake' finds the directory containing the file
  automatically but it is confused if qualified names are given.

* Executables are obtained by linking from object modules (extensions given by
  statement `!objext') and libraries.  All executables and object modules are
  placed in the working directory, irrespective of the location of the
  corresponding sources.  (This arrangement allows to use the same source code
  with different options for one project and with another options for another
  project residing in another directory).


metamake
^^^^^^^^
`makemake' reads file `metamake.mmk' and processes it, generating `makefile'.
It recognizes lines starting with `!' in the 1st column as control codes or
statements.  All other lines are copied unchanged to the output.
`metamake' should contain definitions of the following make macros (variables):
  CC = <name of the C compiler>
  CPLUS = <name of the C++ compiler (if needed)>; the C++ extension is .cc
  OPT = <compile options (used for compiling *.c to object files)>
  LIBOPT = <libraries to link (e.g., `-lm')>

Metamake commands:
^^^^^^^^^^^^^^^^^^
* `list' means items separated by spaces
* lines may be continued in the makefile-like style by trailing `\'
* spaces around `=' are optional

General commands:

!! <any comment>
!home = <fully qualified directory name>
!dir = <list of subdirectory names> max 9, must end with the working directory
!include <include file> "" around name optional, missing=error, nesting allowed
!private <include file> as above, missing=silently skipped
!objext = <extension for object files> default = `.o'
!ignore = <list of files not to include into dependencies> without paths
!treat = <list of files to undo 1 previous !ignore>
!<exe1> : [MOD:] [MOD:] <list of object file names without extension>
!<exe2> : [MOD:] [MOD:] <list of object file names without extension>
etc.

The modifier MOD serves to replace option $(OPT) and $(LIBOPT).
MOD = OPT=<compile option to replace $(OPT) in compiling the target>
      LIBOPT=<linker option to replace $(LIBOPT)>
NOTE: `:' delimiting MODs must be in a single line, spaces around are allowed
          (list of names may then span more lines)
Example (test will be compiled with -O0 and linked with -lm;
         ground is compiled with the OPT set in the file):
  test : OPT=-O0 : LIBOPT=-lm : test ground

Conditionals:

!if NAME     (max nesting level=8)
!ifdef NAME  (= !if NAME)
!else
!endif
!define NAME
!def NAME    (= !define NAME)
!undef NAME

These conditionals have similar meaning as in the C preprocessor except:
* NAMEs cannot have value, they are only defined or undefined
* Limited expressions only are allowed:
!if NAME1 | NAME2 | NAME3 ... (= !if NAME1 || NAME2 || NAME3 ... )
!if NAME1 & NAME2 & NAME3 ... (= !if NAME1 && NAME2 && NAME3 ... )
* `|' and `&' cannot be combined in one statement, parentheses not supported
* `!' in front of variable name means `not' (no space after !)

Example:
!if high & !emulated

Processing:
^^^^^^^^^^^
* The directory names (statement `!dir') are relative to `!home', the value of
  `home' (with `/' added) is prepended before the directory names to obtain
  fully qualified names.

* Dependencies and building (=linking) rules for exe's are created.  The
  modules listed accept extension given by `!objext'.

* Dependencies and building rules for all object files are created.  Each
  object file encountered in the list is assumed to depend on the corresponding
  *.c file and consequently to all files that are included by `#include "NAME"'
  (with # in the 1st column) from it (system files included by `#include
  <NAME>' are not considered).  The `#include' statements may be nested.
  BUG/FEATURE: All `#include "NAME"' statements found are used even though they
  are within comments or not active because of #if clauses.  This strategy is
  robust against compiler options -D and changed #defines.  To hide #include
  from `makemake' prepend it by a space.

* Output make-file has lines at most 80 characters long, longer are
  continued by `\`

Diagnostics:
^^^^^^^^^^^^
* If file (*.c or #included) is not found, $(D?) instead of the
  correct macro for the directory name is written and a message is
  printed.

* If a module is found in several `!dir', the last one is active and a warning
  is printed.

Bugs/caveats/features:
^^^^^^^^^^^^^^^^^^^^^^
* `#include "NAME"' within a comment in source files or within false
  conditionals is still treated as valid dependency (provided that #
  is in 1st column).

* The output listing with object files is in reverse order than one
  would expect.
***/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

char *myfgets(char *s, int size, FILE *stream) /******************** myfgets */
/* legacy/compatibility: version of fgets with converting DOS -> UNIX */
{
  char *ret=fgets(s,size,stream),*c;

  if (ret)
    for (c=s; *c; c++)
      if (*c=='\r') {
        if (c[1]=='\n')
          if (c[2]) *c=' ',c[1]=' '; /* (should not happen...) */
          else *c='\n',c[1]=0; /* CR LF -> LF */
        else
          if (c[1]) *c=' '; /* (should not happen...) */
          else *c=0; /* CR -> LF */ }

  return ret;
}

#define MAXFILELIST 2048
#define BUFLEN 4096
#define FNLEN 256
#define NDIRS (9+1) /* max 9 directories handled */
#define NINCL 6 /* # of nest levels for !include statements */

#define LINELEN 77 /* for output */

char *targetname;

struct namelist_s {
  char name[FNLEN];
  struct namelist_s *next; };

struct project_s {
  char name[FNLEN];
  struct namelist_s *namelist;
  struct project_s *next;
} *project0;

struct objlist_s {
  char name[FNLEN];
  struct objlist_s *next;
} *objlist0;

FILE *in[NINCL],*out;
char *inname[NINCL];
char inpathname[BUFLEN];

char key;
char home[FNLEN],eschome[FNLEN];
int len=0;
int notfound=0;

int ndirs=1;
char dirs[NDIRS][FNLEN];
char objext[FNLEN];

struct file_s {
  int type;
  int cplus; /* set if C++ */
  int todo;
  char name[FNLEN];
  char *xopt;
  struct list_s *dep;
  struct file_s *next;
} *f0;

struct list_s {
  struct file_s *f;
  struct list_s *next; };

struct ignore_s {
  char name[FNLEN];
  struct ignore_s *next;
} *ignore0;

void escape(char *to,char *from) /*********************************** escape */
/* https://stackoverflow.com/questions/7654386/how-do-i-properly-escape-data-for-a-makefile
   space and # are escaped
   $ is doubled
   other characters are unchanged
   backslash and newline are likely treated incorrectly
   many special and control characters are discouraged
*/
{
  char *c,*t;

  for (c=from,t=to; *c; c++ ) {
    if (t-to>=FNLEN-2) {
      fprintf(stderr,"makemake: too long home\n");
      exit(4); }
    switch (*c) {
      case '$': *t++=*c; break;
      case ' ':
      case '\\':
      case '#': *t++='\\'; }
    *t++=*c; }
}

int nottoignore(char *x) /************************************* nottoignore */
{
  struct ignore_s *ig;

  for (ig=ignore0; ig; ig=ig->next)
    if (!strcmp(ig->name,x)) {
      fprintf(stdout,"%s IGNORED\n",x);
      return 0; /* found in the list of ignored */ }

  return 1; /* not found */
}

void scopy_l(char *a,char *b,int l) /****************************** scopy_l */
{
  if (strlen(b)+5>=FNLEN) {
    fprintf(stderr,"makemake: too long name: \"%s\" in line %d\n",b,l);
    exit(4); }
  strcpy(a,b);
}

#define scopy(A,B) scopy_l(A,B,__LINE__)

void *myalloc(unsigned size) /************************************* myalloc */
{
  void *v=malloc(size);

  if (v!=NULL) return v;
  fprintf(stderr,"makemake: cannot malloc\n");
  exit(4);
}

char *mystrdup(const char *s) /************************************ mystrdup */
{
  char *n=strdup(s);

  if (n) return n;
  fprintf(stderr,"makemake: cannot malloc\n");
  exit(4);
}

/*
   Scan one white-separated item with possible other separators around.
   Spaces not allowed inside the item.
   Pointer to item returned, end of parsing returned in *C.
*/
char *scan(char **C,char sep[2]) /************************************ scan */
{
  char *c=*C,*t=c;

  key=' ';
  while (*c==' ' || *c=='\t') c++;
  if (*c==0) return NULL;
  t=c;
  while (*c!=' ' && *c!='\t' && *c!=sep[0] && *c!=sep[1] && *c!=0) c++;
  while (*c==' ' || *c=='\t') *c++=0;
  if (*c==sep[0] || *c==sep[1]) { key=*c; *c++=0; }
  *C=c;

  return t;
} /* scan */

/*
   Scan one white-separated item.
   The item may be enclosed in ""; then, whites etc. allowed inside.
   Pointer to item returned, end of parsing returned in *C.
*/
char *qscan(char **C) /*********************************************** qscan */
{
  char *c=*C,*b=c;

  while (*c==' ' || *c=='\t') c++;
  if (*c==0) return NULL;
  b=c;
  if (*c=='\"') {
    b=c+1;
    c=strchr(b,'\"');
    if (!c) {
      fprintf(stderr,"makemake: missing right \"\n");
      exit(4); }
    *c++=0; }
  else {
    while (*c!=' ' && *c!='\t' && *c!=0) c++; }

  while (*c && *c<=' ') *c++=0;
  *C=c;

  return b;
} /* scan */

void putobj(char *name) /******************************************* putobj */
/* prints name of 1 obj file */
{
  int l=strlen(name)+strlen(objext)+1;
  struct objlist_s *ol;

  if ( (len+=l) > LINELEN ) {
    fprintf(out," \\\n "); len=l+1; }
  fprintf(out," %s%s",name,objext);
}

int omitlastpath;
int depnames; /* to replace `targetname ()' in dependency by `targetname
                 (targetname.c)' because Turbo C does not like empty list */
char *depname=""; /* stored targetname */

char **filelist;
int nfilelist;

void includefilelist(char *di,char*name) /****************** includefilelist */
{
  char c[2*FNLEN];
  int i;

  sprintf(c,"%s/%s",di,name);
  for (i=0; i<nfilelist; i++) if (!strcmp(c,filelist[i])) return;

  if (nfilelist>=MAXFILELIST)
    fprintf(stderr,"makemake: too many files in filelist\n");
  filelist[nfilelist++]=mystrdup(c);
}

int cmpfilelist(const void *a,const void *b) /****************** cmpfilelist */
{
  return strcmp(*(char**)(a),*(char**)(b));
}

void putitem(char *name,int type) /********************************* putitem */
/* prints 1 item (typically file name) */
{
  int l,last=(type==ndirs-1)&omitlastpath;
  if (!type) return;

  if (filelist) includefilelist(dirs[type],name);

  l=strlen(name)+7;
  if ( (len+=l) > LINELEN ) {
    fprintf(out," \\\n "); len=l+1; }
  if (type>0)
    fprintf(out," $(D%d)/%s",type,name);
  else {
    fprintf(stderr,"makemake: >>> $(D?)/%s\n",name);
    fprintf(out," $(D?)/%s",name); }
}

int containscplus=0;
int warnings=0;

/* reentrant */
struct file_s *dependency(char *name,int obj,char *xopt) /******* dependency */
{
  struct file_s *f;
  struct list_s *l;
  char fn[FNLEN];

  for (f=f0; f; f=f->next)
    if (!strcmp(f->name,name)) return f; /* found in the list */

  /* not found */
  fprintf(stdout,"%s\n",name);
  f=myalloc(sizeof(struct file_s));
  f->type=-1;
  f->cplus=0;
  if (xopt) f->xopt=strdup(xopt);
  else f->xopt="$(OPT)";
  strcpy(f->name,name);
  f->next=f0; f0=f;
  f->dep=NULL;

  if (obj) {
    /* *.o file that has *.c direct dependency only (.o is not in the name) */
    scopy(fn,name); strcat(fn,".c");
    l=myalloc(sizeof(struct list_s));
    f->type=0;
    f->dep=l;
    l->next=NULL;
    l->f=dependency(fn,0,NULL);
    if (obj==2) {
      fprintf(stdout,"C++ ");
      containscplus++;
      l->f->cplus=2; }
    else if (l->f->type==-1) {
      strcat(fn,"c");
      fprintf(stdout,"C++ ");
      l->f=dependency(fn,0,NULL);
      containscplus++;
      l->f->cplus=1; } }

  else {
    /* search for file in all directories
       new 2/2017: check multiple directories, select the last one */
    FILE *fi=NULL,*fi1;
    char fullfn[FNLEN*2+128];
    char *x;
    int i,ii;

    for (i=1; i<ndirs; i++) {
      sprintf(fullfn,"%s%s/%s",home,dirs[i],name);
      fi1=fopen(fullfn,"rt");
      if (fi1) {
        ii=i;
        if (fi!=NULL) {
          fprintf(stdout,"%s MULTIPLY DEFINED (previous ignored)\n",name);
          warnings++;
          fclose(fi); }
        fi=fi1; } }

    if (fi!=NULL) {
      char *r=myalloc(BUFLEN);

      f->type=ii;
      while (myfgets(r,BUFLEN,fi)) {
        /* # in C-code in other than 1st column ignored */
        if (r[0]=='#') {
          char *c=r+1, *token=scan(&c,"$$");

          if (!strcmp(token,"include") && *c=='\"') {
            x=qscan(&c);
            if (x && nottoignore(x)) {
              strcpy(fn,x);
              l=myalloc(sizeof(struct list_s));
              l->next=f->dep;
              f->dep=l;
              l->f=dependency(fn,0,NULL); } } } }
      free(r);
      fclose(fi);
      return f; } }

  return f;
} /* dependency */

void alldep(struct file_s *f) /************************************** alldep */
{
  struct list_s *l;

  if (f) if (f->todo) {
    putitem(f->name,f->type);
    f->todo=0;
    for (l=f->dep; l; l=l->next) alldep(l->f); }
} /* alldep */

static char buf[BUFLEN]; /* for line incl. concatenation */

char *getlongline(void) /*************************************** getlongline */
/* reads one line incl. continuation lines ended by backslash */
{
  char *c=buf;

  if (!in[0]) return NULL;

  again:
    if (myfgets(c,(buf+BUFLEN)-c,in[0])==NULL) return NULL;
    while (*c!='\n' && *c!=0 && c<buf+BUFLEN) c++;
    if (c>=buf+BUFLEN) {
      fprintf(stderr,"makemake: concatenated line too long\n");
      exit(6); }
    if (c>buf) c--;
    while (*c==' ') {
      if (c==buf) { *c=0; goto retbuf; }
      c--; }
    if (*c=='\\') { *c++=' '; goto again; }

  *++c=0;

  retbuf:
  if (buf[0]=='\n' && buf[1]==0) buf[0]=0; /* oops ... empty line "\n" -> "" */

  return buf;
} /* getlongline */

struct def_s {
  struct def_s *next;
  char name[1]; } *def0,*def1;

char *name0="?",*name1="?";

int isdef(char *name) /*********************************************** isdef */
{
  struct def_s *def;

  for (def=def0; def; def=def->next) if (!strcmp(def->name,name)) {
    name1=def->name;
    def1=def;
    return 1; }

  return 0;
}

int define(char *name) /********************************************* define */
{
  struct def_s *def;

  if (isdef(name)) return 1;

  def=myalloc(sizeof(struct def_s)+strlen(name));
  def->next=def0; def0=def;
  strcpy(def->name,name);

  return 0;
}

int undef(char *name) /*********************************************** undef */
{
  if (isdef(name)) {
    def1->name[0]=1; /* very cheap... */
    return 1; }

  return 0;
}

char activestat[8]={1}; /* 1 if active (=defined) */
int ifnest;

int getactive(void) /********************************************* getactive */
{
  int i;

  for (i=0; i<=ifnest; i++) if (!activestat[i]) return 0;
  return 1;
}

void getxopt(char **opt,char *target,char *s) /********************* getxopt */
{
  char *c;

  if (!memcmp(target,s,strlen(target))) {
    c=strchr(s,'=');
    if (!c) {
      fprintf(stderr,"makemake: %s not followed by =\n",target);
      exit(12); }
    c++;
    while (*c && *c<=' ') c++;
    *opt=c; }
}

void makepath(void) /********************************************** makepath */
{
  int i;
  char *x,*sl;

  inpathname[0]=0;
  for (i=NINCL-1; i>=0; i--) if (inname[i]) {
    if (inname[i][0]=='/')
      strcpy(inpathname,inname[i]);
    else {
      if (strlen(inpathname)+strlen(inname[i])>1023) {
        fprintf(stderr,"makemake: too long !include arguments\n");
        exit(4); }
      if (i==0 || strchr(inname[i],'/')) strcat(inpathname,inname[i]); }
    if (i && strchr(inname[i],'/')) {
      sl=inpathname;
      for (x=inpathname; *x; x++) if (*x=='/') sl=x;
      if (sl!=inpathname) sl[1]=0; } }
}

int main(int narg,char **arg) /**************************************** main */
{
  char *c,*token,*x,*makefile="";
  struct file_s *ff,*f;
  int i,nobj=0,active=1,ifop,l;
  /* ifnest= !if nesting level, `active' shadows activestat[ifnest] */
  char ifstat[8]; /* 1=!if, 0=!else */
  struct def_s *def;

  struct project_s *project;
  struct namelist_s *namelist;

  fprintf(stderr,"\
MAKEMAKE 02/2024: MACSIMUS makefile generator\n");

  if (narg<2) {
    fprintf(stderr,"\
Generate `makefile' from `metamake'. Call by:\n\
  makemake DEF [DEF..]\n\
where DEFs are keywords (conditionals) handled in `metamake'.\n\
Special DEFs:\n\
  filelist = generate `filelist' with all needed files\n\
  linux = (legacy default)\n\
  makefile = all modules will depend also on the makefile\n\
             => recompile if the makefile changes\n\
Example:\n\
  makemake lj polar gcc\n"
      );
    exit(4); }
  else
    for (i=1; i<narg; i++) define(arg[i]);

  fprintf(stderr,"makemake: PWD=%s\n",getenv("PWD"));

  if (isdef("filelist")) filelist=malloc(MAXFILELIST*sizeof(filelist[0]));
  if (isdef("makefile")) makefile=" makefile";

  if (in[0]==NULL) in[0]=fopen(METAMAKE,"rt");
  if (in[0]==NULL) {
    fprintf(stderr,"makemake: file \""METAMAKE"\" not found\n");
    exit(4); }

  inname[0]=mystrdup("metamake"); /* irrelevant -- only (none) subdirectory needed */

  out=fopen("makefile","wt");
  if (out==NULL) {
    fprintf(stderr,"makemake: cannot write out\n");
    exit(4); }

  strcpy(objext,".o");

  fputs("\
### this makefile was generated from " METAMAKE "using makemake ###\n\
#defined: ",out);
  for (def=def0; def; def=def->next) fprintf(out," %s",def->name);
  fprintf(out,"\n");

  for (;;) {
    while (!(c=getlongline())) {
      int i;

      if (in[0]==NULL) break;
      free(inname[0]);
      for (i=0; i<NINCL-1; i++) {
        in[i]=in[i+1];
        inname[i]=inname[i+1]; }
      in[NINCL-1]=NULL;
      inname[NINCL-1]=NULL; }

    if (!c) break; /* EOF on all !include files */

    if (*c!='!') {
      /* not !-line just copied */
      if (active) fprintf(out,"%s\n",c); }
    else {
      /* !-line */
      char *sem,*xopt,*xlibopt;

      c++;
      if (*c!='!' && (token=scan(&c,":=")) )
        switch (key) {

          case ' ': /* keyword not followed by : nor = */

            if (!strcmp(token,"include") || !strcmp(token,"private")) {
              if (in[NINCL-1]) {
                fprintf(stderr,"makemake: !include %s: too many nested levels\n",c);
                exit(3); }

              for (i=NINCL-1; i>0; i--) {
                in[i]=in[i-1];
                inname[i]=inname[i-1]; }
              x=qscan(&c);
              inname[0]=mystrdup(x);
              makepath();
              if ( (in[0]=fopen(inpathname,"rt")) )
                fprintf(stdout,"including \"%s\"\n",inpathname);
              else if (token[0]=='i') {
                fprintf(stderr,"makemake: cannot open !include file \"%s\"\n",inpathname);
                exit(3); }
              else 
                fprintf(stderr,"makemake: NOTE: no private \"%s\" to include\n",inpathname); }
      
            else if (!strcmp(token,"if") || !strcmp(token,"ifdef")) {
              if (ifnest>=7) {
                fprintf(stderr,"makemake: too many nested !if's\n");
                exit(12); }
              ifnest++;
              ifstat[ifnest]=1;
              ifop=0;
              while ( (x=scan(&c,"&|")) ) if (*x) {
                  if (ifop && key!=' ' && key!=ifop) {
                    fprintf(stderr,"makemake: !if: illegal combination of | and &\n");
                    exit(12); }
                  i=*x=='!';
                  if (i) {
                    x++;
                    if (!*x) {
                      fprintf(stderr,"makemake: !if: '!' misplaced (or followed by space)\n");
                      exit(12); } }
                  i^=isdef(x);
                  if (!ifop) active=i;
                  if (ifop=='&' && active) active=i;
                  if (ifop=='|' && !active) active=i;
                  ifop=key; }
              activestat[ifnest]=active;
              active=getactive(); }

            else if (!strcmp(token,"else")) {
              if (ifnest==0 || !ifstat[ifnest]) {
                fprintf(stderr,"makemake: misplaced !else\n");
                exit(12); }
              ifstat[ifnest]=0;
              activestat[ifnest]^=1;
              active=getactive(); }

            else if (!strcmp(token,"endif")) {
              if (ifnest<=0) {
                fprintf(stderr,"makemake: !endif without matching !if\n");
                exit(12); }
              ifnest--;
              active=getactive(); }

            else if (!strcmp(token,"def") || !strcmp(token,"define")) {
              if (active) if ((x=scan(&c,"$$"))) if (*x) {
                    define(x);
                    fprintf(out,"# %s defined\n",x); } }

            else if (!strcmp(token,"undef")) {
              if (active) if ((x=scan(&c,"$$"))) if (*x) {
                    undef(x);
                    fprintf(out,"# %s undefined\n",x); } }

            else if (!strcmp(token,"error")) {
              if (active) { fprintf(stderr,"makemake: !error %s\n",c);
                exit(12); } }

            break;

          case '=': /* keyword followed by = */

            if (!active) break;

            if (!strcmp(token,"home")) {
              if (home[0]) fprintf(stderr,"makemake: \aWARNING !home redefined\n");
              if (ndirs>1) {
                fprintf(stderr,"makemake: !home after !dir illegal\n");
                exit(4); }
              if ( (x=qscan(&c)) ) {
                strcpy(home,x); strcat(home,"/");
                escape(eschome,home);
                if (strcmp(home,eschome))
                  fprintf(stderr,"makemake: WARNING: special characters in !home will be \\escaped\n");
              } }

            if (!strcmp(token,"ignore"))
              while ( (x=qscan(&c)) ) {
                struct ignore_s *ig=myalloc(sizeof(struct ignore_s));

                if (strlen(x)>=FNLEN) {
                  fprintf(stderr,"makemake: too long name in !ignore : %s\n",x);
                  exit(4); }
                strcpy(ig->name,x);
                ig->next=ignore0; ignore0=ig; }

            if (!strcmp(token,"treat"))
              while ( (x=qscan(&c)) ) {
                struct ignore_s *ig;
                /* masked, not freed (memory leak) */

                for (ig=ignore0; ig; ig=ig->next)
                  if (!strcmp(ig->name,x)) {
                    strcpy(ig->name,"-/void/");
                    break; } }

            if (!strcmp(token,"objext")) {
            if ( (x=qscan(&c)) ) {
              scopy(objext+(x[0]!='.'),x); } }

            if (!strcmp(token,"dir")) {
              if (ndirs>1) fprintf(stderr,"makemake: more than one !dir statement (appended, %d dirs)\n",ndirs);
              if (home[0]==0) {
                fprintf(stderr,"makemake: !dir without !home\n");
                exit(4); }
              while ( (x=qscan(&c)) ) {
                scopy(dirs[ndirs],x);
                fprintf(out,"D%d = %s%s\n",ndirs,eschome,x);
                ndirs++;
                if (ndirs>NDIRS) {
                  fprintf(stderr,"makemake: too many dirs\n");
                  exit(4); } } }
            break;

          case ':': /* keyword followed by : */

            if (!active) break;

            /* exe from obj command */
            containscplus=0;
            x=token;
            nobj++;

            /* list of exe's (used by make all,clean) */
            project=myalloc(sizeof(*project));
            scopy(project->name,x);
            project->namelist=NULL;
            project->next=project0;
            project0=project;

            fprintf(out,"OBJ%d =",nobj); len=6+(nobj>9);

            xlibopt="$(LIBOPT)";
            xopt=NULL; /* it means the default = $(OPT) */

            while ( (sem=strchr(c,':')) ) {
              while (*c<=' ') c++;
              getxopt(&xopt,"OPT",c);
              getxopt(&xlibopt,"LIBOPT",c);
              *sem=0; c=sem+1;
              while (sem[-1] && sem[-1]<=' ') { sem--; *sem=0; } }

            while ( (token=qscan(&c)) ) {
              char *plus=strchr(token,'+');

              if (plus) *plus=0;
              putobj(token);
              dependency(token,1+(plus!=NULL),strcmp(x,token)?NULL:xopt); }

            /* print exe:obj dependencies and linking rule */
            fprintf(out,"\n%s :%s $(OBJ%d)\n",
                    x,makefile,nobj);
            fprintf(out,"\t$(%s) -o %s $(OBJ%d) %s\n\n",
                    /* linking requires CPLUS if at least one module is C++ */
                    containscplus?"CPLUS":"CC",
                    x,nobj,
                    xlibopt);
            break;
        } /* switch (key) */
    }
  } /* for (;;) */

  if (ifnest) {
    fprintf(stderr,"makemake: missing !endif\n");
    exit(12); }

  /* print dependencies and compiling rules for obj-files */
  for (f=f0; f; f=f->next) if (f->type==0) {
      fprintf(out,"%s%s :%s",f->name,objext,makefile);
      len=2+strlen(f->name)+strlen(objext);
      for (ff=f0; ff; ff=ff->next) ff->todo=1;
      alldep(f);
      fprintf(out,"\n\t$(%s)",f->dep->f->cplus?"CPLUS":"CC");
      fprintf(out," -o %s%s",f->name,objext);
      /*        for (i=1; i<ndirs; i++) fprintf(out," -I$(D%d)",i);*/
      for (i=ndirs-1; i; i--) fprintf(out," -I$(D%d)",i);
      fprintf(out," %s",f->xopt);
      fprintf(out,
              " -c $(D%d)/%s.%s\n\n",
              f->dep->f->type,f->name,f->dep->f->cplus&1?"cc":"c"); }

  fprintf(out,"all :");
  len=6;
  for (project=project0; project; project=project->next) {
    l=strlen(project->name)+1;
    if ( (len+=l) > LINELEN ) {
      fprintf(out," \\\n "); len=l+1; }
    fprintf(out," %s",project->name); }
  fprintf(out,"\n\n"); /* NOTHING after TAB ? */

  fprintf(out,"clean :\n\trm");
  len=10;
  for (project=project0; project; project=project->next) {
    l=strlen(project->name)+1;
    if ( (len+=l) > LINELEN ) {
      fprintf(out," \\\n "); len=l+1; }
    fprintf(out," %s",project->name); }
  for (f=f0; f; f=f->next) if (f->type==0) {
      l=strlen(f->name)+strlen(objext)+1;
      if ( (len+=l) > LINELEN ) {
        fprintf(out," \\\n "); len=l+1; }
      fprintf(out," %s%s",f->name,objext); }

  fprintf(out,"\n");
  fclose(out);

  if (filelist) {
    qsort(filelist,nfilelist,sizeof(filelist[0]),cmpfilelist);
    out=fopen("filelist","wt");
    for (i=0; i<nfilelist; i++) fprintf(out,"%s\n",filelist[i]);
    fclose(out); }

  if (warnings)
    fprintf(stderr,"makemake: WARNING: %d multiple file%s (the last one active)\n",warnings,"s"+(warnings==1));

  return notfound;
}
