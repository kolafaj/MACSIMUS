/* cc -O2 -o start start.c

 Start application associated with file extension.
 (c) Jiri Kolafa 2001, GNU General Public Licence

 2011: exts made case-insensitive (.jpg=.JPG)

 Motivation:
   This is a command-prompt analogue of the Windows (or Mac, Gnome, KDE,
   OS/2, ...) mechanism of starting applications according to the type
   of the associated file.  In (some versions of) Windows there is
   command `start' of the same function.  Such a style of work may be
   considered strange by Windows folk.  But my brain is not able to find
   a file of interest among more than ten icons in a graphical folder
   (or Norton/Midnight Commander list) and I consider typing the file
   name much faster (especially with file completition and wildcards).

 Simple examples (unfortunately, xv has died...):
   start pig.jpg             --> xv pig.jpg
   start archive.zip         --> unzip -v archive.zip | less

 Extended examples:
   start -expand 2.5 pig.jpg --> xv -expand 2.5 pig.jpg
   start pig.jpg -expand 2.5 --> xv pig.jpg -expand 2.5 (number .5 not ext)
   start mol cfg.plb -YS     --> show mol cfg -YS

 Details:
 1/ Scans the arguments and determines the extension of the last file
    argument:
    * File argument is an argument not starting with - nor + and with
      an extension.
    * The extension is a suffix after the last . in the file name provided
      that the name before this . is not empty.
    * The extensions are case-insensitive.
    * Decimal digits only without a letter after . are not recognized as
      file extension.
    * Spaces allowed in file names.
 2/ Looks up the associated application.
 3/ Starts this application with arguments, using system().
 4/ If the executable is called `ss' or `starts', and there is more than one
    associated executable, it asks for selection, otherwise the first one is
    started.

File registration:
  Use file ~/.startdata.  Associated applications are separated by TABs.
  These defaults may be overriden by optional .startdata in the working
  directory.

Bugs:
  Some special characters may be treated incorrectly.
  Args are enclosed in "" which mey not be what you wants.
*/

#include "../gen/include.h"

#define N 10240 /* max command line length */

struct reg_s {
  struct reg_s *next;
  char *ext; /* file extension, without . */
  int n; /* # of commands */
  int pass; /* 0=./ 1=~/ */
  char **cmd; /* commands to execute
      "%s.%s" is replaced by NAME.EXT
      "%s" is replaced by NAME (without extension)
      missing "%s" means that "%s.%s" & is appended */
} *head;

int nstr(char *s)
/* counts %s in string s */
{
  int n=0;

  while ( (s=strstr(s,"%s")) ) s++,n++;

  return n;
}

void strtac(char *where,char *what)
/* where:=what||where
   similas as strcat, which is where:=where||what */
{
  memmove(where+strlen(what),where,strlen(where)+1);
  memcpy(where,what,strlen(what));
}


int EQUAL(char *a,char *b)
/* returns 1 if strings identical, case-insensitive */
{
  while (*a && *b) {
    if (toupper(*a)!=toupper(*b)) return 0;
    a++; b++; }

  return *a==*b;
}

int main(int narg,char **arg)
{
  struct reg_s *r,*rr,*r0=NULL;
  int iarg;
  int ii,n,ask=1,pass;
  char *dot,*ext,*fn,*x;
  char fmt[N],sys[N];
  FILE *in;

  getsbufsize=N;
  
  if (narg<2) {
    fprintf(stderr,"\
Starts application associated with file extension. Call by:\n\
  start[s] [application -OPTIONS or other ARGS] FILE.EXT [more -OPTIONS]\n\
  start[s] .[EXT] : list association(s)\n\
If there are several associated commands:\n\
  start takes the first (default) command (also shortcut s)\n\
  starts offers all commands and asks which one to start (also shortcut ss)\n\
Extensions and commands are listed in ~/.startdata (read first) and may be\n\
  modified by ./.startdata (working directory local changes)\n\
See also:\n\
  mc.ext (MACSIMUS bindings for Midnight Commander)\n");
    exit(0); }

  if (!strcmp("s",arg[0])) ask=0;
  if (!strcmp("ss",arg[0])) ask=1;
  if (strstr(arg[0],"start")) ask=0;
  if (strstr(arg[0],"starts")) ask=1;

  /* reading ~/.startdata */
  n=0;
  loop (pass,0,2) {
    if (pass) {
      x=getenv("HOME");
      if (!x) continue; }
    else x=".";
    sprintf(sys,"%s/.startdata",x);
    in=fopen(sys,"rt");
    if (!in) continue;
    n++;
    while (fgets(sys,N,in)) if (sys[0]!='#') if (strchr(sys,'\t')) {
      if ( (x=strchr(sys,'\n')) ) *x=0;
      alloconezero(r);
      if (head) r0->next=r,r0=r;
      else head=r0=r;
      for (x=sys; *x; x++) r->n+=*x=='\t';
      if (!r->n) {
        fprintf(stderr,"%s\n",sys);
        Error(".startdata: no TAB"); }
      r->pass=pass;
      allocarray(r->cmd,r->n);
      x=strtok(sys,"\t");
      looplist (rr,head) if (rr->ext) if (!strcmp(x,rr->ext)) if (rr->pass) {
        fprintf(stderr,"%s\n",sys);
        Error("extension repeats in ~/.startdata"); }
      r->ext=strdup(x);
      ii=0;
      while ( (x=strtok(NULL,"\t")) ) {
        if (ii>=r->n) Error("internal");
        r->cmd[ii]=strdup(x);
        if (nstr(r->cmd[ii])>2) {
          fprintf(stderr,"%s\n",sys);
          Error("too many %%s in the cmd"); }
        ii++; } }
    fclose(in); }

  if (!n) Error("none of ~/.startdata ./.startdata found");

  /* check args */
  for (iarg=1,n=2; iarg<narg; iarg++) n+=strlen(arg[iarg]+1);
  if (n>N) {
    fprintf(stderr,"too long/many arguments");
    exit(0); }

  /* finding file argument */
  dot=NULL;
  for (iarg=narg-1; iarg>0; iarg--) if (!strchr("-+",arg[iarg][0])) {
    char *a;
    int flag;

    /* finding the last `.' ignoring `.' just after start or `/' */
    for (a=arg[iarg],flag=1; *a; a++)
      switch (*a) {
        case '.': if (!flag) dot=a;
        default: flag=0; break;
        case '/': dot=NULL,flag++; }

    if (!dot) continue; /* no extension - no `file arg' */

    /* if there are only dec. digits after ., this is NOT extension */
    flag=1;
    for (a=dot+1; *a; a++) if (!isdigit(*a)) flag=0;
    if (flag) continue;

    fn=strdup(arg[iarg]);
    dot+=fn-arg[iarg];
    ext=dot+1;
    *dot=0;

    /* finding the command */
    looplist (r,head) if (EQUAL(ext,r->ext)) {
      char *cmde;
      int i;

      if (r->n>1 && ask) {
        loop (ii,0,r->n) printf("%d: %s\n",ii+1,r->cmd[ii]);
        gets(fmt);
        ii=atoi(fmt)-1;
        if (ii<0 || ii>=r->n) {
          ii=0;
          fprintf(stderr,"out of range, default used\n"); } }
      else ii=0;

      n=nstr(r->cmd[ii]); /* # of %s in the cmd */

      strcpy(fmt,r->cmd[ii]);

      if (n==0)
        /* no format ==> %s.%s & assumed (& postponed) */
        strcat(fmt," \"%s.%s\"");

      for (cmde=fmt; *cmde; cmde++) if (*cmde==' ') break;
      if (*cmde!=' ') Error("missing space in command");
      cmde++;

      /* fmt += "options" which should precede FILE[.EXT] arg */
      loop (i,1,iarg) {
        strtac(cmde,"\" ");
        strtac(cmde,arg[i]);
        strtac(cmde,"\""); }

      /* fmt += "options" which should follow FILE[.EXT] arg */
      loop (i,iarg+1,narg) {
        strcat(fmt," \"");
        strcat(fmt,arg[i]);
        strcat(fmt,"\""); }

      /* replace %s[.%s] by FILE[.EXT] to create the command */
      if (n==0) strcat(fmt," &");
      sprintf(sys,fmt,fn,ext); /* ext not used if n=1 */

      /* run the command */
      fprintf(stderr,"%s\n",sys);
      return system(sys); }

    Error("unknown ext"); }

  /* list association(s) of command(s) to extension(s) */
  if (arg[1][0]=='.') {
    looplist (r,head)
      if (!arg[1][1] || !strcmp(arg[1]+1,r->ext)) {
        loop (ii,0,r->n)
          printf("%c%-4s\t%s\n",ii?' ':'.',ii?"":r->ext,r->cmd[ii]); } }
  else
    Error("no file");

  return -1;
}
