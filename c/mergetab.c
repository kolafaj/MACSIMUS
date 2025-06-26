/* cc -O2 -Wall -o mergetab mergetab.c -lm
 */
#include "../gen/include.h"

#define N 512 /* max # of files, max # of columns in one file */
#define LINE 8192 /* max line length */
typedef short int shortint;

char *fn;
FILE *f[N];
int ncol[N]; /* ncol[file index], # of columns requested from file */
shortint col[N][N]; /* col[file index][upto ncol[N]], column */
int omit[N]; /* omit line counter */
int stride[N]; /* data stride */
shortint eq[N]; /* eq[file index], column of f[file index]
                      to match column eq0[file index] of file 1 */
shortint eq0[N];
char line[N][LINE];
char *toks[N][N];  /* tokens (pointing into line[][]) found */
int na[N]; /* synchronization flag:
              if set, does not read next line of f[file index] */

char *mystrdup(char *s) /****************************************** mystrdup */
{
  char *n=strdup(s);

  if (!n) Error("mergetab: no heap");

  return n;
}

int stdinused;

void scanarg(char *arg,int i) /************************************* scanarg */
{
  char *tok,*tokcont=NULL,*t;
  static char *lastfn="-";
  int j;

  fn=mystrdup(arg);
  if (fn[0]==':') {
    tokcont=fn+1;
    tok=lastfn; }
  else
    tok=strtok(fn,":");
  if (strcmp(tok,"-"))
    f[i]=fopen(tok,"rt");
  else {
    if (stdinused++) Error("mergetab: cannot repeat - (stdin)");
    f[i]=stdin; }
  if (!f[i]) Error(fn);
  lastfn=mystrdup(tok);
  ncol[i]=0;
  stride[i]=1;
  loop (j,0,N) col[i][j]=-1;

  while ((tok=strtok(tokcont,":")) && ncol[i]<N) {
    j=atoi(tok);
    tokcont=NULL;
    if (tok[0]=='-')
      omit[i]=-j;
    else if (tok[0]=='/')
      stride[i]=atoi(tok+1);
    else {
      if ( (t=strchr(tok,'=')) ) eq[i]=j-1,eq0[i]=atoi(t+1)-1;
      else if ( (t=strchr(tok,'-')) ) {
        int ii,to=atoi(t+1);
        loop (ii,j-1,to) col[i][ncol[i]++]=ii; }
      else col[i][ncol[i]++]=j-1; } }

  if (!stride[i]) Error("mergetab: zero stride");
  if (stride[i]<0) {
    stride[i]=-stride[i];
    omit[i]+=stride[i]-1; }
}

int anyopened(void) /********************************************* anyopened */
{
  int i;

  loop (i,0,N) if (f[i]) return 1;

  return 0;
}

char *sep=" \t\n";

int mygetline(int iarg) /******************************************* getline */
{
  int i=0;
  char *tok;

 again:
  memset(line[iarg],0,LINE);
  loop (i,0,N) toks[iarg][i]=NULL;
  if (!f[iarg])
    return 0;
  else {
    if (!fgets(line[iarg],LINE,f[iarg])) {
      if (f[iarg]!=stdin) fclose(f[iarg]);
      f[iarg]=NULL;
      goto again; }
    if (strchr("!#",line[iarg][0])) goto again; }

  for (i=0,tok=strtok(line[iarg],sep); tok && i<N; tok=strtok(NULL,sep))
    toks[iarg][i++]=tok;

  if (!i) goto again;

  return i;
}

int main(int narg,char **arg) /**************************************** main */
{
  int i,j;
  static int usestride=0;
  char *FMT="%-7s ";
  double err=0;

  if (narg<2) {
    fprintf(stderr,"\
Manipulate white-separated ASCII tables. Call by:\n\
  mergetab [FILE1]:[-OMIT:][/STRIDE:]COL[-COL]:COL[-COL]... \n\
    [FILE2]:[-OMIT:][/[-]STRIDE:]{COL[-COL]|COL=COL1}:{COL[-COL]|COL=COL1}...\n\
    ...\n\
FILEi\t file name, or - for stdin (max once allowed)\n\
\t missing FILEi repeats previous, missing FILE1=stdin\n\
OMIT\t # of noncomment lines omitted from the top of file\n\
STRIDE\t every STRIDE-th noncomment line\n\
-STRIDE\t as above, starting with (STRIDE-1)-th line\n\
\t OMIT and -STRIDE are combined, e.g., \":-2:/-10\" = \":-11:/10\"\n\
COL\t column to print\n\
COL-COL\t range of columns to print, e.g., \":2-4\" = \"2:3:4\"\n\
COL=COL1 synchronizes with column COL1 of FILE1\n\
\t data in columns to synchronize should be in increasing order\n\
Environment:\n\
  MERGEFMT = format (string); default=\"%s\"\n\
  MERGESEP = separators; default=\" \\t\\n\"\n\
  MERGEERR = permitted error of COL=COL, negative: rel.error (data>0 only)\n\
Example:\n\
  mergetab st2.g:1:2 tip4p.g:1=1:2 > waters.g\n\
will print columns 1 and 2 of st2.g and column 2 of tip4p.g shifted so that\n\
it matches column 1. If tip4p.g has finer grid, extra data are omitted,\n\
if st2.g has finer grid, missing data in tip4p.g are n.a.\n\
See also:\n\
  datablock hcat tabproc transtab field shifttab maketab selcol mergeg\n",FMT);
  exit(0); }

  if (getenv("MERGEFMT")) FMT=getenv("MERGEFMT");
  if (!strchr(FMT,'%')) Error("string format (%s) expected");

  if (getenv("MERGESEP")) sep=getenv("MERGESEP");
  if (getenv("MERGEERR")) err=atof(getenv("MERGEERR"));

  for (i=1; i<N; i++) eq[i]=-1;
  for (i=1; i<N && i<narg; i++) scanarg(arg[i],i);
  if (eq[1]>=0) Error("mergetab: COL=COL for 1st file");

  for (i=1; i<N && i<narg; i++) while(omit[i]--) mygetline(i);

  if (!anyopened()) Error("mergetab: no valid data");

  for (;;) {

    for (i=1; i<N && i<narg; i++) if (f[i] && !na[i]) {
      if (usestride) loop (j,0,stride[i]) mygetline(i);
      else mygetline(i); }

    usestride=1;

    if (!anyopened()) break;

    for (i=1; i<N && i<narg; i++) {
      int valid=0;

      na[i]=0;
    again:
      if (eq[i]>=0) {
        char *x0=toks[1][eq0[i]];
        char *x=toks[i][eq[i]];

        if (x0 && x) {
          if (err==0) {
            if (atof(x)<atof(x0)) {
              loop (j,1,stride[i]) mygetline(i);
              if (mygetline(i)) goto again; /* read next line NOW */ }
            else if (atof(x)>atof(x0)) na[i]++; /* skip reading next line */
            else valid++; }
          else if (err>0) {
            if (atof(x)<atof(x0)-err) {
              loop (j,1,stride[i]) mygetline(i);
              if (mygetline(i)) goto again; /* read next line NOW */ }
            else if (atof(x)>atof(x0)+err) na[i]++; /* skip reading next line */
            else valid++; }
          else  {
            if (atof(x)<atof(x0)*exp(-err)) {
              loop (j,1,stride[i]) mygetline(i);
              if (mygetline(i)) goto again; /* read next line NOW */ }
            else if (atof(x)>atof(x0)/exp(err)) na[i]++; /* skip reading next line */
            else valid++; }
        }
      }
      else valid++;

      loop (j,0,ncol[i])
        if (valid && toks[i][col[i][j]]) printf(FMT,toks[i][col[i][j]]);
        else printf(FMT,"n.a."); }
    printf("\n"); }

  return 0;
}
