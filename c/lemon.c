/* cc -O2 -o lemon lemon.c
  Similar to fgrep (no \ before special characters) plus:
  - can list more lines after searched string
  - string replacement in multiple files
  - case, words, ...
*/

#include "../gen/include.h"

#define LLEN 1048578

char *str,*repl;
int match=0;
int pluslines=0;
int plusrepl;
enum { NONE,ALNUM,ALNUMPUNCT,NOWHITE } words=NONE;
int binchar=0;
int yes=0;
int Case=1;
int filenames=1;
char *tmpfn;

char *casestr(const char *where,const char *what) /***************** casestr */
{
  const char *a,*b,*c;

  for (a=where; *a; a++) {
    for (b=what,c=a; *b; b++,c++) if (toupper(*b)!=toupper(*c)) goto cont;
    return (char*)a;
   cont:; }
  return NULL;
}

int word(char x) /***************************************************** word */
{
  switch (words) {
    case ALNUM: return isalnum(x); break;
    case ALNUMPUNCT: return isalnum(x) || strchr("+-.",x); break;
    case NOWHITE: return x>' '; break;
    default: fprintf(stderr,"internal: bad switch %d\n",words); exit(1);
  }
}

char *mystrstr(const char *where,const char *what) /*************** mystrstr */
/* as strstr but possibly as whole words */
{
  char *p=(char*)where,*e;

 nextstrstr:
  if (Case) p=strstr(p,what);
  else p=casestr(p,what);

  if (!p) return NULL;
  if (words!=NONE) {
    if (where<p) if (word(p[-1])) { p++; goto nextstrstr; }
    e=p+strlen(what);
    if (*e) if (word(*e)) { p++; goto nextstrstr; } }

  return p;
}

void search(char *fn) /********************************************** search */
{
  char l[LLEN],*m,*ll;
  FILE *f=stdin,*ff=stdout;
  int changed=0;
  int iplus=0;
  long unsigned bit7=0;

  if (fn[0]) {
    f=fopen(fn,"rt");
    
    if (repl) {
      tmpfn=string("%s~l",fn);
      ff=fopen(tmpfn,"wt"); }
    if (!ff) {
      fprintf(stderr,"lemon: cannot write to %s\n",tmpfn);
      exit (1); } }

  if (!f) {
    fprintf(stderr,"lemon: cannot open %s\n",fn[0]?fn:"(stdin)");
    return; }

  while (fgets(l,LLEN,f)) {

    if (strlen(l)>LLEN-2) {
#if 1
      if (fn[0]) printf("%12s: ? line longer than %d\n",fn,LLEN-2);
      else printf("? line longer than %d\n",LLEN-2);
#endif
      return; }

    for (m=l; *m; m++) bit7 += (unsigned char)*m >= (unsigned char)128;

    if ((m=mystrstr(l,str))) {
      match++;
      if (fn[0]) { 
        if (plusrepl) {
          if (filenames==1) printf("%12s: %s",fn,l); 
          else if (filenames==2) {
            char *x=strdup(l),*nl=strchr(x,'\n');
            if (nl) *nl=0;
            printf("%s [%s]\n",x,fn);
            free(x); }
          else printf("%s",l); }
      }
      else if (!repl) printf("%s",l);
      if (!repl) iplus=pluslines+!!pluslines;
      ll=l;

      if (repl) {
	char answ[16];
	int torepl=!fn[0] || yes;

	if (!torepl && plusrepl) {
	  fprintf(stderr,"replace (y/N) ? ");
	  if (fgets(answ,16,stdin)) torepl = answ[0]=='y' || answ[0]=='Y';
          else torepl=0; }

	if (plusrepl==0) torepl=0;
	else plusrepl--;

	if (torepl) {
	  changed=1;
	  do {
	    char zero=*m;
	    *m=0;
	    fputs(ll,ff);
	    fputs(repl,ff);
	    ll=m+strlen(str);
	    if (strlen(str)==0) { *m=zero; break; }
          } while ((m=mystrstr(ll,str))); }

	fputs(ll,ff); } }

    else if (repl) fputs(l,ff);

    if (iplus) {
      if (iplus<=pluslines) {
	if (fn[0] && filenames==1) printf("%12s: %s","",l); 
	else printf("%s",l); }
      /* if (!fn[0]) */
      if (iplus==1) printf("\n");
      iplus--; }
  }

  if (binchar && bit7) {
    if (fn[0]) printf("%12s! %ld chars >= 128\n",fn,bit7);
    else printf("! %ld chars>=128\n",bit7); }

  fclose(f);
  if (repl) {
    fclose(ff);
    if (fn[0]) {
      if (changed) {
        if (rename(fn,string("%s~",fn))) {
          fprintf(stderr,"lemon: cannot rename %s to %s\n",fn,fn);
          exit(1); }
        if (rename(tmpfn,fn)) {
          fprintf(stderr,"lemon: cannot rename %s to %s\n",tmpfn,fn);
          exit(1); } }
      else
        remove(tmpfn); } }
}

int main(int narg,char **arg) /**************************************** main */
{
  int i,from=2;

  if (narg<2) {
    fprintf(stderr,"\
String (no regexp!) search/replace in files. Call by:\n\
  lemon STRING [-REPLACE] [OPTIONs] [FILE [FILE ..]]\n\
STRING:   string to search for\n\
-REPLACE: all occurrences of STRING within line will be replaced by REPLACE\n\
FILE:     file of lines max 1048578 long\n\
          with -REPLACE, asks for confirmation and makes backup FILE~\n\
          no FILE: filter (does not ask for confirmation with -REPLACE)\n\
OPTIONs:  +# = print additional # lines after each line containing STRING\n\
          with -REPLACE: replace STRING in at most # lines\n\
          +a = whole words {alpha digits}; e.g., 123 matches -123.4\n\
          +w = whole words+numbers {alpha digits +-.}\n\
          +f = fields (all but whites)\n\
          +b = print info on chars>=128\n\
          +i = ignore case\n\
          +nC = character C in REPLACE becomes newline on output\n\
          +y = does not ask for replace confirmation\n\
          +h = do not print filenames\n\
          +hh = print filenames after string (default=before:)\n\
See also:\n\
  grep fgrep head tail unix2dos dos2unix mac2unix\n\
  orange myfgrep itail llemon lgrep strcount count lemonn charhist liat\n\
  repl replalt binrepl charrepl text line maxline oneline revlines padcutrev\n\
  rndlines dellines excl eqfield extract eddata\n\
  kostnice kostnice2 utf8txt\n");
    exit(0); }

  str=arg[1];
  while (from<narg) {
    if (arg[from][0]=='-') {
      repl=arg[from]+1; from++; continue; }
    if (arg[from][0]=='+') {
      if (arg[from][1]=='a') words=ALNUM;
      else if (arg[from][1]=='w') words=ALNUMPUNCT;
      else if (arg[from][1]=='f') words=NOWHITE;
      else if (arg[from][1]=='b') binchar=1;
      else if (arg[from][1]=='i') Case=0;
      else if (arg[from][1]=='y') yes=1;
      else if (arg[from][1]=='h') {
        filenames=0;
        if (arg[from][2]=='h') filenames=2; }
      else if (arg[from][1]=='n') {
	if (repl) {
	  char *c;
	  for (c=repl; *c; c++) if (*c==arg[from][2]) *c='\n'; } }
      else pluslines=atoi(arg[from]+1);
      from++; continue; }
    break; }

  if (narg==from) {
    plusrepl=pluslines?pluslines:0x7fffffff;
    search(""); }
  else for (i=from; i<narg; i++) {

    plusrepl=pluslines?pluslines:0x7fffffff;

    search(arg[i]); }

  return match;
}
