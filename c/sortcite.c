/* cc -o sortcite sortcite.c

To reorded \bibitem{} list in LaTeX files
Original version renum.c
Updated (1995) as sortcite.c
Updated (1998) to handle \Cite and \Hide
Updated (2008) to ignore items containing # (assumed macro definition)
Updated (2009) to handle \citen, \Cite and \Hide removed
*/
#include "../gen/include.h"

struct cite_s {
  struct cite_s *next;
  char ref[64]; } *citehead,*cite,*cc;

struct bib_s {
  struct bib_s *next;
  int printed;
  char ref[64];
  char text[1024]; } *bibhead,*bib;

int citepos;

char *nextcite(char *l)
{
  char *cite=strstr(l,"\\cite{");
  char *citen=strstr(l,"\\citen{");

  citepos=6;
  if (!cite) cite=citen,citepos=7; 
  else if (citen) {
    if (citen<cite) cite=citen,citepos=7; }

  if (cite && cite[citepos]=='#') cite=NULL; /* \citeX{#... assumed to be macro definition */

  return cite;
}

int main(int narg,char **arg)
{
  FILE *fi,*fo;
  char fn[80];
  char line[1024];
  char *l,*e,*t,*c;
  static char keyword[24];
  int nl;
  int found;

  if (narg<2) {
    printf("\
Sort \\bibitem{XXX} list according to the order of arguments in the\n\
\\cite{} and \\citen{} statements in the .tex file\n\
Usage:\n\
  %s FILE [-]\n\
FILE.tex is input, FILE.ren output\n\
optional - means that not referenced bibitems will be removed from output\n",
           arg[0]);
  exit(0); }

  sprintf(fn,"%s.tex",arg[1]);
  fi=fopen(fn,"rt");
  if (!fi) {
    strcat(fn," not found");
    Error(fn); exit(1); }

  sprintf(fn,"%s.ren",arg[1]);
  fo=fopen(fn,"wt");

  nl=0;
  while (fgets(line,1024,fi)) {
    nl++;
    l=line;
    while ( (l=nextcite(l)) ) {
      /* all keywords of the same length !! */
      copy(keyword,l,citepos); keyword[citepos]=0;
      l+=citepos;
      e=strchr(l,'}');
      if (!e) {
        printf("%d:%s",nl,line);
        strcat(keyword," not ended by }");
        Error(keyword);
        break; }
      t=l;
      while (t<e) {
        cc=malloc(sizeof(struct cite_s));
        if (!cc) Error("no memory");
        c=cc->ref;
        while (*t<=' ') t++;
        while (*t!=',' && *t!='}') *c++=*t++;
        *c=0; t++;
        printf("%s%s}\n",keyword,cc->ref);
        cc->next=NULL;
        if (citehead) { cite->next=cc; cite=cc; }
        else cite=citehead=cc; }
      l=e+1; } }

  rewind(fi); nl=0;
  while (fgets(line,1024,fi)) {
    nl++;
    if (strstr(line,"\\begin{thebibliography}")) {
      while (fgets(line,1024,fi)) {
        nl++;
        if (line[0]!='%' && strstr(line,"\\bibitem{")) goto OK; } } }

  Error("no \\begin{thebibliography}");

  OK: nl--;
  l=strstr(line,"\\bibitem{");
  do {
    if (!l) goto END;
    *l='@';
    e=strchr(l+9,'}');
    if (!e) {
      printf("%d:%s",nl,line);
      Error("\\bibitem{ not ended by }"); }
    *e++=0;
    bib=malloc(sizeof(struct bib_s));
    if (!bib) Error("no memory");
    strcpy(bib->ref,l+9);
    bib->printed=0;
    printf("\\bibitem{%s}\n",bib->ref);
    bib->next=bibhead;
    bibhead=bib;
    strcpy(bib->text,e);
    while (fgets(line,1024,fi)) {
      while (line[0]=='%') { fgets(line,1024,fi); nl++; }
      nl++;
      l=strstr(line,"\\bibitem{");
      if (l) {
        *l=0;
        strcat(bib->text,line);
        break; }
      l=strstr(line,"\\end{thebibliography}");
      if (l) {
        *l=0;
        strcat(bib->text,line);
        goto END; }
      strcat(bib->text,line); }
    } while (l);

  END:

  /*.....for (bib=bibhead; bib; bib=bib->next) printf("<%s>\n%s\n",bib->ref,bib->text);*/

  rewind(fi); nl=0;
  while (fgets(line,1024,fi)) {
    nl++; fputs(line,fo);
    if (strstr(line,"\\begin{thebibliography}")) goto OK2; }

  Error("no \\begin{thebibliography}");

  OK2:

  for (cite=citehead; cite; cite=cite->next) {
    again:
    found=0;
    for (bib=bibhead; bib; bib=bib->next) {
      if (!strcmp(cite->ref,bib->ref)) {
        if (!bib->printed) {
          fprintf(fo,"\\bibitem{%s}%s",bib->ref,bib->text);
          bib->printed++; }
        found++; } }
    if (!found) {
      bib=malloc(sizeof(struct bib_s));
      if (!bib) Error("no memory");
      strcpy(bib->ref,cite->ref);
      bib->printed=0;
      printf("missing \\bibitem{%s}\n",bib->ref);
      bib->next=bibhead;
      bibhead=bib;
      strcpy(bib->text,"{\\it unknown}\n");
      goto again; } 
    }

  for (bib=bibhead; bib; bib=bib->next)
    if (!bib->printed) {
      if (narg<3 || strcmp(arg[2],"-"))
        fprintf(fo,"%% not referenced:\n\\bibitem{%s}%s",bib->ref,bib->text);
      printf("\\bibitem{%s} not referenced\n",bib->ref); }

  while (fgets(line,1024,fi)) {
    if (line[0]=='%') fputs(line,fo);
    if ( (l=strstr(line,"\\end{thebibliography}")) ) break; }
  fputs(l,fo);

  while (fgets(line,1024,fi)) fputs(line,fo);

  fclose(fo);
  fclose(fi);

  return 0;
}
