/* cc -O2 -o grouptab grouptab.c -lm
*/
#include "../gen/include.h"

#define NL 32000
char line[NL];

int ncol;
char **data;

int nextline(void)
{
  char *c;
  int i;

  do
    if (fgets(line,NL,stdin)==NULL) return 0;
    while ( strchr("!#\n",line[0]) || !(c=strtok(line," \t\n")) );

  loop (i,0,ncol)
    if (c) {
      data[i]=c;
      c=strtok(NULL," \t\n"); }
    else
      data[i]="n.a.";

  return 1;
}

int main(int narg, char **arg)
{
  int group,i,j;

  getsbufsize=NL;
  
  if (narg<3) {
    fprintf(stderr,"\
Grouping table data (making longer lines).  Call by:\n\
  grouptab BLOCK NCOL < INPUT > OUTPUT\n\
BLOCK\thow many lines will compose one output line\n\
NCOL\t# of columns in a line considered; 0=use whole lines\n\
Example:\n\
  grouptab 2 3\n\
with input\n\
  a b c x\n\
  A B C X\n\
will produce\n\
  a b c A B C\n\
See also:\n\
  ungrouptab blocktab splittab septab smooth mergetab tabproc sum\n");
   exit(0); }

  group=atoi(arg[1]);
  if (group<1) {
    fprintf(stderr,"invalid group size\n");
    exit(1); }

  ncol=atoi(arg[2]);

  if (ncol) {

    allocarray(data,ncol);
  
    j=0;
    while (nextline()) {
      j++;
      loop (i,0,ncol)
        printf("%s%c",data[i],i==ncol-1 && j==group ? '\n' : ' ');
      if (j==group) j=0; }
  
    if (j) {
      fprintf(stderr,"WARNING: incomplete group (%d line%s) was ignored\n",
  	    j,"s"+(j==1));
      printf("\n"); } }
  else {
    j=0;
    while (gets(line)) {
      j++;
      printf("%s%c",line,j%group?' ':'\n'); }
    printf("\n"); }

  return 0;
}
