/* cc -O2 -o eqfield eqfield.c
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
  char line[32000];
  char *key;

  if (narg<2) {
    fprintf(stderr,"Select fields after given keys, tab-separate. Call by:\n\
  eqfield +KEY < INFILE > OUTFILE\n\
  eqfield KEY [KEY..] < INFILE > OUTFILE\n\
Examples - let us filter line:\n\
  a=1 b=2 a=3 d=4\n\
The results are:\n\
  eqfield a= b=  => \"1 3 2\" (all a= first)\n\
  eqfield =      => \"1 2 3 4\"\n\
  eqfield +b=    => \"2 a=3 d=4\"\n\
See also:\n\
  line field data2tab eddata numbers extract\n");
   exit(0); }

  key=arg[1];
  if (key[0]=='+') {
    key++;
    while (fgets(line,32000,stdin)) {
      char *l=line;

      if ( (l=strstr(l,key)) ) {
	l+=strlen(key);
	fputs(l,stdout); } } }
  else
    while (fgets(line,32000,stdin)) {
      char *l,*e;
      int yes=0;
      int iarg;

      for (iarg=1; iarg<narg; iarg++) {
        l=line;
        key=arg[iarg];
        while ( (l=strstr(l,key)) ) {
          l+=strlen(key);
          yes++;
          while (*l==' ') l++;
          for (e=l; isalnum(*e)||strchr(".+-()/*^",*e); e++) putchar(*e);
          putchar('\t'); } }
      if (yes) putchar('\n'); }

 return 0;
}
