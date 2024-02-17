/* cc -Wall -O2 -o cpacol cpacol.c
 */
#include "../gen/include.h"

char line[65536];

int main(int narg,char **arg)
{
  int n,i;
  char *c,*s;

  FILE *in;

  if (narg<3) {
    fprintf(stderr,"\
Convert column/symbol of cpa files.\n\
Note that with -a2, time column is added to .cpa and columns are +1\n\
Print symbol of column number COL:\n\
  cpcol.sh CPA-FILE COLNUM\n\
Print column number of mnemonic NAME:\n\
  cpcol.sh CPA-FILE NAME\n\
Example:\n\
  plot sim.cpa:n:`cpacol sim.cpa Jz`\n\
See also:\n\
  showcp\n\
");
    exit(0); }

  in=fopen(arg[1],"rt");
  if (!in) {
    puts("ERROR no file");
    return 4; }

  n=atoi(arg[2]);

  while (fgets(line,65536,in))
    if (!memcmp(line, "#__",3)) {
      s=line;
      if (*s=='#') s++;
      if (*s=='_') s++;
      if (*s=='_') s++;
      if (*s=='_') s++;

      s=strtok(s," \n");

      for (;;) {
        if (!s) {
          puts("ERROR not found");
          return 2; }
        i=atoi(s);
        c=strchr(s,':');
        if (!c) {
          puts("ERROR format");
          return 3; }
        c++;
        if (n) {
          if (i==n) {
            puts(c);
            return 0; } }
        else {
          if (!strcmp(c,arg[2])) {
            printf("%d\n",i);
            return 0; } }

        s=strtok(NULL," \n"); }
    }
  puts("ERROR EOF");
  return 2;
}
