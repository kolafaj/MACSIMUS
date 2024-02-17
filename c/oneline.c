/* cc -O2 -o oneline oneline.c
  Make 1 line of input data (filter)
  see also: onelineetc
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
  char line[1024];
  int next=0;

  if (narg>1) {
    fprintf(stderr,"\
Make 1 line of input data (filter).  Call by:\n\
  oneline < INPUT > OUTPUT\n\
See also:\n\
  onelineetc ungrouptab\n\
");
    exit(0); }

  while (gets(line)) {
    if (next++) printf(" %s",line);
    else printf("%s",line); }

  _n

  return 0;
}
