/* cc -O2 -o ungrouptab ungrouptab.c -lm
*/
#include "../gen/include.h"

#define NL 1048576

char line[NL];

int main(int narg, char **arg)
{
  int group,ngroups=1048576,omit,i,n,nlin=0,nlout=0,silent;

  if (narg<3) {
    fprintf(stderr,"\
Un-grouping table data (fewer items per line).  Call by:\n\
  ungrouptab [-]OMIT GROUP [NGROUPS] < INPUT > OUTPUT\n\
Parameters:\n\
  OMIT: columns omitted from line beginning\n\
  -OMIT: as above, quiet (do not print # of lines to stderr)\n\
  GROUP: columns in one group (will be printed to one line)\n\
  NGROUPS: number of groups (default=until EOL)\n\
columns are white-separated\n\
!#-comments are omitted\n\
Example:\n\
  ungrouptab 1 2 2\n\
Will make from INPUT:\n\
  a b c d e f\n\
  A B C D E F\n\
The following file:\n\
  b c\n\
  d e\n\
  B C\n\
  D E\n\
See also:\n\
  oneline grouptab blocktab smooth mergetab tabproc sum sumetc\n");
   exit(0); }

  omit=abs(atoi(arg[1]));
  silent=arg[1][0]=='-';
  group=atoi(arg[2]);
  if (narg>3) ngroups=atoi(arg[3]);

  getsbufsize=1048576;
  
  while (gets(line)) if (!strchr("#!",line[0])) {
    char *tok=strtok(line," \t\n\r");

    nlin++;

    loop (i,0,omit) tok=strtok(NULL," \t\n\r");

    loop (n,0,ngroups) {
      loop (i,0,group) {
	if (!tok) break;
	if (i) printf(" ");
	printf("%s",tok);
	tok=strtok(NULL," \t\n\r"); }
      nlout++; 
      printf("\n"); 
      if (!tok) break; } }

  if (!silent) fprintf(stderr,"ungrouptab: %d non-comment input lines, %d output lines\n",nlin,nlout);

  return 0;
}
