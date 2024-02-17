/* cc -O2 -o line line.c
 */
#include "../gen/include.h"
#include <values.h>

int main(int narg,char **arg)
{
  int n=1,n0=1,N,iarg,last=1,minlen=-1,Cstyle=0;
  int buf=1048576;
  char *line;
  char *lf="\n";
  char *comments=NULL;

  if (narg<2) {
    fprintf(stderr,"\
Select line(s) from a text-file:\n\
  line [ARGS|OPTIONS] < INFILE > OUTFILE\n\
ARGS (any order, output ordered as input w. optional repeats,\n\
empty/impossible silently ignored):\n\
  N       extract line number N\n\
  N1-N2   extract range of lines (N2 included unless -C)\n\
  N1-     extract all lines from N1 until EOF\n\
OPTIONS:\n\
  -cCHARS ignore comments (lines starting with any CHAR) [NULL]\n\
  -C      C-style: line N2 is not included in range N1-N2\n\
          consider also -0 to number lines from 0\n\
  -i      ignore empty (only LF) lines\n\
  -iLEN   ignore lines shorter than LEN (LF does not count)\n\
  -lSIZE  max line length [default=1048576]; undefined (but safe) if exceeded\n\
  -n      do not print LF at EOL\n\
  -nCHAR  replace LF by CHAR (1 character)\n\
  -z      ignore last incomplete (not ending by LF) line\n\
  -NUMBER number lines from NUMBER [default=1]\n\
Examples:\n\
  line -i -c 1-10 < a.asc > b.asc # first 10 noncomment not empty lines\n\
  line -z 1- < a.asc > b.asc      # all except data after the last LF\n\
See also:\n\
  tail head  liat revlines maxline oneline tab2line\n\
");
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'c': comments=arg[iarg]+2;
                  if (!comments[0]) comments="!#";
                  break;
        case 'C': Cstyle=1; break;
        case 'l': buf=atoi(arg[iarg]+2); break;
        case 'i': minlen=atoi(arg[iarg]+2); break;
        case 'n': lf=arg[iarg]+2; break;
        case 'z': last=0; break;
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9': n0=n=atoi(arg[iarg]+1); break;
        default: Error("line: bad option"); }

  if (buf<1 || buf>1048576*1024) Error("line: buffer size out of bounds");
  allocarray(line,buf);
  if (!buf) Error("line: no heap");

  while (fgets(line,buf,stdin)) if (comments==NULL || !strchr(comments,line[0])) {
    loop (iarg,1,narg) if (arg[iarg][0]!='-') {
      int from=atoi(arg[iarg]),to;
      char *dash=strchr(arg[iarg]+1,'-');

      //      if (from<n0) Error("line: (starting) line number precedes the first line");

      if (!dash)
        to=from;
      else if (dash[1]) {
        to=-atoi(dash);
        if (Cstyle) to--; }
      else
        to=MAXINT;
      //      fprintf(stderr,"to=%d n0=%d\n",to,n0);
      //      if (to<n0) Error("line: final line number precedes the first line");

      if (n>=from && n<=to)
        if (last || strchr(line,'\n'))
	  if (minlen<0 || strlen(line)>minlen+1) {
	    if (lf) {
	      char *n=strchr(line,'\n');
	      if (n) strcpy(n,lf); }
	    fputs(line,stdout); } }
    n++; }

  return 0;
}
