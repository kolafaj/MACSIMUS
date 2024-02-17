/* cc -O2 -o comment comment.c
 */
#include "../gen/include.h"

#define N 32000

int main(int narg,char **arg)
{
  static char line[N];
  int iarg;

  if (narg<2) {
    fprintf(stderr,"\
Omit or extract comment lines. Call by:\n\
  comment OPTIONS < INPUT > OUTPUT\n\
OPTIONS (can repeat, options are processed sequentially for each line):\n\
  e  extract all lines (the default)\n\
  o  omit all lines (useful before +$)\n\
  s  remove trailing spaces\n\
  t  remove trailing tabs\n\
  x  remove trailing whites (chars<=' ')\n\
  -  remove empty (only LF) lines\n\
  +  extract empty (only LF) lines\n\
  -$ remove lines beginning by string $\n\
  +$ extract lines beginning by string $\n\
Examples:\n\
  comment -\\! -\\# - < file.in > file.out  # omit empty, !-, and #-lines\n\
  comment x - < file.in > file.out        # omit white lines\n\
  comment o +\\# < file.gnu                # show all #-comments\n\
See also:\n\
  dellines numeric\n\
");
    exit(0); }

  getsbufsize=N;

  //  extract0=0;
  //  loop (iarg,1,narg) if (arg[iarg][0]=='-') extract0=1;
  
  while (gets(line)) {
    int extract=1;
    char *c;

    loop (iarg,1,narg)
      switch (arg[iarg][0]) {
        case 'e':
          extract=1;
          break;
        case 'o':
          extract=0;
          break;
        case 'x':
          for (c=line+strlen(line)-1; c>=line && *c<=' '; c-- ) *c=0;
          break;
        case 't':
          for (c=line+strlen(line)-1; c>=line && *c==8; c-- ) *c=0;
          break;
        case 's':
          for (c=line+strlen(line)-1; c>=line && *c==' '; c-- ) *c=0;
          break;
        case '-':
          if (arg[iarg][1]) {
            if (!memcmp(line,arg[iarg]+1,strlen(arg[iarg]+1))) extract=0; }
          else {
            if (strlen(line)==0) extract=0; }
          break;
        case '+':
          if (arg[iarg][1]) {
            if (!memcmp(line,arg[iarg]+1,strlen(arg[iarg]+1))) extract=1; }
          else {
            if (strlen(line)==0) extract=1; }
          break;
        default:
          Error("comment: bad option"); }

    if (extract) puts(line); }

  return 0;
}
