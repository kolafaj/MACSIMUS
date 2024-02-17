/* cc -O2 -o repl repl.c
 */
#include "../gen/include.h"

#define N 65536

int main(int narg,char **arg)
{
  static char line[N];
  int nrepl=N,from=1,totrepl=0;
  enum opt_e { ANY,BOL,EOL } opt=ANY;

  if (narg<3) {
    fprintf(stderr,"Replace strings in a stream. Call by:\n\
  repl [OPTION] FROM TO [FROM TO ...] < SOURCE > RESULT\n\
if the number of args is odd, the 1st arg is OPTION:\n\
  -NUMBER = max. number of replacements of one type (FROM->TO) per line\n\
            (default = all; use -1 in case of possible recursion)\n\
  -b = replace only if FROM is at the beginning of line\n\
  -e = replace only if FROM is at the end of line\n\
total number of replacements is returned\n\
Examples:\n\
  echo ABAA | repl A a  B bb      # -> abbaa\n\
  echo ABAA | repl -2  A a        # -> aBaA\n\
  echo ABAA | repl -e  A \"\"       # -> ABA\n\
See also:\n\
  ifrepl binrepl replace lemon texrepl\n");
    exit(0); }

  if (!(narg&1)) {
    /* odd number of parameters (excl. arg[0]): 1st is OPTION */
    from=2;
    if (arg[1][0]!='-') Error("repl: wrong OPTION: must be one of -NUMBER -b -e");
    if (arg[1][1]=='b') opt=BOL,nrepl=1;
    else if (arg[1][1]=='e') opt=EOL,nrepl=1;  
    else nrepl=atoi(arg[1]+1); }
  if (nrepl<=0) Error("repl: wrong OPTION or args not paired");

  while (fgets(line,N,stdin)) {
    int i;
    
    for (i=from; i<narg; i+=2) {
      int irepl,replaced=0;

      loop (irepl,0,nrepl) {
        char *c=strstr(line,arg[i]);

        switch (opt) {
          case BOL:
            if (c!=line) c=NULL; break;
          case EOL: 
            c=strlast(line);
            c+=(*c!='\n');
            c-=strlen(arg[i]);
            if (c>=line) c=strstr(c,arg[i]);
            else c=NULL; 
            break;
        }
      
        if (c) {
          if (strlen(line)+strlen(arg[i+1])-strlen(arg[i])>=N) 
            Error("repl: buffer overflow or recursion (consider repl -1)");
          memmove(c+strlen(arg[i+1]),
                  c+strlen(arg[i]),
                  strlen(c)-strlen(arg[i])+1);
          memcpy(c,arg[i+1],strlen(arg[i+1]));      
          replaced++;
          totrepl++; } 
        if (!replaced) break; } }
    fputs(line,stdout); }
    
  return totrepl;
}

