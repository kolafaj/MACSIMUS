/* cc -O2 -o sparsetab sparsetab.c
*/
#include "../gen/include.h"

#define NL 32000
char line[NL],stored[NL];

int main(int narg, char **arg) /*************************************** main */
{
  char *fmt="%.15g";
  int COL,init=1,printed=0;
  double DX,x0,xstored;

  if (narg<3) {
   fprintf(stderr,"\
Make a table sparse by omitting data. Call by:\n\
  sparsetab [-]COLX DX [FORMAT] < INPUT > OUTPUT\n\
Arguments:\n\
  COLX = control column of (sorted increasingly) X\n\
  - in front of COL: empty, #-, and !-lines are omitted (default=copied)\n\
  DX = stride, lines will be omitted so that X[i+1]-X[i] < DX, but the\n\
       differences X[i+1]-X[i] are maximal\n\
  FORMAT = optional format, default=%s\n\
See also:\n\
  gblock smooth mergetab tabproc sum redtab datablock blocktab\n",fmt);
   exit(0); }

  COL=atoi(arg[1]);
  DX=atof(arg[2]);
  if (narg>3) fmt=arg[3];

  getsbufsize=NL;

  while (gets(line)) {
    if (strlen(line)==0 || strchr("!#",line[0])) {
      printed=1;
      if (COL>0) puts(line); }
    else {
      int i;
      char *tok=line;
      double x;

      if (tok) while (strchr(" \t\r\n",*tok)) tok++;

      x=atof(tok);

      loop (i,1,abs(COL)) {
        x=atof(tok);
        while (!strchr(" \t\r\n",*tok)) {
          if (!tok) Error("no such column");
          tok++; } }

      if (init) {
        puts(line);
        printed=1;
        x0=x;
        init=0; }
      else {
      again:
        if (x>x0+DX) {
          if (stored[0]) {
            puts(stored);
            printed=0;
            stored[0]=0;
            x0=xstored;
            goto again; }
          else {
            puts(line);
            printed=1;
            x0=x; } }
        else {
          strcpy(stored,line);
          xstored=x;
          printed=0; }
      } } }

  if (!printed) puts(line[0]?line:stored);
  
  return 0;
}
