/* cc -Wall -O0 -o filttab filttab.c -lm
 should be compiled with -O0 to prevent long double in registers
 */
#include "../gen/include.h"

char line[16000],line2[16000];

int main(int narg, char **arg)
{
  double x;
  int iarg,icol,cmp,iarg0;
  unsigned minlen=0; /* NB: must be unsigned because compared to strlen() */
  char *tok,*comment="\r";
  int blank=0;
  struct arg_s {
    int col, not;
    double from,to,by,by0; } *a;

  if (narg<2) {
    fprintf(stderr,"\
Select range of values from a table. Call by:\n\
  filttab [+OPT] [-]COLUMN:FROM[:TO[:BY[:BY0]]] [...] < DATA > SELECTED\n\
  filttab [+OPT] [-]COLUMN:VALUE:+ERROR[:BY[:BY0]]] [...] < DATA > SEL.\n\
Lines with numerical value in given COLUMN in interval [FROM,TO] are selected\n\
- in front of COLUMN selects data outside the range\n\
OPT=[length][ ][comment]:\n\
  length: lines shorter than length are always copied (default=0=none)\n\
  space: blank lines printed instead of omitted lines\n\
  comment: any char from comment as 1st char causes the line to be copied\n\
More conditions are OR-ed (to AND them, use a pipe)\n\
Invalid lines and comment lines are omitted\n\
Default FROM=0, default TO=FROM (exact match)\n\
BY selects data which are an integer multiple of BY, or (datum-BY0)/BY is int\n\
Example (select DATA with COLUMN1 outside [-3,5] and COLUMN2 = 0):\n\
  cat DATA | filttab -1:-3:5 | filttab 2:0\n\
Example (copy lines with 0 in the 1st column + empty lines + !-comments):\n\
  filttab +1\\! 1:0 < DATA\n\
See also:\n\
  key septab splittab mergetab tabproc shifttab field transtab\n");
    exit(0); }

  getsbufsize=16000;
  
  allocarray(a,narg);
  if (arg[1][0]=='+') {
    comment=arg[1]+1;
    minlen=atoi(comment);
    while (isdigit(*comment)) comment++;
    if (*comment==' ') blank++,comment++;
    iarg0=2; }
  else
    iarg0=1;
  loop (iarg,iarg0,narg) {
    a[iarg].col=abs(atoi(arg[iarg]));
    a[iarg].not=atoi(arg[iarg])<0;
    a[iarg].from=a[iarg].to=a[iarg].by=a[iarg].by0=0;
    if ( (tok=strchr(arg[iarg],':')) ) {
      tok++;
      a[iarg].from=a[iarg].to=atof(tok);
      if ( (tok=strchr(tok,':')) ) {
        tok++;
        if (tok[0]=='+') {
          x=atof(tok+1);
          a[iarg].to+=x;
          a[iarg].from-=x; }
        else
          a[iarg].to=atof(tok);
        if ( (tok=strchr(tok,':')) ) {
          tok++;
          a[iarg].by=atof(tok);
          if ( (tok=strchr(tok,':')) ) {
            tok++;
            a[iarg].by0=atof(tok); } } } } }

  while (gets(line)) {
    if (strlen(line)<minlen || (line[0] && strchr(comment,line[0])))
      puts(line);
    else {
      strcpy(line2,line);
      tok=strtok(line,"\t \r,");
      icol=1;
      while (tok) {
        loop (iarg,iarg0,narg) if (a[iarg].col==icol) {
          char *end;
          x=strtod(tok,&end);
          if (end<=tok) continue;
          cmp = x>=a[iarg].from && x<=a[iarg].to;
          if (a[iarg].not) cmp=!cmp;
          if (cmp && a[iarg].by!=0) {
            x=(x-a[iarg].by0)/a[iarg].by;
            cmp = fabs(floor(x+0.5)-x)<1e-3; }
          if (cmp) {
            puts(line2);
            goto nextline; } }
        tok=strtok(NULL,"\t \r,");
        icol++; }
      if (blank) puts("");
    nextline: ; } }

  return 0;
}
