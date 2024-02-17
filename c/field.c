/* cc -O2 -o field field.c

update 3/2021, default COMMENT chaged to "#!"
 */
#include "../gen/include.h"

#define LINE 260000
#define N 10002

int main(int narg,char **arg) /**************************************** main */
{
  char line[LINE];
  char *tok[N],*minus,nr[16];
  int iarg;
  int nf,f,t,i,nl=0;
  int numbers=0;
  unsigned nprt=0;
  char *SEP=" \t\r\n";
  char *COMMENT="#!";

  if (narg<2) {
    fprintf(stderr,"\
Select white-separated field(s) except comments. Call by:\n\
  field ARGUMENTs < INFILE > OUTFILE\n\
ARGUMENT:\n\
  FIELD = dec. number of the field (fields are numbered from 1)\n\
  FIELD-FIELD = range of fields\n\
  FIELD- = the same as FIELD-BIG\n\
  + = next field is advanced by 1\n\
  +RELPOS = next field is advanced by RELPOS\n\
  TEXT = any string not beginning with digit/+ is unchanged copied to output\n\
  -TEXT = TEXT is any text (may begin with a digit)\n\
The ARGUMENTs may appear in any order (except those for the next argument)\n\
and may repeat or overlap. THE OUTPUT IS IN THE ORIGINAL ORDER OF COLUMNS!!!\n\
Returns 0 on success, 1 if no field has been printed (all lines too short)\n\
Environment:\n\
  SEP = list of separators, default = \" \\t\\r\\n\"\n\
  empty SEP = extract numbers (e.g., a12T1e4 -> 12 1e4)\n\
  COMMENT = string, omit lines beginning by characters given, default=\"%s\"\n\
Example (plot all *.dat files, one by one, only TAB as separator):\n\
  export SEP=$\'\\t\'\n\
  ls -1 *.dat | field plot 1 | sh\n\
See also:\n\
  mergetab filttab tabproc deleol numbers\n",COMMENT);
    exit(0); }

  if (getenv("SEP")) SEP=getenv("SEP");
  if (!strlen(SEP)) numbers++,SEP=" \t\r\n";
  if (getenv("COMMENT")) COMMENT=getenv("COMMENT");

  tok[0]=nr;

  while (fgets(line,LINE,stdin)) if (!strchr(COMMENT,line[0])) {

    if (numbers) {
      /* extract only numbers */
      char *c;
      char *wasnr=NULL,*isnr;

      for (c=line; *c; c++) {
        isnr=strchr("0123456789.+-",*c);
        if (!isnr) {
          if (wasnr && c[-1]=='.') c[-1]=' ';
          if (!strchr("eE",*c) || !wasnr) *c=' '; } 
        wasnr=isnr; } }

    sprintf(nr,"%d",nl++);

    nf=1;
    tok[nf]=strtok(line,SEP);

    while (tok[nf++]) {
      if (nf>=N) Error("too many columns");
      tok[nf]=strtok(NULL,SEP); }
    nf--;

    loop (iarg,1,narg) {
      if (iarg>1) printf(" ");
      if (isdigit(arg[iarg][0])) {
        f=t=atoi(arg[iarg]);
        if ( (minus=strchr(arg[iarg],'-')) ) {
          if (minus[1]) t=atoi(minus+1);
          else t=nf; }
        if (t>=nf) t=nf-1;
        if (f<1) f=0;
        loopto (i,f,t) {
          if (i>f) printf(" ");
          nprt++;
          printf("%s",tok[i]); } }
      else if (arg[iarg][0]=='+') {
        if (isdigit(arg[iarg][1])) f+=atoi(arg[iarg]+1); 
        else f++;
        if (t>=nf) t=nf-1;
        t=f; /* just one... */
        loopto (i,f,t) {
          if (i>f) printf(" ");
          nprt++;
          printf("%s",tok[i]); } }
      else
        printf("%s",arg[iarg]+(arg[iarg][0]=='-')); }
      printf("\n"); }

  return !nprt;
}
