/* cc -O2 -o extract extract.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char Line[1024];

char toksep[32]=" \t\n";
char *strtofind="";

void edittoksep(char *arg) /************************************* edittoksep */
{
 if (isalnum(arg[1])) return;
 if (arg[0]=='+') {
   if (strlen(toksep)+strlen(arg+1)>30) {
     fprintf(stderr,"too many field separators (option +<SEP>...)\n");
     exit(-5); }
   strcat(toksep,arg+1); }
 else {
   char *del,*t,*d,ch;

   for (del=arg+1; *del; del++)
     for (t=toksep; *t; )
       if (*t==*del)
         for (ch=1,d=t; ch; d++) ch=*d=d[1];
       else
         t++;
   }
}

int token(char *l,int i) /******************************************** token */
{
  char *tok;

  if (i<=0) return 0;

  strcpy(Line,l);
  tok=strtok(Line,toksep);
  while (tok) {
    if (--i<=0) return tok-Line;
    tok=strtok(NULL,toksep); }

  return strlen(l);
}

int main(int narg,char **arg) /**************************************** main */
{
  int i,j,i1,i2;
  FILE *in;
  char line[1024],*rg;
  char *fmt=NULL,*ch;
  double x=0,xx;
  char op=0;
  int ifrom=-1,count;

  if (narg<2) {
    fprintf(stderr,"\
Extracts parts from a text file, (c) J.Kolafa 1997. Call by:\n\
  %s OPTION [ OPTION ... ] { FILE | - }\n\
where OPTION is\n\
  /STRING      find next line with STRING, position at BOL\n\
  =STRING      find next line with STRING and extract one item after STRING\n\
  /  =         repeat previous /STRING or =STRING\n\
  +NUMBER      go to line NUMBER relatively from the current position\n\
  FIELD:FIELD  extract fields in given range (default separators=SPACE TAB NL)\n\
  :FIELD       extract one field (FIELD=decimal number)\n\
  :FIELD%%fmt   extract one field, convert to double, print by double format\n\
  :            extract whole line\n\
  (NUMBER .. ) loop: do NUMBER times commands enclosed in (..)\n\
  _TEXT        print TEXT\n\
  _            print newline\n\
  @NUMBER      print argument[NUMBER]\n\
  @            print filename (last argument)\n\
  +s +m +M     start sum, min, Max (only 1st field from range treated)\n\
  -s -m -M     print sum, min, Max; also -s%%fmt etc. with format (df.=%%g)\n\
  +STRING      add field separator(s) (no alpha nor numeric)\n\
  -STRING      remove field separator(s) (no alpha nor numeric)\n\
  COL-COL      extract columns in given range\n\
FILE  file to process (- instead of FILE = standard input)\n\
Example (to find \"Tin  \" and extract block of 24 lines):\n\
  extract \"/Tin  \" \\(24 : +1 \\) test.prt\n\
See also:\n\
  eqfield field\n\
",arg[0]);
  exit(0); }

  narg--;
  if (arg[narg][0]=='-') in=stdin;
  else in=fopen(arg[narg],"rt");
  if (!in) {
    fprintf(stderr,"%s: no such file\n",arg[narg]);
    exit(-4); }

  line[0]=0;

  for (i=1; i<narg; i++) {

    if (line[strlen(line)-1]=='\n') line[strlen(line)-1]=0;

    switch (arg[i][0]) {

      case '(':
        ifrom=i;
        count=atoi(arg[i]+1);
        break;

      case ')':
	if (ifrom<0) {
	  fprintf(stderr,") without opening (COUNT\n");
	  exit(1); }
        if (--count) i=ifrom;
        break;

      case '@': {
        int j=narg;
        if (arg[i][1]) j=atoi(arg[i]+1);
        if (j>=0 && j<=narg) printf("%s ",arg[j]); }
        break;

      case '/':
        if (arg[i][1]) strtofind=arg[i]+1;
        if (!*strtofind) {
          fprintf(stderr,"nothing to find (empty option / or =)\n");
          exit(-6); }
        while (fgets(line,1024,in)) if (strstr(line,strtofind)) goto found;
          goto eof;
        found: break;

      case '=':
        if (arg[i][1]) strtofind=arg[i]+1;
        if (!*strtofind) {
          fprintf(stderr,"nothing to find (empty option = or /)\n");
          exit(-6); }
        while (fgets(line,1024,in)) if (ch=strstr(line,strtofind)) {
            i1=ch-line+strlen(strtofind);
            while (*++ch) if (strchr(toksep,*ch)) break;
            i2=ch-line;
            goto columns; }
        goto eof;

      case '-':
        //        if (arg[i][1]=='n') sep=' '; else
        if (toupper(arg[i][1])=='S' || toupper(arg[i][1])=='M') {
          printf("%g ",x);
          op=0; }
        edittoksep(arg[i]);
        break;

      case '+':
        edittoksep(arg[i]);
        //        if (arg[i][1]=='n') sep='\n'; else
        if (arg[i][1]=='s' || arg[i][1]=='S') op='s',x=0;
        else if (arg[i][1]=='m') op='m',x=9e99;
        else if (arg[i][1]=='M') op='M',x=-9e99;
        else
          for (j=0; j<atoi(arg[i]+1); j++)
            if (!fgets(line,1024,in)) goto eof;
        break;

      case '_':
        if (arg[i][1])
          printf("%s ",arg[i]+1);
        else
          printf("\n");
        break;

      case ':':
        if (arg[i][1]==0) {
          /* print whole line */
          i1=0; i2=strlen(line);
          goto columns; }
        else {
          i1=i2=atoi(arg[i]+1);
          fmt=strchr(arg[i]+2,'%');
          goto fields; }

      default:
        if ( (rg=strchr(arg[i],':')) ) {
          i1=atoi(arg[i]);
          i2=atoi(rg+1);
        fields:
          i1=token(line,i1);
          i2=token(line,i2+1)-1;
          while (i2>0 && strchr(toksep,line[i2])) i2--;
          i2++;
          goto columns; }

        else if ( (rg=strchr(arg[i],'-')) ) {
          i1=atoi(arg[i])-1;
          i2=atoi(rg+1);
        columns:
          if (i2<=i1 || i1<0 || i2<0) break;
          strcpy(Line,line);
          Line[i2]=0;
          xx=atof(Line+i1);
          switch (op) {
            case 's': case'S': x+=xx; break;
            case 'm': if (xx<x) x=xx; break;
            case 'M': if (xx>x) x=xx; break;
            case 0:
              if (fmt)
                printf(fmt,atof(Line+i1));
              else
                printf("%s",Line+i1);
              break;
            default: printf("?"); }
          printf(" ");
          fmt=NULL; }

        else {
          fprintf(stderr,"%s: bad option\n",arg[i]);
          exit(-2); }

      } /* switch */
    }

 eof: fclose(in);

  return 0;
}
