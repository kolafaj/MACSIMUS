/* cc -O2 -o cutprt cutprt.c
  7/2018 new version - DOS support removed, code simplified, longer CP chunk kept
*/
#include "../gen/include.h"

#define LINE 1024

int main(int narg,char **arg)
{
  FILE *in,*out;
  char line[LINE],last[LINE],fn[256];
  char *tmp="shortened~.prt";
  int iarg,l,skipped,wasshortened;
  enum stat_e { BODY,STARTED,HDR,RED } stat;
  int NPLUS=3; /* default # of lines kept */

  if (narg<2) {
    fprintf(stderr,"\
Shorten MACSIMUS prt-files by removing most of the convergence profiles\n\
Call by:\n\
  [CUTPRT=N] cutprt FILE.prt [FILE.prt ..]\n\
N =  number of lines shown after convergence profile header or ERROR/WARNING\n\
     [default=2]\n\
The shortened file is FILE.prt, time stamp is copied\n\
The original file becomes FILE.prt~\n\
File beginning by '... shortened by cutprt ...' is left unchanged\n\
");
    exit(0); }

  if (getenv("CUTPRT")) NPLUS=atoi(getenv("CUTPRT"));

  loop (iarg,1,narg) {

    fprintf(stderr,"cutprt: opening %s .. ",arg[iarg]);

    in=fopen(arg[iarg],"rt");
    if (!in) {
      fprintf(stderr,"not found\n");
      goto cont; }

    if (!fgets(line,LINE,in)) {
      fprintf(stderr,"empty\n");
      fclose(in);
      goto cont; }

    if (!strcmp("... shortened by cutprt ...\n",line)) {
      fprintf(stderr,"already shortened\n");
      fclose(in);
      goto cont; }

    out=fopen(tmp,"wt");
    if (!out) Error("cannot open tmp file");
    fputs("... shortened by cutprt ...\n",out);
    strcpy(last,line);

    stat=BODY;
    l=skipped=wasshortened=0;

    while (fgets(line,LINE,in)) {

      if (stat==BODY && strstr(line," started at ") && strstr(line,"cycle"))
        stat=STARTED;

      if (!memcmp("===============================================================================",line,75)) {
        l=NPLUS+2;
        if (stat==STARTED) stat=HDR;
        else stat=BODY; }

      if (!memcmp("-------------------------------------------------------------------------------",line,75)) {
        if (stat==HDR) stat=RED,l=NPLUS+2;
        else stat=BODY; }

      if (stat==RED) {
        if (!strcmp("...\n",line)) {
          fprintf(stderr,"already shortened\n");
          fclose(in);
          fclose(out);
          goto cont; }
        if (atof(line)==0) l=NPLUS+2; }

      l--;

      if (stat<RED || l>=0) {
        if (skipped) fputs("...\n",out);
        skipped=0;
        fputs(last,out); }
      else 
        skipped=1,wasshortened=1;

      strcpy(last,line); 

    } /* while(fgets) */

    fputs(last,out);

    if (!wasshortened) {
      fprintf(stderr,"nothing to shorten\n");
      fclose(in);
      fclose(out);
      goto cont; }

    if (stat!=BODY) {
      fprintf(stderr,"unexpected EOF or not MACSIMUS prt-file\n");
      fclose(in);
      fclose(out);
      goto cont; }

    fprintf(stderr,"%s written\n",tmp);

    fclose(in);
    fclose(out);

    strcpy(fn,arg[iarg]);
    strcat(fn,"~");

    if (!remove(fn)) fprintf(stderr,"cutprt: %s removed\n",fn);
    fprintf(stderr,"cutprt: renaming %s to %s\n",arg[iarg],fn);
    if (rename(arg[iarg],fn)) Error("could not rename");
    fprintf(stderr,"cutprt: renaming %s to %s\n",tmp,arg[iarg]);
    if (rename(tmp,arg[iarg])) Error("could not rename temporary file\n");
    fprintf(stderr,"cutprt: copying time stamp\n");
    if (system(string("touch -r \"%s\" \"%s\"",fn,arg[iarg])))
      fprintf(stderr,"cutprt: %s: cannot copy time stamp\n",arg[iarg]);

  cont:; }

  return 0;
}
