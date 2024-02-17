/* cc -O2 -o eddata eddata.c -lm

   07/2022 - tmpnam no longer used, FN renamed to FN~ and new FN made

 */
#include "../gen/include.h"

#define LEN 1024

int main(int narg,char **arg)
{
  FILE *in,*out;
  char line[LEN],*fn,*fnbak;
  int i,iarg,ifile=0;
  
  if (narg<2) {
    fprintf(stderr,"\
Edit data file (.def, .get) in the MACSIMUS `get data' format. Call by:\n\
  eddata ID=VALUE [ID=VALUE ...] [-]FILE [[-]FILE ...]\n\
Description:\n\
  Any ID=OLDVALUE field in FILE is replaced by ID=VALUE.\n\
  If (the first) FILE contains `=', it must be prefixed by `-'.\n\
  VALUE may be generally any string.\n\
  No spaces are allowed around `='.\n\
  If VALUE is missing (ID=), whole ID=OLDVALUE is removed from FILE.\n\
  ID may be qualified (STRUCT.ID) and subscripted (ARRAY[1]).\n\
  The old FILEs are backed up as FILE~.\n\
Example:\n\
  eddata init=2 T=300 \"tau.rho= tau.P=10\" \"L[0]=1\" sim-[abc].def\n\
See also:\n\
  eqfield repl replace lemon\n");
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-' || !strchr(arg[iarg],'=')) { 
      ifile=iarg;
      break; }
       
  if (!ifile) Error("no file arg");
  if (ifile==1) Error("1st arg must be ID=[VALUE]");

  loop (iarg,ifile,narg) {
    fn=arg[iarg]+(arg[iarg][0]=='-');
    fnbak=string("%s~",fn);
    rename(fnbak,string("%s~",fnbak));

    if (rename(fn,fnbak)) {
      fprintf(stderr,"eddata: cannot rename %s to %s, file skipped\n",fn,fnbak);
      continue; }
    else {
      in=fopen(fnbak,"rt");
      if (!in) {
        fprintf(stderr,"eddata: cannot open %s, exit (rename back to %s?)\n",fnbak,fn);
        exit(1); }
      out=fopen(fn,"wt");
      if (!out) {
        fprintf(stderr,"eddata: cannot write to %s, exit\n",fn);
        exit(1); }
       
      while (fgets(line,LEN,in)) {
        loop (i,1,ifile) {
          char needle[1024],*n,*c,*repl,*exclm;
        
          for (c=arg[i],n=needle;*c;c++,n++) {
            *n=*c;
            if (*c=='=') break; }
          *++n=0;
          repl=strstr(line,needle);
          exclm=strchr(line,'!');
          if (repl && repl!=line)
            if (isalnum(repl[-1]) || strchr("[].",repl[-1]) || (exclm && repl>exclm)) repl=NULL;
          if (repl) {
            char *end=repl,*eq=strchr(arg[i],'=')+1;
  
            while (*end>' ') {
              putchar(*end);
              end++; }
            if (*eq>' ') eq=arg[i];
            printf(" --> %s\n",eq);
            memmove(repl+strlen(eq),end,strlen(end)+1);
            memcpy(repl,eq,strlen(eq)); } }
        fputs(line,out); }
  
      fclose(in);
      fclose(out); }
  }
  
  return 0;
}
