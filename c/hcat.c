/* cc -O2 -o hcat hcat.c
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
  char line[4100];
  FILE **f;
  char *nl;
  int iarg,nopen;
  char *SEP=getenv("HCAT");

  if (narg<2) {
    fprintf(stderr,"\
Concatenate text files horizontally. Call by:\n\
  [HCAT=] hcat FILE [FILE..] > OUTFILE\n\
Where\n\
  HCAT = separator string (default=" ")\n\
See also:\n\
  cat mergetab\n\
");
    exit(0); }

  if (SEP[0]==0) SEP=" ";
  
  allocarray(f,narg);
  
  loop (iarg,1,narg) 
    if (! (f[iarg]=fopen(arg[iarg],"rt")) ) 
      fprintf(stderr,"%s: cannot open\n",arg[iarg]);

  do {
    nopen=0;
    loop (iarg,1,narg)
      if (f[iarg] && fgets(line,4099,f[iarg])) {
        nopen++;
        if (iarg!=narg-1) {
          nl=strchr(line,'\n');
          if (nl) *nl=0;
          strcat(line,SEP); }
        printf("%s",line); }
    if (!strchr(line,'\n')) printf("\n");  
  } while (nopen);

  return 0;
}
