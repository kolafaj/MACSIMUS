/* cc -O2 -o binrepl binrepl.c
 */
#include "../gen/include.h"

#define N 0x10000000
unsigned char fff[N+16];
unsigned char *file=fff+16;

int yeskey;
char *extarg;

int yes(int pos) /****************************************************** yes */
{
  int i,from,to,l=100/*??*/;
  char line[1024];

  if (yeskey) return 1;

  for (from=pos-5; from>pos-77 && from>=0; from--) if (file[from]=='\n') break;
  for (to=pos+l; file[to]; to++) if (file[to]=='\n') break;
  
  loopto (i,from,to) putchar(file[i]);
  for (;;) {
    printf("\nREPLACE in %s? (n)o (y)es (a)ll_to_EOF (A)ll_files (q)uit\n",extarg);
    gets(line);
    switch (toupper(line[0])) {
      case 'Q': exit(0);
      case 'N': return 0;
      case 'Y': return 1;
      case '!':
      case 'A': { yeskey=line[0]; return 1; } } }

}

int excl;

char *str(char *coded,int *len) /*************************************** str */
{
  char *s=strdup(coded),*a,*b;
  
  *len=0;
  for (a=s,b=coded; *b; b++,a++) {
    if (excl && *b==excl) {
      b++;
      if (*b==excl) *a=*b;
      else if (*b=='n') *a='\n';
      else if (*b=='r') *a='\r';
      else if (*b=='t') *a='\t';
      else if (isdigit(*b)) {
        if (isdigit(b[1])) { b++; 
          if (isdigit(b[1])) { b++; 
            *a=(b[-2]-'0')*64+(b[-1]-'0')*8+(b[0]-'0'); }
          else {
            *a=(b[-1]-'0')*8+(b[0]-'0'); } }
        else 
          *a=(b[0]-'0'); } }
    else 
      *a=*b; 
    (*len)++; }
  a[1]=0;
  
  return s;
}

int main(int narg, char **arg) /*************************************** main */
{
  int i,n=0,iarg;
  int c,filter;
  char *from,*to;
  int lfrom,lto,argfrom;

  if (narg<3) {
    fprintf(stderr,"\
Replace strings in files, binary. Call by:\n\
  binrepl FROM TO [OPTIONS] [-] < INPUT > OUTPUT\n\
  binrepl FROM TO [OPTIONS] FILE [FILE ...]\n\
FROM,TO = STRING of characters\n\
OPTIONS:\n\
  -i  asks for confirmation:\n\
      - not if filter\n\
      - not good if control characters changed\n\
  -eX treat X as escape character in STRING (-e = none); coding:\n\
      XNNN = octal (max. 3 digits)\n\
      Xn = LF, Xr = CR, Xt = TAB, XX = X\n\
Example: (the same as unix to dos but the last ^z)\n\
  binrepl /n /r/n -e/ file.ext\n\
Example: replace all LF to LFLF\n\
  binrepl /n /n/n -e/ file.ext\n\
See also:\n\
  repl texrepl lemon\n");
    exit(0); }

  argfrom=0;

  yeskey='A';

  loop (iarg,3,narg) 
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'i': yeskey=0; break;
        case 'e': excl=arg[iarg][2]; break; 
        case 0: break;
        default: Error("binrepl: bad option"); }
    else 
      argfrom=iarg; 

  from=str(arg[1],&lfrom);
  to=str(arg[2],&lto);
//  printf("%s %d  %s %d\n",from,lfrom, to,lto);
  
  filter=argfrom==0;
  if (filter) yeskey='A',argfrom=narg-1;
  loop (iarg,argfrom,narg) {
    FILE *in =filter?stdin: fopen(arg[iarg],"rt");
    FILE *out=filter?stdout:fopen("binrepl.aux","wt");
    int repl=0,j;

    if (yeskey!='A') yeskey=0;
    extarg=arg[iarg];

    if (!filter) fprintf(stderr,"%s: ",arg[iarg]);

    n=0;
    while ( (c=getc(in))>0 ) {
      if (n>=N) Error("binrepl: buffer overflow - increase N and recompile");
      file[n++]=c; }

    for (i=0; i<n; )
      if (!memcmp(file+i,from,lfrom) && i<=n-lfrom && yes(i)) {
	repl++;
	for (j=0; j<lto; j++) putc(to[j],out); 
	i+=lfrom; }
      else {
	putc(file[i],out);
	i++; }
    fclose(in);
    fclose(out);
    if (!filter) {
      if (repl) {
        fprintf(stderr,"<<< %d replaced >>>\n",repl);
        rename(arg[iarg],string("%s.b~",arg[iarg]));
        rename("binrepl.aux",arg[iarg]); }
      else 
        fprintf(stderr,"unchanged\n"); } }

  return 0;
}
