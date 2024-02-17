/* cc -o plbscale -O2 plbscale.c -lm
*/

#include "../gen/include.h"

void OP(float *x,char op,double nr)
{
  switch (op) {
    case '=': *x=nr; break;
    case '+': *x+=nr; break;
    case '-': *x-=nr; break;
    case 'x':
    case '*': *x*=nr; break;
    case 'e': *x*=exp(nr); break;
    case '/': *x/=nr; break;
    case '<': if (*x>nr) *x=nr; break;
    case '>': if (*x<nr) *x=nr; break; }
}

typedef float vector[3];

int main(int narg,char **arg)
{
  FILE *f1,*f2;
  int ns,i,k,n=0;
  static float header[2];
  vector r;
  vector L,oldL;
  char op[3]="...";
  double nr[3]={0,0,0};

  if (narg<4) {
    fprintf(stderr,"\
Rescale plb-file (all coordinates and L). Call by:\n\
  plbscale {INPUT_FILE|-} {OUTPUT_FILE|-} LxKEY [LyKEY [LzKEY]]\n\
where\n\
  LxKEY,LyKEY,LzKEY = [=+-*/<>.]NUMBER\n\
    = L is replaced by NUMBER\n\
    + NUMBER added to L\n\
    - NUMBER subtracted from L\n\
    * L multiplied by NUMBER (also x)\n\
    e L multiplied by exp(NUMBER)\n\
    / L divided by NUMBER\n\
    > if L>NUMBER then unchanged, otherwise replaced by NUMBER\n\
    < if L<NUMBER then unchanged, otherwise replaced by NUMBER\n\
    . no change (NUMBER may be missing)\n\
  the configuration is rescaled accordingly\n\
  if LyKEY(LzKEY) is missing, it is the same as LxKEY(LyKEY)\n\
  not suitable for free b.c.\n\
Example:\n\
  plbscale small.plb larger.plb '*1.1'\n\
See also:\n\
  plbbox plbstack plbtran plbrot\n");
    exit(0); }

  loop (i,3,6) {
    int ii=min(i,narg-1);

    op[i-3]=arg[ii][0];
    if (!strchr("e+-*x=.<>",op[i-3])) Error("plbscale: bad operator, use one of e+-*x/=.<>");
    nr[i-3]=atof(arg[ii]+1); }

  if (strcmp(arg[1],"-"))
    if (!strcmp(arg[1],arg[2])) Error("plbscale: INPUT_FILE = OUTPUT_FILE");


  if (strcmp(arg[1],"-")) f1=fopen(arg[1],"rb");
  else f1=stdin;
  if (f1==NULL) Error("plbscale: no INPUT_FILE");

  if (fread(header,sizeof(header),1,f1)!=1) Error("file too short");
  ns=header[0];

  if (header[1]!=-3)
    Error("plbscale: wrong format: only variable-L plb-file supported (use plbconv)");

  printf("plbscale: %s: of sites=%d  size of 1 frame=%d\n",
	 arg[1],ns,(int)(ns*sizeof(vector)));
  if ((header[0]-ns)!=0 || ns<1 || header[0]>16777216.) Error("plbscale: bad header");

  if (strcmp(arg[2],"-")) f2=fopen(arg[2],"wb");
  else f2=stdout;
  if (f2==NULL) Error("plbscale: cannot open OUTPUT_FILE");

  fwrite(header,sizeof(header),1,f2);

  for (;;) {
    if (1!=fread(L,sizeof(L),1,f1)) break;
    if (L[0]*L[1]*L[2]==0) Error("plbscale: box size must not be zero");
    n++;
    loop (i,0,3) {
      oldL[i]=L[i];
      OP(L+i,op[i],nr[i]); }
    fwrite(L,sizeof(L),1,f2);
    printf("frame %d: scaling %g %g %g times\n",
               n,         L[0]/oldL[0],L[1]/oldL[1],L[2]/oldL[2]);
    loop (i,0,ns) {
      if (1!=fread(r,sizeof(r),1,f1)) Error("plbscale: frame truncated");
      loop (k,0,3) r[k]*=L[k]/oldL[k];
      fwrite(r,sizeof(r),1,f2); } }

  fclose(f2);
  fclose(f1);

  return 0;
}
