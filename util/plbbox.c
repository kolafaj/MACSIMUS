/* cc -o plbbox -O2 plbbox.c
*/

#include "../gen/include.h"

void OP(float *x,char op,double nr) /************************************ OP */
{
  switch (op) {
    case '=': *x=nr; break;
    case '+': *x+=nr; break;
    case '-': *x-=nr; break;
    case 'x':
    case '*': *x*=nr; break;
    case '/': *x/=nr; break;
    case '<': if (*x>nr) *x=nr; break;
    case '>': if (*x<nr) *x=nr; break; }
}

int main(int narg,char **arg) /**************************************** main */
{
  FILE *f;
  float header[2];
  float L[3];
  char op[3]="...";
  double nr[3]={0,0,0};
  int i;

  if (narg<2) {
    fprintf(stderr,"\
Change box size L of a plb-file. Call by:\n\
  plbbox PLB-FILE { KNUMBER KNUMBER KNUMBER | KNUMBER }\n\
where\n\
  K = \n\
    = L is replaced by NUMBER\n\
    + NUMBER added to L\n\
    - NUMBER subtracted from L\n\
    * L multiplied by NUMBER (also x)\n\
    / L divided by NUMBER\n\
    > if L>NUMBER then unchanged, otherwise replaced by NUMBER\n\
    < if L<NUMBER then unchanged, otherwise replaced by NUMBER\n\
    . no change (NUMBER may be missing)\n\
Three parameters KNUMBER: Lx,Ly,Lz are changed\n\
One parameter KNUMBER: all box sides are changed in the same way\n\
Only the first frame is changed for new (variable L) format\n\
The file is changed in place (NO BACKUP!)\n\
Example (add 1 AA to Lz):\n\
  plbbox slab.plb . . +1\n\
See also:\n\
  plbscale plbstack plbinfo plbcut plb2asc asc2plb plb2plb plbconv\n");
    exit(0); }

  loop (i,0,3) if (i+2<narg) {
    op[i]=arg[i+2][0];
    if (!strchr("+-*x/=.<>",op[i])) Error("plbbox: bad operator, use one of +-*x/=.<>");
    nr[i]=atof(arg[i+2]+1); }

  f=fopen(arg[1],"r+");
  if (!f) Error("plbbox: cannot open plb-file");
  
  if (fread(header,4,2,f)!=2) Error("plbbox: too short plb-file");
  
  if (header[1]<0) {
    printf("%s: ns=%.0f L=",arg[1],header[0]);
    if (fread(L,4,3,f)!=3) Error("plbbox too short plb-file (no L)");
    printf("%.8g %.8g %.8g\n",L[0],L[1],L[2]);
    if (narg==4) Error("plbbox: either one LcubeKEY argument or three LxKEY,LyKEY,LzKEY allowed");
    else if (narg==3) op[1]=op[2]=op[0],nr[1]=nr[2]=nr[0];
    loop (i,0,3) OP(L+i,op[i],nr[i]);

    printf("   --> %.8g %.8g %.8g\n",L[0],L[1],L[2]);
    rewind(f);
    if (fwrite(header,4,2,f)!=2) Error("plbbox: cannot write header");
    if (fwrite(L,4,3,f)!=3) Error("plbbox: cannot write L");
    fclose(f); }
  else {
    printf("%s: ns=%.0f old format L=%.8g",arg[1],header[0],header[1]);
    if (narg>3) Error("plbbox: only one LcubeKEY argument expected");
    OP(header+1,op[0],nr[0]);
    printf(" --> %.8g\n",header[1]);
    rewind(f);
    if (fwrite(header,4,2,f)!=2) Error("plbbox: cannot write");
    fclose(f); }

  return 0;
}

