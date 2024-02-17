/* cc -Wall -O2 -o pdb2atm pdb2atm.c
 */
#include "/home/jiri/macsimus/gen/include.h"
#include "/home/jiri/macsimus/gen/forscanf.c"

int main(int narg,char **arg)
{
  int i,n;
  char id[256],line[256];
  char *name,*c;
  char L[3][8],x[8],y[8],z[8],atm[8],dummy[16];
  FILE *in,*out;

  L[0][0]=0;

  if (narg<2) {
    fprintf(stderr,"Convert PDB to ATM (simple). Call by:\n\
  pdb2atm NAME[.pdb]\n\
NAME.pdb will be converted to NAME.atm\n\
");
    exit(0); }

  name=arg[1];
  if (strlen(name)>4 && !strcmp(strend(name)-4,".pdb")) {
    name=strdup(name);
    strend(name)[-4]=0; }

  //  if (narg>2) =atoi(arg[2]);

  c=string("%s.pdb",name);
  in=fopen(c,"rt");
  if (!in) Error(c);

  n=0;
  while (fgets(line,256,in)) {
    if (!memcmp(line,"CRYST",5)) { // box size, only rectangular
      /*
CRYST1   47.096   47.096  450.000  90.00  90.00  90.00 P 1           1*/
      if (4!=sforscanf(line,"%6s%9s%9s%9s",x,L[0],L[1],L[2])) Error("CRYST"); }
    if (!memcmp(line,"ATOM ",5) || !memcmp(line,"HETATM",6)) n++; }

  rewind(in);

  c=string("%s.atm",name);
  out=fopen(c,"wt");
  if (!out) Error(c);

  fprintf(out,"%d\n",n);
  if (L[0][0])
    fprintf(out,"  %s %s %s box\n",L[0],L[1],L[2]);
  else
    fprintf(out,"\n");

  i=0;
  while (fgets(line,256,in)) {
    if (!memcmp(line,"ATOM ",5) || !memcmp(line,"HETATM",6)) {
      i++;
      if (i>n) Error("PDB changed as I read it");
      /*
Dummy_Dummy_At__DummDummy_____X_______Y_______Z_______
ATOM  31910  HW2 SOL X7002      18.690  24.440 237.220  0.00  0.00
ATOM    324  CG  ASN    46      12.538   4.304  14.922  1.00  7.98      1CRN 393
HETATM    1  O           1      -2.411   6.250  -4.773
      */

      sforscanf(line,"%6a%6a %4s%4a%10a%8s%8s%8s",
                      dummy,
                         dummy,
                             atm,
                                dummy,
                                   dummy,
                                       x, y, z);
      strcpy(id,atm);
      atm[1]=tolower(atm[1]);
      for (c=atm+1; *c; c++) if (!isalpha(*c)) *c=0;
      fprintf(out,"%-2s %7s %7s %7s\n",atm,x,y,z);
    }
  }

  return 0;
}
