/* cc -O2 -o pdb2pdb pdb2pdb.c
  this filter reads a pdb file and writes it in a more standard order 
*/

#include "../gen/include.h"

int foratoi(char *c)
{
while (*c==' ') c++;
return atoi(c);
}

char (*line)[82];
int maxn=10000,n,i;
int head=1; /* header lines */
int tail=0; /* tail lines */
static int start=-1,end,resn=-1,newi=0,i0;
/* warning: extern end is sometimes prohibited symbol */

void selectid(char *id)
{
int j;
char si[6];

loop (j,i0,i) if ((id[0] && !memcmp(line[j]+13,id,3)) || (!id[0] && line[j][0]!='~')) {
  newi++;
  sprintf(si,"%5i",newi);
  memcpy(line[j]+6,si,5);
  fputs(line[j],stdout);
  line[j][0]='~'; }
}

void renameid(char *id,char *newid)
{
int j;

loop (j,i0,i) if (!memcmp(line[j]+13,id,3)) memcpy(line[j]+13,newid,3);
}

int main(int narg,char **arg)
{

fprintf(stderr,"Rearrange PDB to N-CA-C-O[-sidechain]. Call by:\n\
  pdb2pdb [MAXLINES] < inputpdb > outputpdb\n");

if (narg>1) maxn=atoi(arg[1]);
line=(void*)malloc(sizeof(line[0])*maxn);
if (!line) Error("no heap");

while (fgets(line[n],82,stdin)) {
  n++;
  if (n>=maxn) Error("too many lines");
  if (strlen(line[n])>81) Error("line >80 chars"); }

fprintf(stderr,"%d lines read\n",n);

loop (i,0,n) if (!memcmp("ATOM  ",line[i],6) || !memcmp("HETATM",line[i],6)) {
  end=i;
  if (start<0) start=i; }

if (start<0) Error("no ATOM or HETATM line");

loop (i,0,head) fputs(line[i],stdout);
  
if (head>start) Error("header (no ATOM) lines expected");

for (i=start; i<end; ) {
/*
ATOM    326  ND2 ASN    46      13.407   3.298  15.015  1.00 10.32      1CRN 395
0         1         2         3         4         5         6         7         8
*/
  i0=i;
  resn=foratoi(line[i]+21);
  do i++; while (resn==foratoi(line[i]+21));
  
  { int j; loop (j,i0,i) if (line[j][13]=='H') line[j][0]='~'; }
  selectid("N  ");
  selectid("CA ");
  selectid("C  ");
  if (end<i) {
    renameid("OC1","O  ");
    renameid("OC2","OXT"); }
  selectid("O  ");
  selectid("");
  }

fprintf(stderr,"%d lines written, last residue=%d\n",newi,resn);
return 0;
}
