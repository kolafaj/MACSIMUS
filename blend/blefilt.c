/* cc -o blefilt blefilt.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int narg,char **arg)
{
char *ch,r[128];
int iarg;
int ic,ncol,col[16];
char field[16][64];

if (narg<2) {
  fprintf(stderr,"Blend-file filter. Call by:\n\
  blefilt TABLE [COL [COL2 ...]] [TABLE2 ... ] < INPUT.ble > OUTPUT\n\
where TABLE is one of { sites bonds angles dihedrals impropers aromatics }\n\
and optional list of columns follows\n\
The default column is the column with the value of bond length or angle\n\
or x y z for sites\n\
The output is in the order of the blend-file\n\
Examples:\n\
  blefilt angles dihedrals < cyto.ble > cyto.int\n\
  blefilt sites < cyto.ble > cyto.xyz\n\
  blefilt bonds 2 4 7 < cyto.ble > cyto.int\n");
  exit(0); }

while (fgets(r,128,stdin)) {
  if (r[0]=='!') continue;
  if ( (ch=strchr(r,'\n')) ) *ch=0;

  for (iarg=1; iarg<narg; iarg++) {
    if (!strchr("0123456789",arg[iarg][0]) && !strcmp(arg[iarg],r)) {
      ncol=0;
      for (ic=iarg+1; ic<narg; ic++)
        if (strchr("0123456789",arg[ic][0])) col[ncol++]=atoi(arg[ic]);
      if (ncol==0) {
        ncol=1;
        col[0]=1;
        if (!strcmp(r,"sites")) { col[0]=4; col[1]=5; col[2]=6; ncol=3; }
        if (!strcmp(r,"bonds")) col[0]=7;
        if (!strcmp(r,"angles")) col[0]=9;
        if (!strcmp(r,"dihedrals")) col[0]=12;
        if (!strcmp(r,"impropers")) col[0]=12;
        if (!strcmp(r,"aromatics")) col[0]=12; }

      while (fgets(r,128,stdin)) {
        if (r[0]=='!') continue;
        if ( (ch=strchr(r,'\n')) ) *ch=0;
        if (r[0]==0) break;
        sscanf(r,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
          field[0], field[1], field[2], field[3],
          field[4], field[5], field[6], field[7],
          field[8], field[9], field[10], field[11],
          field[12], field[13], field[14], field[15]);
        for (ic=0; ic<ncol; ic++) printf("%s%c",field[col[ic]-1]," \n"[ic==ncol-1]); }
      }
    }
  }
return 0;
}
