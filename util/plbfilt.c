/* cc -Wall -O2 -o plbfilt plbfilt.c -lm
*/
#include "../gen/include.h"
typedef double real;
#include "../gen/vector3d.h"

struct site_s {
  int inc;
  char *id; } *mol1,*mol2;

int ns1,ns2;
int rev;

#define Endian(X) endian((char*)((void*)(X)))

void endian(char *a) /********************************************** endian */
{
  char x;
  if (!rev) return;

  x=a[0]; a[0]=a[3]; a[3]=x;
  x=a[1]; a[1]=a[2]; a[2]=x;
}

int readMOL(char *fn,struct site_s **mol) /************************* readMOL */
{
  char line[256];
  int ignorenum=fn[0]=='-';
  FILE *file=fopen(fn+ignorenum,"rt");
  int ns=0,i;

  char *tok;
  if (!file) Error(fn);

  printf("reading %s\n",fn);

  for (;;) {
    if (!fgets(line,256,file)) Error("plbfilt: open");
    if ( (tok=strtok(line," =\t\n")) )
      if (!strcmp(tok,"number_of_atoms")) {
	tok=strtok(NULL," =\t\n");
	ns=atoi(tok);
	printf("%s: number_of_atoms = %i\n",fn,ns);
	break; } }

  if (ns<=0) Error("plbfilt: undefined or bad number_of_atoms");
  alloczero(*mol,ns*sizeof(struct site_s));

  for (;;) {
    if (!fgets(line,256,file)) Error("plbfilt: format");
    if ( (tok=strtok(line," =\t\n")) )
      if (!strcmp(tok,"atoms")) break; }

  loop (i,0,ns) {
    for (;;) {
      if (!fgets(line,256,file)) Error("plbfilt: too short");
      if (line[0]!='!') break; }
    tok=strtok(line," =\t\n");
    if (!tok) Error("plbfilt: format");
    if (!ignorenum && atoi(tok)!=i) Error("plbfilt: numbering");
    (*mol)[i].id=strdup(strtok(NULL," \n\t")); }

  return ns;
}

int main(int narg,char **arg) /**************************************** main */
{
  int i,j,f=0;
  float r[3];
  FILE *plb1,*plb2;
  int varL=0;

  if (narg<4) {
    fprintf(stderr,"Playback filter. Call by:\n\
  plbfilt [-]RICH.mol [-]RICH.plb [-]POOR.mol POOR.plb\n\
Reads playback file RICH.plb which corresponds to molecule RICH.mol,\n\
selects only sites contained in POOR.mol and creates POOR.plb\n\
POOR.mol must be a subset of RICH.mol!\n\
\"-\" in front of RICH.plb means reversed endian\n\
\"-\" in front of mol-files ignores wrong numbering\n\
    (so that files with manually removed lines can be used)\n\
Output plb format is the same as input (cf. plbconv)\n\
See also: plbmerge plbcut plbinfo plbasc plbbox plb2plb\n\
WARNING: this is unchecked variable-L version of filtplb\n");
    exit(0); }

  ns1=readMOL(arg[1],&mol1);
  ns2=readMOL(arg[3],&mol2);
  rev=arg[2][0]=='-';
  plb1=fopen(arg[2]+rev,"rb");
  if (!plb1) Error(arg[2]);
  plb2=fopen(arg[4],"wb");

  loop (j,0,ns2) {
    loop (i,0,ns1) {
      if (!strcmp(mol1[i].id,mol2[j].id)) {
	mol1[i].inc++;
	goto done; } }
    Error(mol2[j].id);
   done:; } 

  if (2!=fread(r,4,2,plb1)) Error("plbfilt: plb no header");
  Endian(r); Endian(r+1);
  if (r[1]<0) {
    if (r[1]<0) varL=1;
    else Error("plbfilt: wrong L in plb file"); }

  if (ns1!=r[0]) Error("plbfilt: wrong ns - endian?");
  r[0]=ns2;
  fwrite(r,4,2,plb2);

  for (f=0;;f++) {
    if (varL) {
      if (3!=fread(r,4,3,plb1)) goto end;
      Endian(r); Endian(r+1); Endian(r+2);
      fwrite(r,4,3,plb2); }
    loop (i,0,ns1) {
      if (3!=fread(r,4,3,plb1)) goto end;
      if (mol1[i].inc) {
	Endian(r); Endian(r+1); Endian(r+2);
	if (3!=fwrite(r,4,3,plb2)) Error("plbfilt: cannot write"); }
      }
    }
  
 end: 
  printf("%d frames\n",f);
  fclose(plb1);
  fclose(plb2);

  return 0;
}
