#include "ground.h"
#include "alloc.h"
typedef float vector[3];

#include "simopt.h"
#include "intrapot.h"

char fn[128];

angleparm_t a={1,1};
torsionparm_t t={1,{0}};
int measure;

int main(int narg,char **arg) /*************************************** main */
{
FILE *plb,*mol,*intc;
float hdr[2];
int from=1,to=-1,by=1,ns=0,i,ii,frame;
char *dot;
char line[128],dummy[64],type[64];
#define ATLEN 3
char (*atoms)[ATLEN];
vector *r;
static vector f[4];

if (narg<2) {
  fprintf(stderr,"\
Internal coordinates from playback files.  Call by:\n\
  intcoord SIMNAME [FROM TO BY]\n\
SIMNAME.plb will be read and files SIMNAME.1, SIMNAME.2 ... generated\n");
  exit(0); }

initscroll(0);

if (narg>2) from=atoi(arg[2]);
if (narg>3) to=atoi(arg[3]);
if (narg>4) by=atoi(arg[4]);

if (from<1) Error("bad FROM");
if (by<1) Error("bad BY");

strcpy(fn,arg[1]);
dot=fn+strlen(fn);

strcpy(dot,".plb");
plb=fopen(fn,"rb");
if (!plb) {
  prt_("WARNING %s not found",fn);
  strcpy(dot,".p00");
  prt(", trying %s",fn); 
  plb=fopen(fn,"rb"); }

if (!plb) Error("no playback file");
prt("playback file %s");
if (fread(hdr,sizeof(hdr),1,plb)!=1) Error("playback file too short");
ns=hdr[0];
printf("ns=%d hdr[1]=%f\n",ns,hdr[1]); 
if (ns<1 || ns>100000) Error("bad ns in playback file");
alloczero(atoms,ATLEN*ns);
alloc(r,sizeof(vector)*ns);

strcpy(dot,".mol");
mol=fopen(fn,"rt");
if (mol) {
  prt("reading %s",fn);
  i=0;
  while (fgets(line,128,mol)) 
    if (strstr(line,"number_of_atoms")) {
      strtok(line,"=\n ");
      i=atoi(strtok(NULL,"=\n "));
      if (i!=ns) Error("ns in plb and mol files inconsistent");
      break; } 
  if (!i) Error("no number_of_atoms in mol-file");
  while (fgets(line,128,mol)) 
    if (!memcmp(line,"atoms",5)) break; 
  while (fgets(line,128,mol)) if (line[0]!='!') break;
  loop (i,0,ns) 
    if (!fgets(line,128,mol)) Error("bad mol-file");
    else {
      sscanf(line,"%d%s%s",&ii,dummy,type);
      atoms[i][0]=type[0];
      if (!strcmp(type,"MFE") || !strcmp(type,"FE")) strcpy(atoms[i],"Fe"); }
  }

if (to<0) to=0x7fff7fff;
measure=2;

loopto (frame,0,to) {
  if (fread(r,sizeof(vector),ns,plb)!=ns) break;
  if (frame<from || (frame-from)%by) continue;
  sprintf(dot,".%d",frame);
  intc=fopen(fn,"wt");
  prt("writing %s",fn);
  loop (i,0,ns) {
    fprintf(intc,"%-2s",atoms[i]);
    if (i) fprintf(intc," %10.8f",
      sqrt(Sqr(r[i][0]-r[i-1][0])
          +Sqr(r[i][1]-r[i-1][1])
          +Sqr(r[i][2]-r[i-1][2])));
    if (i>1) {
      anglepot(r[i],r[i-1],r[i-2],f[0],f[1],f[2],&a);
      fprintf(intc," %10.6f",phi*(180/PI)); }
    if (i>2) {
      improperpot(r[i],r[i-1],r[i-2],r[i-3],f[0],f[1],f[2],f[3],&t);
      if (phi<0) phi+=PI*2;
      fprintf(intc," %10.6f",phi*(180/PI)); }
    fprintf(intc,"\n"); }
  fclose(intc); }

return 0;
}
